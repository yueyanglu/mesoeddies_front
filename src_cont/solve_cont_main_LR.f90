program solve_cont_main
    !
    ! Compiling and linking in one step:
    ! ifort -o exe test.f90
    !
  use netcdf
  use mod_nc_wr_rd
  use mod_grid_my, only : ocean_grid_type, build_grid_cartesian
  use mod_continuity

  implicit none

  ! 
  character*200 :: dir_exp, dir_inith, flx_dir, h_dir, flnm_uv, flnm_h, flnmuv_ful, IC_h_fnm
  character*10  :: funm, fvnm ! name of flux components in the .nc file
  character*200 :: dir_save, flnm_o, flnm_ful_o
  logical       :: ifexist, ifsave, ifsave_last
  integer       :: ncid

  ! times 
  integer :: time_s(3), time_e(3) ! Array [yr, day of year, hr],
                                ! consistent w/ how nc files are saved
  real    :: dt_arch   ! intervals between archive input [hr]
  real    :: dt_save   ! intervals between the output [hr]
  real    :: dt_in    ! actual time step for h [hr]
  real    :: dt        ! time step for continuity [s]
  real    :: dt_tend ! [sec] for tendency
  real    :: hr_elapsed, dy_elapsed ! hr/dy that elapsed wrt the initial time
  integer :: narch     ! # of archieve to be read
  integer :: nintrp     ! # of time steps of interp within 1 archive interval 
  integer :: iar, istep, i, j, k
  integer :: iyear1, iday1, ihr1, iyear15, iday15, ihr15, iyear2, iday2, ihr2 &
             ,iyear_last, iday_last, ihr_last
  real    :: q1, q2 ! weight of interp
  real, parameter :: yr2dy = 365., dy2hr = 24., hr2sc = 3600.
  real :: missing = 1.0E20 ! missing values in nc file

  ! dimensions
  integer, parameter :: NDIMS = 4 ! # of dim of vars in file
  integer :: nih, njh, niq, njq, nrecs
  integer, parameter :: nihalo = 64, njhalo = 64, nz = 3

  ! ------ model geometry & grid -------
  character*200 :: geom_fnm, geom_fnm_ful, vert_fnm, vert_fnm_ful
  type(ocean_grid_type) :: g !< The horizontal grid type

  ! ------ Key variables (data domain w/ halos) -------
  ! layer thickness solution
  real, allocatable, dimension(:,:,:) :: h, h_init
  real, allocatable, dimension(:,:,:) :: h1, h2, h_n1_last, dhdt ! for calc tendency term
  ! snapshot read from nc file
  real, allocatable, dimension(:,:,:) :: u, v, h_mod, uh, vh
  real, allocatable, dimension(:,:,:) :: u1, v1 ! t1
  real, allocatable, dimension(:,:,:) :: u2, v2 ! t2

  ! ------ filenames -------
  !dir_exp = '/glade/scratch/yueyanglu/smooth_33/'
  IC_h_fnm = '/glade/work/yueyanglu/MOM6_OUT/forc_uvh_64/sol_h_180d/h_snap__0021_001_00.nc'
  !dir_exp = '/glade/work/yueyanglu/MOM6_OUT/tr_off_64/'
  flnm_h = 'h_snap__YEAR_DAY_HR.nc' !h_snap__YEAR_DAY_HR
  ! 
  flx_dir = '/glade/work/yueyanglu/MOM6_OUT/forc_uvh_64/uvhm_CS_decomp_180d/'
  flnm_uv = 'uvh_mean__YEAR_DAY_HR.nc'
  funm = "uh";  fvnm = "vh"
  ! 
  geom_fnm = 'ocean_geometry.nc' 
  vert_fnm = 'Vertical_coordinate.nc'
  ! 
  dir_save = '/glade/work/yueyanglu/MOM6_OUT/forc_uvh_64/sol_h_180d/'
  flnm_o = 'h_snap__YEAR_DAY_HR.nc'

! check save dir
  inquire(directory = dir_save, exist = ifexist)
  if ( .NOT. ifexist ) then
    call system('mkdir ' // trim(dir_save)) 
    print*, 'OUTPUT will be saved to: ', trim(dir_save)
  endif

  !============================================================
  !                  read model geometry 
  !============================================================
  ! geom_fnm_ful = trim(dir_exp)//trim(geom_fnm)
  ! vert_fnm_ful = trim(dir_exp)//trim(vert_fnm)
  ! get dimensions (not necessary)
  ! call get_dims_geom(geom_fnm_ful, nih, njh, niq, njq)
  ! call get_dims_vert(vert_fnm_ful, nz)
  ! get model grids
  call build_grid_cartesian(g, 3840., 3840., nihalo, njhalo, nz)
  ! call build_grid(g, geom_fnm_ful, vert_fnm_ful)
  ! print*, 'Read model geom from:', trim(geom_fnm_ful)
  print*, 'g extents:      ', g%isg,g%ieg,g%jsg,g%jeg
  print*, 'gB-sym extents: ', g%isgb,g%iegb,g%jsgb,g%jegb
  print*, 'c extents:      ', g%isc,g%iec,g%jsc,g%jec
  print*, 'cB-sym extents: ', g%iscb,g%iecb,g%jscb,g%jecb
  print*, 'd extents:      ', g%isd,g%ied,g%jsd,g%jed
  print*, 'dB-sym extents: ', g%isdb,g%iedb,g%jsdb,g%jedb
  print*, 'halos (i,j):  ', g%nihalo, g%njhalo
  print*, 'g%idg_offset=', g%idg_offset, 'g%jdg_offset=',g%jdg_offset
  print*, 'g%ke:  ', g%ke
  print*, 'g%bathyT:', g%bathyT(1:6,20)
  print*, 'g%dy_Cu(0:6,)', g%dy_Cu(0:6,20)
  print*, 'g%dyCu(0:6,)', g%dyCu(0:6,20)
  print*, 'g%areaCu(0:6,)', g%areaCu(0:6,20)

  ! allocate main vars
  allocate( h(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
  allocate( u(g%isdb:g%iedb, g%jsd:g%jed,  g%ke) )
  allocate( uh(g%isdb:g%iedb, g%jsd:g%jed,  g%ke) )
  allocate( v(g%isd:g%ied,   g%jsdb:g%jedb, g%ke) )
  allocate( vh(g%isd:g%ied,   g%jsdb:g%jedb, g%ke) )
  allocate( h_mod(g%isd:g%ied, g%jsd:g%jed, g%ke) )
  allocate( h_init(g%isd:g%ied, g%jsd:g%jed, g%ke) )
  allocate( h1(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
  allocate( h_n1_last(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
  allocate( h2(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
  allocate( dhdt(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
  ! 
  allocate( u1(g%isdb:g%iedb, g%jsd:g%jed,  g%ke) )
  allocate( v1(g%isd:g%ied,   g%jsdb:g%jedb, g%ke) )
  allocate( u2(g%isdb:g%iedb, g%jsd:g%jed,  g%ke) )
  allocate( v2(g%isd:g%ied,   g%jsdb:g%jedb, g%ke) )
    ! 
  print*, 'Shape of h: ', shape(h)
  print*, 'Shape of u: ', shape(u)
  print*, 'Shape of v: ', shape(v)

  !============================================================
  !                  set the times  
  !============================================================
  time_s = (/ 21, 1, 0 /)
  time_e = (/ 30, 365, 0 /)
  dt_arch = 6.  ! [hr]
  dt_save = 6. ! [hr]
  dt_in = 300. ! [s] 300 or 50
! # of archives to be read
  narch = nint( (time_e(1) - time_s(1))*yr2dy*(dy2hr/dt_arch) + &
      (time_e(2) - time_s(2))*(dy2hr/dt_arch) + &
      (time_e(3) - time_s(3))/dt_arch + 1. );
! # of time steps of interp within one archive interval 
  nintrp = nint(dt_arch*hr2sc / dt_in)
  ! 
  print*, 'time_s=',time_s,', time_e=',time_e, &
      ', dt_arch=',dt_arch,', narch=', narch,'dt_in [s] =', dt_in, 'nintrp=',nintrp

  !============================================================
  !            initial layer thickness; read
  !============================================================
  ! NaN in nc is read as extreme values, e.g., 1e20
  print*, NEW_LINE('A')//"Read initial fields from:", trim(IC_h_fnm)
  call get_ncvar(IC_h_fnm, "h", h, g%nihalo, g%njhalo)
  ! print*, u(0+g%nihalo,555+g%njhalo,1), u(1+g%nihalo,555+g%njhalo,1), &
  !  u(1024+g%nihalo,555+g%njhalo,1), u(1025+g%nihalo,555+g%njhalo,1) ! matlab: u(6+1,7,1)
  ! print*, "h_sol_max=",  maxval(h(:,:,1)),  minval(h(:,:,1))

  !============================================================
  !            main loop  
  !============================================================
  ifsave = 0

  do iar = 1, narch ! t1

    ifsave_last = ifsave ! last archv step
    !
    ifsave = 0
    if (mod(iar*dt_arch, dt_save) .eq. 0.0)  ifsave = 1

    !----------  read u & v averaged over t1~t2 ----------
    hr_elapsed  = (iar-1)*dt_arch 
    call get_timeint(iyear1, iday1, ihr1, time_s, hr_elapsed, yr2dy, dy2hr)
    hr_elapsed  = iar*dt_arch 
    call get_timeint(iyear2, iday2, ihr2, time_s, hr_elapsed, yr2dy, dy2hr)
    hr_elapsed  = (iar-1)*dt_arch + .5*dt_arch
    call get_timeint(iyear15, iday15, ihr15, time_s, hr_elapsed, yr2dy, dy2hr)

    ! 
    call wr_flnm_time(flnm_uv, iyear15, iday15, ihr15)
    flnmuv_ful = trim(flx_dir)//trim(flnm_uv)
    print*, NEW_LINE('A')//"Read uv-t1 from:", trim(flnmuv_ful)
    ! missing values will be fill as 0's
    call get_ncvar(flnmuv_ful, funm, uh, g%nihalo, g%njhalo)
    call get_ncvar(flnmuv_ful, fvnm, vh, g%nihalo, g%njhalo)
    ! print*, 'maxval(u) = ', maxval(u), 'minval(u) = ', minval(u)
    ! print*, u(0+g%nihalo,33+g%njhalo,1), u(1+g%nihalo,33+g%njhalo,1)
    ! print*, 'maxval(v) = ', maxval(v), 'minval(v) = ', minval(v)

    !----------  time step h ----------
    do istep = 1, nintrp

      print*, "Updating layer thickness at istep=", istep 
      call continuity_useuh(h, h, uh, vh, dt_in, g)

      if (ifsave .and. istep == nintrp-1) then
        h1 = h
        print*, 'record h1'
      endif

      if (ifsave_last .and. istep == 1) then
        h2 = h
        print*, 'record h2' 
        ! save (h2-h1)/2dt 
        dt_tend = 2.0*dt_in
        dhdt = (h2 - h_n1_last) / dt_tend
        call wr_flnm_time(flnm_o, iyear1, iday1, ihr1)
        flnm_ful_o = trim(dir_save)//trim(flnm_o)
        call wr_ncfile(flnm_ful_o, "dhdt", dhdt, g%nihalo, g%njhalo)
        print*, "dhdt saved to: ", trim(flnm_ful_o)
      endif

    enddo

    h_n1_last = h1
    print*, "h_sol_max=", maxval(h(g%isc:g%iec,g%jsc:g%jec,1)),  minval(h(g%isc:g%iec,g%jsc:g%jec,1))

    !----------  save h at t2 ----------
    if (ifsave) then 
      ! 
      call wr_flnm_time(flnm_o, iyear2, iday2, ihr2)
      flnm_ful_o = trim(dir_save)//trim(flnm_o)
      !
      call wr_ncfile(flnm_ful_o, "h", h, g%nihalo, g%njhalo)
      print*, "h saved to: ", trim(flnm_ful_o)
    endif

  enddo ! iar

contains

  ! update tracer conc by the point-wise forcing 
  !only works for POSITIVE 'hr_elap' 
  subroutine get_timeint(iyear, iday, ihr, time_s, hr_elap, yr2dy, dy2hr)
    integer, intent(out) :: iyear, iday, ihr
    integer, intent(in)  :: time_s(3)
    real,    intent(in)  :: hr_elap
    real,    optional    :: yr2dy, dy2hr
    !
    real :: dy_elap
    !
    if(.not. present(yr2dy))  yr2dy = 365. 
    if(.not. present(dy2hr))  dy2hr = 24. 

    !------- 1. check input 'time_s', make sure it is within 365 dys and 24 hr
    if (time_s(3) > dy2hr .OR. time_s(2) > yr2dy) then 
      print*, 'ERROR: time_s must be within hr/dy or dy/yr limits !!!'
      RETURN
    endif

    !------- 2. 
    dy_elap = hr_elap / dy2hr
    iyear = time_s(1) + int( dy_elap/yr2dy )
    iday  = time_s(2) + mod( dy_elap, yr2dy ) 
    ihr   = time_s(3) + mod( hr_elap, dy2hr )
    ! if hr >= 24 (ihr will not exceed 48 due to eqn above)
    if (ihr >= dy2hr) then 
      ihr = ihr - dy2hr
      iday = iday + 1 
    endif
    ! if day > 365 (iday will not exceed 730 due to eqn above)
    if (iday > yr2dy) then
      iday = iday - yr2dy
      iyear = iyear + 1 
    endif
  end subroutine get_timeint

  subroutine wr_flnm_time(flnm, iyear, iday, ihr)
    character*200, intent(inout) :: flnm
    integer, intent(in) :: iyear, iday, ihr
    integer :: l

    l = len_trim(flnm)
    write(flnm(l-4:l-3),'(i2.2)') ihr
    write(flnm(l-8:l-6),'(i3.3)') iday
    write(flnm(l-13:l-10),'(i4.4)') iyear
  end subroutine wr_flnm_time 

end ! program


