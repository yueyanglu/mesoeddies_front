program solve_trac_main
    !
    ! Compiling and linking in one step:
    ! ifort -o exe test.f90
    ! 
    ! Two relaxations added;
    ! 1. toward reference (optional): after ADV but before params
    ! 2. toward passive temperature: after all steps
    ! 
  use netcdf
  use mod_nc_wr_rd
  use mod_grid_my, only : ocean_grid_type, build_grid_cartesian
  use mom_tracer_registry, only : tracer_registry_type, call_tracer_register, initialize_tracer
  use mom_offline_main, only: offline_advection_layer
  use mom_tracer_hor_diff
  use mom_tracer_aux

  implicit none

  ! 
  character*200 :: dir_exp, IC_tr_fnm, flx_dir, h_dir, eforc_dir, rel_dir, &
      flnm_uv, flnm_h, flnm_forc, flnmuv_ful, flnmh_ful, flnmforc_ful, &
      flnm_rel, flnmrel_ful, flnmrel_pt_ful, flnm, &
      advStr, forcStr, relStr, parmStr, flnm_params, &
      cStr, k_dir, flnm_k, flnm_k_ful
  character*20  :: funm, fvnm, forcnm, relnm ! name of flux components in the .nc file
  character*200 :: dir_save, flnm_o, flnm_ful_o
  character*200 :: reldir_save, fnm_relsv, fnm_relsv_ful, relforcname
  character*200 :: dir_kforc, fnm_kforc, fnm_kforc_ful, forcname_k, forcname_ka, forcname_chi
  character*5   :: trname = 'tr#' ! e.g. "tr1"
  logical       :: ifexist
  integer       :: ncid, l
  integer       :: clock_rate, clock_start, clock_stop ! running time

  ! times 
  integer :: time_s(3), time_e(3) ! Array [yr, day of year, hr],
                                ! consistent w/ how nc files are saved
  real    :: dt_arch   ! intervals between archive input [hr]
  real    :: dt_save   ! intervals between the output [hr]
  real    :: dt_in    ! actual time step for h [hr]
  real    :: dt_rel   ! intervals for enforcing the relaxation [hr]
  real    :: dt        ! time step for continuity [s]
  real    :: hr_elapsed, dy_elapsed ! hr/dy that elapsed wrt the initial time
  integer :: narch     ! # of archieve to be read
  integer :: nintrp     ! # of time steps of interp within 1 archive interval 
  integer :: iyear1, iday1, ihr1, iyear15, iday15, ihr15, iyear2, iday2, ihr2 
  integer :: iar, istep, i, j, k, m
  real, parameter :: yr2dy = 365., dy2hr = 24., hr2sc = 3600.
  real    :: missing = 1.0E20 ! missing values in nc file
  logical :: ifadv, ifforc, iffct, ifrelax ! if update by pre-diagnosed tracer forcing
  real    :: T_relx, T_relx_pt, r_relx_pt ! relaxation time scale [s]
  integer :: paramflg 

  ! dimensions
  integer, parameter :: NDIMS = 4 ! # of dim of vars in file
  integer :: nih, njh, niq, njq, nrecs
  integer, parameter :: nihalo = 64, njhalo = 64, nz = 3
  real, parameter :: khtr = 100.0 ! m2/s
  real :: fac ! factors upon readed flux 

  ! tracer initial
  integer, parameter :: NTR = 1
  integer, parameter :: carries(NTR) = [1] ! which tracers, eg. 'tr2, tr4,...'
  ! [1, 2, 3, 4, 5, 6, 7, 8]
  character(len=5)   :: cnames(NTR) ! names of tracers to be solved

  real, dimension(nz, NTR)  :: cm_prev, cm_afte, cm1, cm2, cm3, q ! tracer-mass in each layer
  
  ! ------ model geometry & grid -------
  character*200 :: geom_fnm, geom_fnm_ful, vert_fnm, vert_fnm_ful
  type(ocean_grid_type) :: g !< The horizontal grid type
  type(tracer_registry_type) :: tr_Reg, tr_Reg0
  type(tracer_hor_diff_cs), pointer :: tr_diff_CSp => null() 

  ! ------ Key variables (data domain w/ halos) -------
  ! layer thickness solution
  real, allocatable, dimension(:,:,:) :: h, h_t2, h_cp &
      , trTemp
  real, allocatable, dimension(:,:,:) :: dp1, dp2
  ! snapshot read from nc file
  real, allocatable, dimension(:,:,:) :: uhtr, vhtr, h_mod
  real, allocatable, dimension(:,:,:) :: u1, v1 ! t1
  real, allocatable, dimension(:,:,:) :: u2, v2 ! t2
  real, allocatable, dimension(:,:,:) :: uh_cp, vh_cp
  real, allocatable, dimension(:,:,:,:) :: forc, relxforc, relx, relx_pt ! [ijk-ntr]
  real, allocatable, dimension(:,:,:) :: kxx, kxy,  kyx, kyy
  real, allocatable, dimension(:,:,:) :: ka, chiu, chiv
  real, allocatable, dimension(:,:,:,:) :: forc_ka, forc_chi, forc_kchi
  real :: q1, q2

  ! ======================== end allocate ========================
  ! ------ disp some tracer
  do m=1,NTR
    write(trname(3:3),'(i1)') carries(m)
    cnames(m) = trname
  enddo
  print*,NEW_LINE('A'),'This run will do TRACERS: ', cnames, NEW_LINE('A')

  ! ------ filenames -------
  ! ----- exp 1
  IC_tr_fnm = '/glade/work/yueyanglu/MOM6_OUT/tr_off_64/ICs/trac_init_pt_cs.nc'
  flx_dir = '/glade/work/yueyanglu/MOM6_OUT/forc_uvh_64/uvhm_CS_decomp/'  
  flnm_uv = 'uvh_mean__YEAR_DAY_HR.nc'; funm = "uh"; fvnm = "vh"; 
  h_dir = '/glade/work/yueyanglu/MOM6_OUT/forc_uvh_64/sol_h/' ! solh
  ! eforc_dir = '/glade/work/yueyanglu/MOM6_OUT/tr_off_64/cforc_pt_uvhm_addreld10/'
  eforc_dir = '/glade/work/yueyanglu/MOM6_OUT/tr_off_64/cforc_pt_fromKL_C0102040507/'
  flnm_forc = 'forc__YEAR_DAY_HR.nc'; forcnm = 'cforc#' ! cforc_m_rel# cforc_m#
  rel_dir = '/glade/work/yueyanglu/MOM6_OUT/tr_off_64/tr_pt_onl_CS/'
  flnm_rel = 'tr__YEAR_DAY_HR.nc'; relnm = 'tr#'
  flnmrel_pt_ful = '/glade/work/yueyanglu/MOM6_OUT/tr_off_64/ICs/trac_init_pt_cs.nc'
  ifadv = 1;  ifforc = 1; ifrelax = 0 ! ifrelax controls the extra relax!
  paramflg = 0
  cStr = '0102030507'
  flnm_k = 'K_C'//trim(cStr)//'_YEAR_DAY_HR.nc'

  if (ifadv) then
    advStr = ''
  else
    advStr = '_noadv'
  endif
  if (ifforc) then;  forcStr = '_forc'; else; forcStr = '_noforc'; endif
  if (ifrelax) then;  relStr = '_addrel'; else; relStr = ''; endif
  if (paramflg == 0) then
    parmStr = '_noparam'
  elseif (paramflg == 1) then
    parmStr = '_paramKt' 
    k_dir = '/glade/work/yueyanglu/MOM6_OUT/tr_analysis/ktens/' &
        //'grid64/C'//trim(cStr)//'/'
    print*, 'K-tens will be read from: ', trim(k_dir)
  elseif (paramflg == 2) then
    parmStr = '_paramkChi' 
    k_dir = '/glade/work/yueyanglu/MOM6_OUT/tr_analysis/klmd_iso/' &
        //'grid64/C'//trim(cStr)//'/'
    print*, 'kChi will be read from: ', trim(k_dir)
  else
    parmStr = '_paramNOTREADY' 
  endif

  if (paramflg == 0) then
    dir_save = '/glade/work/yueyanglu/MOM6_OUT/tr_off_64/sols_pt/tr'//trim(advStr)&
                //trim(forcStr)//'_ptrel'//trim(parmStr)//trim(relStr)//'/'
  else
    dir_save = '/glade/work/yueyanglu/MOM6_OUT/tr_off_64/sols_pt/tr'//trim(advStr)&
                //trim(forcStr)//'_ptrel'//''//trim(parmStr)//'_'&
                //trim(cStr)//trim(relStr)//'/'
  endif
  reldir_save = trim(dir_save)//'relforc/' ! save the extra relaxation forc
  flnm_params = trim(dir_save)//'A_README.txt'

  ! for passive temp
  r_relx_pt = 2.9e-8
  T_relx_pt = 1. / r_relx_pt

  flnm_params = trim(dir_save)//'A_README.txt'
  dir_kforc = trim(dir_save)//'paramforc/'

  ! for artificial relax! 
  T_relx = 10. * dy2hr*hr2sc ! 1day in [sec]
  dt_rel = 6. ! [hr]

  !
  flnm_h = 'h_snap__YEAR_DAY_HR.nc'
  flnm_o = 'tr__YEAR_DAY_HR.nc'
  fnm_relsv = 'forc__YEAR_DAY_HR.nc';   relforcname = 'rel#'
  fnm_kforc = 'forc__YEAR_DAY_HR.nc'; 
  forcname_k = 'cforc#'; forcname_ka = 'cforc_ka#'; forcname_chi = 'cforc_chi#'

!  check save dir
  inquire(directory = dir_save, exist = ifexist)
  if ( .NOT. ifexist ) then
    call system('mkdir ' // trim(dir_save)) 
  endif
  print*, NEW_LINE('A')//'Tr sol will be saved to: ', trim(dir_save)
! 
  if (paramflg /= 0) then
    inquire(directory = dir_kforc, exist = ifexist)
    if ( .NOT. ifexist ) call system('mkdir ' // trim(dir_kforc)) 
  endif

  if(ifrelax) then
    inquire(directory = reldir_save, exist = ifexist)
    if ( .NOT. ifexist ) call system('mkdir ' // trim(reldir_save)) 
  endif

  geom_fnm = 'ocean_geometry.nc';  vert_fnm = 'Vertical_coordinate.nc'

  !============================================================
  !                  read model geometry 
  !============================================================
  ! get model grids
  call build_grid_cartesian(g, 3840., 3840., nihalo, njhalo, nz)
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
  ! print*, 'g%bathyT:', g%bathyT(1:6,20)
  ! print*, 'g%dy_Cu(0:6,)', g%dy_Cu(0:6,20)
  ! print*, 'g%dyCu(0:6,)', g%dyCu(0:6,20)
  print*, 'g%mask2dT(0:6,)', g%mask2dT(0:6,20)

  ! allocate main vars
  allocate( h(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
  allocate( h_t2(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
  allocate( h_cp(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
  allocate( dp1(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
  allocate( dp2(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
  ! 
  allocate( uhtr(g%isdb:g%iedb, g%jsd:g%jed,  g%ke) )
  allocate( vhtr(g%isd:g%ied,   g%jsdb:g%jedb, g%ke) )
  allocate( uh_cp(g%isdb:g%iedb, g%jsd:g%jed,  g%ke) )
  allocate( vh_cp(g%isd:g%ied,   g%jsdb:g%jedb, g%ke) )
  ! 
  allocate( trTemp(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
  !
  allocate( forc    (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
  allocate( relxforc(g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
  allocate( relx    (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
  allocate( relx_pt (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
  allocate( forc_ka    (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
  allocate( forc_chi    (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
  allocate( forc_kchi    (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
  ! 
  allocate( ka(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
  allocate( chiu(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
  allocate( chiv(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
  ! 
  allocate( kxx(g%isdb:g%iedb, g%jsd:g%jed,  g%ke) )
  allocate( kxy(g%isdb:g%iedb, g%jsd:g%jed,  g%ke) )
  allocate( kyx(g%isd:g%ied,   g%jsdb:g%jedb,  g%ke) )
  allocate( kyy(g%isd:g%ied,   g%jsdb:g%jedb,  g%ke) )

  print*, 'Shape of h: ', shape(h)
  print*, 'Shape of uhtr: ', shape(uhtr)
  print*, 'Shape of vhtr: ', shape(vhtr)

  !============================================================
  !            initialize CS & tracers
  !============================================================
  ! register tracer CS
  call call_tracer_register(G, tr_Reg, NTR, cnames)
  print*, NEW_LINE('A')//'Shape of tr_Reg%tr(1)%t: ', shape(tr_Reg%tr(1)%t) 
  print*, 'tr_Reg%ntr: ', tr_Reg%ntr 
  print*, 'tr_reg%tr(1:NTR)%name= ', tr_Reg%tr(1:NTR)%name
  print*, 'call_tracer_register done !!!'
  call call_tracer_register(G, tr_Reg0, NTR, cnames)
  print*, 'ctr_Reg0 registered !!!', NEW_LINE('A')

  ! initialize tracer distribution
  call initialize_tracer(G, tr_Reg, IC_tr_fnm)
  print*, 'Initialize tracers done !!!', NEW_LINE('A')

  ! initialize tracer diffusion CS
  call tracer_hor_diff_init(G, tr_diff_CSp, khtr)
  ! print*, 'Diffusivity tr_diff_CSp%khtr = ', tr_diff_CSp%khtr

  !============================================================
  !                  set the times  
  !============================================================
  time_s = (/ 21, 1, 0 /)
  time_e = (/ 22, 365, 0 /)
  dt_arch = 6.  ! [hr]
  dt_save = 24. ! [hr]
  dt_in = dt_arch*hr2sc ! dt_arch*hr2sc or 300.
  fac = dt_in ! fac to convert uh_mean [m3/s] to uhtr_sum [m3], dt_arch*hr2sc
  print*, 'fac = ', fac
  ! # of time steps of interp within one archive interval 
  nintrp = nint(dt_arch*hr2sc / dt_in)

! # of archives to be read
  narch = nint( (time_e(1) - time_s(1))*yr2dy*(dy2hr/dt_arch) + &
      (time_e(2) - time_s(2))*(dy2hr/dt_arch) + &
      (time_e(3) - time_s(3))/dt_arch + 1. );
  ! 
  print*,NEW_LINE('A')//'time_s=',time_s,', time_e=',time_e, &
      ', dt_arch=',dt_arch,', narch=', narch,'dt_in [s] =', dt_in, 'nintrp=',nintrp

  !============================================================
  !    record params in a .txt file 
  !============================================================
  inquire(FILE = flnm_params, exist = ifexist)
  if(.NOT. ifexist ) then
    open(unit=10, file=flnm_params, status='new')  

    write(10,'(a,a)') 'IC_tr_fnm = ',IC_tr_fnm
    write(10,'(a,a)') 'flx_dir = ',flx_dir
    write(10,'(a,a)') 'h_dir = ',h_dir

    write(10,'(a, f7.1, a)') 'khtr = ', khtr, ' [m2/s]'
    write(10,'(a, f7.1, a)') 'dt_arch = ', dt_arch, ' [h]'
    write(10,'(a, f7.1, a)') 'dt_in = ', dt_in, ' [s]'
    write(10,'(a, f7.1, a)') 'fac = ', fac, ' [s]'
    write(10,'(a)') 'cnames = '; write(10,'(a)') cnames
    write(10,'(a, i2)') 'NTR = ', NTR
    write(10,'(a, i2)') 'paramflg = ', paramflg

    if (ifforc) then
      write(10,'(a,a)') 'eforc_dir = ',eforc_dir
      write(10,'(a,a)') 'forcnm = ', forcnm
    endif
    if (ifrelax) then
      write(10,'(a,a)') 'relnm = ', relnm
      write(10,'(a,a)') 'relforcname = ', relforcname
      write(10,'(a, f5.2, a)') 'T_relx = ', T_relx /dy2hr/hr2sc, ' [d]'
      write(10,'(a, f7.1, a)') 'dt_rel = ', dt_rel, ' [h]'
    endif
    if (paramflg /= 0) then
      write(10,'(a,a)') 'k_dir = ', k_dir
      write(10,'(a,a)') 'Params diagnosed from tracers: ', cStr
    endif
    write(10,'(a,a)') 'dir_save = ',dir_save

    close(unit=10)
  endif

  !============================================================
  !            main loop: time step tracer at t1 --> t2
  !============================================================
  call system_clock(count_rate=clock_rate) !Find the time rate
  call system_clock(count=clock_start)     !Start Timer

  do iar = 1, narch ! t1

    print*,NEW_LINE('A')//'========== DO iar = ', iar, ' of narch = ', narch, '=========='
    !---------- times: t1, t2, and t1.5 ----------
    hr_elapsed  = (iar-1)*dt_arch 
    call get_timeint(iyear1, iday1, ihr1, time_s, hr_elapsed, yr2dy, dy2hr)
    hr_elapsed  = iar*dt_arch 
    call get_timeint(iyear2, iday2, ihr2, time_s, hr_elapsed, yr2dy, dy2hr)
    hr_elapsed  = (iar-1)*dt_arch + .5*dt_arch
    call get_timeint(iyear15, iday15, ihr15, time_s, hr_elapsed, yr2dy, dy2hr)
    
    !----------  read h at t1&t2, and u/v averaged over t1~t2 ----------
    !----- h at t1
    call wr_flnm_time(flnm_h, iyear1, iday1, ihr1)
    flnmh_ful = trim(h_dir)//trim(flnm_h)
    print*, NEW_LINE('A')//"Read h at t1 from:", trim(flnmh_ful)
    call get_ncvar(flnmh_ful, "h", h, g%nihalo, g%njhalo)

    !----- h at t2
    call wr_flnm_time(flnm_h, iyear2, iday2, ihr2)
    flnmh_ful = trim(h_dir)//trim(flnm_h)
    print*, "Read h at t2 from:", trim(flnmh_ful)
    call get_ncvar(flnmh_ful, "h", h_t2, g%nihalo, g%njhalo)
    ! copy the pre-calcultaed "h" (cs%h in "step_offline") to "h_end" 

    !----- accumulated uh t1~t2; its .nc name is 't1.5'
    call wr_flnm_time(flnm_uv, iyear15, iday15, ihr15)
    flnmuv_ful = trim(flx_dir)//trim(flnm_uv)
    print*, "Read uvh(tr_sum or mean)-t1~t2 from:",trim(flnmuv_ful)//NEW_LINE('A')
    call get_ncvar(flnmuv_ful, funm, uhtr, g%nihalo, g%njhalo)
    call get_ncvar(flnmuv_ful, fvnm, vhtr, g%nihalo, g%njhalo)
    ! apply factor: uh_mean [m3/s] --> uhtr_sum [m3]
    uhtr = uhtr * fac; vhtr = vhtr * fac; 
    uh_cp = uhtr; vh_cp = vhtr ! uvh_cp [m3]

    !----- save the tracer before adv
    h_cp = h
    tr_Reg0 = tr_Reg

    !----------  time step tracer ----------
    do istep = 1, nintrp
      !----- interp  ----
      q1 = float(istep-1) / float(nintrp)
      q2 = float(istep) / float(nintrp)
      do i = g%isc, g%iec
        do j = g%jsc, g%jec
          do k = 1, g%ke
            dp1(i,j,k) = h(i,j,k)*(1.-q1) + h_t2(i,j,k)*q1
            dp2(i,j,k) = h(i,j,k)*(1.-q2) + h_t2(i,j,k)*q2
          enddo
        enddo
      enddo

      !---------- advevtion  
      if (ifadv) then
        call tracer_mass(cm1,tr_Reg,dp1,G)
        ! "dp1" will be updated; "uhtr" will become 0 after "offline_advection_layer"!
        call offline_advection_layer(dt_in, tr_Reg, G, dp1, dp2, uhtr, vhtr)
        call tracer_mass(cm2,tr_Reg,dp2,G)

        print*, '--- CMASS (ADV) for tr-1 cm1=', cm1(:,1), 'cm2=', cm2(:,1), & 
            'Err=', (cm2(:,1)-cm1(:,1)) / cm1(:,1)
        q = cm1 / cm2; call corr_tracer(q,tr_Reg,G)
      else ! no advection
        do m = 1, NTR
          do k = 1, nz
            do j = g%jsc, g%jec
              do i = g%isc, g%iec
                ! c_new = (c_old*h_old + delta_ch) / h_new
                tr_Reg%tr(m)%t(i,j,k) = (tr_Reg%tr(m)%t(i,j,k) * dp1(i,j,k) &
                  + 0.) / dp2(i,j,k)
              enddo               
            enddo                  
          enddo !k
        enddo !ntr

      endif  ! adv

      !---------- diffusion 
      print*, '------ DIFF ------'
      call tracer_hordiff(h_cp, dt_in, g, tr_diff_CSp, tr_Reg)
    enddo ! interp

    !---------- update tracer by forcing [m/s*c] pre-saved at t1 !----------
    if (ifforc) then
      print*, NEW_LINE('A')//'----- FORCING -----'
      call wr_flnm_time(flnm_forc, iyear1, iday1, ihr1)
      flnmforc_ful = trim(eforc_dir)//trim(flnm_forc)
      l = len_trim(forcnm)
      do m = 1, NTR
        write(forcnm(l:l),'(i1)') m
        call get_ncvar(flnmforc_ful, forcnm, forc(:,:,:,m), g%nihalo, g%njhalo)
        print*, trim(forcnm)," readed"
      enddo
      print*, "from:"//trim(flnmforc_ful)

      ! update tracer by forcing
      call tracer_forc(tr_Reg, dt_arch*hr2sc, forc, h_t2, h_t2, g)
    endif

    !---------- relax tracer by ATM forcing ----------
    print*, NEW_LINE('A')//'----- RELAX passive temp tr -----'
    l = len_trim(relnm)
    do m = 1, NTR
      write(relnm(l:l),'(i1)') carries(m)
      call get_ncvar(flnmrel_ful, relnm, relx(:,:,:,m), g%nihalo, g%njhalo)
      print*, trim(relnm)," readed (relx profile)"
    enddo
    print*, "from:"//trim(flnmrel_ful) //NEW_LINE('A') &
        //'with T_relx [d]=',T_relx/dy2hr/hr2sc
    ! relax tracer 
    call tracer_relx(tr_Reg, tr_Reg0, dt_arch*hr2sc, relx, T_relx, g, h_t2, h_t2, relxforc)

    !---------- update tracer by parameters (use K saved at t1 times c before adv (at t1))-----
    if (paramflg.eq.1) then
      print*, NEW_LINE('A')//'----- PARAMS (', paramflg, ') -----'
      call wr_flnm_time(flnm_k, iyear1, iday1, ihr1)
      flnm_k_ful = trim(k_dir)//trim(flnm_k)
      ! read
      call get_ncvar(flnm_k_ful, 'Kxx', kxx, g%nihalo, g%njhalo)
      call get_ncvar(flnm_k_ful, 'Kxy', kxy, g%nihalo, g%njhalo)
      call get_ncvar(flnm_k_ful, 'Kyx', kyx, g%nihalo, g%njhalo)
      call get_ncvar(flnm_k_ful, 'Kyy', kyy, g%nihalo, g%njhalo)
      print*, "K-tens read from:"//trim(flnm_k_ful)
      ! kxx = 100.; kxy = 0.; kyx = 0.; kyy = 100.
      call tracer_ktensor(tr_Reg, tr_Reg0, dt_arch*hr2sc, kxx, kxy, &
             kyx, kyy, h_cp, h_t2, h_t2, G) ! or h_cp, h_t2, h_t2

    elseif (paramflg.eq.2) then
      print*, NEW_LINE('A')//'----- PARAMS (', paramflg, ') -----'
      call wr_flnm_time(flnm_k, iyear1, iday1, ihr1)
      flnm_k_ful = trim(k_dir)//trim(flnm_k)
      ! read
      call get_ncvar(flnm_k_ful, 'K11', ka, g%nihalo, g%njhalo)
      call get_ncvar(flnm_k_ful, 'lmdu', chiu, g%nihalo, g%njhalo)
      call get_ncvar(flnm_k_ful, 'lmdv', chiv, g%nihalo, g%njhalo)
      print*, "kChi read from:"//trim(flnm_k_ful)
      ! print*, "max kappa Z1 ", maxval(ka(g%isc:g%iec,g%jsc:g%jec,1)), &
      !  minval(ka(g%isc:g%iec,g%jsc:g%jec,1))
      ! ka = 100.; chiu = 0.; chiv = 0.
      call tracer_kachi(tr_Reg, tr_Reg0, dt_arch*hr2sc, ka, chiu, chiv, & 
                h_cp, h_t2, h_t2, G, forc_ka, forc_chi, forc_kchi)
      
      ! save parameterized forcing to t1.nc
      print*, '----- Save k-Chi forcing -----'
      call wr_flnm_time(fnm_kforc, iyear1, iday1, ihr1)
      fnm_kforc_ful = trim(dir_kforc)//trim(fnm_kforc)
      do m = 1, NTR
        l = len_trim(forcname_ka); write(forcname_ka(l:l),'(i1)') carries(m)
        l = len_trim(forcname_chi); write(forcname_chi(l:l),'(i1)') carries(m)
        l = len_trim(forcname_k); write(forcname_k(l:l),'(i1)') carries(m)

        call wr_ncfile(fnm_kforc_ful, forcname_ka, forc_ka(:,:,:,m), g%nihalo, g%njhalo)
        call wr_ncfile(fnm_kforc_ful, forcname_chi, forc_chi(:,:,:,m), g%nihalo, g%njhalo)
        call wr_ncfile(fnm_kforc_ful, forcname_k, forc_kchi(:,:,:,m), g%nihalo, g%njhalo)
        print*, trim(forcname_ka),' & ',trim(forcname_chi),' & ',trim(forcname_k)," saved"
      enddo
      print*, "to: ", trim(fnm_kforc_ful)
    endif

    !---------- (OPT) relax tracer towards the ref at t2, & save RF to t1.nc ----------
    if (ifrelax) then
      if (mod(iar*dt_arch, dt_rel)==0.0 .OR. iar==1) then 

        print*, NEW_LINE('A')//'----- RELAXATION (artificial) -----'
        call wr_flnm_time(flnm_rel, iyear2, iday2, ihr2)
        flnmrel_ful = trim(rel_dir)//trim(flnm_rel)
        l = len_trim(relnm)
        do m = 1, NTR
          write(relnm(l:l),'(i1)') carries(m)
          call get_ncvar(flnmrel_ful, relnm, relx(:,:,:,m), g%nihalo, g%njhalo)
          print*, trim(relnm)," readed (relx profile)"
        enddo
        print*, "from:"//trim(flnmrel_ful) //NEW_LINE('A') &
            //'with T_relx [d]=',T_relx/dy2hr/hr2sc
        ! relax tracer 
        call tracer_relx(tr_Reg, tr_Reg0, dt_arch*hr2sc, relx, T_relx, g, h_t2, h_t2, relxforc)

        ! save rel forcing at each relaxation to t1.nc
        print*, '----- SAVE relax forc -----'
        call wr_flnm_time(fnm_relsv, iyear1, iday1, ihr1)
        fnm_relsv_ful = trim(reldir_save)//trim(fnm_relsv)
        l = len_trim(relforcname)
        do m = 1, NTR
          write(relforcname(l:l),'(i1)') carries(m)
          call wr_ncfile(fnm_relsv_ful, relforcname, relxforc(:,:,:,m), g%nihalo, g%njhalo)
          print*, trim(relforcname)," saved"
        enddo
        print*, "to: ", trim(fnm_relsv_ful)
      endif
    endif

    !----------  save tracer at t2 ----------
    if (mod(iar*dt_arch, dt_save) == 0.0) then 
      print*, NEW_LINE('A')//'----- SAVE tr -----'
      call wr_flnm_time(flnm_o, iyear2, iday2, ihr2)
      flnm_ful_o = trim(dir_save)//trim(flnm_o)
      l = len_trim(trname)
      do m = 1, NTR
        trname = tr_reg%tr(m)%name
        call wr_ncfile(flnm_ful_o, trname, tr_Reg%tr(m)%t, g%nihalo, g%njhalo)
        print*, trname," saved"
      enddo
      print*, "to: ", trim(flnm_ful_o)
    endif

  enddo ! iar

  call system_clock(count=clock_stop)      ! Stop Timer
  write(*,*) 'Job took ',real(clock_stop-clock_start)/real(clock_rate),' seconds'

contains

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
