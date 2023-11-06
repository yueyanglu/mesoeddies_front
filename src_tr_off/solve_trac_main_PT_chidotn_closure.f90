program solve_trac_main
    !
    ! Compiling and linking in one step:
    ! ifort -o exe test.f90
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
      flnm_rel, flnmrel_ful, flnm, &
      advStr, forcStr, relStr, parmStr, flnm_params, &
      cStr, k_dir, flnm_k, flnm_k_ful, tmstr, exp_dir, forc_dir, &
      chidotn_dir, flnm_chidotn, flnm_chidotn_ful
  character*20  :: funm, fvnm, forcnm, relnm ! name of flux components in the .nc file
  character*200 :: dir_save, flnm_o, flnm_ful_o
  character*200 :: reldir_save, fnm_relsv, fnm_relsv_ful, relforcname
  character*200 :: dir_kforc, fnm_kforc, fnm_kforc_ful, forcname_k, forcname_ka, forcname_chi
  character*20   :: trname = 'tr#' ! e.g. "tr1"
  character*20  :: trtendname = 'tr#_tendency' ! e.g. "tr1"
  logical       :: ifexist, ifsave, ifsave_last
  integer       :: ncid, l
  integer       :: clock_rate, clock_start, clock_stop ! running time

  ! times 
  integer :: time_s(3), time_e(3) ! Array [yr, day of year, hr],
                                ! consistent w/ how nc files are saved
  real    :: dt_arch   ! intervals between archive input [hr]
  real    :: dt_save   ! intervals between the output [hr]
  real    :: dt_in    ! actual time step for h [hr]
  real    :: dt        ! time step for continuity [s]
  real    :: hr_elapsed, dy_elapsed ! hr/dy that elapsed wrt the initial time
  integer :: narch     ! # of archieve to be read
  integer :: nintrp     ! # of time steps of interp within 1 archive interval 
  integer :: iyear1, iday1, ihr1, iyear15, iday15, ihr15, iyear2, iday2, ihr2 
  integer :: iar, istep, i, j, k, m
  real, parameter :: yr2dy = 365., dy2hr = 24., hr2sc = 3600.
  real    :: missing = 1.0E20 ! missing values in nc file
  logical :: ifadv, ifforc, iffct, ifrelax ! if update by pre-diagnosed tracer forcing
  real    :: T_relx, r_relx ! relaxation time scale [s]
  integer :: paramflg 
  real    :: ka_hold

  ! dimensions
  integer, parameter :: NDIMS = 4 ! # of dim of vars in file
  integer :: nih, njh, niq, njq, nrecs
  integer, parameter :: ni = 64, nj = 64, nz = 3
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
  tmstr = '_180d'
  exp_dir = '/glade/work/yueyanglu/MOM6_OUT/tr_off_64'//trim(tmstr)//'/'
  forc_dir = '/glade/work/yueyanglu/MOM6_OUT/forc_uvh_64/'
  !
  IC_tr_fnm = trim(forc_dir)//'ICs/trac_init_pt_cs.nc'
  flx_dir = trim(forc_dir)//'uvhm_CS_decomp'//trim(tmstr)//'/'  
  h_dir = trim(forc_dir)//'sol_h'//trim(tmstr)//'/' 
  eforc_dir = trim(exp_dir)//'cforc_pt_prog_addreld10/'
  flnmrel_ful = trim(forc_dir)//'ICs/trac_init_pt_cs.nc'; relnm = 'tr1'
  !
  flnm_uv = 'uvh_mean__YEAR_DAY_HR.nc'; funm = "uh"; fvnm = "vh"; 
  flnm_forc = 'forc__YEAR_DAY_HR.nc'; forcnm = 'cforc_rel#' ! cforc_m_rel# cforc_m#

  ifadv = 1;  iffct = 0;  ifforc = 0
  paramflg = 1
  flnm_k = 'K__YEAR_DAY_HR.nc'

  if (ifadv) then
    if (iffct) then; advStr = '_fctadv'; else; advStr = ''; endif
  else
    advStr = '_noadv'
  endif
  if (ifforc) then;  forcStr = '_forc'; else; forcStr = '_noforc'; endif


  chidotn_dir = trim(exp_dir)//'/params/ideal_prof/'
  flnm_chidotn_ful = trim(chidotn_dir)//'chinxyz_closure.nc'
  print*, 'Time-indepent Chi*n will be read from: ', trim(flnm_chidotn_ful)

  dir_save = trim(exp_dir)//'sols_pt/tr'//trim(advStr)//trim(forcStr)&
              //'_ptrel_paramka400Chidotn_chinclosure/'

  reldir_save = trim(dir_save)//'relforc/'
  r_relx = 2.9e-8
  T_relx = 1. / r_relx

  flnm_params = trim(dir_save)//'A_README.txt'
  dir_kforc = trim(dir_save)//'paramforc/'

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
  
  if(ifrelax .and. ifforc) then
    inquire(directory = reldir_save, exist = ifexist)
    if ( .NOT. ifexist ) call system('mkdir ' // trim(reldir_save)) 
  endif

  geom_fnm = 'ocean_geometry.nc';  vert_fnm = 'Vertical_coordinate.nc'

  !============================================================
  !                  model geometry & vars alloc
  !============================================================
  ! call get_dims_geom(geom_fnm_ful, nih, njh, niq, njq)
  ! call get_dims_vert(vert_fnm_ful, nz)
  ! get model grids
  call build_grid_cartesian(g, 3840., 3840., ni, nj, nz)
  ! print*, 'Read model geom from:', trim(geom_fnm_ful)
  print*, 'g extents:      ', g%isg,g%ieg,g%jsg,g%jeg
  ! print*, 'gB-sym extents: ', g%isgb,g%iegb,g%jsgb,g%jegb
  print*, 'c extents:      ', g%isc,g%iec,g%jsc,g%jec
  ! print*, 'cB-sym extents: ', g%iscb,g%iecb,g%jscb,g%jecb
  print*, 'd extents:      ', g%isd,g%ied,g%jsd,g%jed
  ! print*, 'dB-sym extents: ', g%isdb,g%iedb,g%jsdb,g%jedb
  ! print*, 'halos (i,j):  ', g%nihalo, g%njhalo
  ! print*, 'g%idg_offset=', g%idg_offset, 'g%jdg_offset=',g%jdg_offset
  print*, 'g%ke:  ', g%ke
  ! print*, 'g%bathyT:', g%bathyT(1:6,20)
  ! print*, 'g%dy_Cu(0:6,)', g%dy_Cu(0:6,20)
  ! print*, 'g%dyCu(0:6,)', g%dyCu(0:6,20)
  print*, 'g%mask2dT(0:6,)', g%mask2dT(0:6,20)

  !===== allocate key vars
  call vars_alloc(g, NTR)
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
  print*, 'tr_Reg0 registered !!!', NEW_LINE('A')

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
  time_e = (/ 23, 1, 0 /)
  dt_arch = 6.  ! [hr]
  dt_save = 24. ! [hr]
  if (iffct) then ! use "fct3d.f" with time interp 
    dt_in = 300. ![s]
  else ! interp done inside the MOM6's adv code
    dt_in = dt_arch*hr2sc/36 ! dt_arch*hr2sc (360min) or 300.
  endif
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

    write(10,'(a,a)') 'flnmrel_ful = ',flnmrel_ful
    write(10,'(a,a)') 'relforcname = ',relforcname
    write(10,'(a, e15.6, a)') 'r_relx = ', r_relx, ' [1/s]'
    write(10,'(a, f5.2, a)') 'T_relx = ', T_relx /dy2hr/hr2sc, ' [d]'

    if (paramflg /= 0) then
      write(10,'(a,a)') 'k_dir = ', k_dir
      write(10,'(a,a)') 'Params diagnosed from tracers: ', cStr
    endif
    write(10,'(a,a)') 'dir_save = ',dir_save

    close(unit=10)
  endif

  !============================================================
  !            read relaxation profile
  !============================================================ 
  l = len_trim(relnm)
  do m = 1, NTR
    write(relnm(l:l),'(i1)') carries(m)
    call get_ncvar(flnmrel_ful, relnm, relx(:,:,:,m), g%nihalo, g%njhalo)
    print*, trim(relnm)," readed (relx profile)"
  enddo
  print*, "from:"//trim(flnmrel_ful) //NEW_LINE('A') &
      //'with T_relx [d]=',T_relx/dy2hr/hr2sc

  !============================================================
  !            main loop: time step tracer at t1 --> t2
  !============================================================
  call system_clock(count_rate=clock_rate) !Find the time rate
  call system_clock(count=clock_start)     !Start Timer

  ifsave = 0
  do iar = 1, narch ! t1

    print*,NEW_LINE('A')//'========== DO iar = ', iar, ' of narch = ', narch, '=========='
    
    ! if the tracer is saved at the last arch t-step 
    ifsave_last = ifsave 
    if (mod(iar*dt_arch, dt_save) .eq. 0.0) then
      ifsave = 1
    else
      ifsave = 0
    endif

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

    !----- accumulated uh t1~t2; its .nc name is 't1.5'
    call wr_flnm_time(flnm_uv, iyear15, iday15, ihr15)
    flnmuv_ful = trim(flx_dir)//trim(flnm_uv)
    print*, "Read uvh(tr_sum or mean)-t1~t2 from:",trim(flnmuv_ful)//NEW_LINE('A')
    call get_ncvar(flnmuv_ful, funm, uhtr, g%nihalo, g%njhalo)
    call get_ncvar(flnmuv_ful, fvnm, vhtr, g%nihalo, g%njhalo)

    !---------- read: eddy forcing & parameters (opt) at t1.nc ----------
    if (ifforc) then
      print*, NEW_LINE('A')//'----- read FORCING -----'
      call wr_flnm_time(flnm_forc, iyear1, iday1, ihr1)
      flnmforc_ful = trim(eforc_dir)//trim(flnm_forc)
      l = len_trim(forcnm)
      do m = 1, NTR
        write(forcnm(l:l),'(i1)') m
        call get_ncvar(flnmforc_ful, forcnm, forc(:,:,:,m), g%nihalo, g%njhalo)
        print*, trim(forcnm)," readed"
      enddo
      print*, "from:"//trim(flnmforc_ful)
    endif

    !---------- read kappa & chi*n -----
    print*, NEW_LINE('A')//'----- read PARAMS (ka & chi*n) -----'
    ka = 400.
    call get_ncvar(flnm_chidotn_ful, 'chidotn', chin, g%nihalo, g%njhalo)
    print*, 'chin=,', chin(g%isc+22,g%jsc+33,1)

    ! 
    ! convert to mass transport in the interp interval [m3]
    ! uh_mean [m3/s] --> uhtr_sum [m3]
    uhtr = uhtr * fac; vhtr = vhtr * fac; 

    !----- save the tracer and layer thick before adv
    h_cp = h
    tr_Reg0 = tr_Reg
    ! c_tm = 0.

    !----------  time step tracer ----------
    do istep = 1, nintrp

      uh_cp = uhtr; vh_cp = vhtr ! uvh_cp [m3]

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

      !---------- advection  
      if (ifadv) then

        print*, '------ ADV at istep = ', istep,' of ', nintrp, '------'
        if (iffct) then
          call fct3d(dt_in, tr_Reg, G, dp1, dp2, uhtr, vhtr)
        else
          call tracer_mass(cm1,tr_Reg,dp1,G)
          ! "dp1" will be updated; "uhtr" will become 0 after "offline_advection_layer"!
          call offline_advection_layer(dt_in, tr_Reg, G, dp1, dp2, uh_cp, vh_cp)
          call tracer_mass(cm2,tr_Reg,dp2,G)

          print*, '--- CMASS (ADV) for tr-1 cm1=', cm1(:,1), 'cm2=', cm2(:,1), & 
              'Err=', (cm2(:,1)-cm1(:,1)) / cm1(:,1)
          q = cm1 / cm2; call corr_tracer(q,tr_Reg,G)
        endif

      else ! no advection
        call tracer_mass(cm1,tr_Reg,dp1,G)
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
        call tracer_mass(cm2,tr_Reg,dp2,G)
        print*, '--- CMASS (ADV) for tr-1 cm1=', cm1(:,1), 'cm2=', cm2(:,1), & 
            'Err=', (cm2(:,1)-cm1(:,1)) / cm1(:,1)
        q = cm1 / cm2; call corr_tracer(q,tr_Reg,G)
      endif  ! adv

      !---------- diffusion 
      print*, '------ DIFF ------'
      call tracer_mass(cm1,tr_Reg,h_cp,G)
      call tracer_hordiff(h_cp, dt_in, g, tr_diff_CSp, tr_Reg)
      call tracer_mass(cm2,tr_Reg,h_cp,G)
      ! q = cm1 / cm2; call corr_tracer(q,tr_Reg,G)
      ! print*, '-----CMASS (DIF) bf=', cm1, 'af=', cm2, 'Err=', (cm2-cm1)/cm1

      !---------- update tracer by forcing [m/s*c] pre-saved at t1 !----------
      if (ifforc) then
        print*, NEW_LINE('A')//'----- apply FORCING -----'
        ! update tracer by forcing [in the interval]
        call tracer_mass(cm1,tr_Reg,dp2,G)
        call tracer_forc(tr_Reg, dt_in, forc, h_t2, h_t2, g) ! dt_arch*hr2sc
        call tracer_mass(cm2,tr_Reg,dp2,G)
        print*, '   Conserve cmass: cm1=', cm1(:,1), 'cm2=', cm2(:,1), &
          'Err=', (cm2(:,1)-cm1(:,1)) / cm1(:,1)
        q = cm1 / cm2;
        call corr_tracer(q,tr_Reg,G)
      endif

      !---------- update tracer by parameters (use K saved at t1 times c before adv (at t1))-----
      if (paramflg) then
        call tracer_mass(cm1,tr_Reg,dp2,G)
        call tracer_kachi_chidotn_sign(tr_Reg, tr_Reg0, dt_in, ka, chin, & 
              dp1, dp2, dp2, G, forc_ka, forc_chi, forc_kchi) ! h_cp, h_t2, h_t2
        call tracer_mass(cm2,tr_Reg,dp2,G)
        !
        print*, '--- Conserve cmass for params: cm1=', cm1(:,1), 'cm2=', cm2(:,1), & 
            'Err=', (cm2(:,1)-cm1(:,1)) / cm1(:,1)
        q = cm1 / cm2; 
        call corr_tracer(q,tr_Reg,G)
      endif

      !---------- relax tracer by ATM forcing ----------
      print*, NEW_LINE('A')//'----- RELAX passive temp tr -----'
      call tracer_relx(tr_Reg, tr_Reg0, dt_in, relx, T_relx, g, dp2, dp2, relxforc)
      !h_t2, h_t2

      !----------  calc dc/dt ----------
      !!!!!!! NOTE that this does NOT work when nintrp=1 !!!!!!!
      ! 1. record the tracer at the (nt-1) step
      if (ifsave .and. istep == nintrp-1) then
        do m = 1, NTR
          c_t1(:,:,:,m) = tr_Reg%tr(m)%t
          ! cm_t1(:,:,:,m) = tr_Reg%tr(m)%t * dp2
        enddo
        print*, 'record c at nt-1'
      endif

      ! 2. calc tendency at the nt step, using snapshots at (nt-1) and
      !    (nt+1) steps
      if (ifsave_last .and. istep == 1) then
        do m = 1, NTR
          c_t2(:,:,:,m) = tr_Reg%tr(m)%t
          ! cm_t2(:,:,:,m) = tr_Reg%tr(m)%t * dp2
        enddo
        dcdt = (c_t2 - c_t1) / (2.0*dt_in)
        ! dcmdt = (cm_t2 - cm_t1) / (2.0*dt_in)
        ! here dcmdt is very close to h*dcdt+c*dhdt calc offline
        print*, 'calc: dc/dt = [c(nt+1)-c(nt-1)]/(2dt_in)'
        !
        call wr_flnm_time(flnm_o, iyear1, iday1, ihr1)
        flnm_ful_o = trim(dir_save)//trim(flnm_o)
        do m = 1, NTR
          trtendname = trim(tr_reg%tr(m)%name)//'_tendency'
          call wr_ncfile(flnm_ful_o, trtendname, dcdt(:,:,:,m), g%nihalo, g%njhalo)
        enddo
        print*, "dcdt saved to: ", trim(flnm_ful_o)
      endif

    enddo ! interp


    !----------  save ----------
    if (ifsave) then 

      ! save parameterized forcing to t1.nc
      if (paramflg /= 0) then
        print*,NEW_LINE('A')//'----- SAVE k-Chi forcing -----'
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

      !----------  save tracer to t2.nc
      print*, NEW_LINE('A')//'----- SAVE tr -----'
      call wr_flnm_time(flnm_o, iyear2, iday2, ihr2)
      flnm_ful_o = trim(dir_save)//trim(flnm_o)
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
