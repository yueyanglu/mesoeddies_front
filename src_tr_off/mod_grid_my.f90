module mod_grid_my
!> Provides the ocean grid type, modified from mom_grid.
use netcdf
use mod_nc_wr_rd

implicit none; private

public build_grid_file, build_grid_cartesian

!> Ocean grid type. See mom_grid for details.
type, public :: ocean_grid_type

  integer :: isc !< The start i-index of cell centers within the computational domain
  integer :: iec !< The end i-index of cell centers within the computational domain
  integer :: jsc !< The start j-index of cell centers within the computational domain
  integer :: jec !< The end j-index of cell centers within the computational domain

  integer :: isd !< The start i-index of cell centers within the data domain
  integer :: ied !< The end i-index of cell centers within the data domain
  integer :: jsd !< The start j-index of cell centers within the data domain
  integer :: jed !< The end j-index of cell centers within the data domain

  integer :: isg !< The start i-index of cell centers within the global domain
  integer :: ieg !< The end i-index of cell centers within the global domain
  integer :: jsg !< The start j-index of cell centers within the global domain
  integer :: jeg !< The end j-index of cell centers within the global domain

  integer :: isdb !< The start i-index of cell vertices within the data domain
  integer :: iedb !< The end i-index of cell vertices within the data domain
  integer :: jsdb !< The start j-index of cell vertices within the data domain
  integer :: jedb !< The end j-index of cell vertices within the data domain

  integer :: iscb !< The start i-index of cell vertices within the computational domain
  integer :: iecb !< The end i-index of cell vertices within the computational domain
  integer :: jscb !< The start j-index of cell vertices within the computational domain
  integer :: jecb !< The end j-index of cell vertices within the computational domain

  integer :: isgb !< The start i-index of cell vertices within the global domain
  integer :: iegb !< The end i-index of cell vertices within the global domain
  integer :: jsgb !< The start j-index of cell vertices within the global domain
  integer :: jegb !< The end j-index of cell vertices within the global domain

  integer :: nihalo, njhalo, ke, idg_offset, jdg_offset

  real :: len_lon, len_lat, south_lat, west_lon

  real, allocatable, dimension(:,:) :: &

    bathyT, &   ! ocean bottome depth [m], given by MAXIMUM_DEPTH
    mask2dT, &   !< 0 for land points and 1 for ocean points on the h-grid [nondim].
    dxT, &       !< dxT is delta x at h points [L ~> m].
    dyT, &       !< dyT is delta y at h points [L ~> m].
    areaT, &     !< The area of an h-cell [L2 ~> m2].
    IareaT, &    !< 1/areaT [L-2 ~> m-2].
    IdxT, &      !< 1/dxT [L-1 ~> m-1].
    IdyT, &      !< IdyT is 1/dyT [L-1 ~> m-1].
    geolatT, &   !< The geographic latitude at h points in degrees of latitude or m.
    geolonT, &   !< The geographic longitude at h points in degrees of longitude or m.
! 
    mask2dCu, &  !< 0 for boundary points and 1 for ocean points on the u grid [nondim].
    dxCu, &      !< dxCu is delta x at u points [L ~> m].
    IdxCu, &     !< 1/dxCu [L-1 ~> m-1].
    dyCu, &      !< dyCu is delta y at u points [L ~> m].
    IdyCu, &     !< 1/dyCu [L-1 ~> m-1].
    dy_Cu, &     !< The unblocked lengths of the u-faces of the h-cell [L ~> m].
    IareaCu, &   !< The masked inverse areas of u-grid cells [L-2 ~> m-2].
    areaCu, &       !< The areas of the u-grid cells [L2 ~> m2].
! 
    mask2dCv, &  !< 0 for boundary points and 1 for ocean points on the v grid [nondim].
    dxCv, &      !< dxCv is delta x at v points [L ~> m].
    IdxCv, &      !< dxCv is delta x at v points [L ~> m].
    dyCv, &      !< dyCv is delta y at v points [L ~> m].
    IdyCv, &      !< dyCv is delta y at v points [L ~> m].
    dx_Cv, &     !< The unblocked lengths of the v-faces of the h-cell [L ~> m].
    areaCv, &    !< The areas of the v-grid cells [L2 ~> m2].
    IareaCv, &    !< The areas of the v-grid cells [L2 ~> m2].
! 
    mask2dBu, &  !< 0 for boundary points and 1 for ocean points on the q grid [nondim].
    dxBu, &      !< dxBu is delta x at q points [L ~> m].
    IdxBu, &      !< dxBu is delta x at q points [L ~> m].
    dyBu, &      !< dyBu is delta y at q points [L ~> m].
    IdyBu, &      !< dyBu is delta y at q points [L ~> m].
    areaBu, &       !< areaBu is the area of a q-cell [L2 ~> m2]
    IareaBu       !< areaBu is the area of a q-cell [L2 ~> m2]

end type ocean_grid_type

contains

!=====================================================================
!  SUBTOUTINE: 
!=====================================================================
! calc horizontal indices (symmetric mode used) from given grid sizes
! THIS NEEDS TO BE CHECKED WITH MOM6 run !!!!!
! 
subroutine calc_index(g, nih, njh)
  type(ocean_grid_type), intent(inout) :: g 
  integer, intent(in) :: nih, njh

  ! local vars
  integer, parameter :: nihalo = 4, njhalo = 4, idg_offset = -4, jdg_offset = -4

  ! indices. see "hor_index_init" w/ symmetric mode

  g%nihalo = nihalo; g%njhalo = njhalo
  g%idg_offset = idg_offset; g%jdg_offset = jdg_offset;


! global extends: sizes are [nih,njh]
  g%isg = 1; g%ieg = nih;
  g%jsg = 1; g%jeg = njh;

! c extends: sizes are [nih,njh]
  g%isc = 1 + nihalo; g%iec = nih + nihalo 
  g%jsc = 1 + njhalo; g%jec = njh + njhalo

! d extents, derived from c-extends: sizes are [nih+2*nihalo,njh+2*njhalo]
  g%isd = g%isc - nihalo; g%ied = g%iec + nihalo
  g%jsd = g%jsc - njhalo; g%jed = g%jec + njhalo

! gB (sym-) extents, derived from g-extends
  g%isgb = g%isg-1;  g%iegb = g%ieg
  g%jsgb = g%jsg-1;  g%jegb = g%jeg
! cB (sym-) extents, derived from c-extends
  g%iscb = g%isc-1;  g%iecb = g%iec
  g%jscb = g%jsc-1;  g%jecb = g%jec
! dB (sym-) extents, derived from d-extends
  g%isdb = g%isd-1;  g%iedb = g%ied
  g%jsdb = g%jsd-1;  g%jedb = g%jed

end subroutine calc_index

!=====================================================================
!  SUBTOUTINE: 
!=====================================================================
! 
! allocate g with known indices!
! 
subroutine allocate_metrics(g)
  type(ocean_grid_type), intent(inout) :: g 
  ! Local variables
  integer :: isd, ied, jsd, jed, isdb, iedb, jsdb, jedb
  
  ! 
  isd = g%isd ; ied = g%ied ; jsd = g%jsd ; jed = g%jed
  isdb = g%isdb ; iedb = g%iedb ; jsdb = g%jsdb ; jedb = g%jedb

  ! allocate 

  allocate(g%bathyT(isd:ied,jsd:jed))    ; g%bathyT(:,:) = 0.0
  allocate(g%dxT(isd:ied,jsd:jed))       ; g%dxT(:,:) = 0.0
  allocate(g%dxCu(isdb:iedb,jsd:jed))    ; g%dxCu(:,:) = 0.0
  allocate(g%dxCv(isd:ied,jsdb:jedb))    ; g%dxCv(:,:) = 0.0
  allocate(g%dxBu(isdb:iedb,jsdb:jedb))  ; g%dxBu(:,:) = 0.0
  allocate(g%IdxT(isd:ied,jsd:jed))      ; g%IdxT(:,:) = 0.0
  allocate(g%IdxCu(isdb:iedb,jsd:jed))   ; g%IdxCu(:,:) = 0.0
  allocate(g%IdxCv(isd:ied,jsdb:jedb))   ; g%IdxCv(:,:) = 0.0
  allocate(g%IdxBu(isdb:iedb,jsdb:jedb)) ; g%IdxBu(:,:) = 0.0

  allocate(g%dyT(isd:ied,jsd:jed))       ; g%dyT(:,:) = 0.0
  allocate(g%dyCu(isdb:iedb,jsd:jed))    ; g%dyCu(:,:) = 0.0
  allocate(g%dyCv(isd:ied,jsdb:jedb))    ; g%dyCv(:,:) = 0.0
  allocate(g%dyBu(isdb:iedb,jsdb:jedb))  ; g%dyBu(:,:) = 0.0
  allocate(g%IdyT(isd:ied,jsd:jed))      ; g%IdyT(:,:) = 0.0
  allocate(g%IdyCu(isdb:iedb,jsd:jed))   ; g%IdyCu(:,:) = 0.0
  allocate(g%IdyCv(isd:ied,jsdb:jedb))   ; g%IdyCv(:,:) = 0.0
  allocate(g%IdyBu(isdb:iedb,jsdb:jedb)) ; g%IdyBu(:,:) = 0.0

  allocate(g%areaT(isd:ied,jsd:jed))       ; g%areaT(:,:) = 0.0
  allocate(g%IareaT(isd:ied,jsd:jed))      ; g%IareaT(:,:) = 0.0
  allocate(g%areaBu(isdb:iedb,jsdb:jedb))  ; g%areaBu(:,:) = 0.0
  allocate(g%IareaBu(isdb:iedb,jsdb:jedb)) ; g%IareaBu(:,:) = 0.0

  allocate(g%mask2dT(isd:ied,jsd:jed))      ; g%mask2dT(:,:) = 0.0
  allocate(g%mask2dCu(isdb:iedb,jsd:jed))   ; g%mask2dCu(:,:) = 0.0
  allocate(g%mask2dCv(isd:ied,jsdb:jedb))   ; g%mask2dCv(:,:) = 0.0
  allocate(g%mask2dBu(isdb:iedb,jsdb:jedb)) ; g%mask2dBu(:,:) = 0.0

  allocate(g%dx_Cv(isd:ied,jsdb:jedb))     ; g%dx_Cv(:,:) = 0.0
  allocate(g%dy_Cu(isdb:iedb,jsd:jed))     ; g%dy_Cu(:,:) = 0.0

  allocate(g%areaCu(isdb:iedb,jsd:jed))  ; g%areaCu(:,:) = 0.0
  allocate(g%areaCv(isd:ied,jsdb:jedb))  ; g%areaCv(:,:) = 0.0
  allocate(g%IareaCu(isdb:iedb,jsd:jed)) ; g%IareaCu(:,:) = 0.0
  allocate(g%IareaCv(isd:ied,jsdb:jedb)) ; g%IareaCv(:,:) = 0.0

  allocate(g%geoLatT(isd:ied,jsd:jed))      ; g%geoLatT(:,:) = 0.0
  allocate(g%geoLonT(isd:ied,jsd:jed))      ; g%geoLonT(:,:) = 0.0

end subroutine allocate_metrics

!=====================================================================
!  SUBTOUTINE:  get_vars_geom
!=====================================================================
! 
! get vars from the geom file (e.g. dxT) and write them into "g"
! 
subroutine get_vars_geom(g, fnm)
  type(ocean_grid_type), intent(inout) :: g         
  character*200, intent(in) :: fnm
  ! vars id
  integer :: ncid
  integer :: dxT_id, dxCu_id, dxCv_id, dxBu_id, &
             dyT_id, dyCu_id, dyCv_id, dyBu_id, mask2dT_id 
  real, allocatable, dimension(:,:) :: dxT, dxCu, dxCv, dxBu, &
                                       dyT, dyCu, dyCv, dyBu, mask2dT
  integer :: isd, ied, jsd, jed, isdb, iedb, jsdb, jedb
  integer :: isc, iec, jsc, jec, iscb, iecb, jscb, jecb

  ! known grid size in "g"
  isc = g%isc ; iec = g%iec ; jsc = g%jsc ; jec = g%jec
  isd = g%isd ; ied = g%ied ; jsd = g%jsd ; jed = g%jed
  iscb = g%iscb ; iecb = g%iecb ; jscb = g%jscb ; jecb = g%jecb
  isdb = g%isdb ; iedb = g%iedb ; jsdb = g%jsdb ; jedb = g%jedb

  ! allocate the readed arrays, no halos in the saved data !!!
  allocate( dxT(isc:iec, jsc:jec) )
  allocate( dxCu(iscb:iecb, jsc:jec) )    
  allocate( dxCv(isc:iec, jscb:jecb) )   
  allocate( dxBu(iscb:iecb, jscb:jecb) )  
  allocate( dyT(isc:iec, jsc:jec) )
  allocate( dyCu(iscb:iecb, jsc:jec) )    
  allocate( dyCv(isc:iec, jscb:jecb) )   
  allocate( dyBu(iscb:iecb, jscb:jecb) )  
  allocate( mask2dT(isc:iec, jsc:jec) )

  ! open
  call check( nf90_open(fnm, nf90_nowrite, ncid) )
  ! IDs
  call check(nf90_inq_varid(ncid, 'dxT', dxT_id))
  call check(nf90_inq_varid(ncid, 'dxCu', dxCu_id))
  call check(nf90_inq_varid(ncid, 'dxCv', dxCv_id))
  call check(nf90_inq_varid(ncid, 'dxBu', dxBu_id))
  call check(nf90_inq_varid(ncid, 'dyT', dyT_id))
  call check(nf90_inq_varid(ncid, 'dyCu', dyCu_id))
  call check(nf90_inq_varid(ncid, 'dyCv', dyCv_id))
  call check(nf90_inq_varid(ncid, 'dyBu', dyBu_id))
  call check(nf90_inq_varid(ncid, 'wet', mask2dT_id))
  ! read
  call check( nf90_get_var(ncid, dxT_id, dxT) )
  call check( nf90_get_var(ncid, dxCu_id, dxCu) )
  call check( nf90_get_var(ncid, dxCv_id, dxCv) )
  call check( nf90_get_var(ncid, dxBu_id, dxBu) )
  call check( nf90_get_var(ncid, dyT_id, dyT) )
  call check( nf90_get_var(ncid, dyCu_id, dyCu) )
  call check( nf90_get_var(ncid, dyCv_id, dyCv) )
  call check( nf90_get_var(ncid, dyBu_id, dyBu) )
  call check( nf90_get_var(ncid, mask2dT_id, mask2dT) )

  ! pass to "g", note to deal with the HALOs !!!!
  g%dxT(isc:iec, jsc:jec)      = dxT(:,:)
  g%dxCu(iscb:iecb, jsc:jec)   = dxCu(:,:)
  g%dxCv(isc:iec, jscb:jecb)   = dxCv(:,:)
  g%dxBu(iscb:iecb, jscb:jecb) = dxBu(:,:)
  g%dyT(isc:iec, jsc:jec)      = dyT(:,:)
  g%dyCu(iscb:iecb, jsc:jec)   = dyCu(:,:)
  g%dyCv(isc:iec, jscb:jecb)   = dyCv(:,:)
  g%dyBu(iscb:iecb, jscb:jecb) = dyBu(:,:)
  g%mask2dT(isc:iec, jsc:jec)  = mask2dT(:,:)

  ! close
  call check( nf90_close(ncid) )

end subroutine get_vars_geom

!=====================================================================
!  SUBTOUTINE: 
!=====================================================================
! Calculate the index and metric terms for a Cartesian grid 
! see "set_grid_metrics_cartesian" 
! 
subroutine build_grid_cartesian(g, len_lon, len_lat, niglobal, njglobal, nz)
  type(ocean_grid_type), intent(inout) :: g   
  real, intent(in) :: len_lon, len_lat      ! [km]
  integer, intent(in) :: niglobal, njglobal, nz ! 
  ! local
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, I1off, J1off
  real, allocatable :: grid_latT(:), grid_lonT(:)
  real :: dx_everywhere, dy_everywhere ! Grid spacings [m].
  real :: I_dx, I_dy                   ! Inverse grid spacings [m-1].
  real :: m_to_L  ! A unit conversion factor [L m-1 ~> nondim]
  real :: L_to_m  ! A unit conversion factor [m L-1 ~> nondim]
  real :: min_depth, dmin    ! The depth for masking in the same units as G%bathyT [m].
  real :: max_depth    ! for "flat" topo,  MAXIMUM_DEPTH

  min_depth = 1. ; max_depth = 4000.

  !------------- 1. set index (halos applied implicitly)  -------------
  call calc_index(g, niglobal, njglobal)
  call allocate_metrics(g)

  ! some aux params
  g%south_lat = 0.
  g%west_lon = 0.
  g%len_lon = len_lon
  g%len_lat = len_lat
  g%ke = nz

  is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec
  isd = g%isd ; ied = g%ied ; jsd = g%jsd ; jed = g%jed
  isdb = g%IsdB ; iedb = g%IedB ; jsdb = g%JsdB ; jedb = g%JedB
  i1off = g%idg_offset ; j1off = g%jdg_offset

  m_to_l = 1.0 ;  l_to_m = 1.0 

  dx_everywhere = 1000.0 * g%len_lon / (real(niglobal))
  dy_everywhere = 1000.0 * g%len_lat / (real(njglobal))
  i_dx = 1.0 / dx_everywhere ; i_dy = 1.0 / dy_everywhere

  !------------- 2. set grid space "set_grid_metrics_cartesian" -------------

  allocate( grid_latt(G%jsd:G%jed) )
  allocate( grid_lonT(G%isd:G%ied) )

  do j = jsd,jed 
    grid_latt(j) = g%south_lat + g%len_lat*(real(j+j1off-g%jsg)+0.5)/real(njglobal)
  enddo
  do i = isd,ied
    grid_lont(i) = g%west_lon + g%len_lon*(real(i+i1off-g%isg)+0.5)/real(niglobal)
  enddo
  print*,'grid_lont(55)', grid_lont(isd:isd+6)


  do j=jsdb,jedb
    do i=isdb,iedb
      g%dxBu(i,j) = m_to_l*dx_everywhere ; g%IdxBu(i,j) = l_to_m*i_dx
      g%dyBu(i,j) = m_to_l*dy_everywhere ; g%IdyBu(i,j) = l_to_m*i_dy
      g%areaBu(i,j) = m_to_l**2*dx_everywhere * dy_everywhere 
      g%IareaBu(i,j) = l_to_m**2*i_dx * i_dy
    enddo
  enddo

  do j=jsd,jed
   	do i=isd,ied
      g%geoLonT(i,j) = grid_lont(i) ; g%geoLatT(i,j) = grid_latt(j)
      g%dxT(i,j) = m_to_l*dx_everywhere ; g%IdxT(i,j) = l_to_m*i_dx
      g%dyT(i,j) = m_to_l*dy_everywhere ; g%IdyT(i,j) = l_to_m*i_dy
      g%areaT(i,j) = m_to_l**2*dx_everywhere * dy_everywhere 
      g%IareaT(i,j) = l_to_m**2*i_dx * i_dy
  	enddo
  enddo

  do j=jsd,jed
   	do i=isdb,iedb
      g%dxCu(i,j) = m_to_l*dx_everywhere ; g%IdxCu(i,j) = l_to_m*i_dx
      g%dyCu(i,j) = m_to_l*dy_everywhere ; g%IdyCu(i,j) = l_to_m*i_dy
  	enddo
  enddo

  do j=jsdb,jedb
   	do i=isd,ied
      g%dxCv(i,j) = m_to_l*dx_everywhere ; g%IdxCv(i,j) = l_to_m*i_dx
      g%dyCv(i,j) = m_to_l*dy_everywhere ; g%IdyCv(i,j) = l_to_m*i_dy
   	enddo
  enddo

  !------------- 3. set derived grid "set_derived_dyn_horgrid" -------------
  do j=jsd,jed
   	do i=isd,ied
   		if (g%dxT(i,j) < 0.0) g%dxT(i,j) = 0.0
     	if (g%dyT(i,j) < 0.0) g%dyT(i,j) = 0.0
        g%IdxT(i,j) = adcroft_reciprocal(g%dxT(i,j))
        g%IdyT(i,j) = adcroft_reciprocal(g%dyT(i,j))
        g%IareaT(i,j) = adcroft_reciprocal(g%areaT(i,j))
  	enddo
  enddo

  do j=jsd,jed
   	do i=isdb,iedb
      if (g%dxCu(i,j) < 0.0) g%dxCu(i,j) = 0.0
      if (g%dyCu(i,j) < 0.0) g%dyCu(i,j) = 0.0
      g%IdxCu(i,j) = adcroft_reciprocal(g%dxCu(i,j))
      g%IdyCu(i,j) = adcroft_reciprocal(g%dyCu(i,j))
  	enddo
  enddo

  do j=jsdb,jedb
   	do i=isd,ied
      if (g%dxCv(i,j) < 0.0) g%dxCv(i,j) = 0.0
      if (g%dyCv(i,j) < 0.0) g%dyCv(i,j) = 0.0
      g%IdxCv(i,j) = adcroft_reciprocal(g%dxCv(i,j))
      g%IdyCv(i,j) = adcroft_reciprocal(g%dyCv(i,j))
   	enddo
  enddo

  do j=jsdb,jedb
    do i=isdb,iedb
      if (g%dxBu(i,j) < 0.0) g%dxBu(i,j) = 0.0
      if (g%dyBu(i,j) < 0.0) g%dyBu(i,j) = 0.0
  
      g%IdxBu(i,j) = adcroft_reciprocal(g%dxBu(i,j))
      g%IdyBu(i,j) = adcroft_reciprocal(g%dyBu(i,j))
      ! areaBu has usually been set to a positive area elsewhere.
      if (g%areaBu(i,j) <= 0.0) g%areaBu(i,j) = g%dxBu(i,j) * g%dyBu(i,j)
      g%IareaBu(i,j) =  adcroft_reciprocal(g%areaBu(i,j))
    enddo
  enddo

  !------------- 4. set bathymetry "initialize_topography_named" (only flat works!) -------------
  do j=js,je
    do i=is,ie
      g%bathyT(i,j) = max_depth 
  	enddo
  enddo
  ! This is here just for safety.  Hopefully it doesn't do anything.
  do j=js,je
    do i=is,ie
      if (g%bathyT(i,j) > max_depth) g%bathyT(i,j) = max_depth
      if (g%bathyT(i,j) < min_depth) g%bathyT(i,j) = 0.5*min_depth
  	enddo
  enddo

  !------------- 5. masks "initialize_masks" -------------
  g%mask2dCu(:,:) = 0.0 ; g%mask2dCv(:,:) = 0.0 ; g%mask2dBu(:,:) = 0.0
  dmin = min_depth

  ! Construct the h-point or T-point mask
  do j=jsd,jed
   	do i=isd,ied
      if (g%bathyT(i,j) <= dmin) then
        g%mask2dT(i,j) = 0.0
      else
        g%mask2dT(i,j) = 1.0
     endif
  	enddo
  enddo

  do j=jsd,jed
    do i=isd,ied-1
     if ((g%bathyT(i,j) <= dmin) .or. (g%bathyT(i+1,j) <= dmin)) then
       g%mask2dCu(i,j) = 0.0
     else
       g%mask2dCu(i,j) = 1.0
     endif
    enddo
  enddo

  do j=jsd,jed-1
    do i=isd,ied
      if ((g%bathyT(i,j) <= dmin) .or. (g%bathyT(i,j+1) <= dmin)) then
        g%mask2dCv(i,j) = 0.0
      else
        g%mask2dCv(i,j) = 1.0
     endif
    enddo
  enddo 

  do j=jsd,jed-1 
  	do i=isd,ied-1
     	if ((g%bathyT(i+1,j) <= dmin) .or. (g%bathyT(i+1,j+1) <= dmin) .or. &
         (g%bathyT(i,j) <= dmin) .or. (g%bathyT(i,j+1) <= dmin)) then
        g%mask2dBu(i,j) = 0.0
     	else
        g%mask2dBu(i,j) = 1.0
     	endif
  	enddo
  enddo

  ! for dyCu and dxCv
  do j=jsd,jed
   	do i=isdb,iedb
      g%dy_Cu(i,j) = g%mask2dCu(i,j) * g%dyCu(i,j)
      g%areaCu(i,j) = g%dxCu(i,j) * g%dy_Cu(i,j)
      g%IareaCu(i,j) = g%mask2dCu(i,j) * adcroft_reciprocal(g%areaCu(i,j))
  	enddo
  enddo

  do j=jsdb,jedb
   	do i=isd,ied
      g%dx_Cv(i,j) = g%mask2dCv(i,j) * g%dxCv(i,j)
      g%areaCv(i,j) = g%dyCv(i,j) * g%dx_Cv(i,j)
      g%IareaCv(i,j) = g%mask2dCv(i,j) * adcroft_reciprocal(g%areaCv(i,j))
   	enddo
  enddo


end subroutine build_grid_cartesian


!=====================================================================
!  SUBTOUTINE: 
!=====================================================================
! 
! build grid from given nc file 
! 
subroutine build_grid_file(g, fnm_geom, fnm_vert)
  type(ocean_grid_type), intent(inout) :: g         
  character*200, intent(in) :: fnm_geom, fnm_vert ! geom/vert file, from which g is read
  ! 
  integer :: nih, njh, niq, njq, nz, ncid, i, j

  ! get dimensions
  call get_dims_geom(fnm_geom, nih, njh, niq, njq)
  call get_dims_vert(fnm_vert, nz)
  ! calc index (from dims) and allocate g
  call calc_index(g, nih, njh)
  call allocate_metrics(g)

  ! read from geom file and write into g
  call get_vars_geom(g, fnm_geom)
  ! write vert into g
  g%ke = nz

  ! some additional fields needs to be calc
  do j=g%jsd,g%jed
   	do i=g%isd,g%ied
        g%areaT(i,j) = g%dxT(i,j) * g%dyT(i,j)
        g%IareaT(i,j) = adcroft_reciprocal(g%areaT(i,j))
        g%IdxT(i,j) = adcroft_reciprocal(g%dxT(i,j))
        g%IdyT(i,j) = adcroft_reciprocal(g%dyT(i,j))
  	enddo
  enddo
  do j=g%jsd,g%jed
   	do i=g%isdb,g%iedb
       g%areaCu(i,j) = g%dxCu(i,j) * g%dyCu(i,j)
  	enddo
  enddo
  do j=g%jsdb,g%jedb
   	do i=g%isd,g%ied
       g%areaCv(i,j) = g%dxCv(i,j) * g%dyCv(i,j)
   	enddo
  enddo
  do j=g%jsdb,g%jedb
    do i=g%isdb,g%iedb
       g%areaBu(i,j) = g%dxBu(i,j) * g%dyBu(i,j)
    enddo
  enddo

end subroutine build_grid_file

!> Adcroft_reciprocal(x) = 1/x for |x|>0 or 0 for x=0.
function adcroft_reciprocal(val) result(I_val)
  real, intent(in) :: val  !< The value being inverted.
  real :: i_val            !< The Adcroft reciprocal of val.
  
  i_val = 0.0 ; if (val /= 0.0) i_val = 1.0/val
end function adcroft_reciprocal

end module mod_grid_my