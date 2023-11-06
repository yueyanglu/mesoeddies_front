module mod_nc_wr_rd
! This module saves subroutines that write and read nc. files, based on
! the 'netcdf' package.

use netcdf

implicit none; public 

contains

!=====================================================================
!  SUBTOUTINE:  get_dims_geom
!=====================================================================
! 
! get the lengths of all dimensions in the GEOMETRY file
! 
subroutine get_dims_geom(fnm, nih, njh, niq, njq)
  character*200, intent(in) :: fnm
  integer, intent(out) :: nih, njh, niq, njq
  ! 
  integer :: ncid, lonh_id, lath_id, lonq_id, latq_id
  ! name of all the dimensions
  character (len = *), parameter :: LONH_NAME = "lonh"
  character (len = *), parameter :: LATH_NAME = "lath"
  character (len = *), parameter :: LONQ_NAME = "lonq"
  character (len = *), parameter :: LATQ_NAME = "latq"

  ! open
  call check( nf90_open(fnm, nf90_nowrite, ncid) )
  ! get dimension IDs from the dimension names
  call check( nf90_inq_dimid(ncid, LONH_NAME, lonh_id) )
  call check( nf90_inq_dimid(ncid, LATH_NAME, lonh_id) )
  call check( nf90_inq_dimid(ncid, LONQ_NAME, latq_id) )
  call check( nf90_inq_dimid(ncid, LATQ_NAME, latq_id) )
  ! get lens of each dimension from dimID
  call check( nf90_inquire_dimension(ncid, lonh_id, len = nih) )
  call check( nf90_inquire_dimension(ncid, lonh_id, len = njh) )
  call check( nf90_inquire_dimension(ncid, latq_id, len = niq) )
  call check( nf90_inquire_dimension(ncid, latq_id, len = njq) )
  ! close
  call check( nf90_close(ncid) )

end subroutine get_dims_geom

!=====================================================================
!  SUBTOUTINE:  get_dims_vert
!=====================================================================
! 
! get the lengths of all dimensions in the GEOMETRY file
! 
subroutine get_dims_vert(fnm, nz)
  character*200, intent(in) :: fnm
  integer, intent(out) :: nz
  ! 
  integer :: ncid, layer_id
  ! name of all the dimensions
  character (len = *), parameter :: LAY_NAME = "Layer"

  ! open
  call check( nf90_open(fnm, nf90_nowrite, ncid) )
  ! get dimension IDs from the dimension names
  call check( nf90_inq_dimid(ncid, LAY_NAME, layer_id) )
  ! get lens of each dimension from dimID
  call check( nf90_inquire_dimension(ncid, layer_id, len = nz) )
  ! close
  call check( nf90_close(ncid) )

end subroutine get_dims_vert

!=====================================================================
!  read var
!=====================================================================
! 
! Get one var from the nc file, using the VAR_NAME and the known sizes
! of the var!!!
! Note that the halos are dealt with!!!!
! 
subroutine get_ncvar(fnm, VAR_NAME, var, nihalo, njhalo)
  character(len = *), intent(in) :: fnm, VAR_NAME
  real, dimension(:,:,:), intent(inout) ::  var ! var must be 3D
  integer, intent(in) :: nihalo, njhalo

  ! local vars
  real, allocatable, dimension(:,:,:) ::  var_read  ! var in the file, size may be 
                                       ! different from "var"
  integer :: ncid, varid, recid, nrecs ! len of time coord in the file
  integer :: ni, nj, nz, ni_sub, nj_sub, nz_sub
  integer :: is, ie, js, je,  is_sub, ie_sub, js_sub, je_sub
  character(len = *), parameter :: REC_NAME = "Time"
  integer :: start(4), count_var(4) ! in nc, [x,y,z,t]

  ! the given var's sizes
  ni = size(var,1)
  nj = size(var,2)
  nz = size(var,3)
  is = lbound(var, 1); ie = ubound(var, 1)
  js = lbound(var, 2); je = ubound(var, 2)
  is_sub = is + nihalo; ie_sub = ie - nihalo
  js_sub = js + njhalo; je_sub = je - njhalo
  ! print*, is, ie, js, je 
  ! print*, is_sub, ie_sub, js_sub, je_sub

  ! size of var in the saved file
  ni_sub = ni - 2*nihalo
  nj_sub = nj - 2*njhalo
  nz_sub = nz
  allocate( var_read(ni_sub, nj_sub, nz_sub) )

  ! open 
  call check( nf90_open(fnm, nf90_nowrite, ncid) )
  ! varID
  call check( nf90_inq_dimid(ncid, REC_NAME, recid) )
  call check( nf90_inq_varid(ncid, VAR_NAME, varid) )
  ! check the time coordinate !!!
  call check( nf90_inquire_dimension(ncid, recid, len = nrecs) )
  if(nrecs /= 1) print*, "WARNING: ", trim(fnm), "has ",nrecs,&
     "record, but only the first is read !!!!!!!!!!!!"

 ! Read 1 record of NX*NY*NZ values, starting at the beginning 
 ! of the record (the (1, 1, 1, rec) element in the netCDF file).
  start = (/ 1, 1, 1, 1 /)
  count_var = (/ ni_sub, nj_sub, nz_sub, 1 /)
  call check( nf90_get_var(ncid, varid, var_read, start = start, &
                              count = count_var) )
  ! 
  call check( nf90_close(ncid) )

  ! write var_read into var, note the halos
  var(is_sub:ie_sub, js_sub:je_sub, 1:nz_sub) = var_read(:,:,:)     

  ! fill the missing values
  call fill_miss(var, 1.0E15, 0.0)

end subroutine get_ncvar

!=====================================================================
! 
!=====================================================================
! Write a 3D data into nc file: ni x nj x nz grid
! Var(1+nihalo:end-nihalo, 1+njhalo:end-njhalo,:) are used!
! 
! 
subroutine wr_ncfile(fnm, VAR_NAME, var, nihalo, njhalo)
  character(len = *), intent(in) :: fnm, VAR_NAME
  real, dimension(:,:,:), intent(in) ::  var ! var must be 4D 
  integer, intent(in) :: nihalo, njhalo
  ! local vars
  integer :: ncid, rec ! len of time coord in the file
  integer :: ni, nj, nz, ni_sub, nj_sub
  integer :: is, ie, js, je,  is_sub, ie_sub, js_sub, je_sub
  character (len = *), parameter :: X_NAME = "x"
  character (len = *), parameter :: Y_NAME = "y"
  character (len = *), parameter :: Z_NAME = "zl"
  character (len = *), parameter :: REC_NAME = "Time"
  integer, parameter :: NDIMS = 4, nrecs = 1

  integer :: x_dimid, y_dimid, z_dimid, rec_dimid

  ! The start and count arrays will tell the netCDF library where to
  ! write our data.
  integer :: start(NDIMS), count(NDIMS)

  ! Info for the netCDF variables
  integer :: varid, dimids(NDIMS)
  real, allocatable, dimension(:,:,:) :: var_save
  logical :: ifexist

  ! the given var's sizes
  ni = size(var,1)
  nj = size(var,2)
  nz = size(var,3)
  is = lbound(var, 1); ie = ubound(var, 1)
  js = lbound(var, 2); je = ubound(var, 2)
  is_sub = is + nihalo; ie_sub = ie - nihalo
  js_sub = js + njhalo; je_sub = je - njhalo
  ! size of var to be saved
  ni_sub = ni - 2*nihalo
  nj_sub = nj - 2*njhalo

  ! allocate( var_save(nz,nj_sub,ni_sub) )
  allocate( var_save(ni_sub,nj_sub,nz) )
  var_save(:,:,:) = var(is_sub:ie_sub,js_sub:je_sub,:)


  count = (/ ni_sub, nj_sub, nz, 1 /)
  start = (/ 1, 1, 1, 1 /)

  ! ------------------ begin save ----------------
  inquire(FILE = fnm, exist = ifexist)
  if( .not. ifexist ) then

    ! Create a new file. 
    call check( nf90_create(path=fnm, cmode=nf90_noclobber, ncid=ncid) )

      ! Define dimensions. 
      call check( nf90_def_dim(ncid, X_NAME, ni_sub, x_dimid) )
      call check( nf90_def_dim(ncid, Y_NAME, nj_sub, y_dimid) )
      call check( nf90_def_dim(ncid, Z_NAME, nz, z_dimid) )
      call check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )

      ! The dimids array is used to pass the dimids of the dimensions of
      ! the netCDF variables. Both of the netCDF variables we are creating
      ! share the same four dimensions. 
      dimids = (/ x_dimid, y_dimid, z_dimid, rec_dimid /)

      ! Define variables.
      call check( nf90_def_var(ncid, VAR_NAME, NF90_DOUBLE, dimids, varid) )

    ! End define mode.
    call check( nf90_enddef(ncid) )  

  else 
    ! if already exists, open file, def and add new vars 
    call check( nf90_open(path=fnm, mode=NF90_WRITE, ncid=ncid) )

    ! put it into define mode
    call check( nf90_redef(ncid) )

      ! inqire existing dims
      call check( nf90_inq_dimid(ncid, X_NAME, x_dimid) )
      call check( nf90_inq_dimid(ncid, Y_NAME, y_dimid) )
      call check( nf90_inq_dimid(ncid, Z_NAME, z_dimid) )
      call check( nf90_inq_dimid(ncid, REC_NAME, rec_dimid) )
      dimids = (/ x_dimid, y_dimid, z_dimid, rec_dimid /)

      ! define additional vars
      call check( nf90_def_var(ncid, VAR_NAME, NF90_DOUBLE, dimids, varid) )

    ! End define mode.
    call check( nf90_enddef(ncid) )  

  endif

  ! data mode -- put data in nc file
  do rec = 1, NRECS
     start(4) = rec
     call check( nf90_put_var(ncid, varid, var_save, start = start, &
                              count = count) )
  enddo

  ! close file
  call check( nf90_close(ncid) )

end subroutine wr_ncfile

!!=====================================================================
!  set values that are >= "miss" to "fill"
!====================================================================
subroutine fill_miss(var,miss, fill)
    real, dimension(:,:,:), intent(inout) ::  var 
    real, intent(in) ::  miss
    real, intent(in) :: fill 
    ! 
    integer :: is, ie, js, je, nz
    integer :: i,j,k

    is = lbound(var, 1); ie = ubound(var, 1)
    js = lbound(var, 2); je = ubound(var, 2)
    nz = size(var,3)

    do k = 1,nz
      do j = js,je
        do i = is,ie
          if (var(i,j,k) >= miss ) then
            var(i,j,k) = fill
          endif
        enddo
      enddo
    enddo
end subroutine fill_miss

!=====================================================================
!  SUBTOUTINE: check
!=====================================================================
subroutine check(status)
   integer, intent ( in) :: status
    
   if(status /= nf90_noerr) then 
     print *, trim(nf90_strerror(status))
     stop "Stopped"
   endif
end subroutine check  

end module mod_nc_wr_rd