!> This module contains the tracer_registry_type and the subroutines
!! that handle registration of tracers and related subroutines.
!! The primary subroutine, register_tracer, is called to indicate the
!! tracers advected and diffused.
! 
! Also does initialization!!
! 

module mom_tracer_registry

use mod_grid_my, only : ocean_grid_type
use mod_nc_wr_rd, only : get_ncvar

implicit none ; private

public call_tracer_register, initialize_tracer

 !> The tracer type
 type, public :: tracer_type
  
   real, allocatable, dimension(:,:,:) :: t     !< tracer concentration array [conc]
   ! real, allocatable, dimension(:,:,:) :: forc     !< tracer forc array [c*m/s]
   real, allocatable, dimension(:,:,:) :: ad_x           
   real, allocatable, dimension(:,:,:) :: ad_y        !< diagnostic array for y-advective tracer flux
                                                               !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
   real, allocatable, dimension(:,:)  :: ad2d_x     ! < diagnostic vertical sum x-advective tracer flux
                                                               !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
   real, allocatable, dimension(:,:) :: ad2d_y     !< diagnostic vertical sum y-advective tracer flux
  
   real, allocatable, dimension(:,:,:) :: df_x       !< diagnostic array for x-diffusive tracer flux
                                                               !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
   real, allocatable, dimension(:,:,:) :: df_y       !< diagnostic array for y-diffusive tracer flux
                                                               !! [conc H L2 T-1 ~> conc m3 s-1 or conc kg s-1]
   real, allocatable, dimension(:,:,:) :: advection_xy   !< convergence of lateral advective tracer fluxes
 
   character(len=5) :: name     !< name for each tracer, e.g. 'tr2'

 end type tracer_type
  
 !> Type to carry basic tracer information
 type, public :: tracer_registry_type
   integer              :: ntr = 0       !< number of registered tracers
   ! character(len=5), allocatable, dimension(:)    :: cnames
   type(tracer_type)    :: tr(10)   !< array of registered tracers
 end type tracer_registry_type

contains

!===========================================================
! call_tracer_register "mom_tracer_flow_control.F90"
!===========================================================
 !> This subroutine determines which tracer packages are to be used and does the calls to
 !! register their tracers to be advected, diffused, and read from restarts.
subroutine call_tracer_register(G, tr_reg, NTR, cnames) ! tr_Reg
   type(ocean_grid_type),        intent(in) :: G         !< A horizontal index type structure.
   type(tracer_registry_type), intent(inout)   :: tr_reg     !< A pointer that is set to point to the
                                                          !! control structure for the tracer
                                                          !! advection and diffusion module.
   integer, intent(in)  :: NTR
   character(len=*), intent(in)  :: cnames(NTR) ! names of tracers to be solved

   ! local vars
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke) :: tr_ptr
   integer :: isd, ied, jsd, jed, nz, m
   isd = g%isd ; ied = g%ied ; jsd = g%jsd ; jed = g%jed ; nz = g%ke

   ! allocate( tr_reg%cnames(NTR) )
   do m = 1,ntr
      tr_reg%ntr = tr_reg%ntr + 1
      allocate( tr_reg%tr(m)%t(g%isd:g%ied, g%jsd:g%jed, g%ke) )
      tr_reg%tr(m)%name = cnames(m)
     ! Register the tracer for horizontal advection, diffusion, and restarts.
     ! call register_tracer(tr_ptr, tr_reg, g)
   enddo
   ! print*, 'tr_reg%tr(2)%name= ', tr_reg%tr(2)%name
  
end subroutine call_tracer_register


!===========================================================
!  modified from  user_initialize_tracer in "user_tracer_example.f90"
!===========================================================
 !> This subroutine initializes the NTR tracer fields in reg
 !! and it sets up the tracer output.
 subroutine initialize_tracer(G, tr_reg, tracer_IC_file)
   type(ocean_grid_type),              intent(in) :: g    !< The ocean's grid structure
   type(tracer_registry_type)         :: tr_reg   !< The control structure returned by a previous
                                                          !! call to USER_register_tracer_example.
   character*100, optional, intent(in) :: tracer_IC_file
 ! Local variables
   character(len=32) :: name     ! A variable's name in a NetCDF file.
   integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
   integer :: isdb, iedb, jsdb, jedb, ntr
  
   ! if (.not.associated(tr_reg)) return

   is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec ; nz = g%ke
   isd = g%isd ; ied = g%ied ; jsd = g%jsd ; jed = g%jed
   isdb = g%IsdB ; iedb = g%IedB ; jsdb = g%JsdB ; jedb = g%JedB
  
   ntr = tr_reg%ntr ! Avoids compile-time warning when NTR<2
   if (present(tracer_IC_file)) then
 !  Read the tracer concentrations from a netcdf file.
      do m=1,ntr
         name = trim(tr_reg%tr(m)%name)
         ! print*, name, '-1 %t(22,22,1): ', tr_reg%tr(m-1)%t(22,22,1)
         call get_ncvar(tracer_IC_file, name, tr_reg%tr(m)%t, g%nihalo, g%njhalo)
         print*, 'Read ', trim(name), ' from ', trim(tracer_IC_file)
      enddo
      ! print*, 'tr5 %t(22,22,1): ', tr_reg%tr(5)%t(22,22,1)
   else
      do m=1,ntr
         do k=1,nz 
            do j=js,je
               do i=is,ie
                  tr_reg%tr(m)%t(i,j,k) = m*1.0e-20 
               enddo
            enddo
         enddo
      enddo
   endif
   
   ! do m=1,ntr
   !    print*, m, '%t(22,22,1): ', tr_reg%tr(m)%t(22,22,1)
   ! enddo

 end subroutine initialize_tracer


end 
