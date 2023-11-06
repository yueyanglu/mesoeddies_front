module mod_continuity
!> Solve the layer continuity equation, based on "mom_continuity.f90".

use mod_grid_my, only : ocean_grid_type

implicit none; private

public continuity_useuh

!> Control structure for mom_continuity_ppm
type, public :: continuity_ppm_cs ; private
  logical :: upwind_1st      
  logical :: monotonic       
  logical :: simple_2nd      
  real :: tol_eta            
  real :: tol_vel            
  real :: tol_eta_aux       
  real :: cfl_limit_adjust   
  logical :: aggress_adjust  
  logical :: vol_cfl         
  logical :: better_iter     
  logical :: use_visc_rem_max 
  logical :: marginal_faces  
end type continuity_ppm_cs

!> A container for loop bounds
type :: loop_bounds_type ; private
  integer :: ish, ieh, jsh, jeh
end type loop_bounds_type
  
contains

!===========================================================
!  continuity_ppm
!===========================================================
!> Time steps layer thicknesses using PPM scheme based on "continuity_ppm" of MOM6.
subroutine continuity_useuh(hin, h, uh, vh, dt, g)
  type(ocean_grid_type),   intent(inout)    :: g  
  real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke),   &
                           intent(in)    :: hin !< Initial layer thickness [H ~> m or kg m-2].
  real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke),   &
                           intent(inout) :: h   !< Final layer thickness.
  real, dimension(g%isdb:g%iedb, g%jsd:g%jed, g%ke), &
                           intent(out)   :: uh  !< Zonal volume flux, u*h*dy [m3 s-1].
  real, dimension(g%isd:g%ied, g%jsdb:g%jedb, g%ke), &
                           intent(out)   :: vh  !< Meridional volume flux, v*h*dx [m3 s-1].
  real,                    intent(in)    :: dt  
  
  ! Local variables
  real, parameter :: h_min = 1.0E-10  ! The minimum layer thickness [gv%Angstrom_H]
  type(loop_bounds_type) :: lb
  integer :: nz, stencil
  integer :: i, j, k, nj, ni
  logical :: x_first = 1

  nz = g%ke
  lb%ish = g%isc ; lb%ieh = g%iec ; lb%jsh = g%jsc ; lb%jeh = g%jec

  !!! ------  update h using uh and vh ------
  !$OMP PARALLEL DO PRIVATE(i,j,k)
  do k = 1, nz
    do j = lb%jsh, lb%jeh
      do i = lb%ish, lb%ieh
        h(i,j,k) = h(i,j,k) - dt * g%IareaT(i,j) * &
            ( uh(i,j,k) - uh(i-1,j,k) + vh(i,j,k) - vh(i,j-1,k) )
        !   Uncomment this line to prevent underflow.
        if (h(i,j,k) < h_min) h(i,j,k) = h_min
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  !!! ------  relaxation of layer thickness ------
  nj = LB%jeh - LB%jsh + 1
  ni = LB%ieh - LB%ish + 1
  call h_relax_user_tilt(h,LB,G,ni,nj, dt)
  print*, 'Relaxation appplied!'

end subroutine continuity_useuh

!===========================================================
!  zonal_mass_flux
!===========================================================
!> Calculates the mass or volume fluxes through the zonal faces, and other related quantities.
subroutine zonal_mass_flux(u, h_in, uh, dt, g, cs, lb)
  type(ocean_grid_type),   intent(inout) :: g    
  real, dimension(g%isdb:g%iedb, g%jsd:g%jed, g%ke), &
                           intent(in)    :: u   
  real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke),   &
                           intent(in)    :: h_in 
  real, dimension(g%isdb:g%iedb, g%jsd:g%jed, g%ke), &
                           intent(out)   :: uh   !< u*h*dy [m3 s-1]
  real,                    intent(in)    :: dt   !< Time increment [T ~> s].
  type(continuity_ppm_cs), pointer    :: CS   !< This module's control structure.
  type(loop_bounds_type),  intent(in)    :: lb   !< Loop bounds structure.

  ! Local variables
  real, dimension(g%isdb:g%iedb, g%jsd:g%jed) :: duhdu ! Partial derivative of uh with u [H L ~> m2 or kg m-1].
  real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke) :: h_L, h_R ! Left and right face thicknesses [H ~> m or kg m-2].
  real, dimension(g%isdb:g%iedb) :: &
   du, &      ! Corrective barotropic change in the velocity [L T-1 ~> m s-1].
   du_min_CFL, & ! Min/max limits on du correction
   du_max_CFL, & ! to avoid CFL violations
   duhdu_tot_0, & ! Summed partial derivative of uh with u [H L ~> m2 or kg m-1].
   uh_tot_0, & ! Summed transport with no barotropic correction [H L2 T-1 ~> m3 s-1 or kg s-1].
   visc_rem_max  ! The column maximum of visc_rem.
  logical, dimension(g%isdb:g%iedb) :: do_I
  real, dimension(g%isdb:g%iedb,g%ke) :: &
   visc_rem      ! A 2-D copy of visc_rem_u or an array of 1's.
  real, dimension(g%isdb:g%iedb) :: FAuI  ! A list of sums of zonal face areas [H L ~> m2 or kg m-1].
  real :: FA_u    ! A sum of zonal face areas [H m ~> m2 or kg m-1].
  real :: I_vrm   ! 1.0 / visc_rem_max, nondim.
  real :: CFL_dt  ! The maximum CFL ratio of the adjusted velocities divided by
                 ! the time step [T-1 ~> s-1].
  real :: I_dt    ! 1.0 / dt [T-1 ~> s-1].
  real :: du_lim  ! The velocity change that give a relative CFL of 1 [L T-1 ~> m s-1].
  real :: dx_E, dx_W ! Effective x-grid spacings to the east and west [L ~> m].
  integer :: i, j, k, ish, ieh, jsh, jeh, n, nz
  logical :: local_specified_BC, use_visc_rem, set_BT_cont, any_simple_OBC
  logical :: local_Flather_OBC, local_open_BC, is_simple
  real, parameter :: h_min = 1.0E-10 ! gv%Angstrom_H

  local_specified_bc = .false. ; set_bt_cont = .false. ; local_flather_obc = .false.
  local_open_bc = .false.; use_visc_rem = .false.

  ish = lb%ish ; ieh = lb%ieh ; jsh = lb%jsh ; jeh = lb%jeh ; nz = g%ke
  cfl_dt = cs%CFL_limit_adjust / dt
  i_dt = 1.0 / dt
  if (cs%aggress_adjust) cfl_dt = i_dt

  ! print*, 'ish, ieh, jsh, jeh = ', ish, ieh, jsh, jeh
  ! print*, 'cs%upwind_1st, monotonic, simple_2nd', cs%upwind_1st, cs%monotonic, cs%simple_2nd

!---------------- Sets h_L and h_R in each layer ----------------
!$OMP parallel do default(none) shared(ish,ieh,jsh,jeh,nz,CS,h_L,h_in,h_R,G,LB,visc_rem)
  do k = 1,nz
    ! print*, cs%upwind_1st
    if (cs%upwind_1st) then
      do j=jsh,jeh
       	do i=ish-1,ieh+1
          h_l(i,j,k) = h_in(i,j,k) ; h_r(i,j,k) = h_in(i,j,k)
      	enddo
      enddo 
    else
      call ppm_reconstruction_x(h_in(:,:,k), h_l(:,:,k), h_r(:,:,k), g, lb, &
                                 2.0*h_min, cs%monotonic, simple_2nd=cs%simple_2nd)
    endif
    do i=ish-1,ieh ; visc_rem(i,k) = 1.0 ; enddo
  enddo ! k

!---------------- calc zonal flux, each j and each k (do a lat band at a time) ----------------
!$OMP parallel do default(none) shared(ish,ieh,jsh,jeh,nz,u,h_in,h_L,h_R,uh,dt,G,CS) &
!$OMP                          private(do_I,duhdu,visc_rem_max) &
!$OMP      firstprivate(visc_rem)
  do j = jsh, jeh
    do i=ish-1,ieh ; do_i(i) = .true. ; visc_rem_max(i) = 0.0 ; enddo
    ! Set uh and duhdu.
    do k = 1,nz

      call zonal_flux_layer(u(:,j,k), h_in(:,j,k), h_l(:,j,k), h_r(:,j,k), &
                            uh(:,j,k), duhdu(:,k), visc_rem(:,k), &
                            dt, g, j, ish, ieh, do_i, cs%vol_CFL)
    enddo !k
  enddo !j-loop

end subroutine zonal_mass_flux

!=====================================================================
!  DONE
!=====================================================================
! > Evaluates the zonal mass or volume fluxes in a layer.
subroutine zonal_flux_layer(u, h, h_L, h_R, uh, duhdu, visc_rem, dt, G, j, &
                             ish, ieh, do_I, vol_CFL)
  type(ocean_grid_type),        intent(inout) :: G        
  real, dimension(g%isdb:g%iedb),    intent(in)   :: u   
  real, dimension(g%isdb:g%iedb),    intent(in)    :: visc_rem    
  real, dimension(g%isd:g%ied),     intent(in)    :: h        
  real, dimension(g%isd:g%ied),     intent(in)    :: h_L 
  real, dimension(g%isd:g%ied),     intent(in)    :: h_R      
  real, dimension(g%isdb:g%iedb),    intent(inout) :: uh       
  real, dimension(g%isdb:g%iedb),    intent(inout) :: duhdu   
  real,                         intent(in)    :: dt  !< Time increment [T ~> s].
  integer,                      intent(in)    :: j        !< Spatial index.
  integer,                      intent(in)    :: ish      !< Start of index range.
  integer,                      intent(in)    :: ieh      !< End of index range.
  logical, dimension(g%isdb:g%iedb), intent(in)    :: do_I     !< Which i values to work on.
  logical,                      intent(in)    :: vol_CFL  
  real :: CFL  ! The CFL number based on the local velocity and grid spacing [nondim]
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                  ! with the same units as h_in.
  real :: h_marg ! The marginal thickness of a flux [H ~> m or kg m-2].
  integer :: i
  
   do i=ish-1,ieh 
    if (do_i(i)) then
     ! Set new values of uh and duhdu.
     if (u(i) > 0.0) then
       if (vol_cfl) then ; cfl = (u(i) * dt) * (g%dy_Cu(i,j) * g%IareaT(i,j))
       else ; cfl = u(i) * dt * g%IdxT(i,j) ; endif
       curv_3 = h_l(i) + h_r(i) - 2.0*h(i)
       uh(i) = g%dy_Cu(i,j) * u(i) * &
           (h_r(i) + cfl * (0.5*(h_l(i) - h_r(i)) + curv_3*(cfl - 1.5)))
       h_marg = h_r(i) + cfl * ((h_l(i) - h_r(i)) + 3.0*curv_3*(cfl - 1.0))
     elseif (u(i) < 0.0) then
       if (vol_cfl) then ; cfl = (-u(i) * dt) * (g%dy_Cu(i,j) * g%IareaT(i+1,j))
       else ; cfl = -u(i) * dt * g%IdxT(i+1,j) ; endif
       curv_3 = h_l(i+1) + h_r(i+1) - 2.0*h(i+1)
       uh(i) = g%dy_Cu(i,j) * u(i) * &
           (h_l(i+1) + cfl * (0.5*(h_r(i+1)-h_l(i+1)) + curv_3*(cfl - 1.5)))
       h_marg = h_l(i+1) + cfl * ((h_r(i+1)-h_l(i+1)) + 3.0*curv_3*(cfl - 1.0))
     else
       uh(i) = 0.0
       h_marg = 0.5 * (h_l(i+1) + h_r(i))
     endif
     duhdu(i) = g%dy_Cu(i,j) * h_marg * visc_rem(i)
   	endif 
   enddo
   
end subroutine zonal_flux_layer
  

!=====================================================================
!  ongoing
!=====================================================================
subroutine meridional_mass_flux(v, h_in, vh, dt, G, CS, LB)
  type(ocean_grid_type),   intent(inout) :: g    
  real, dimension(g%isd:g%ied, g%jsdb:g%jedb, g%ke), &
                           intent(in)    :: v   
  real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke),   &
                           intent(in)    :: h_in 
  real, dimension(g%isd:g%ied, g%jsdb:g%jedb, g%ke), &
                           intent(out)   :: vh   !< u*h*dy [m3 s-1]
  real,                    intent(in)    :: dt   !< Time increment [T ~> s].
  type(continuity_ppm_cs), pointer    :: CS   !< This module's control structure.
  type(loop_bounds_type),  intent(in)    :: lb   !< Loop bounds structure.


  ! Local variables
  real, dimension(g%isd:g%ied, g%ke) :: &
     dvhdv      ! Partial derivative of vh with v [H L ~> m2 or kg m-1].
  real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke) :: &
     h_L, h_R   ! Left and right face thicknesses [H ~> m or kg m-2].
  real, dimension(g%isd:g%ied) :: &
     dv, &      ! Corrective barotropic change in the velocity [L T-1 ~> m s-1].
     dv_min_CFL, & ! Min/max limits on dv correction
     dv_max_CFL, & ! to avoid CFL violations
     dvhdv_tot_0, & ! Summed partial derivative of vh with v [H L ~> m2 or kg m-1].
     vh_tot_0, &   ! Summed transport with no barotropic correction [H L2 T-1 ~> m3 s-1 or kg s-1].
     visc_rem_max  ! The column maximum of visc_rem.
  logical, dimension(g%isd:g%ied) :: do_I
  real, dimension(g%isd:g%ied) :: FAvi  ! A list of sums of meridional face areas [H L ~> m2 or kg m-1].
  real :: FA_v    ! A sum of meridional face areas [H m ~> m2 or kg m-1].
  real, dimension(g%isd:g%ied,g%ke) :: &
     visc_rem      ! A 2-D copy of visc_rem_v or an array of 1's.
  real :: I_vrm   ! 1.0 / visc_rem_max, nondim.
  real :: CFL_dt  ! The maximum CFL ratio of the adjusted velocities divided by
                   ! the time step [T-1 ~> s-1].
  real :: I_dt    ! 1.0 / dt [T-1 ~> s-1].
  real :: dv_lim  ! The velocity change that give a relative CFL of 1 [L T-1 ~> m s-1].
  real :: dy_N, dy_S ! Effective y-grid spacings to the north and south [L ~> m].
  integer :: i, j, k, ish, ieh, jsh, jeh, n, nz
  logical :: local_specified_BC, use_visc_rem, set_BT_cont, any_simple_OBC
  logical :: local_Flather_OBC, is_simple, local_open_BC
  real, parameter :: h_min = 1.0E-10 ! gv%Angstrom_H

  local_specified_bc = .false. ; set_bt_cont = .false. ; local_flather_obc = .false.
  local_open_bc = .false.; use_visc_rem = .false.

  ish = lb%ish ; ieh = lb%ieh ; jsh = lb%jsh ; jeh = lb%jeh ; nz = g%ke
  cfl_dt = cs%CFL_limit_adjust / dt
  i_dt = 1.0 / dt
  if (cs%aggress_adjust) cfl_dt = i_dt

!---------------- Sets h_L and h_R in each layer ----------------
  do k = 1,nz
    if (cs%upwind_1st) then
      do j=jsh-1,jeh+1
       	do i=ish,ieh
          h_l(i,j,k) = h_in(i,j,k) ; h_r(i,j,k) = h_in(i,j,k)
      	enddo
      enddo 
    else
      call ppm_reconstruction_y(h_in(:,:,k), h_l(:,:,k), h_r(:,:,k), g, lb, &
                                 2.0*h_min, cs%monotonic, simple_2nd=cs%simple_2nd)
    endif
    do i=ish,ieh ; visc_rem(i,k) = 1.0 ; enddo
  enddo ! k

!---------------- calc meird flux, each j and each k (do a lat band at a time) ----------------
  do j=jsh-1,jeh
    do i=ish,ieh ; do_i(i) = .true. ; visc_rem_max(i) = 0.0 ; enddo
    ! set vh and dvhdv.
    do k=1,nz
    	if (use_visc_rem) then 
      ! .... 
      endif
      call merid_flux_layer(v(:,j,k), h_in(:,:,k), h_l(:,:,k), h_r(:,:,k), &
                             vh(:,j,k), dvhdv(:,k), visc_rem(:,k), &
                             dt, g, j, ish, ieh, do_i, cs%vol_CFL)
    enddo !k

  enddo !j-loop

end subroutine meridional_mass_flux


!=====================================================================
!  DONE
!=====================================================================
!> Evaluates the meridional mass or volume fluxes in a layer.
subroutine merid_flux_layer(v, h, h_L, h_R, vh, dvhdv, visc_rem, dt, G, J, &
                             ish, ieh, do_I, vol_CFL)
  type(ocean_grid_type),        intent(inout) :: G        
  real, dimension(g%isd:g%ied),    intent(in)   :: v   
  real, dimension(g%isd:g%ied),     intent(in)    :: visc_rem 
  real, dimension(g%isd:g%ied, g%jsd:g%jed),     intent(in)    :: h        
  real, dimension(g%isd:g%ied, g%jsd:g%jed),     intent(in)    :: h_L 
  real, dimension(g%isd:g%ied, g%jsd:g%jed),     intent(in)    :: h_R 
  real, dimension(g%isd:g%ied),    intent(inout) :: vh       
  real, dimension(g%isd:g%ied),    intent(inout) :: dvhdv  
  real,                         intent(in)    :: dt       !< Time increment [T ~> s].
  integer,                      intent(in)    :: j        !< Spatial index.
  integer,                      intent(in)    :: ish      !< Start of index range.
  integer,                      intent(in)    :: ieh      !< End of index range.
  logical, dimension(g%isdb:g%iedb), intent(in)    :: do_I     !< Which i values to work on.
  logical,                      intent(in)    :: vol_CFL  
  ! 
  real :: CFL  ! The CFL number based on the local velocity and grid spacing [nondim]
  real :: curv_3 ! A measure of the thickness curvature over a grid length,
                  ! with the same units as h_in.
  real :: h_marg ! The marginal thickness of a flux [H ~> m or kg m-2].
  integer :: i

  do i=ish,ieh 
    if (do_i(i)) then
    	if (v(i) > 0.0) then
       	if (vol_cfl) then ; cfl = (v(i) * dt) * (g%dx_Cv(i,j) * g%IareaT(i,j))
       	else ; cfl = v(i) * dt * g%IdyT(i,j) ; endif
        curv_3 = h_l(i,j) + h_r(i,j) - 2.0*h(i,j)
        vh(i) = g%dx_Cv(i,j) * v(i) * ( h_r(i,j) + cfl * &
           (0.5*(h_l(i,j) - h_r(i,j)) + curv_3*(cfl - 1.5)) )
        h_marg = h_r(i,j) + cfl * ((h_l(i,j) - h_r(i,j)) + &
                                   3.0*curv_3*(cfl - 1.0))
     	elseif (v(i) < 0.0) then
        if (vol_cfl) then ; cfl = (-v(i) * dt) * (g%dx_Cv(i,j) * g%IareaT(i,j+1))
        else ; cfl = -v(i) * dt * g%IdyT(i,j+1) ; endif
        curv_3 = h_l(i,j+1) + h_r(i,j+1) - 2.0*h(i,j+1)
        vh(i) = g%dx_Cv(i,j) * v(i) * ( h_l(i,j+1) + cfl * &
           (0.5*(h_r(i,j+1)-h_l(i,j+1)) + curv_3*(cfl - 1.5)) )
        h_marg = h_l(i,j+1) + cfl * ((h_r(i,j+1)-h_l(i,j+1)) + &
                                     3.0*curv_3*(cfl - 1.0))
     	else
        vh(i) = 0.0
        h_marg = 0.5 * (h_l(i,j+1) + h_r(i,j))
     	endif ! (v(i) > 0.0)
     dvhdv(i) = g%dx_Cv(i,j) * h_marg * visc_rem(i)
   	endif
 	enddo

end subroutine merid_flux_layer


!=====================================================================
!  DONE
!=====================================================================
!> Calculates left/right edge values for PPM reconstruction.
subroutine ppm_reconstruction_x(h_in, h_L, h_R, G, LB, h_min, monotonic, simple_2nd)
 type(ocean_grid_type),             intent(in)  :: G    !< Ocean's grid structure.
 real, dimension(g%isd:g%ied, g%jsd:g%jed),  intent(in)  :: h_in 
 real, dimension(g%isd:g%ied, g%jsd:g%jed),  intent(out) :: h_L  !< Left thickness in the reconstruction,
                                                        !! [H ~> m or kg m-2].
 real, dimension(g%isd:g%ied, g%jsd:g%jed),  intent(out) :: h_R  !< Right thickness in the reconstruction,
                                                        !! [H ~> m or kg m-2].
 type(loop_bounds_type),            intent(in)  :: LB   !< Active loop bounds structure.
 real,                              intent(in)  :: h_min !< The minimum thickness
                   !! that can be obtained by a concave parabolic fit.
 logical, optional,                 intent(in)  :: monotonic 
 logical, optional,                 intent(in)  :: simple_2nd 

 ! Local variables with useful mnemonic names.
 real, dimension(g%isd:g%ied, g%jsd:g%jed)  :: slp ! The slopes.
 real, parameter :: oneSixth = 1./6.
 real :: h_ip1, h_im1
 real :: dMx, dMn
 logical :: use_CW84, use_2nd
 character(len=256) :: mesg
 integer :: i, j, isl, iel, jsl, jel, n, stencil

 use_cw84 = .false. ; if (present(monotonic)) use_cw84 = monotonic
 use_2nd = .false. ; if (present(simple_2nd)) use_2nd = simple_2nd


 isl = lb%ish-1 ; iel = lb%ieh+1 ; jsl = lb%jsh ; jel = lb%jeh

 ! This is the stencil of the reconstruction, not the scheme overall.
 stencil = 2 ; if (use_2nd) stencil = 1

 if ((isl-stencil < g%isd) .or. (iel+stencil > g%ied)) then
   write(*,*) "ERROR stencil in ppm_reconstruction_x !!!!!"
 endif
 if ((jsl < g%jsd) .or. (jel > g%jed)) then
   write(*,*) "ERROR stencil in ppm_reconstruction_x !!!!!"
 endif

 if (use_2nd) then
   write(*,*) "ERROR: simple_2nd code not repared!!!!!"
   return
 else
   do j=jsl,jel 
    do i=isl-1,iel+1
     if ((g%mask2dT(i-1,j) * g%mask2dT(i,j) * g%mask2dT(i+1,j)) == 0.0) then
       slp(i,j) = 0.0
     else
       ! This uses a simple 2nd order slope.
       slp(i,j) = 0.5 * (h_in(i+1,j) - h_in(i-1,j))
       ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
       dmx = max(h_in(i+1,j), h_in(i-1,j), h_in(i,j)) - h_in(i,j)
       dmn = h_in(i,j) - min(h_in(i+1,j), h_in(i-1,j), h_in(i,j))
       slp(i,j) = sign(1.,slp(i,j)) * min(abs(slp(i,j)), 2. * min(dmx, dmn))
     endif
   	enddo 
   enddo

   do j=jsl,jel 
    do i=isl,iel
     ! Neighboring values should take into account any boundaries.  The 3
     ! following sets of expressions are equivalent.
     h_im1 = g%mask2dT(i-1,j) * h_in(i-1,j) + (1.0-g%mask2dT(i-1,j)) * h_in(i,j)
     h_ip1 = g%mask2dT(i+1,j) * h_in(i+1,j) + (1.0-g%mask2dT(i+1,j)) * h_in(i,j)
     ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
     h_l(i,j) = 0.5*( h_im1 + h_in(i,j) ) + onesixth*( slp(i-1,j) - slp(i,j) )
     h_r(i,j) = 0.5*( h_ip1 + h_in(i,j) ) + onesixth*( slp(i,j) - slp(i+1,j) )
   	enddo 
   enddo
 endif ! not use_2nd

 if (use_cw84) then
   write(*,*) "ERROR: use_cw84 code not repared!!!!!"
   return
 else
   call ppm_limit_pos(h_in, h_l, h_r, h_min, g, isl, iel, jsl, jel)
 endif

 return
end subroutine ppm_reconstruction_x

!=====================================================================
!  DONE
!=====================================================================
!> Calculates left/right edge values for PPM reconstruction.
subroutine ppm_reconstruction_y(h_in, h_L, h_R, G, LB, h_min, monotonic, simple_2nd)
  type(ocean_grid_type),             intent(in)  :: G    !< Ocean's grid structure.
  real, dimension(g%isd:g%ied, g%jsd:g%jed),  intent(in)  :: h_in 
  real, dimension(g%isd:g%ied, g%jsd:g%jed),  intent(out) :: h_L  !< Left thickness in the reconstruction,
                                                        !! [H ~> m or kg m-2].
  real, dimension(g%isd:g%ied, g%jsd:g%jed),  intent(out) :: h_R  !< Right thickness in the reconstruction,
                                                        !! [H ~> m or kg m-2].
  type(loop_bounds_type),            intent(in)  :: LB   !< Active loop bounds structure.
  real,                              intent(in)  :: h_min !< The minimum thickness
                   !! that can be obtained by a concave parabolic fit.
  logical, optional,                 intent(in)  :: monotonic 
  logical, optional,                 intent(in)  :: simple_2nd 

! Local variables with useful mnemonic names.
  real, dimension(g%isd:g%ied, g%jsd:g%jed)  :: slp ! The slopes.
  real, parameter :: oneSixth = 1./6.
  real :: h_jp1, h_jm1
  real :: dMx, dMn
  logical :: use_CW84, use_2nd
  integer :: i, j, isl, iel, jsl, jel, n, stencil

  use_cw84 = .false. ; if (present(monotonic)) use_cw84 = monotonic
  use_2nd = .false. ; if (present(simple_2nd)) use_2nd = simple_2nd
  
  isl = lb%ish ; iel = lb%ieh ; jsl = lb%jsh-1 ; jel = lb%jeh+1
  
  ! This is the stencil of the reconstruction, not the scheme overall.
  stencil = 2 ; if (use_2nd) stencil = 1
  
  if ((isl < g%isd) .or. (iel > g%ied)) then
    write(*,*) "ERROR stencil in ppm_reconstruction_x !!!!!"
  endif
  if ((jsl-stencil < g%jsd) .or. (jel+stencil > g%jed)) then
    write(*,*) "ERROR stencil in ppm_reconstruction_x !!!!!"
  endif
  
	if (use_2nd) then
    write(*,*) "ERROR: simple_2nd code not repared!!!!!"
    return
 	else
    do j=jsl-1,jel+1
      do i=isl,iel
        if ((g%mask2dT(i,j-1) * g%mask2dT(i,j) * g%mask2dT(i,j+1)) == 0.0) then
          slp(i,j) = 0.0
        else
          ! This uses a simple 2nd order slope.
          slp(i,j) = 0.5 * (h_in(i,j+1) - h_in(i,j-1))
         ! Monotonic constraint, see Eq. B2 in Lin 1994, MWR (132)
          dmx = max(h_in(i,j+1), h_in(i,j-1), h_in(i,j)) - h_in(i,j)
          dmn = h_in(i,j) - min(h_in(i,j+1), h_in(i,j-1), h_in(i,j))
          slp(i,j) = sign(1.,slp(i,j)) * min(abs(slp(i,j)), 2. * min(dmx, dmn))
        endif
      enddo 
    enddo
  
    do j=jsl,jel
      do i=isl,iel
       ! Neighboring values should take into account any boundaries.  The 3
       ! following sets of expressions are equivalent.
       h_jm1 = g%mask2dT(i,j-1) * h_in(i,j-1) + (1.0-g%mask2dT(i,j-1)) * h_in(i,j)
       h_jp1 = g%mask2dT(i,j+1) * h_in(i,j+1) + (1.0-g%mask2dT(i,j+1)) * h_in(i,j)
       ! Left/right values following Eq. B2 in Lin 1994, MWR (132)
       h_l(i,j) = 0.5*( h_jm1 + h_in(i,j) ) + onesixth*( slp(i,j-1) - slp(i,j) )
       h_r(i,j) = 0.5*( h_jp1 + h_in(i,j) ) + onesixth*( slp(i,j) - slp(i,j+1) )
     	enddo
    enddo 

  endif ! use_2nd
  
  if (use_cw84) then
     write(*,*) "ERROR: use_cw84 code not repared!!!!!"
     return
  else
     call ppm_limit_pos(h_in, h_l, h_r, h_min, g, isl, iel, jsl, jel)
  endif
  
  return
end subroutine ppm_reconstruction_y

!=====================================================================
!  DONE
!=====================================================================
!
!> This subroutine limits the left/right edge values of the PPM reconstruction
!! to give a reconstruction that is positive-definite.  Here this is
!! reinterpreted as giving a constant thickness if the mean thickness is less
!! than h_min, with a minimum of h_min otherwise.
subroutine ppm_limit_pos(h_in, h_L, h_R, h_min, G, iis, iie, jis, jie)
  type(ocean_grid_type),             intent(in)  :: G    !< Ocean's grid structure.
  real, dimension(g%isd:g%ied, g%jsd:g%jed),  intent(in)  :: h_in 
  real, dimension(g%isd:g%ied, g%jsd:g%jed),  intent(inout) :: h_L 
  real, dimension(g%isd:g%ied, g%jsd:g%jed),  intent(inout) :: h_R 
  real,                              intent(in)  :: h_min 
  integer,                           intent(in)  :: iis      !< Start of i index range.
  integer,                           intent(in)  :: iie      !< End of i index range.
  integer,                           intent(in)  :: jis      !< Start of j index range.
  integer,                           intent(in)  :: jie      !< End of j index range.

! Local variables
  real    :: curv, dh, scale
  integer :: i,j

	do j=jis,jie 
  do i=iis,iie
    ! This limiter prevents undershooting minima within the domain with
    ! values less than h_min.
    curv = 3.0*(h_l(i,j) + h_r(i,j) - 2.0*h_in(i,j))
   	if (curv > 0.0) then ! Only minima are limited.
      dh = h_r(i,j) - h_l(i,j)
      if (abs(dh) < curv) then ! The parabola's minimum is within the cell.
       if (h_in(i,j) <= h_min) then
         h_l(i,j) = h_in(i,j) ; h_r(i,j) = h_in(i,j)
       elseif (12.0*curv*(h_in(i,j) - h_min) < (curv**2 + 3.0*dh**2)) then
         ! The minimum value is h_in - (curv^2 + 3*dh^2)/(12*curv), and must
         ! be limited in this case.  0 < scale < 1.
         scale = 12.0*curv*(h_in(i,j) - h_min) / (curv**2 + 3.0*dh**2)
         h_l(i,j) = h_in(i,j) + scale*(h_l(i,j) - h_in(i,j))
         h_r(i,j) = h_in(i,j) + scale*(h_r(i,j) - h_in(i,j))
       endif
     endif
   	endif
 	enddo 
	enddo

end subroutine ppm_limit_pos

!=====================================================================
!--- relaxation forcing on the upper-layer thickness
!=====================================================================
subroutine h_relax_user_tilt(h,LB,G,ni,nj,dt)

  type(ocean_grid_type),                     intent(inout) :: G   !< The ocean's grid structure
  real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke),  intent(inout)  :: h
  type(loop_bounds_type),              intent(in)     :: LB
  integer,                             intent(in)     :: nj,ni
  real, dimension(g%isd:g%ied, g%jsd:g%jed)           :: h_relax
  real, dimension(g%isd:g%ied, g%jsd:g%jed)           :: h1_restore    !! USER-added variable, used to store the profile towards which the upper layer thickness is restored. 
  real,                             intent(in)     :: dt ! time step [s] for correcting relax rate

  character(len=40)  :: mod = "GOLD_continuity_PPM"

  real :: relax_const  !! The rate of relaxation.
  real :: dh      !! 1/2 the SSH drop

  integer :: i, j, je, js, ie, is, ig_50, ih_50
  integer :: j_12, j_25, j_37, j_50
  integer :: j_62, j_75, j_87
  integer :: j_shift, jg_50

  real :: x, y, xy
  real :: PI

  ! These should match with values in user_surface_forcing.F90 if using both relaxation and USER wind.
  real :: m, yf, y_shift
  ! one wave length (-pi to pi) of the tilted sinusoidal reference profile
  real :: L_pro

  real :: a, b, c

  !==

  a = 162.0
  b = -80.0 
  c = -80.0

  PI = 4.0*atan(1.0)

  js = LB%jsh ; je = LB%jeh
  is = LB%ish ; ie = LB%ieh

  !call get_param(param_file, mod, "relax_const", relax_const, &
  !               "Rate of relaxation", default=5e-5)
  !call get_param(param_file, mod, "delta_h", delta_h, &
  !               "1/2 size of SSH sine curve", default=80.0)
  relax_const = 2.0e-8 * dt !1.0e-5 or 0.7e-5 for dt=300s  
  ! This value is dimensionless, actual rate of relaxation is relax_const/dt in [1/s]. 1.0 for 256, 0.5 for 512, 2.0 for 128
  ! For N = 513, use relaxation rate 0.5e-5.
  dh = 150.0 ![m]
  ! For h1 = 500.0, suitable dh = 120.
  ! For h1 = 300.0, suitable dh = ???.
  
  ! Shift is caused by the asymmetry term. With asym = 1.0: use 28 (N=1024); use 14 (N=512); use 2 (N=50).
  j_shift = 0 !13 !24
  y_shift = G%geoLatT(is,j_shift+1) - G%geoLatT(is,1)
  !print *, y_shift ! has been tested

  ! related to the geometry of USER-defined wind
  m = 0.1
  yf = G%len_lat/2.0 !G%geoLatT(is,j_50);
  ! relaxation all over the domain
  L_pro = G%len_lat/1.0 ! (G%geoLatT(is,j_62) - G%geoLatT(is,j_37))

  h_relax(:,:) = 0.0
  do j=js,je 
  	do i=is,ie
    x = G%geoLonT(i,j) - G%west_lon
    y = G%geoLatT(i,j) - G%south_lat
    if (m * x + yf - 0.6 * L_pro < y - y_shift .AND. y - y_shift < m * x + yf + 0.6 * L_pro) then
      xy = 2.0 * PI * (m * x - y + y_shift + yf) / L_pro 
      !SIN(xy) = xy - xy**3 / 6.0 + xy**5 / 120.0 - xy**7 / 5040.0 + xy**9 / 362880.0 - xy**11 / 39916800.0 + xy**13 / 6227020800.0 - xy**15 / 1307674368000 + xy**17 / 3.55687428e14 - xy**19 / 1.216451e17
      h1_restore(i,j) = 300.0 + dh * SIN(xy)
      ! Apply the relaxation in the jet region.
      h_relax(i,j) = relax_const * (h1_restore(i,j) - h(i,j,1))
      h(i,j,1) = h(i,j,1) + h_relax(i,j)

    endif 
  	enddo
	enddo

end subroutine h_relax_user_tilt

!=====================================================================
! 		initialize CS of continuity_PPM
!=====================================================================
!> Initializes continuity_ppm_cs
subroutine continuity_ppm_init(CS)
  type(continuity_ppm_cs), pointer    :: cs   !< Module's control structure.

  allocate(cs)
  ! see MOM_parameter_doc.all, some of which are different for numerical reason
  cs%monotonic = .false.
  cs%simple_2nd = .false.
  cs%upwind_1st = .false.
  cs%tol_eta = 1.0E-12
  cs%tol_vel = 3.0E+08
  cs%tol_eta_aux = 1.0E-12
  cs%cfl_limit_adjust = 0.5
  cs%aggress_adjust = .false.
  cs%vol_cfl = .false.
  ! the following 3 params in the online model are true, for offline consider change?
  cs%better_iter = .false. 
  cs%use_visc_rem_max = .false.
  cs%marginal_faces = .false.

end subroutine continuity_ppm_init

end module mod_continuity
