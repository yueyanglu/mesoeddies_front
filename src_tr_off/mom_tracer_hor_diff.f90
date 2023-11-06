module mom_tracer_hor_diff

use mom_tracer_registry,      only : tracer_registry_type, tracer_type
use mod_grid_my, only : ocean_grid_type

implicit none ; private

public tracer_hor_diff_init, tracer_hordiff

!===========================================================
! 
!===========================================================
!> The ocntrol structure for along-layer and epineutral tracer diffusion
type, public :: tracer_hor_diff_cs ; private
  real    :: khtr           !< The along-isopycnal tracer diffusivity [L2 T-1 ~> m2 s-1].
  real    :: khtr_slope_cff !< The non-dimensional coefficient in KhTr formula [nondim]
  real    :: khtr_min       !< Minimum along-isopycnal tracer diffusivity [L2 T-1 ~> m2 s-1].
  real    :: khtr_max       !< Maximum along-isopycnal tracer diffusivity [L2 T-1 ~> m2 s-1].
  real    :: khtr_passivity_coeff !< Passivity coefficient that scales Rd/dx (default = 0)
                                 !! where passivity is the ratio between along-isopycnal
                                 !! tracer mixing and thickness mixing [nondim]
  real    :: khtr_passivity_min   !< Passivity minimum (default = 1/2) [nondim]
  real    :: ml_khtr_scale        !< With Diffuse_ML_interior, the ratio of the
                                 !! truly horizontal diffusivity in the mixed
                                 !! layer to the epipycnal diffusivity [nondim].
  real    :: max_diff_cfl         !< If positive, locally limit the along-isopycnal
                                 !! tracer diffusivity to keep the diffusive CFL
                                 !! locally at or below this value [nondim].
  logical :: diffuse_ml_interior  !< If true, diffuse along isopycnals between
                                 !! the mixed layer and the interior.
  logical :: check_diffusive_cfl  !< If true, automatically iterate the diffusion
                                 !! to ensure that the diffusive equivalent of
                                 !! the CFL limit is not violated.
  logical :: use_neutral_diffusion !< If true, use the neutral_diffusion module from within
                                  !! tracer_hor_diff.
  logical :: use_lateral_boundary_diffusion !< If true, use the lateral_boundary_diffusion module from within
                                        !! tracer_hor_diff.
  logical :: recalc_neutral_surf   !< If true, recalculate the neutral surfaces if CFL has been
                                  !! exceeded
  logical :: debug                 !< If true, write verbose checksums for debugging purposes.
  logical :: show_call_tree        !< Display the call tree while running. Set by VERBOSITY level.
  logical :: first_call = .true.   !< This is true until after the first call
  !>@{ Diagnostic IDs
  integer :: id_khtr_u  = -1
  integer :: id_khtr_v  = -1
  integer :: id_khtr_h  = -1
  integer :: id_cfl     = -1
  integer :: id_khdt_x  = -1
  integer :: id_khdt_y  = -1
  !!@}

end type tracer_hor_diff_cs

contains

!===========================================================
! 
!===========================================================
!> Compute along-coordinate diffusion of all tracers
!! using the diffusivity in CS%KhTr, or using space-dependent diffusivity.
!! Multiple iterations are used (if necessary) so that there is no limit
!! on the acceptable time increment.
subroutine tracer_hordiff(h, dt, G, CS, Reg, do_online_flag, read_khdt_x, read_khdt_y)
  type(ocean_grid_type),      intent(in) :: g       !< Grid type
  real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                             intent(in)    :: h       !< Layer thickness [H ~> m or kg m-2]
  real,                       intent(in)    :: dt      !< time step [T ~> s]
  type(tracer_hor_diff_cs),   pointer       :: cs      !< module control structure
  type(tracer_registry_type), intent(inout)       :: reg     !< registered tracers

  ! Optional inputs for offline tracer transport
  logical,          optional, intent(in)    :: do_online_flag !< If present and true, do online
                                                      !! tracer transport with stored velcities.
  real, dimension(g%isdb:g%iedb, g%jsd:g%jed), &
                   optional, intent(in)    :: read_khdt_x !< If present, these are the zonal
                                                      !! diffusivities from previous run.
  real, dimension(g%isd:g%ied,   g%jsdb:g%jedb), &
                   optional, intent(in)    :: read_khdt_y !< If present, these are the meridional
                                                      !! diffusivities from previous run.

  real, dimension(g%isd:g%ied, g%jsd:g%jed) :: &
   ihdxdy, &     ! The inverse of the volume or mass of fluid in a layer in a
                 ! grid cell [H-1 L-2 ~> m-3 or kg-1].
   kh_h, &       ! The tracer diffusivity averaged to tracer points [L2 T-1 ~> m2 s-1].
   cfl, &        ! A diffusive CFL number for each cell [nondim].
   dtr           ! The change in a tracer's concentration, in units of concentration [Conc].

  real, dimension(g%isdb:g%iedb, g%jsd:g%jed) :: &
   khdt_x, &     ! The value of Khtr*dt times the open face width divided by
                 ! the distance between adjacent tracer points [L2 ~> m2].
   coef_x, &     ! The coefficients relating zonal tracer differences
                 ! to time-integrated fluxes [H L2 ~> m3 or kg].
   kh_u          ! Tracer mixing coefficient at u-points [L2 T-1 ~> m2 s-1].
  real, dimension(g%isd:g%ied,   g%jsdb:g%jedb) :: &
   khdt_y, &     ! The value of Khtr*dt times the open face width divided by
                 ! the distance between adjacent tracer points [L2 ~> m2].
   coef_y, &     ! The coefficients relating meridional tracer differences
                 ! to time-integrated fluxes [H L2 ~> m3 or kg].
   kh_v          ! Tracer mixing coefficient at u-points [L2 T-1 ~> m2 s-1].

  real :: khdt_max ! The local limiting value of khdt_x or khdt_y [L2 ~> m2].
  real :: max_cfl ! The global maximum of the diffusive CFL number.
  logical :: use_varmix, resoln_scaled, do_online, use_eady
  integer :: s_idx, t_idx ! Indices for temperature and salinity if needed
  integer :: i, j, k, m, is, ie, js, je, nz, ntr, itt, num_itts
  real :: i_numitts  ! The inverse of the number of iterations, num_itts.
  real :: scale      ! The fraction of khdt_x or khdt_y that is applied in this
                    ! layer for this iteration [nondim].
  real :: idt        ! The inverse of the time step [T-1 ~> s-1].
  real :: h_neglect  ! A thickness that is so small it is usually lost
                    ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: kh_loc     ! The local value of Kh [L2 T-1 ~> m2 s-1].
  real :: res_fn     ! The local value of the resolution function [nondim].
  real :: rd_dx      ! The local value of deformation radius over grid-spacing [nondim].
  real :: normalize  ! normalization used for diagnostic Kh_h; diffusivity averaged to h-points.
  real, parameter :: H_subroundoff = 1.0e-20 * max(1.0e-10, 1.0e-17)

  is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec ; nz = g%ke
  do_online = .true.
  if (present(do_online_flag)) do_online = do_online_flag

  if (.not. associated(cs))  print*, "ERROR: MOM_tracer_hor_diff CS !!!"

  ntr = reg%ntr
  idt = 1.0 / dt
  h_neglect = H_subroundoff

  use_varmix = .false. ; resoln_scaled = .false. ; use_eady = .false.

  ! print*, 'cs%khtr = ', cs%khtr

  if (do_online) then

    !--------------- calc the diffusivities --------------
    if (use_varmix) then
      print*, "ERROR: NOT finished yet: use_varmix"
    elseif (resoln_scaled) then
      print*, "ERROR: NOT finished yet: resoln_scaled"
    else  ! Use a simple constant diffusivity.
      ! khdt_x
      if (cs%id_KhTr_u > 0) then
        !$OMP parallel do default(shared)
        do j=js,je
          do i=is-1,ie
            kh_u(i,j) = cs%KhTr
            khdt_x(i,j) = dt*(cs%KhTr*(g%dy_Cu(i,j)*g%IdxCu(i,j)))
          enddo
        enddo
        !$OMP END PARALLEL DO
      else
        !$OMP parallel do default(shared)
        do j=js,je
          do i=is-1,ie
            khdt_x(i,j) = dt*(cs%KhTr*(g%dy_Cu(i,j)*g%IdxCu(i,j)))
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif
      ! khdt_y
      if (cs%id_KhTr_v > 0) then
        !$OMP parallel do default(shared)
        do j=js-1,je
          do i=is,ie
            kh_v(i,j) = cs%KhTr
            khdt_y(i,j) = dt*(cs%KhTr*(g%dx_Cv(i,j)*g%IdyCv(i,j)))
          enddo
        enddo
        !$OMP END PARALLEL DO
      else
        !$OMP parallel do default(shared)
        do j=js-1,je
          do i=is,ie
            khdt_y(i,j) = dt*(cs%KhTr*(g%dx_Cv(i,j)*g%IdyCv(i,j)))
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif 

    endif ! VarMix
  
    if (cs%max_diff_CFL > 0.0) then

      if ((cs%id_KhTr_u > 0) .or. (cs%id_KhTr_h > 0)) then
        !$OMP parallel do default(shared) private(khdt_max)
        do j=js,je
          do i=is-1,ie
            khdt_max = 0.125*cs%max_diff_CFL * min(g%areaT(i,j), g%areaT(i+1,j))
            if (khdt_x(i,j) > khdt_max) then
              khdt_x(i,j) = khdt_max
              if (dt*(g%dy_Cu(i,j)*g%IdxCu(i,j)) > 0.0) &
                kh_u(i,j) = khdt_x(i,j) / (dt*(g%dy_Cu(i,j)*g%IdxCu(i,j)))
            endif
          enddo
        enddo
        !$OMP END PARALLEL DO
      else
        !$OMP parallel do default(shared) private(khdt_max)
        do j=js,je
          do i=is-1,ie
            khdt_max = 0.125*cs%max_diff_CFL * min(g%areaT(i,j), g%areaT(i+1,j))
            khdt_x(i,j) = min(khdt_x(i,j), khdt_max)
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif
      if ((cs%id_KhTr_v > 0) .or. (cs%id_KhTr_h > 0)) then
        !$OMP parallel do default(shared) private(khdt_max)
        do j=js-1,je
          do i=is,ie
            khdt_max = 0.125*cs%max_diff_CFL * min(g%areaT(i,j), g%areaT(i,j+1))
            if (khdt_y(i,j) > khdt_max) then
              khdt_y(i,j) = khdt_max
              if (dt*(g%dx_Cv(i,j)*g%IdyCv(i,j)) > 0.0) &
              kh_v(i,j) = khdt_y(i,j) / (dt*(g%dx_Cv(i,j)*g%IdyCv(i,j)))
            endif
          enddo
        enddo
        !$OMP END PARALLEL DO
      else
        !$OMP parallel do default(shared) private(khdt_max)
        do j=js-1,je
          do i=is,ie
            khdt_max = 0.125*cs%max_diff_CFL * min(g%areaT(i,j), g%areaT(i,j+1))
            khdt_y(i,j) = min(khdt_y(i,j), khdt_max)
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif

    endif ! cs%max_diff_CFL > 0.0
  
  else ! .not. do_online
    !$OMP parallel do default(shared)
    do j=js,je
      do i=is-1,ie
        khdt_x(i,j) = read_khdt_x(i,j)
      enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP parallel do default(shared)
    do j=js-1,je
      do i=is,ie
        khdt_y(i,j) = read_khdt_y(i,j)
      enddo
    enddo
    !$OMP END PARALLEL DO
  endif ! do_online
  

  if (cs%check_diffusive_CFL) then ! False
  elseif (cs%max_diff_CFL > 0.0) then ! -1.0
    num_itts = max(1, ceiling(cs%max_diff_CFL - 4.0*epsilon(cs%max_diff_CFL)))
    i_numitts = 1.0 / (real(num_itts))
  else
    num_itts = 1 ; i_numitts = 1.0
  endif
  
  
  if (cs%use_neutral_diffusion) then 
  else    ! following if not using neutral diffusion, but instead along-surface diffusion
    do itt=1,num_itts
      !$OMP parallel do default(shared) private(scale,Coef_y,Coef_x,Ihdxdy,dTr)
      do k=1,nz
        scale = i_numitts
        do j=js-1,je
          do i=is,ie
            coef_y(i,j) = ((scale * khdt_y(i,j))*2.0*(h(i,j,k)*h(i,j+1,k))) / &
                                                    (h(i,j,k)+h(i,j+1,k)+h_neglect)
          enddo
        enddo
  
        do j=js,je
          do i=is-1,ie
            coef_x(i,j) = ((scale * khdt_x(i,j))*2.0*(h(i,j,k)*h(i+1,j,k))) / &
                                                      (h(i,j,k)+h(i+1,j,k)+h_neglect)
          enddo
          do i=is,ie
            ihdxdy(i,j) = g%IareaT(i,j) / (h(i,j,k)+h_neglect)
          enddo
        enddo

        ! ------------ update tracers ------------
        do m=1,ntr
          do j=js,je
            do i=is,ie
              dtr(i,j) = ihdxdy(i,j) * &
                ((coef_x(i-1,j) * (reg%Tr(m)%t(i-1,j,k) - reg%Tr(m)%t(i,j,k)) - &
                 coef_x(i,j) * (reg%Tr(m)%t(i,j,k) - reg%Tr(m)%t(i+1,j,k))) + &
                (coef_y(i,j-1) * (reg%Tr(m)%t(i,j-1,k) - reg%Tr(m)%t(i,j,k)) - &
                 coef_y(i,j) * (reg%Tr(m)%t(i,j,k) - reg%Tr(m)%t(i,j+1,k))))
            enddo
          enddo

          do j=js,je
            do i=is,ie
              reg%Tr(m)%t(i,j,k) = reg%Tr(m)%t(i,j,k) + dtr(i,j)
            enddo
          enddo
        enddo ! m
  
      enddo ! End of k loop.
      !$OMP END PARALLEL DO
    enddo ! End of num_itts
  
  endif   ! endif for CS%use_neutral_diffusion  

end subroutine tracer_hordiff

!===========================================================
!   tracer_hor_diff_init
!===========================================================
!> Initialize lateral tracer diffusion module
subroutine tracer_hor_diff_init(G, CS, khtr)
  type(ocean_grid_type),      intent(in)    :: g          !< ocean grid structure
  type(tracer_hor_diff_cs),   pointer       :: cs         !< horz diffusion control structure
  real,      intent(in)    :: khtr

  if (associated(cs)) then
    print*, 'ERROR: tracer_hor_diff_init called with associated control structure'
   return
  endif
  allocate(cs)

  cs%khtr = khtr
  cs%use_neutral_diffusion = .false.
  cs%check_diffusive_CFL = .false.
  cs%max_diff_CFL = -1.0
  cs%id_KhTr_u = -1
  cs%id_KhTr_v = -1
  cs%id_KhTr_h = -1
  cs%id_CFL    = -1
end subroutine tracer_hor_diff_init



end module mom_tracer_hor_diff

