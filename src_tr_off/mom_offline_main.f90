module mom_offline_main

use mom_tracer_registry,      only : tracer_registry_type, tracer_type
use mod_grid_my, only : ocean_grid_type

implicit none ; private


public offline_advection_layer

contains

!===========================================================
!  offline_advection_layer
!===========================================================
!> When in layer mode, 3D horizontal advection using stored mass fluxes must be used. Horizontal advection is
!! done via tracer_advect, whereas the vertical component is actually handled by vertdiff in tracer_column_fns
subroutine offline_advection_layer(time_interval, tracer_Reg, &
    G, h_pre, h_end, uhtr, vhtr, alpha_input)

real,                  intent(in)        :: time_interval !< Offline transport time interval
type(ocean_grid_type), intent(in)  :: g  
type(tracer_registry_type), intent(inout)          :: tracer_Reg            !< Control structure for offline module
real, dimension(g%isd:g%ied, g%jsd:g%jed,   g%ke), intent(inout) :: h_pre !< h at t1 before adv (from onl model)
real, dimension(g%isd:g%ied, g%jsd:g%jed,   g%ke), intent(in)    :: h_end !< h at t2 after adv (added by YUEYANG)
real, dimension(g%isdb:g%iedb, g%jsd:g%jed, g%ke), intent(inout) :: uhtr  !< Zonal mass transport
real, dimension(g%isd:g%ied, g%jsdb:g%jedb, g%ke), intent(inout) :: vhtr  !< Meridional mass transport
real, dimension(g%isd:g%ied, g%jsd:g%jed,   g%ke), optional, intent(in) &
                                :: alpha_input ! for (1-alpha)* div{Uc} 

! Local pointers
! Remaining zonal mass transports
real, dimension(g%isdb:g%iedb, g%jsd:g%jed, g%ke)   :: uhtr_sub
! Remaining meridional mass transports
real, dimension(g%isd:g%ied, g%jsdb:g%jedb, g%ke)   :: vhtr_sub
!
real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke)   :: alpha_inside ! 

real :: sum_abs_fluxes, sum_u, sum_v  ! Used to keep track of how close to convergence we are
real :: dt_offline

! Variables used to keep track of layer thicknesses at various points in the code
real, dimension(g%isd:g%ied, g%jsd:g%jed,   g%ke) :: &
       h_new, &
       h_vol, &
       h_end_vol, h_aft_iter

integer :: niter, iter, num_off_iter
real    :: inum_iter
real    :: dt_iter  ! The timestep of each iteration [T ~> s]
logical :: converged
character(len=160) :: mesg  ! The text of an error message
integer :: i, j, k, m, is, ie, js, je, isd, ied, jsd, jed, nz
integer :: isv, iev, jsv, jev ! The valid range of the indices.
integer :: isdb, iedb, jsdb, jedb
logical :: z_first, x_before_y

if (present(alpha_input)) then
   alpha_inside = alpha_input
 else
   alpha_inside = 0.
endif

is  = g%isc ; ie  = g%iec ; js  = g%jsc ; je  = g%jec ; nz = g%ke
isd = g%isd ; ied = g%ied ; jsd = g%jsd ; jed = g%jed
isdb = g%IsdB ; iedb = g%IedB ; jsdb = g%JsdB ; jedb = g%JedB
  
x_before_y = .true.
num_off_iter = 60 ! offline_transport_cs%num_off_iter

dt_iter = time_interval / real(max(1, num_off_iter))
print*,'offline_advection_layer: iteration dt_iter [s] = ',dt_iter

! --------- begin iteration
do iter=1,num_off_iter

  do k = 1, nz
    do j=js-1,je+1
      do i=is-2,ie+1
        uhtr_sub(i,j,k) = uhtr(i,j,k)
      enddo
    enddo
  enddo

  do k = 1, nz
    do j=js-2,je+1
      do i=is-1,ie+1
        vhtr_sub(i,j,k) = vhtr(i,j,k)
      enddo
    enddo
  enddo

  ! Calculate 3d mass transports to be used in this iteration
  ! "h_pre" unchaged! 
  call limit_mass_flux_3d(g, uhtr_sub, vhtr_sub, h_pre)

  ! print*, 'h_pre(22,22,:) = ', h_pre(22,22,:)

  ! --------- vertical advection ---------
  ! --- call update_h_vertical_flux(g, gv, eatr_sub, ebtr_sub, h_pre, h_new)
  ! --- call call_tracer_column_fns(h_pre, h_new, eatr_sub, ebtr_sub,...)
  h_new = h_pre ! h_new after vertical adv (= h_pre because vertical flux is zero!)
  ! We are now done with the vertical mass transports, so now h_new is h_sub
  do k = 1, nz
    do j=js-1,je+1
      do i=is-1,ie+1
        h_pre(i,j,k) = h_new(i,j,k)
      enddo
    enddo
  enddo

  ! print*, 'AFTER vert ADV h_pre(22,22,:) = ', h_pre(22,22,:)

  ! --------- horizontal advection ---------
  ! "h_pre" not changed; "h_new" is layer thick after mass flux adv (output) !
  call update_h_horizontal_flux(g, uhtr_sub, vhtr_sub, h_pre, h_new)
  do k = 1, nz
    do i=is-1,ie+1
      do j=js-1,je+1
        h_vol(i,j,k) = h_pre(i,j,k)*g%areaT(i,j)
      enddo
    enddo
  enddo

  ! advection occurs here, iteratively !!!
  ! "h_end" not changed, h_vol is the updated one, ideally h_vol = h_end but not 
  call advect_tracer(h_end, uhtr_sub, vhtr_sub, dt_iter, g, tracer_Reg, &
    x_first_in=x_before_y, vol_prev=h_vol, max_iter_in=30, update_vol_prev=.true., &
    alpha = alpha_inside) ! alpha added by Yueyang

 ! Done with horizontal so now h_pre should be h_new
  do k = 1, nz
    do i=is-1,ie+1
      do j=js-1,je+1
        h_pre(i,j,k) = h_new(i,j,k)
      enddo
    enddo
  enddo

  ! print*, 'AFTER hori ADV h_pre(22,22,:) = ', h_pre(22,22,:)

  ! Update remaining transports
  do k = 1, nz
    do j=js-1,je+1
      do i=is-2,ie+1
        uhtr(i,j,k) = uhtr(i,j,k) - uhtr_sub(i,j,k)
      enddo
    enddo
  enddo
  do k = 1, nz
    do j=js-2,je+1
      do i=is-1,ie+1
        vhtr(i,j,k) = vhtr(i,j,k) - vhtr_sub(i,j,k)
      enddo
    enddo
  enddo

  ! Calculate how close we are to converging by summing the remaining fluxes at each point
  sum_abs_fluxes = 0.0
  sum_u = 0.0
  sum_v = 0.0
  do k = 1, nz
    do i=is,ie
      do j=js,je
        sum_u = sum_u + abs(uhtr(i-1,j,k))+abs(uhtr(i,j,k))
        sum_v = sum_v + abs(vhtr(i,j-1,k))+abs(vhtr(i,j,k))
        sum_abs_fluxes = sum_abs_fluxes + abs(uhtr(i-1,j,k)) + &
            abs(uhtr(i,j,k)) + abs(vhtr(i,j-1,k)) + abs(vhtr(i,j,k))
      enddo
    enddo
  enddo
  print*, "offline_advection_layer: Remaining u-flux, v-flux:", sum_u, sum_v
  if (sum_abs_fluxes==0) then
    print*, 'offline_advection_layer: Converged after iteration', iter
    
    !----- added by YUEYANG --- at the end of iteration, update tracer by h_t2 !
    ! c_new = c * h_vol / h_end
    do k = 1, nz
      do j = js,je
        do i = is,ie
          ! h_end_vol(i,j,k) = h_end(i,j,k)*g%areaT(i,j)
          h_aft_iter(i,j,k) = h_vol(i,j,k) * g%IareaT(i,j)
        enddo
      enddo
    enddo
    ! print*, 'h_vol(22,22,1)=', h_vol(22,22,1), 'h_end_vol(22,22,1)=', h_end_vol(22,22,1)
    ! call update_tracer_by_h(tracer_Reg, G, h_aft_iter, h_end)

    exit
  endif

  ! Switch order of Strang split every iteration
  x_before_y = .not. x_before_y

enddo ! iter

end subroutine offline_advection_layer

! ====================================================================
!           SUBTOUTINE: update tracer at each point by a h_t1 and h_t2
! =====================================================================
subroutine update_tracer_by_h(reg, G, h_pre, h_end)
  type(ocean_grid_type),    intent(in) :: g   
  real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                            intent(in) :: h_pre, h_end
  type(tracer_registry_type), intent(inout)  :: reg   !< pointer to tracer registry

  ! local
  integer :: i, j, k, m, ntr, is, ie, js, je, nz
  is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec ; nz = g%ke

  ntr = reg%ntr

  do m = 1,ntr
    do k = 1, nz
      do j = js,je
        do i = is,ie
          if (g%mask2dT(i,j) /= 0.) reg%Tr(m)%t(i,j,k) = & 
            reg%Tr(m)%t(i,j,k) * h_pre(i,j,k) / h_end(i,j,k)
        enddo
      enddo
    enddo
  enddo

end subroutine update_tracer_by_h

!===========================================================
!  limit_mass_flux_3d from "mom_offline_aux"
!===========================================================
!> This routine limits the mass fluxes so that the a layer cannot be completely depleted.
!! NOTE: Only used in non-ALE mode
subroutine limit_mass_flux_3d(G, uh, vh, h_pre)
  type(ocean_grid_type),   intent(in)       :: g     !< ocean grid structure
  real, dimension(g%isdb:g%iedb, g%jsd:g%jed,  g%ke), &
                          intent(inout) :: uh    !< Mass flux through zonal face [kg]
  real, dimension(g%isd:g%ied,   g%jsdb:g%jedb, g%ke), &
                          intent(inout) :: vh    !< Mass flux through meridional face [kg]
  real, dimension(g%isd:g%ied,   g%jsd:g%jed,  g%ke), &
                          intent(in)    :: h_pre !< Layer thicknesses at the end of the previous
                                                 !! step [kg m-2].

  ! Local variables
  integer :: i, j, k, m, is, ie, js, je, nz
  real :: pos_flux, hvol, h_neglect, scale_factor, max_off_cfl

  max_off_cfl = 0.5

  ! In this subroutine, fluxes out of the box are scaled away if they deplete
  ! the layer, note that we define the positive direction as flux out of the box.
  ! Hence, uh(I-1) is multipled by negative one, but uh(I) is not

  ! Set index-related variables for fields on T-grid
  is  = g%isc ; ie  = g%iec ; js  = g%jsc ; je  = g%jec ; nz = g%ke

  ! Calculate sum of positive fluxes (negatives applied to enforce convention)
  ! in a given cell and scale it back if it would deplete a layer
  do k = 1, nz
    do j=js-1,je+1
      do i=is-1,ie+1

        hvol = h_pre(i,j,k)*g%areaT(i,j)
        pos_flux  = max(0.0,-uh(i-1,j,k)) + max(0.0, -vh(i,j-1,k)) + &
            max(0.0, uh(i,j,k)) + max(0.0, vh(i,j,k))

        ! if (i==20 .and. j==20 .and. k==1)   print*, pos_flux, hvol

        if (pos_flux>hvol .and. pos_flux>0.0) then
          scale_factor = ( hvol )/pos_flux*max_off_cfl
        else ! Don't scale
          scale_factor = 1.0
        endif

      ! Scale horizontal fluxes
        if (-uh(i-1,j,k)>0) uh(i-1,j,k) = uh(i-1,j,k)*scale_factor
        if (uh(i,j,k)>0)    uh(i,j,k)   = uh(i,j,k)*scale_factor
        if (-vh(i,j-1,k)>0) vh(i,j-1,k) = vh(i,j-1,k)*scale_factor
        if (vh(i,j,k)>0)    vh(i,j,k)   = vh(i,j,k)*scale_factor
      enddo
    enddo
  enddo

end subroutine limit_mass_flux_3d


!===========================================================
!  update_h_horizontal_flux from "mom_offline_aux"
!===========================================================
!> This updates thickness based on the convergence of horizontal mass fluxes
!! NOTE: Only used in non-ALE mode
subroutine update_h_horizontal_flux(G, uhtr, vhtr, h_pre, h_new)
  type(ocean_grid_type),   intent(in)       :: g     
  real, dimension(g%isdb:g%iedb, g%jsd:g%jed,  g%ke), &
                          intent(in)    :: uhtr  !< Accumulated mass flux through zonal face [kg]
  real, dimension(g%isd:g%ied,   g%jsdb:g%jedb, g%ke), &
                          intent(in)    :: vhtr  !< Accumulated mass flux through meridional face [kg]
  real, dimension(g%isd:g%ied,   g%jsd:g%jed,  g%ke), &
                          intent(in)    :: h_pre !< Previous layer thicknesses [kg m-2].
  real, dimension(g%isd:g%ied,   g%jsd:g%jed,  g%ke), &
                          intent(inout) :: h_new !< Updated layer thicknesses [kg m-2].
  ! Local variables
  real, parameter :: Angstrom_H = 1.0E-10
  integer :: i, j, k, m, is, ie, js, je, nz
  ! Set index-related variables for fields on T-grid
  is  = g%isc ; ie  = g%iec ; js  = g%jsc ; je  = g%jec ; nz = g%ke

  do k = 1, nz
    do i=is-1,ie+1
      do j=js-1,je+1
        h_new(i,j,k) = max( 0.0, g%areaT(i,j)*h_pre(i,j,k) + &
         ((uhtr(i-1,j,k) - uhtr(i,j,k)) + (vhtr(i,j-1,k) - vhtr(i,j,k))) )

        ! In the case that the layer is now dramatically thinner than it was previously,
        ! add a bit of mass to avoid truncation errors.  This will lead to
        ! non-conservation of tracers
        h_new(i,j,k) = h_new(i,j,k) + &
         max(Angstrom_H, 1.0e-13*h_new(i,j,k) - g%areaT(i,j)*h_pre(i,j,k))

        ! Convert back to thickness
        h_new(i,j,k) = h_new(i,j,k) / g%areaT(i,j)
      enddo
    enddo 
  enddo
end subroutine update_h_horizontal_flux

!===========================================================
!  advect_tracer from "mom_tracer_advect"
!===========================================================
 !> This routine time steps the tracer concentration using a
 !! monotonic, conservative, weakly diffusive scheme.
subroutine advect_tracer(h_end, uhtr, vhtr, dt, G, Reg, x_first_in, &
      vol_prev, max_iter_in, update_vol_prev, uhr_out, vhr_out, alpha)
  type(ocean_grid_type),   intent(in) :: g     !< ocean grid structure
  real, dimension(g%isd:g%ied,   g%jsd:g%jed,  g%ke), &
                          intent(in)    :: h_end !< layer thickness after advection [H ~> m or kg m-2]
  real, dimension(g%isdb:g%iedb, g%jsd:g%jed,  g%ke), &
                          intent(in)    :: uhtr  !< accumulated volume/mass flux through zonal face [H L2 ~> m3 or kg]
  real, dimension(g%isd:g%ied,   g%jsdb:g%jedb, g%ke), &
                          intent(in)    :: vhtr  !< accumulated volume/mass flux through merid face [H L2 ~> m3 or kg]
  real,                    intent(in)    :: dt    !< time increment [T ~> s]
  type(tracer_registry_type), intent(inout)    :: reg   !< pointer to tracer registry
  logical,       optional, intent(in)    :: x_first_in !< If present, indicate whether to update
                                                 !! first in the x- or y-direction.

  ! The remaining optional arguments are only used in offline tracer mode.
  real, dimension(g%isd:g%ied,   g%jsd:g%jed,  g%ke), &
                 optional, intent(inout) :: vol_prev !< Cell volume before advection [H L2 ~> m3 or kg].
                                                  !! If update_vol_prev is true, the returned value is
                                                  !! the cell volume after the transport that was done
                                                  !! by this call, and if all the transport could be
                                                  !! accommodated it should be close to h_end*G%areaT.
  integer,       optional, intent(in)    :: max_iter_in !< The maximum number of iterations
  logical,       optional, intent(in)    :: update_vol_prev !< If present and true, update vol_prev to
                                                  !! return its value after the tracer have been updated.
  real, dimension(g%isdb:g%iedb, g%jsd:g%jed,  g%ke), &
                optional, intent(out)    :: uhr_out  !< accumulated volume/mass flux through zonal face
                                                 !! [H L2 ~> m3 or kg]
  real, dimension(g%isd:g%ied,   g%jsdb:g%jedb, g%ke), &
                optional, intent(out)    :: vhr_out  !< accumulated volume/mass flux through merid face
                                                 !! [H L2 ~> m3 or kg]
  real, dimension(g%isd:g%ied,   g%jsd:g%jed,  g%ke), &
                           intent(in)    :: alpha 
  ! local vars
  type(tracer_type) :: tr(10) ! The array of registered tracers
  real, dimension(g%isd:g%ied,   g%jsd:g%jed,  g%ke) ::  &
            hprev           ! cell volume at the end of previous tracer change [H L2 ~> m3 or kg]
  real, dimension(g%isdb:g%iedb, g%jsd:g%jed,  g%ke) ::  &
            uhr             ! The remaining zonal thickness flux [H L2 ~> m3 or kg]
  real, dimension(g%isd:g%ied,   g%jsdb:g%jedb, g%ke)::  &
            vhr             ! The remaining meridional thickness fluxes [H L2 ~> m3 or kg]
  real :: uh_neglect(g%isdb:g%iedb, g%jsd:g%jed) ! uh_neglect and vh_neglect are the
  real :: vh_neglect(g%isd:g%ied,   g%jsdb:g%jedb) ! magnitude of remaining transports that
                                      ! can be simply discarded [H L2 ~> m3 or kg].

  real :: landvolfill                   ! An arbitrary? nonzero cell volume [H L2 ~> m3 or kg].
  real :: idt                           ! 1/dt [T-1 ~> s-1].
  logical :: domore_u(g%jsd:g%jed,g%ke)  ! domore__ indicate whether there is more
  logical :: domore_v(g%jsdb:g%jedb,g%ke) ! advection to be done in the corresponding
                                       ! row or column.
  logical :: x_first            ! If true, advect in the x-direction first.
  integer :: max_iter           ! maximum number of iterations in each layer
  integer :: domore_k(g%ke)
  integer :: stencil            ! stencil of the advection scheme
  integer :: nsten_halo         ! number of stencils that fit in the halos
  integer :: i, j, k, m, is, ie, js, je, isd, ied, jsd, jed, nz, itt, ntr, do_any
  integer :: isv, iev, jsv, jev ! The valid range of the indices.
  integer :: isdb, iedb, jsdb, jedb
  ! the following parameters were in "tracer_advect_cs"
  integer :: first_direction = 0 ! FIRST_DIRECTION
  logical :: usePPM = .false.,  useHuynh = .false.
  real, parameter :: H_subroundoff = 1.0e-20 * max(1.0e-10, 1.0e-17)
  real, parameter :: cs_dt = 300. ! dt from "tracer_advect_cs"

  domore_u(:,:) = .false.
  domore_v(:,:) = .false.
  is  = g%isc ; ie  = g%iec ; js  = g%jsc ; je  = g%jec ; nz = g%ke
  isd = g%isd ; ied = g%ied ; jsd = g%jsd ; jed = g%jed
  isdb = g%IsdB ; iedb = g%IedB ; jsdb = g%JsdB ; jedb = g%JedB
  landvolfill = 1.0e-20         ! This is arbitrary, but must be positive.
  stencil = 2                   ! The scheme's stencil; 2 for PLM and PPM:H3

  if (reg%ntr==0) return

  x_first = (mod(first_direction,2) == 0)

  ! increase stencil size for Colella & Woodward PPM
  if (usePPM .and. .not. useHuynh) stencil = 3

  ntr = reg%ntr
  do m=1,ntr
    tr(m) = reg%Tr(m)
  enddo
  idt = 1.0 / dt

  max_iter = 2*int(ceiling(dt/cs_dt)) + 1
  print*, '(advect_tracer) dt, cs_dt, max_iter = ', dt, cs_dt, max_iter

  if (present(max_iter_in)) max_iter = max_iter_in
  if (present(x_first_in))  x_first = x_first_in

  !$OMP parallel default(shared)

  ! This initializes the halos of uhr and vhr because pass_vector might do
  ! calculations on them, even though they are never used.
  !$OMP do
  do k=1,nz
    ! --- initialize
    do j=jsd,jed
      do i=isdb,iedb
        uhr(i,j,k) = 0.0
      enddo
    enddo
    do j=jsdb,jedb
      do i=isd,ied
        vhr(i,j,k) = 0.0
      enddo
    enddo
    do j=jsd,jed
      do i=isd,ied
        hprev(i,j,k) = 0.0
      enddo
    enddo
    domore_k(k)=1
    !  Put the remaining (total) thickness fluxes into uhr and vhr.
    do j=js,je
      do i=is-1,ie
        uhr(i,j,k) = uhtr(i,j,k)
      enddo
    enddo
    do j=js-1,je
      do i=is,ie
        vhr(i,j,k) = vhtr(i,j,k)
      enddo
    enddo

    if (.not. present(vol_prev)) then
   ! This loop reconstructs the thickness field the last time that the
   ! tracers were updated, probably just after the diabatic forcing.  A useful
   ! diagnostic could be to compare this reconstruction with that older value.
      do i=is,ie
        do j=js,je
          hprev(i,j,k) = max(0.0, g%areaT(i,j)*h_end(i,j,k) + &
            ((uhr(i,j,k) - uhr(i-1,j,k)) + (vhr(i,j,k) - vhr(i,j-1,k))))
   ! In the case that the layer is now dramatically thinner than it was previously,
   ! add a bit of mass to avoid truncation errors.  This will lead to
   ! non-conservation of tracers
          hprev(i,j,k) = hprev(i,j,k) + &
                      max(0.0, 1.0e-13*hprev(i,j,k) - g%areaT(i,j)*h_end(i,j,k))
        enddo
      enddo
    else
      do i=is,ie
        do j=js,je
          hprev(i,j,k) = vol_prev(i,j,k) ! [m3]
        enddo
      enddo
    endif
  enddo ! k


  !$OMP do
  do j=jsd,jed
    do i=isd,ied-1
      uh_neglect(i,j) = H_subroundoff*min(g%areaT(i,j),g%areaT(i+1,j))
    enddo
  enddo

  !$OMP do
  do j=jsd,jed-1
    do i=isd,ied
      vh_neglect(i,j) = H_subroundoff*min(g%areaT(i,j),g%areaT(i,j+1))
    enddo
  enddo
  !$OMP end parallel

  isv = is ; iev = ie ; jsv = js ; jev = je

  print*, '(advect_tracer) max_iter=', max_iter, 'x_first=', x_first

  do itt=1,max_iter
  
    if (isv > is-stencil) then
  
      nsten_halo = min(is-isd,ied-ie,js-jsd,jed-je)/stencil
      isv = is-nsten_halo*stencil ; jsv = js-nsten_halo*stencil
      iev = ie+nsten_halo*stencil ; jev = je+nsten_halo*stencil

      ! Reevaluate domore_u & domore_v unless the valid range is the same size as
      ! before.  Also, do this if there is Strang splitting.
      if ((nsten_halo > 1) .or. (itt==1)) then
        !$OMP parallel do default(shared)
        do k=1,nz
          if (domore_k(k) > 0) then

            ! do u & v
            do j=jsv,jev
              if (.not.domore_u(j,k)) then
                do i=isv+stencil-1,iev-stencil
                  if (uhr(i,j,k) /= 0.0) then
                    domore_u(j,k) = .true.
                    exit
                  endif
                enddo ! i-loop
              endif
            enddo ! j-loop
            do j=jsv+stencil-1,jev-stencil
              if (.not.domore_v(j,k)) then
                do i=isv+stencil,iev-stencil
                  if (vhr(i,j,k) /= 0.0) then
                    domore_v(j,k) = .true.
                    exit
                  endif
                enddo ! i-loop
              endif
            enddo ! j-loop
  
            !   At this point, domore_k is global.  Change it so that it indicates
            ! whether any work is needed on a layer on this processor.
            domore_k(k) = 0
            do j=jsv,jev ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
            do j=jsv+stencil-1,jev-stencil ; if (domore_v(j,k)) domore_k(k) = 1 ; enddo
  
          endif ! (domore_k(k) > 0)
        enddo ! k-loop
       endif ! ((nsten_halo > 1) .or. (itt==1))
    endif ! (isv > is-stencil)
  

    ! Set the range of valid points after this iteration.
    isv = isv + stencil ; iev = iev - stencil
    jsv = jsv + stencil ; jev = jev - stencil
  
    !$OMP parallel default(shared)
  
    !  To ensure positive definiteness of the thickness at each iteration, the
    !  mass fluxes out of each layer are checked each step, and limited to keep
    !  the thicknesses positive.  This means that several iterations may be required
    !  for all the transport to happen.  The sum over domore_k keeps the processors
    !  synchronized.  This may not be very efficient, but it should be reliable.

    ! print*, 'hprev(22,22,:)/areaT = ', hprev(22,22,:) / g%areaT(22,22)
    ! print*, 'tr(1)%t(22,22,:) = ', tr(1)%t(22,22,:)

    do k=1,nz
      if (domore_k(k) > 0) then
  
        if (x_first) then
  
         ! First, advect zonally. "hprev" will be updated (no iteration in advect_x/y) !!!
         call advect_x(tr, hprev, h_end, uhr, uh_neglect, domore_u, ntr, idt, &
                       isv, iev, jsv-stencil, jev+stencil, k, g, usePPM, useHuynh, &
                       alpha)
  
         !  Next, advect meridionally.
         call advect_y(tr, hprev, h_end, vhr, vh_neglect, domore_v, ntr, idt, &
                       isv, iev, jsv, jev, k, g, usePPM, useHuynh, &
                       alpha)
  
          domore_k(k) = 0
          do j=jsv-stencil,jev+stencil
            if (domore_u(j,k)) domore_k(k) = 1
          enddo
          do j=jsv-1,jev
            if (domore_v(j,k)) domore_k(k) = 1
          enddo
        else
  
         ! First, advect meridionally.
         call advect_y(tr, hprev, h_end, vhr, vh_neglect, domore_v, ntr, idt, &
                       isv-stencil, iev+stencil, jsv, jev, k, g, usePPM, useHuynh, &
                       alpha)

         ! Next, advect zonally.
         call advect_x(tr, hprev, h_end, uhr, uh_neglect, domore_u, ntr, idt, &
                       isv, iev, jsv, jev, k, g, usePPM, useHuynh, &
                       alpha)

          domore_k(k) = 0
          do j=jsv,jev
            if (domore_u(j,k)) domore_k(k) = 1
          enddo
          do j=jsv-1,jev
            if (domore_v(j,k)) domore_k(k) = 1
          enddo
  
        endif ! x_first
      endif ! (domore_k(k) > 0)
    enddo ! End of k-loop
    !$OMP END PARALLEL DO

    ! print*, 'hprev(22,22,:)/areaT = ', hprev(22,22,:) / g%areaT(22,22)
    ! print*, 'tr(1)%t(22,22,:) = ', tr(1)%t(22,22,:)
    ! 

    ! If the advection just isn't finishing after max_iter, move on.
    if (itt >= max_iter) then
      print*, '(advect_tracer) exit w/ itt=', itt
      exit
    endif

    ! Exit if there are no layers that need more iterations.
    if (isv > is-stencil) then
      do_any = 0
      ! call sum_across_pes(domore_k(:), nz)
      do k=1,nz
        do_any = do_any + domore_k(k)
      enddo
      if (do_any == 0) then
        print*, '(advect_tracer) exit w/ itt=', itt
        exit
      endif
    endif
  
  enddo ! Iterations loop
  
  if (present(uhr_out)) uhr_out(:,:,:) = uhr(:,:,:)
  if (present(vhr_out)) vhr_out(:,:,:) = vhr(:,:,:)
  if (present(vol_prev) .and. present(update_vol_prev)) then
    if (update_vol_prev) vol_prev(:,:,:) = hprev(:,:,:)
  endif
  
  !------ This is necessary, if 'reg' is not a pointer (by YY)------
  do m=1,ntr
    reg%Tr(m) = tr(m) 
  enddo

end subroutine advect_tracer

!===========================================================
!  advect_x from "mom_tracer_advect"
!===========================================================
!> This subroutine does 1-d flux-form advection in the zonal direction using
!! a monotonic piecewise linear scheme.
subroutine advect_x(Tr, hprev, h_end, uhr, uh_neglect, domore_u, ntr, Idt, &
                     is, ie, js, je, k, G, usePPM, useHuynh, alpha)
  type(ocean_grid_type),                     intent(in) :: G    !< The ocean's grid structure
  type(tracer_type), dimension(ntr),         intent(inout) :: Tr   !< The array of registered tracers to work on
  real, dimension(g%isd:g%ied,   g%jsd:g%jed,  g%ke),  intent(inout) :: hprev !< cell volume at the end of previous
                                                                 !! tracer change [H L2 ~> m3 or kg]
  real, dimension(g%isd:g%ied,   g%jsd:g%jed,  g%ke),  intent(in) :: h_end !< [m]

  real, dimension(g%isdb:g%iedb, g%jsd:g%jed,  g%ke), intent(inout) :: uhr !< accumulated volume/mass flux through
                                                                 !! the zonal face [H L2 ~> m3 or kg]
  real, dimension(g%isdb:g%iedb, g%jsd:g%jed),      intent(inout) :: uh_neglect !< A tiny zonal mass flux that can
                                                                 !! be neglected [H L2 ~> m3 or kg]
  logical, dimension(g%isdb:g%iedb,g%ke),       intent(inout) :: domore_u !< If true, there is more advection to be
                                                                 !! done in this u-row
  real,                                      intent(in)    :: Idt !< The inverse of dt [T-1 ~> s-1]
  integer,                                   intent(in)    :: ntr !< The number of tracers
  integer,                                   intent(in)    :: is  !< The starting tracer i-index to work on
  integer,                                   intent(in)    :: ie  !< The ending tracer i-index to work on
  integer,                                   intent(in)    :: js  !< The starting tracer j-index to work on
  integer,                                   intent(in)    :: je  !< The ending tracer j-index to work on
  integer,                                   intent(in)    :: k   !< The k-level to work on
  logical,                                   intent(in)    :: usePPM !< If true, use PPM instead of PLM
  logical,                                   intent(in)    :: useHuynh !< If true, use the Huynh scheme
                                                                    !! for PPM interface values
  real, dimension(g%isd:g%ied,   g%jsd:g%jed,  g%ke), &
                 intent(in)    :: alpha 

  real, dimension(g%isd:g%ied,ntr) :: &
   slope_x             ! The concentration slope per grid point [conc].
  real, dimension(g%isdb:g%iedb,ntr) :: &
   flux_x              ! The tracer flux across a boundary [H L2 conc ~> m3 conc or kg conc].
  real, dimension(g%isd:g%ied,ntr) :: &
   T_tmp               ! The copy of the tracer concentration at constant i,k [H m2 conc ~> m3 conc or kg conc].

  real :: maxslope      ! The maximum concentration slope per grid point
                       ! consistent with monotonicity [conc].
  real :: hup, hlos     ! hup is the upwind volume, hlos is the
                       ! part of that volume that might be lost
                       ! due to advection out the other side of
                       ! the grid box, both in [H L2 ~> m3 or kg].
  real :: uhh(g%isdb:g%iedb) ! The zonal flux that occurs during the
                       ! current iteration [H L2 ~> m3 or kg].
  real, dimension(g%isdb:g%iedb) :: &
   hlst, &             ! Work variable [H L2 ~> m3 or kg].
   Ihnew, &            ! Work variable [H-1 L-2 ~> m-3 or kg-1].
   CFL                 ! A nondimensional work variable [nondim].
  real :: min_h         ! The minimum thickness that can be realized during
                       ! any of the passes [H ~> m or kg m-2].
  real :: h_neglect     ! A thickness that is so small it is usually lost
                       ! in roundoff and can be neglected [H ~> m or kg m-2].
  logical :: do_i(g%isdb:g%iedb)     ! If true, work on given points.
  logical :: do_any_i
  integer :: i, j, m, n, i_up, stencil
  real :: aR, aL, dMx, dMn, Tp, Tc, Tm, dA, mA, a6
  real :: fac1,u_L_in,u_L_out  ! terms used for time-stepping OBC reservoirs
  logical :: usePLMslope, use_h_end
  real, parameter :: Angstrom_H = 1.0E-10, H_subroundoff = 1.0e-20 * max(Angstrom_H, 1.0e-17)
  
  useplmslope = .not. (useppm .and. usehuynh)

  use_h_end = .false.

  ! stencil for calculating slope values
  stencil = 1
  if (useppm .and. .not. usehuynh) stencil = 2

  min_h = 0.1*Angstrom_H
  h_neglect = H_subroundoff

  do i=is-1,ie
    cfl(i) = 0.0
  enddo
  
  ! --------- begin j-loop
  do j=js,je
    if (domore_u(j,k)) then
    domore_u(j,k) = .false.
  
    ! Calculate the i-direction profiles (slopes) of each tracer that
    ! is being advected.
    if (useplmslope) then
      do m=1,ntr
        do i=is-stencil,ie+stencil
          tp = tr(m)%t(i+1,j,k) ; tc = tr(m)%t(i,j,k) ; tm = tr(m)%t(i-1,j,k)
          dmx = max( tp, tc, tm ) - tc
          dmn= tc - min( tp, tc, tm )
          slope_x(i,m) = g%mask2dCu(i,j)*g%mask2dCu(i-1,j) * &
            sign( min(0.5*abs(tp-tm), 2.0*dmx, 2.0*dmn), tp-tm )
        enddo
      enddo
    endif ! usePLMslope
  
    ! make a copy of the tracers in case values need to be overridden for OBCs
    do m = 1,ntr
      do i=g%isd,g%ied
        t_tmp(i,m) = tr(m)%t(i,j,k)
      enddo
    enddo
    ! loop through open boundaries and recalculate flux terms
  
    ! Calculate the i-direction fluxes of each tracer, using as much
    ! the minimum of the remaining mass flux (uhr) and the half the mass
    ! in the cell plus whatever part of its half of the mass flux that
    ! the flux through the other side does not require.
    do i=is-1,ie
      if (uhr(i,j,k) == 0.0) then
        uhh(i) = 0.0
        cfl(i) = 0.0
      elseif (uhr(i,j,k) < 0.0) then
        hup = hprev(i+1,j,k) - g%areaT(i+1,j)*min_h
        hlos = max(0.0,uhr(i+1,j,k))
        if ((((hup - hlos) + uhr(i,j,k)) < 0.0) .and. &
            ((0.5*hup + uhr(i,j,k)) < 0.0)) then
           uhh(i) = min(-0.5*hup,-hup+hlos,0.0)
           domore_u(j,k) = .true.
        else
           uhh(i) = uhr(i,j,k)
        endif
        cfl(i) = - uhh(i) / (hprev(i+1,j,k) + h_neglect*g%areaT(i+1,j)) ! CFL is positive
      else
         hup = hprev(i,j,k) - g%areaT(i,j)*min_h
         hlos = max(0.0,-uhr(i-1,j,k))
         if ((((hup - hlos) - uhr(i,j,k)) < 0.0) .and. &
             ((0.5*hup - uhr(i,j,k)) < 0.0)) then
           uhh(i) = max(0.5*hup,hup-hlos,0.0)
           domore_u(j,k) = .true.
         else
           uhh(i) = uhr(i,j,k)
         endif
         cfl(i) = uhh(i) / (hprev(i,j,k) + h_neglect*g%areaT(i,j)) ! CFL is positive
      endif !(uhr(i,j,k) == 0.0) 
    enddo ! i
  
    if (useppm) then
      do m=1,ntr
        do i=is-1,ie
          ! centre cell depending on upstream direction
          if (uhh(i) >= 0.0) then
            i_up = i
          else
            i_up = i+1
          endif
  
          ! Implementation of PPM-H3
          tp = t_tmp(i_up+1,m) ; tc = t_tmp(i_up,m) ; tm = t_tmp(i_up-1,m)
  
          if (usehuynh) then
            al = ( 5.*tc + ( 2.*tm - tp ) )/6. ! H3 estimate
            al = max( min(tc,tm), al) ; al = min( max(tc,tm), al) ! Bound
            ar = ( 5.*tc + ( 2.*tp - tm ) )/6. ! H3 estimate
            ar = max( min(tc,tp), ar) ; ar = min( max(tc,tp), ar) ! Bound
          else
            al = 0.5 * ((tm + tc) + (slope_x(i_up-1,m) - slope_x(i_up,m)) / 3.)
            ar = 0.5 * ((tc + tp) + (slope_x(i_up,m) - slope_x(i_up+1,m)) / 3.)
          endif
  
          da = ar - al ; ma = 0.5*( ar + al )
          if (g%mask2dCu(i_up,j)*g%mask2dCu(i_up-1,j)*(tp-tc)*(tc-tm) <= 0.) then
            al = tc ; ar = tc ! PCM for local extremum and bounadry cells
          elseif ( da*(tc-ma) > (da*da)/6. ) then
            al = 3.*tc - 2.*ar
          elseif ( da*(tc-ma) < - (da*da)/6. ) then
            ar = 3.*tc - 2.*al
          endif
          a6 = 6.*tc - 3. * (ar + al) ! Curvature
  
          if (uhh(i) >= 0.0) then
            flux_x(i,m) = uhh(i)*( ar - 0.5 * cfl(i) * ( &
                ( ar - al ) - a6 * ( 1. - 2./3. * cfl(i) ) ) )
          else
            flux_x(i,m) = uhh(i)*( al + 0.5 * cfl(i) * ( &
                ( ar - al ) + a6 * ( 1. - 2./3. * cfl(i) ) ) )
          endif
        enddo !i
      enddo ! ntr
    else ! PLM
      do m=1,ntr
        do i=is-1,ie
          if (uhh(i) >= 0.0) then
            ! Alternative implementation of PLM
            tc = t_tmp(i,m)
            flux_x(i,m) = uhh(i)*( tc + 0.5 * slope_x(i,m) * ( 1. - cfl(i) ) )
          else
            ! Alternative implementation of PLM
            tc = t_tmp(i+1,m)
            flux_x(i,m) = uhh(i)*( tc - 0.5 * slope_x(i+1,m) * ( 1. - cfl(i) ) )
          endif
        enddo !i 
      enddo ! ntr
    endif ! usePPM
  
    ! Calculate new tracer concentration in each cell after accounting
    ! for the i-direction fluxes.
    do i=is-1,ie
      uhr(i,j,k) = uhr(i,j,k) - uhh(i)
      if (abs(uhr(i,j,k)) < uh_neglect(i,j)) uhr(i,j,k) = 0.0
    enddo
    do i=is,ie
      if ((uhh(i) /= 0.0) .or. (uhh(i-1) /= 0.0)) then
        do_i(i) = .true.
        ! vol [h*area] at t1 (from onl sol)
        hlst(i) = hprev(i,j,k)
        ! vol after this iteration (not necessarily t2)
        if (use_h_end) then
          hprev(i,j,k) = h_end(i,j,k) * g%areaT(i,j)
        else
          hprev(i,j,k) = hprev(i,j,k) - (uhh(i) - uhh(i-1)) 
        endif


        if (hprev(i,j,k) <= 0.0) then
          do_i(i) = .false.
        elseif (hprev(i,j,k) < h_neglect*g%areaT(i,j)) then
          hlst(i) = hlst(i) + (h_neglect*g%areaT(i,j) - hprev(i,j,k))
          ihnew(i) = 1.0 / (h_neglect*g%areaT(i,j))
        else 
          ihnew(i) = 1.0 / hprev(i,j,k)
        endif
      else
        do_i(i) = .false.
      endif
    enddo
  
    ! update tracer concentration from i-flux and save some diagnostics
    do m=1,ntr
  
      ! update tracer
      do i=is,ie
       if (do_i(i)) then
         if (ihnew(i) > 0.0) then
           tr(m)%t(i,j,k) = ( tr(m)%t(i,j,k) * hlst(i) - &
                        (1-alpha(i,j,k)) * (flux_x(i,m) - flux_x(i-1,m)) ) * ihnew(i)
         endif
       endif
      enddo
    
    enddo !ntr
  
  endif ! domore_u
  enddo ! End of j-loop.
  
end subroutine advect_x

!===========================================================
!  advect_y from "mom_tracer_advect"
!===========================================================
!> This subroutine does 1-d flux-form advection using a monotonic piecewise
!! linear scheme.
subroutine advect_y(Tr, hprev, h_end, vhr, vh_neglect, domore_v, ntr, Idt, &
                     is, ie, js, je, k, G, usePPM, useHuynh, alpha)
  type(ocean_grid_type),                     intent(in) :: G    !< The ocean's grid structure
  type(tracer_type), dimension(ntr),         intent(inout) :: Tr   !< The array of registered tracers to work on
  real, dimension(g%isd:g%ied,   g%jsd:g%jed,  g%ke),  intent(inout) :: hprev !< cell volume at the end of previous
                                                                 !! tracer change [H L2 ~> m3 or kg]
  real, dimension(g%isd:g%ied,   g%jsd:g%jed,  g%ke),  intent(in) :: h_end !< [m]
  real, dimension(g%isd:g%ied,   g%jsdb:g%jedb, g%ke), intent(inout) :: vhr !< accumulated volume/mass flux through
                                                                 !! the meridional face [H L2 ~> m3 or kg]
  real, dimension(g%isd:g%ied,   g%jsdb:g%jedb),         intent(inout) :: vh_neglect !< A tiny meridional mass flux that can
                                                                 !! be neglected [H L2 ~> m3 or kg]
  logical, dimension(g%jsdb:g%jedb,g%ke),      intent(inout) :: domore_v !< If true, there is more advection to be
                                                                 !! done in this v-row
  real,                                      intent(in)    :: Idt !< The inverse of dt [T-1 ~> s-1]
  integer,                                   intent(in)    :: ntr !< The number of tracers
  integer,                                   intent(in)    :: is  !< The starting tracer i-index to work on
  integer,                                   intent(in)    :: ie  !< The ending tracer i-index to work on
  integer,                                   intent(in)    :: js  !< The starting tracer j-index to work on
  integer,                                   intent(in)    :: je  !< The ending tracer j-index to work on
  integer,                                   intent(in)    :: k   !< The k-level to work on
  logical,                                   intent(in)    :: usePPM !< If true, use PPM instead of PLM
  logical,                                   intent(in)    :: useHuynh !< If true, use the Huynh scheme
                                                                    !! for PPM interface values
  real, dimension(g%isd:g%ied,   g%jsd:g%jed,  g%ke), &
                                             intent(in)    :: alpha 


  real, dimension(g%isd:g%ied,ntr,g%jsd:g%jed) :: &
   slope_y                     ! The concentration slope per grid point [conc].
  real, dimension(g%isd:g%ied,ntr,g%jsdb:g%jedb) :: &
    flux_y                      ! The tracer flux across a boundary [H m2 conc ~> m3 conc or kg conc].
  real, dimension(g%isd:g%ied,ntr,g%jsdb:g%jedb) :: &
   T_tmp               ! The copy of the tracer concentration at constant i,k [H m2 conc ~> m3 conc or kg conc].
  real :: maxslope              ! The maximum concentration slope per grid point
                               ! consistent with monotonicity [conc].
  real :: vhh(g%isd:g%ied, g%jsdb:g%jedb) ! The meridional flux that occurs during the
                                 ! current iteration [H L2 ~> m3 or kg].
  real :: hup, hlos             ! hup is the upwind volume, hlos is the
                               ! part of that volume that might be lost
                               ! due to advection out the other side of
                               ! the grid box, both in  [H L2 ~> m3 or kg].
  real, dimension(g%isdb:g%iedb) :: &
   hlst, &             ! Work variable [H L2 ~> m3 or kg].
   Ihnew, &            ! Work variable [H-1 L-2 ~> m-3 or kg-1].
   CFL                 ! A nondimensional work variable.
  real :: min_h         ! The minimum thickness that can be realized during
                       ! any of the passes [H ~> m or kg m-2].
  real :: h_neglect     ! A thickness that is so small it is usually lost
                       ! in roundoff and can be neglected [H ~> m or kg m-2].
  logical :: do_j_tr(g%jsd:g%jed)   ! If true, calculate the tracer profiles.
  logical :: do_i(g%isdb:g%iedb)     ! If true, work on given points.
  logical :: do_any_i
  integer :: i, j, j2, m, n, j_up, stencil
  real :: aR, aL, dMx, dMn, Tp, Tc, Tm, dA, mA, a6
  real :: fac1,v_L_in,v_L_out  ! terms used for time-stepping OBC reservoirs
  logical :: usePLMslope, use_h_end
  real, parameter :: Angstrom_H = 1.0E-10, H_subroundoff = 1.0e-20 * max(Angstrom_H, 1.0e-17)

  useplmslope = .not. (useppm .and. usehuynh)
  use_h_end = .false.

  ! stencil for calculating slope values
  stencil = 1
  if (useppm .and. .not. usehuynh) stencil = 2

  min_h = 0.1*Angstrom_H
  h_neglect = H_subroundoff
  
  ! We conditionally perform work on tracer points: calculating the PLM slope,
  ! and updating tracer concentration within a cell
  ! this depends on whether there is a flux which would affect this tracer point,
  ! as indicated by domore_v. In the case of PPM reconstruction, a flux requires
  ! slope calculations at the two tracer points on either side (as indicated by
  ! the stencil variable), so we account for this with the do_j_tr flag array
  !
  ! Note: this does lead to unnecessary work in updating tracer concentrations,
  ! since that doesn't need a wider stencil with the PPM advection scheme, but
  ! this would require an additional loop, etc.
  do_j_tr(:) = .false.

  do j=js-1,je
    if (domore_v(j,k)) then
      do j2=1-stencil,stencil
        do_j_tr(j+j2) = .true.
      enddo
    endif
  enddo

  ! Calculate the j-direction profiles (slopes) of each tracer that
  ! is being advected.
  if (useplmslope) then
    do j=js-stencil,je+stencil
      if (do_j_tr(j)) then
        do m=1,ntr
          do i=is,ie
            tp = tr(m)%t(i,j+1,k) ; tc = tr(m)%t(i,j,k) ; tm = tr(m)%t(i,j-1,k)
            dmx = max( tp, tc, tm ) - tc
            dmn= tc - min( tp, tc, tm )
            slope_y(i,m,j) = g%mask2dCv(i,j)*g%mask2dCv(i,j-1) * &
            sign( min(0.5*abs(tp-tm), 2.0*dmx, 2.0*dmn), tp-tm )
          enddo
        enddo
      endif
    enddo ! End of i-, m-, & j- loops.
  endif ! usePLMslope
  
  ! make a copy of the tracers in case values need to be overridden for OBCs
  
  do j=g%jsd,g%jed
    do m=1,ntr
      do i=g%isd,g%ied
        t_tmp(i,m,j) = tr(m)%t(i,j,k)
      enddo
    enddo
  enddo
  ! loop through open boundaries and recalculate flux terms

  
  ! Calculate the j-direction fluxes of each tracer, using as much
  ! the minimum of the remaining mass flux (vhr) and the half the mass
  ! in the cell plus whatever part of its half of the mass flux that
  ! the flux through the other side does not require.
  do j=js-1,je
    if (domore_v(j,k)) then
      domore_v(j,k) = .false.
  
      do i=is,ie
        if (vhr(i,j,k) == 0.0) then
          vhh(i,j) = 0.0
          cfl(i) = 0.0
        elseif (vhr(i,j,k) < 0.0) then
          hup = hprev(i,j+1,k) - g%areaT(i,j+1)*min_h
          hlos = max(0.0,vhr(i,j+1,k))
          if ((((hup - hlos) + vhr(i,j,k)) < 0.0) .and. &
             ((0.5*hup + vhr(i,j,k)) < 0.0)) then
            vhh(i,j) = min(-0.5*hup,-hup+hlos,0.0)
            domore_v(j,k) = .true.
          else
            vhh(i,j) = vhr(i,j,k)
          endif
          cfl(i) = - vhh(i,j) / (hprev(i,j+1,k) + h_neglect*g%areaT(i,j+1)) ! CFL is positive
        else
         hup = hprev(i,j,k) - g%areaT(i,j)*min_h
         hlos = max(0.0,-vhr(i,j-1,k))
         if ((((hup - hlos) - vhr(i,j,k)) < 0.0) .and. &
             ((0.5*hup - vhr(i,j,k)) < 0.0)) then
           vhh(i,j) = max(0.5*hup,hup-hlos,0.0)
           domore_v(j,k) = .true.
         else
           vhh(i,j) = vhr(i,j,k)
        endif
        cfl(i) = vhh(i,j) / (hprev(i,j,k) + h_neglect*g%areaT(i,j)) ! CFL is positive
        endif
      enddo
  
      if (useppm) then
        do m=1,ntr
          do i=is,ie
            ! centre cell depending on upstream direction
            if (vhh(i,j) >= 0.0) then
              j_up = j
            else
              j_up = j + 1
            endif
    
            ! Implementation of PPM-H3
            tp = tr(m)%t(i,j_up+1,k) ; tc = tr(m)%t(i,j_up,k) ; tm = tr(m)%t(i,j_up-1,k)
    
            if (usehuynh) then
              al = ( 5.*tc + ( 2.*tm - tp ) )/6. ! H3 estimate
              al = max( min(tc,tm), al) ; al = min( max(tc,tm), al) ! Bound
              ar = ( 5.*tc + ( 2.*tp - tm ) )/6. ! H3 estimate
              ar = max( min(tc,tp), ar) ; ar = min( max(tc,tp), ar) ! Bound
            else
              al = 0.5 * ((tm + tc) + (slope_y(i,m,j_up-1) - slope_y(i,m,j_up)) / 3.)
              ar = 0.5 * ((tc + tp) + (slope_y(i,m,j_up) - slope_y(i,m,j_up+1)) / 3.)
            endif
    
            da = ar - al ; ma = 0.5*( ar + al )
            if (g%mask2dCv(i,j_up)*g%mask2dCv(i,j_up-1)*(tp-tc)*(tc-tm) <= 0.) then
              al = tc ; ar = tc ! PCM for local extremum and bounadry cells
            elseif ( da*(tc-ma) > (da*da)/6. ) then
              al = 3.*tc - 2.*ar
            elseif ( da*(tc-ma) < - (da*da)/6. ) then
              ar = 3.*tc - 2.*al
            endif
    
            a6 = 6.*tc - 3. * (ar + al) ! Curvature
    
            if (vhh(i,j) >= 0.0) then
             flux_y(i,m,j) = vhh(i,j)*( ar - 0.5 * cfl(i) * ( &
                  ( ar - al ) - a6 * ( 1. - 2./3. * cfl(i) ) ) )
            else
             flux_y(i,m,j) = vhh(i,j)*( al + 0.5 * cfl(i) * ( &
                  ( ar - al ) + a6 * ( 1. - 2./3. * cfl(i) ) ) )
            endif
          enddo
        enddo ! ntr & i
      else ! PLM
        do m=1,ntr
          do i=is,ie
            if (vhh(i,j) >= 0.0) then
              tc = tr(m)%t(i,j,k)
              flux_y(i,m,j) = vhh(i,j)*( tc + 0.5 * slope_y(i,m,j) * ( 1. - cfl(i) ) )
            else
              tc = tr(m)%t(i,j+1,k)
              flux_y(i,m,j) = vhh(i,j)*( tc - 0.5 * slope_y(i,m,j+1) * ( 1. - cfl(i) ) )
            endif
          enddo
        enddo
      endif ! usePPM
  
    else ! not domore_v.
      do i=is,ie
        vhh(i,j) = 0.0
      enddo
      do m=1,ntr
        do i=is,ie
          flux_y(i,m,j) = 0.0
        enddo
      enddo
    endif ! domore_v.
  enddo ! End of j-loop
  
  do j=js-1,je
    do i=is,ie
      vhr(i,j,k) = vhr(i,j,k) - vhh(i,j)
      if (abs(vhr(i,j,k)) < vh_neglect(i,j)) vhr(i,j,k) = 0.0
    enddo
  enddo
  
  ! Calculate new tracer concentration in each cell after accounting
  ! for the j-direction fluxes.
  do j=js,je
    if (do_j_tr(j)) then

      do i=is,ie
        if ((vhh(i,j) /= 0.0) .or. (vhh(i,j-1) /= 0.0)) then
          do_i(i) = .true.
          hlst(i) = hprev(i,j,k)

          ! vol after this iteration (not necessarily t2)
          if (use_h_end) then
            hprev(i,j,k) = h_end(i,j,k) * g%areaT(i,j)
          else
            hprev(i,j,k) = max(hprev(i,j,k) - (vhh(i,j) - vhh(i,j-1)), 0.0)
          endif

          if (hprev(i,j,k) <= 0.0) then
            do_i(i) = .false.
          elseif (hprev(i,j,k) < h_neglect*g%areaT(i,j)) then
            hlst(i) = hlst(i) + (h_neglect*g%areaT(i,j) - hprev(i,j,k))
            ihnew(i) = 1.0 / (h_neglect*g%areaT(i,j))
          else
            ihnew(i) = 1.0 / hprev(i,j,k)
          endif
        else
          do_i(i) = .false.
        endif
      enddo ! i
  
      ! update tracer and save some diagnostics
      do m=1,ntr

        do i=is,ie
          if (do_i(i)) then
            tr(m)%t(i,j,k) = ( tr(m)%t(i,j,k) * hlst(i) - &
                  (1-alpha(i,j,k))* (flux_y(i,m,j) - flux_y(i,m,j-1)) ) * ihnew(i)
          endif
        enddo
  
        
      enddo ! m

    endif !(do_j_tr(j))
  enddo ! End of j-loop.
  
 end subroutine advect_y

!===========================================================
!  apply tracer eddy forcing
!===========================================================



end module mom_offline_main
