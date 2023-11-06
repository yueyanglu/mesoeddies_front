
module mom_tracer_aux

use mom_tracer_registry,      only : tracer_registry_type, tracer_type
use mod_grid_my, only : ocean_grid_type

implicit none ; public
real, save, allocatable, dimension (:,:,:) :: var_test, alpha, gam
real, save, allocatable, dimension(:,:,:) :: h, h_t2, h_cp &
     , dp1, dp2, trTemp  &
     , uhtr, vhtr, h_mod &
     , u1, v1, u2, v2, uh_cp, vh_cp &
     , kxx, kxy,  kyx, kyy, ka, chiu, chiv &
     , chin, ka_cp

real, save, allocatable, dimension(:,:,:,:) :: forc, relxforc, relx, &
       relx_pt, forc_ka, forc_chi, forc_kchi, forc_divuc, forc_cdivu, &
      forc_udelc, forc_temp, dcdt, c_t1, c_t2, dcmdt, cm_t1, cm_t2, &
      c_tm, trsign


contains

! ====================================================================
!           allocate vars
! ====================================================================
subroutine vars_alloc(g, NTR)
   implicit none
   type(ocean_grid_type),      intent(in) :: g  
   integer,                    intent(in) :: NTR 

   !
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
   allocate( ka(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
   allocate( ka_cp(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
   allocate( chiu(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
   allocate( chiv(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
   ! 
   allocate( kxx(g%isdb:g%iedb, g%jsd:g%jed,  g%ke) )
   allocate( kxy(g%isdb:g%iedb, g%jsd:g%jed,  g%ke) )
   allocate( kyx(g%isd:g%ied,   g%jsdb:g%jedb,  g%ke) )
   allocate( kyy(g%isd:g%ied,   g%jsdb:g%jedb,  g%ke) )
   !
   allocate( forc    (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   allocate( relxforc(g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   allocate( relx    (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   allocate( relx_pt (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   allocate( forc_ka     (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   allocate( forc_chi    (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   allocate( forc_kchi   (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   !
   allocate( forc_divuc  (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   allocate( forc_cdivu  (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   allocate( forc_udelc  (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   allocate( forc_temp   (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   allocate( dcdt   (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   allocate( c_t1   (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   allocate( c_t2   (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   allocate( cm_t1  (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   allocate( cm_t2  (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   allocate( dcmdt  (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   allocate( c_tm   (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )
   allocate( trsign  (g%isd:g%ied,   g%jsd:g%jed,  g%ke, NTR) )

   !
   allocate( alpha(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
   allocate( gam(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )
   allocate( chin(g%isd:g%ied,   g%jsd:g%jed,  g%ke) )

end subroutine vars_alloc

! ====================================================================
!           update tracer conc by the point-wise forcing 
! =====================================================================
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

! ====================================================================
!           update tracer conc by the point-wise forcing 
! =====================================================================
! solve d(ch)/dt = forc [m/s*c]
! c2 = (c1*h1 + forc*dt) / h2
! 
subroutine tracer_forc(Reg, dt, forc, h1, h2, G)
   type(ocean_grid_type),      intent(in) :: g       !< Grid type
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                              intent(in)    :: h1, h2   !< t1 and t2
   real,                       intent(in)    :: dt      !< time step [T ~> s]
   type(tracer_registry_type), intent(inout)       :: reg     !< registered tracers
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke, reg%ntr), &
                              intent(in)    :: forc   ! [m/s*c]  

   ! 
   integer :: i, j, k, m, is, ie, js, je, nz, ntr
   real, parameter :: h_neglect = 1.0e-20 * max(1.0e-10, 1.0e-17)

   is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec ; nz = g%ke
   ntr = reg%ntr

   ! ------------ update tracers ------------
   do m = 1,ntr
      do k = 1,nz
         do j = js,je
            do i = is,ie
               reg%Tr(m)%t(i,j,k) = ( reg%Tr(m)%t(i,j,k) * h1(i,j,k) &
                  + dt*forc(i,j,k,m) ) / (h2(i,j,k)+h_neglect)          
            enddo
         enddo
      enddo ! k
   enddo ! ntr

end subroutine tracer_forc

! ====================================================================
!           relax tracer toward a profile, and output relax forcing
! =====================================================================
subroutine tracer_relx(Reg, Reg0, dt, cprof, T, g, h1, h2, relxforc)
   ! d(ch)/dt + div{Uc} = ... + 1/T*h*(c_prof - c)
   ! c2 = ( c1*h1 + dt * h2/T*(c_prof- c) ) / h2
   type(ocean_grid_type),      intent(in)    :: g       !< Grid type
   type(tracer_registry_type), intent(inout)       :: reg     !< registered tracers
   type(tracer_registry_type), intent(in)       :: reg0     !used to calc RF
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke, reg%ntr), &
                               intent(in)    :: cprof   ! rel profile [c] 
   real,                       intent(in)    :: dt      !< time step [T ~> s]
   real,                       intent(in)    :: T       ! time scale of the relax [s]
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                               intent(in)    :: h1, h2   ! lay thck 
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke, reg%ntr), &
                               intent(out)   :: relxforc ! forc [m/s*c] 
                                                         ! = h*r*(c_prof - c)
   ! local vars
   real :: rate
   integer :: i, j, k, m, is, ie, js, je, nz, ntr
   real, parameter :: h_neglect = 1.0e-20 * max(1.0e-10, 1.0e-17)

   is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec ; nz = g%ke
   ntr = reg%ntr

   rate = 1. / T ! rate (strength) of relax [s]

   relxforc = 0.
   ! ------------ update tracers ------------
   do m = 1,ntr
      do k = 1,1
         do j = js,je
            do i = is,ie
               relxforc(i,j,k,m) = h2(i,j,k) * rate * &
                     ( cprof(i,j,k,m) - reg0%Tr(m)%t(i,j,k) ) ! m/s*c
               !
               reg%Tr(m)%t(i,j,k) = ( reg%Tr(m)%t(i,j,k) * h1(i,j,k) &
                  + dt * relxforc(i,j,k,m) ) / (h2(i,j,k)+h_neglect)  
            enddo
         enddo
      enddo
   enddo ! ntr

end subroutine tracer_relx

! ====================================================================
!           calc the div{U} from given U 
! =====================================================================
subroutine calc_div(uh, vh, div, g)
   type(ocean_grid_type),      intent(in)    :: g       !< Grid type
   real, dimension(g%isdb:g%iedb, g%jsd:g%jed, g%ke), &
                               intent(in)    :: uh      ! uhL [m3/s] or [m3]
   real, dimension(g%isd:g%ied, g%jsdb:g%jedb, g%ke), &
                               intent(in)    :: vh      ! vhL 
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                               intent(out)   :: div   ! [uh / m2] 

   ! local vars
   integer :: i, j, k, is, ie, js, je, nz
   is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec ; nz = g%ke

   do k = 1,nz
      do j = js,je
         do i = is,ie
            div(i,j,k) = g%IareaT(i,j) * &
               ( uh(i,j,k) - uh(i-1,j,k) + vh(i,j,k) - vh(i,j-1,k) )          
         enddo
      enddo
  enddo

end subroutine calc_div

! ====================================================================
!           calc the forcing = div{Uc} from given U and c
! =====================================================================
subroutine calc_divforc(Reg, uh, vh, divuc, cdivu, udelc, g)
   type(ocean_grid_type),      intent(in)    :: g       !< Grid type
   type(tracer_registry_type), intent(in)       :: reg     !< registered tracers
   real, dimension(g%isdb:g%iedb, g%jsd:g%jed, g%ke), &
                               intent(in)    :: uh      ! uhL [m3/s]
   real, dimension(g%isd:g%ied, g%jsdb:g%jedb, g%ke), &
                               intent(in)    :: vh      ! vhL [m3/s]
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke, reg%ntr), &
                               intent(out)   :: divuc   ! [m/s*c] 
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke, reg%ntr), &
                     optional, intent(out)   :: cdivu, udelc  

   ! local vars
   integer :: i, j, k, m, is, ie, js, je, nz, ntr
   real, dimension(g%isdb:g%iedb, g%jsd:g%jed) :: uhc ! uhLc [m3/s*c]
   real, dimension(g%isd:g%ied, g%jsdb:g%jedb) :: vhc

   is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec ; nz = g%ke
   ntr = reg%ntr

   divuc = 0.; cdivu = 0.; udelc = 0.

   do m = 1,ntr
      do k = 1,nz

         ! see "tracer_hordiff"
         do j = js,je
            do i = is-1,ie
               ! m3/s*c
               uhc(i,j) = uh(i,j,k) * 0.5 * ( reg%Tr(m)%t(i,j,k) + &
                  reg%Tr(m)%t(i+1,j,k) ) 
            enddo
         enddo
         do j = js-1,je
            do i = is,ie
               vhc(i,j) = vh(i,j,k) * 0.5 * ( reg%Tr(m)%t(i,j,k) + &
                  reg%Tr(m)%t(i,j+1,k) ) 
            enddo
         enddo

         do j = js,je
            do i = is,ie
               ! div{UC}
               divuc(i,j,k,m) = g%IareaT(i,j) * &
                  ( uhc(i,j) - uhc(i-1,j) + vhc(i,j) - vhc(i,j-1) )      
               ! c*divU
               cdivu(i,j,k,m) = reg%Tr(m)%t(i,j,k) * g%IareaT(i,j) * &
               ( uh(i,j,k) - uh(i-1,j,k) + vh(i,j,k) - vh(i,j-1,k) ) 
               ! U*delc =  div{UC} - c*divU
               udelc(i,j,k,m) = divuc(i,j,k,m) - cdivu(i,j,k,m)    
            enddo
         enddo
      enddo
   enddo ! ntr

end subroutine calc_divforc

! ====================================================================
!           calc the forcing = div{Uc} from given U and c
! =====================================================================
subroutine calc_divforc_trmat(tr, ntr, uh, vh, divuc, g)
   type(ocean_grid_type),      intent(in)    :: g       !< Grid type
   real, dimension(g%isdb:g%iedb, g%jsd:g%jed, g%ke), &
                               intent(in)    :: uh      ! uhL [m3/s]
   real, dimension(g%isd:g%ied, g%jsdb:g%jedb, g%ke), &
                               intent(in)    :: vh      ! vhL [m3/s]
   integer, intent(in)                       :: ntr
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke, ntr), &
                               intent(in)    :: tr
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke, ntr), &
                               intent(out)   :: divuc   ! [m/s*c] 

   ! local vars
   integer :: i, j, k, m, is, ie, js, je, nz
   real, dimension(g%isdb:g%iedb, g%jsd:g%jed) :: uhc ! uhLc [m3/s*c]
   real, dimension(g%isd:g%ied, g%jsdb:g%jedb) :: vhc

   is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec ; nz = g%ke

   divuc = 0.

   do m = 1,ntr
      do k = 1,nz

         ! see "tracer_hordiff"
         do j = js,je
            do i = is-1,ie
               ! m3/s*c
               uhc(i,j) = uh(i,j,k) * 0.5 * ( tr(i,j,k,m) + &
                  tr(i+1,j,k,m) ) 
            enddo
         enddo
         do j = js-1,je
            do i = is,ie
               vhc(i,j) = vh(i,j,k) * 0.5 * ( tr(i,j,k,m) + &
                  tr(i,j+1,k,m) ) 
            enddo
         enddo

         do j = js,je
            do i = is,ie
               ! div{UC}
               divuc(i,j,k,m) = g%IareaT(i,j) * &
                  ( uhc(i,j) - uhc(i-1,j) + vhc(i,j) - vhc(i,j-1) ) 
            enddo
         enddo
      enddo
   enddo ! ntr

end subroutine calc_divforc_trmat

! ====================================================================
!           Parameterization: K-tensor model
! =====================================================================
! update tracer w/ the parameterized forc: 
!   d(ch)/dt = div{ Fu/v } [m/s*c],  where Fu/v = K*hdelC [m2/s*c]
! c2 = ( c1*h1 + dt*EF ) / h2
! 
! Fu = Kxx * h*dCdx + Kxy h*dCdy
! Fv = Kyx * h*dCdx + Kyy h*dCdy
! 
subroutine tracer_ktensor(Reg, Reg0, dt, kxx, kxy, kyx, kyy, h0, h1, h2, G)
   type(tracer_registry_type), intent(inout)    :: reg  
   type(tracer_registry_type), intent(in)       :: reg0 ! 
   type(ocean_grid_type),      intent(in)       :: g      
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                               intent(in)       :: h1, h2, h0 ! before adv  
   real,                       intent(in)       :: dt   ![s]  
   real, dimension(g%isdb:g%iedb, g%jsd:g%jed, g%ke), &
                               intent(in)       :: kxx, kxy  ! [m2/s]
   real, dimension(g%isd:g%ied, g%jsdb:g%jedb, g%ke), &
                               intent(in)       :: kyx, kyy
! local
   real, parameter :: h_neglect = 1.0e-20 * max(1.0e-10, 1.0e-17)
   integer :: i, j, k, m, is, ie, js, je, nz, ntr  &
            , ia, ib, ja, jb
   real, dimension(g%isdb:g%iedb, g%jsd:g%jed)   :: &
               coef_x, dcdx, fu
               ! dt times the open face width divided by
               ! the distance between adjacent tracer points [s].
   real, dimension(g%isd:g%ied,   g%jsdb:g%jedb) :: &
               coef_y, dcdy, fv
   real, dimension(g%isd:g%ied, g%jsd:g%jed) :: &
               dtr           &  ! [c*m]
              ,fmn, fmx, fld, fld0
   real :: q_new  ! updated tracer mass at one point [cm]

   is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec ; nz = g%ke
   ntr = reg%ntr

   ! this is in case that Fortran assign non-zero values to initial fld
   ! Note the 4-point average will use the points that were not touched
   ! in the "dcdx" loop
   dcdx = 0.; dcdy = 0.

   ! ------------ each layer
   do k = 1,nz

      ! --------- dt/dx*h [s]
      do j = js,je
         do i = is-1,ie
            coef_x(i,j) = dt*g%IdxCu(i,j) *2.0*h0(i,j,k)*h0(i+1,j,k) / &
                          ( h0(i,j,k)+h0(i+1,j,k)+h_neglect ) 
         enddo
      enddo
      do j = js-1,je
         do i = is,ie
            coef_y(i,j) = dt*g%IdyCv(i,j) *2.0*h0(i,j,k)*h0(i,j+1,k) / &
                          ( h0(i,j,k)+h0(i,j+1,k)+h_neglect )
         enddo
      enddo

      ! --------- update tracers
      do m = 1,ntr

         ! --------- assign
         fld(:,:) = reg%Tr(m)%t(:,:,k)
         fld0(:,:) = reg0%Tr(m)%t(:,:,k)

         ! --------- min & max values
         !$OMP parallel do private(i,j,ja,jb,ia,ib)
         do j = js, je
            ja = j - 1
            jb = j + 1
            do i = is,ie
               if (g%mask2dT(i,j) /= 0) then
                  ia = i - 1
                  if (g%mask2dT(ia,j)==0) ia = i
                  ib = i + 1
                  if (g%mask2dT(ib,j)==0) ib = i
                  ja = j - 1
                  if (g%mask2dT(i,ja)==0) ja = j
                  jb = j + 1
                  if (g%mask2dT(i,jb)==0) jb = j
                  fmx(i,j) = max( fld(i,j),     &
                      fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
                  fmn(i,j) = min( fld(i,j),     &
                      fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
               endif
            enddo
          enddo
         !$OMP end parallel do

         ! --------- trac grad: dt/dx*h *dc [s*c] and diag eddy flux [m2*c]
         do j = js,je
            do i = is-1,ie
               dcdx(i,j) = g%mask2dCu(i,j)* coef_x(i,j) * &
                  ( fld0(i+1,j) - fld0(i,j) )
               ! diagonal: kxx * dcdx
               fu(i,j) = kxx(i,j,k) * dcdx(i,j)
            enddo
         enddo
         do j = js-1,je
            do i = is,ie
               dcdy(i,j) = g%mask2dCv(i,j)* coef_y(i,j) * &
                  ( fld0(i,j+1) - fld0(i,j) )
               fv(i,j) = kyy(i,j,k) * dcdy(i,j)
            enddo
         enddo

         ! --------- off-diagonal eddy flux [m2*c]
         ! kxy * dcdy (v aver onto u-, clockwise)
         do j = js,je
            do i = is-1,ie
               fu(i,j) = fu(i,j) + kxy(i,j,k) * 0.25 * & 
                     (dcdy(i,j) + dcdy(i+1,j) + dcdy(i+1,j-1) + dcdy(i,j-1))
            enddo
         enddo
         ! kyx * dcdx (u aver onto v-, clockwise)
         do j = js-1,je
            do i = is,ie
               ! off-diagonal:  kyx * dcdx(u aver onto v-, clockwise)
               fv(i,j) = fv(i,j) + kyx(i,j,k) * 0.25 * & 
                  (dcdx(i,j) + dcdx(i-1,j) + dcdx(i-1,j+1) + dcdx(i,j+1))
            enddo
         enddo

         ! --------- flux div [cm] ~ delta(c*h)
         do j = js,je
            do i = is,ie
               dtr(i,j) = g%IareaT(i,j) * &
                  ( g%dy_Cu(i,j)*fu(i,j) - g%dy_Cu(i-1,j)*fu(i-1,j) &
                  + g%dx_Cv(i,j)*fv(i,j) - g%dx_Cv(i,j-1)*fv(i,j-1) )
            enddo
         enddo

         ! --------- update tr
         do j = js,je
            do i = is,ie
               q_new = fld(i,j)*h1(i,j,k) + dtr(i,j)
              !  q_new = max( 0.0, q_new)
               q_new = max( fmn(i,j)*h2(i,j,k),       &
                            min(q_new, fmx(i,j)*h2(i,j,k)) )
               reg%Tr(m)%t(i,j,k) = q_new / (h2(i,j,k)+h_neglect)  
            enddo
         enddo
      enddo ! ntr
      
   enddo ! k
end subroutine tracer_ktensor

! ====================================================================
!           Parameterization: k-Chi model
! =====================================================================
! update tracer ('Reg') by the parameterized forc: 
!   d(ch)/dt = ka*del2C - Chi* h*delC [m/s*c], 
!    where ka [m2/s], Chi [m/s], 
! Note the forcing is calculated using the tracer before adv (Reg0) !
! c2 = ( c1*h1 + dt*EF ) / h2
! 
subroutine tracer_kachi(Reg, Reg0, dt, ka, chiu, chiv, h0, h1, h2, G &
            , forc_ka, forc_chi, forc_kchi)
   type(tracer_registry_type), intent(inout)    :: reg 
   type(tracer_registry_type), intent(in)       :: reg0 ! 
   type(ocean_grid_type),      intent(in)       :: g      
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                               intent(in)       :: h1, h2, h0  
   real,                       intent(in)       :: dt   ![s]  
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                               intent(in)       :: ka, chiu, chiv   
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke, reg%ntr), &
                     optional, intent(out)   :: forc_ka, &
                                                forc_chi, &
                                                forc_kchi ! forc [m/s*c] 

! local
   real, parameter :: h_neglect = 1.0e-20 * max(1.0e-10, 1.0e-17)
   integer :: i, j, k, m, is, ie, js, je, nz, ntr &
            , ia, ib, ja, jb
   real, dimension(g%isdb:g%iedb, g%jsd:g%jed)   :: &
               dt_x, coef_x, dcdx
               ! dt times the open face width divided by
               ! the distance between adjacent tracer points [s].
   real, dimension(g%isd:g%ied,   g%jsdb:g%jedb) :: &
               dt_y, coef_y, dcdy
   real, dimension(g%isd:g%ied,   g%jsd:g%jed) :: &
               dtr_ka, dtr_chi, dtr         & ! [cm]
              ,fmn, fmx, fld, fld0
   real :: q_new  ! updated tracer mass at one point [cm]

   is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec ; nz = g%ke
   ntr = reg%ntr

   ! ------------ dt/dx [s/m] 
   do j = js,je
      do i = is-1,ie
         dt_x(i,j) = dt * g%IdxCu(i,j)
      enddo
   enddo
   do j = js-1,je
      do i = is,ie
         dt_y(i,j) = dt * g%IdyCv(i,j)
      enddo
   enddo

   do k = 1,nz

      ! --------- dt/dx*h [s]
      do j = js,je
         do i = is-1,ie
            coef_x(i,j) = dt_x(i,j) *2.0*h0(i,j,k)*h0(i+1,j,k) / &
                          ( h0(i,j,k)+h0(i+1,j,k)+h_neglect ) 
         enddo
      enddo
      do j = js-1,je
         do i = is,ie
            coef_y(i,j) = dt_y(i,j) *2.0*h0(i,j,k)*h0(i,j+1,k) / &
                          ( h0(i,j,k)+h0(i,j+1,k)+h_neglect )
         enddo
      enddo

      ! --------- update tracers
      do m = 1,ntr
         ! assign
         fld(:,:) = reg%Tr(m)%t(:,:,k)
         fld0(:,:) = reg0%Tr(m)%t(:,:,k)

         ! --------- min & max values
         !$OMP parallel do private(i,j,ja,jb,ia,ib)
         do j = js, je
            ja = j - 1
            jb = j + 1
            do i = is,ie
               if (g%mask2dT(i,j) /= 0) then
                  ia = i - 1
                  if (g%mask2dT(ia,j)==0) ia = i
                  ib = i + 1
                  if (g%mask2dT(ib,j)==0) ib = i
                  ja = j - 1
                  if (g%mask2dT(i,ja)==0) ja = j
                  jb = j + 1
                  if (g%mask2dT(i,jb)==0) jb = j
                  fmx(i,j) = max( fld(i,j),     &
                      fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
                  fmn(i,j) = min( fld(i,j),     &
                      fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
               endif
            enddo
          enddo
         !$OMP end parallel do

         ! --------- trac grad: dt/dx*h *dc [s*c], using un-adv TRACER!!! 
         ! loop over all sym-grid points; 0 for boundary points !!
         do j = js,je
            do i = is-1,ie
               dcdx(i,j) = g%mask2dCu(i,j)* coef_x(i,j) * &
                  ( fld0(i+1,j) - fld0(i,j) )
            enddo
         enddo
         do j = js-1,je
            do i = is,ie
               dcdy(i,j) = g%mask2dCv(i,j)* coef_y(i,j) * &
                  ( fld0(i,j+1) - fld0(i,j) )
            enddo
         enddo

        ! --------- param EF [cm] ~ delta(c*h)
         do j = js,je
            do i = is,ie
               ! ka * del2C [cm] = m2/s * 1/m2 * m*cs 
               dtr_ka(i,j) = ka(i,j,k) * g%IareaT(i,j) * &
                  ( g%dy_Cu(i,j)*dcdx(i,j) - g%dy_Cu(i-1,j)*dcdx(i-1,j) &
                  + g%dx_Cv(i,j)*dcdy(i,j) - g%dx_Cv(i,j-1)*dcdy(i,j-1) )
               ! chi*delC [cm] = m/s * cs
               dtr_chi(i,j) =  &
                  - chiu(i,j,k) * 0.5*(dcdx(i,j)+dcdx(i-1,j)) &
                  - chiv(i,j,k) * 0.5*(dcdy(i,j)+dcdy(i,j-1))
               dtr(i,j) = dtr_ka(i,j) + dtr_chi(i,j)
               
               if(present(forc_ka) .and. present(forc_chi) .and. present(forc_kchi)) then
                  forc_ka(i,j,k,m) = dtr_ka(i,j) / dt
                  forc_chi(i,j,k,m) = dtr_chi(i,j) / dt
                  forc_kchi(i,j,k,m) = dtr(i,j) / dt
               endif
            enddo
         enddo
         !5:68 ~ 1:64
         if (k == 1) then
           !print*, 'Tr:', fld0(22,44:48)
           !print*, 'params:', ka(22,44:48,k), chiu(22,44:48,k), chiv(22,44:48,k)
           !print*, 'EF kchi:', forc_kchi(22,44:48,1,1)
           !print*, 'dcdx:', dcdx(22,44:48) / dt
           !print*, 'dcdy:', dcdy(22,44:48) / dt
         endif
         ! print*, "max dtr_ka(1)", maxval(dtr_ka(is:ie,js:je)), minval(dtr_ka(is:ie,js:je))
         ! print*, "max dtr_chi(1)", maxval(dtr_chi(is:ie,js:je)), minval(dtr_chi(is:ie,js:je))

        ! --------- update tr
         do j = js,je
            do i = is,ie
               q_new = fld(i,j)*h1(i,j,k) + dtr(i,j)
               ! q_new = max( 0.0, q_new)
               q_new = max( fmn(i,j)*h2(i,j,k),       &
                            min(q_new, fmx(i,j)*h2(i,j,k)) )
               reg%Tr(m)%t(i,j,k) = q_new / (h2(i,j,k)+h_neglect) 
               !  Eliminate very small values far from tracer patch
               if (reg%Tr(m)%t(i,j,k) < 1.e-7)  reg%Tr(m)%t(i,j,k) = 0.0
            enddo
         enddo
      enddo !m

   enddo ! k

end subroutine tracer_kachi

! ====================================================================
!           Parameterization: feed chi*n (\chi_{\parallel})
! =====================================================================
! update tracer by: 
!   d(ch)/dt = ka*del2C - Chi*  h*delC [m/s*c], 
!            = ka*del2C - chidotn |h*delC|, 
!   where chidotn = |chi*n| is a scalar input
! 
! Note the forcing is calculated using the tracer before adv (Reg0) !
! c2 = ( c1*h1 + dt*EF ) / h2
subroutine tracer_kachi_chidotn(Reg, Reg0, dt, ka, chidotn, &
            h0, h1, h2, G, forc_ka, forc_chi, forc_kchi)
   type(tracer_registry_type), intent(inout)    :: reg 
   type(tracer_registry_type), intent(in)       :: reg0 ! 
   type(ocean_grid_type),      intent(in)       :: g      
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                               intent(in)       :: h1, h2, h0  
   real,                       intent(in)       :: dt   ![s]  
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                               intent(in)       :: ka, chidotn 
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke, reg%ntr), &
                     optional, intent(out)   :: forc_ka, &
                                                forc_chi, &
                                                forc_kchi  

   ! local
   real, parameter :: h_neglect = 1.0e-20 * max(1.0e-10, 1.0e-17)
   integer :: i, j, k, m, is, ie, js, je, nz, ntr &
            , ia, ib, ja, jb
   real, dimension(g%isdb:g%iedb, g%jsd:g%jed)   :: &
               dt_x, coef_x, dcdx
               ! dt times the open face width divided by
               ! the distance between adjacent tracer points [s].
   real, dimension(g%isd:g%ied,   g%jsdb:g%jedb) :: &
               dt_y, coef_y, dcdy
   real, dimension(g%isd:g%ied,   g%jsd:g%jed) :: &
               dtr_ka, dtr_chi, dtr         & ! [cm]
              ,fmn, fmx, fld, fld0, delc2   &
              , hdelc_norm 
   real :: q_new  ! updated tracer mass at one point [cm]


   is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec ; nz = g%ke
   ntr = reg%ntr

   ! ------------ dt/dx [s/m] 
   do j = js,je
      do i = is-1,ie
         dt_x(i,j) = dt * g%IdxCu(i,j)
      enddo
   enddo
   do j = js-1,je
      do i = is,ie
         dt_y(i,j) = dt * g%IdyCv(i,j)
      enddo
   enddo

   do k = 1,nz

      ! --------- dt/dx*h [s]
      do j = js,je
         do i = is-1,ie
            coef_x(i,j) = dt_x(i,j) *2.0*h0(i,j,k)*h0(i+1,j,k) / &
                          ( h0(i,j,k)+h0(i+1,j,k)+h_neglect ) 
         enddo
      enddo
      do j = js-1,je
         do i = is,ie
            coef_y(i,j) = dt_y(i,j) *2.0*h0(i,j,k)*h0(i,j+1,k) / &
                          ( h0(i,j,k)+h0(i,j+1,k)+h_neglect )
         enddo
      enddo

      ! --------- update tracers
      do m = 1,ntr
         ! assign
         fld(:,:) = reg%Tr(m)%t(:,:,k)
         fld0(:,:) = reg0%Tr(m)%t(:,:,k)

         ! --------- min & max values
         !$OMP parallel do private(i,j,ja,jb,ia,ib)
         do j = js, je
            ja = j - 1
            jb = j + 1
            do i = is,ie
               if (g%mask2dT(i,j) /= 0) then
                  ia = i - 1
                  if (g%mask2dT(ia,j)==0) ia = i
                  ib = i + 1
                  if (g%mask2dT(ib,j)==0) ib = i
                  ja = j - 1
                  if (g%mask2dT(i,ja)==0) ja = j
                  jb = j + 1
                  if (g%mask2dT(i,jb)==0) jb = j
                  fmx(i,j) = max( fld(i,j),     &
                      fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
                  fmn(i,j) = min( fld(i,j),     &
                      fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
               endif
            enddo
          enddo
         !$OMP end parallel do

         ! --------- trac grad: dt/dx*h *dc [s*c], using un-adv TRACER!!! 
         ! loop over all sym-grid points; 0 for boundary points !!
         do j = js,je
            do i = is-1,ie
               dcdx(i,j) = g%mask2dCu(i,j)* coef_x(i,j) * &
                  ( fld0(i+1,j) - fld0(i,j) )
            enddo
         enddo
         do j = js-1,je
            do i = is,ie
               dcdy(i,j) = g%mask2dCv(i,j)* coef_y(i,j) * &
                  ( fld0(i,j+1) - fld0(i,j) )
            enddo
         enddo

        ! --------- param EF [cm] ~ delta(c*h)
         do j = js,je
            do i = is,ie
               ! ka * del2C [cm] = m2/s * 1/m2 * m*c*s 
               dtr_ka(i,j) = ka(i,j,k) * g%IareaT(i,j) * &
                  ( g%dy_Cu(i,j)*dcdx(i,j) - g%dy_Cu(i-1,j)*dcdx(i-1,j) &
                  + g%dx_Cv(i,j)*dcdy(i,j) - g%dx_Cv(i,j-1)*dcdy(i,j-1) )
               ! -chi*delC [cm] = - chidotn  |h*delC| = [m/s] [c*s]
               hdelc_norm(i,j) = (  (0.5*(dcdx(i,j)+dcdx(i-1,j)))**2. + &
                     (0.5*(dcdy(i,j)+dcdy(i,j-1)))**2.  ) **.5
               dtr_chi(i,j) = - chidotn(i,j,k) * hdelc_norm(i,j)
               !
               dtr(i,j) = dtr_ka(i,j) + dtr_chi(i,j)
               
               if(present(forc_ka) .and. present(forc_chi) .and. present(forc_kchi)) then
                  forc_ka(i,j,k,m) = dtr_ka(i,j) / dt
                  forc_chi(i,j,k,m) = dtr_chi(i,j) / dt
                  forc_kchi(i,j,k,m) = dtr(i,j) / dt
               endif
            enddo
         enddo

        ! --------- update tr
         do j = js,je
            do i = is,ie
               q_new = fld(i,j)*h1(i,j,k) + dtr(i,j)
               ! q_new = max( 0.0, q_new)
               q_new = max( fmn(i,j)*h2(i,j,k),       &
                            min(q_new, fmx(i,j)*h2(i,j,k)) )
               reg%Tr(m)%t(i,j,k) = q_new / (h2(i,j,k)+h_neglect) 
               !  Eliminate very small values far from tracer patch
               if (reg%Tr(m)%t(i,j,k) < 1.e-7)  reg%Tr(m)%t(i,j,k) = 0.0
            enddo
         enddo
      enddo !m
   enddo ! k

end subroutine tracer_kachi_chidotn

! ====================================================================
!           Parameterization: feed chi*n (\chi_{\parallel})
! =====================================================================
! update tracer by: 
!   d(ch)/dt = ka*del2C - Chi*  h*delC [m/s*c], 
!            = ka*del2C - chidotn |h*delC|, 
!   where chidotn = |chi*n| is a scalar input
! 
! Note the forcing is calculated using the tracer before adv (Reg0) !
! c2 = ( c1*h1 + dt*EF ) / h2
subroutine tracer_kachi_chidotn_sign(Reg, Reg0, dt, ka, chidotn, &
            h0, h1, h2, G, forc_ka, forc_chi, forc_kchi)
   type(tracer_registry_type), intent(inout)    :: reg 
   type(tracer_registry_type), intent(in)       :: reg0 ! 
   type(ocean_grid_type),      intent(in)       :: g      
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                               intent(in)       :: h1, h2, h0  
   real,                       intent(in)       :: dt   ![s]  
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                               intent(in)       :: ka, chidotn 
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke, reg%ntr), &
                     optional, intent(out)   :: forc_ka, &
                                                forc_chi, &
                                                forc_kchi  

   ! local
   real, parameter :: h_neglect = 1.0e-20 * max(1.0e-10, 1.0e-17)
   integer :: i, j, k, m, is, ie, js, je, nz, ntr &
            , ia, ib, ja, jb
   real, dimension(g%isdb:g%iedb, g%jsd:g%jed)   :: &
               dt_x, coef_x, dcdx
               ! dt times the open face width divided by
               ! the distance between adjacent tracer points [s].
   real, dimension(g%isd:g%ied,   g%jsdb:g%jedb) :: &
               dt_y, coef_y, dcdy
   real, dimension(g%isd:g%ied,   g%jsd:g%jed) :: &
               dtr_ka, dtr_chi, dtr         & ! [cm]
              ,fmn, fmx, fld, fld0, delc2   &
              , hdelc_norm                  &
              , dcdy_p

   real, dimension(g%jsd:g%jed) :: dcdy_xmean, sign1d

   real :: q_new  ! updated tracer mass at one point [cm]
   real :: dtr_sum, dtr_abs_sum ! [D], [|D|]


   is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec ; nz = g%ke
   ntr = reg%ntr

   ! ------------ dt/dx [s/m] 
   do j = js,je
      do i = is-1,ie
         dt_x(i,j) = dt * g%IdxCu(i,j)
      enddo
   enddo
   do j = js-1,je
      do i = is,ie
         dt_y(i,j) = dt * g%IdyCv(i,j)
      enddo
   enddo

   do k = 1,nz

      ! --------- dt/dx*h [s]
      do j = js,je
         do i = is-1,ie
            coef_x(i,j) = dt_x(i,j) *2.0*h0(i,j,k)*h0(i+1,j,k) / &
                          ( h0(i,j,k)+h0(i+1,j,k)+h_neglect ) 
         enddo
      enddo
      do j = js-1,je
         do i = is,ie
            coef_y(i,j) = dt_y(i,j) *2.0*h0(i,j,k)*h0(i,j+1,k) / &
                          ( h0(i,j,k)+h0(i,j+1,k)+h_neglect )
         enddo
      enddo

      ! --------- update tracers
      do m = 1,ntr
         ! assign
         fld(:,:) = reg%Tr(m)%t(:,:,k)
         fld0(:,:) = reg0%Tr(m)%t(:,:,k)

         ! --------- min & max values
         !$OMP parallel do private(i,j,ja,jb,ia,ib)
         do j = js, je
            ja = j - 1
            jb = j + 1
            do i = is,ie
               if (g%mask2dT(i,j) /= 0) then
                  ia = i - 1
                  if (g%mask2dT(ia,j)==0) ia = i
                  ib = i + 1
                  if (g%mask2dT(ib,j)==0) ib = i
                  ja = j - 1
                  if (g%mask2dT(i,ja)==0) ja = j
                  jb = j + 1
                  if (g%mask2dT(i,jb)==0) jb = j
                  fmx(i,j) = max( fld(i,j),     &
                      fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
                  fmn(i,j) = min( fld(i,j),     &
                      fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
               endif
            enddo
          enddo
         !$OMP end parallel do

         ! --------- trac grad: dt/dx*h *dc [s*c], using un-adv TRACER!!! 
         ! loop over all sym-grid points; 0 for boundary points !!
         do j = js,je
            do i = is-1,ie
               dcdx(i,j) = g%mask2dCu(i,j)* coef_x(i,j) * &
                  ( fld0(i+1,j) - fld0(i,j) )
            enddo
         enddo
         do j = js-1,je
            do i = is,ie
               dcdy(i,j) = g%mask2dCv(i,j)* coef_y(i,j) * &
                  ( fld0(i,j+1) - fld0(i,j) )
            enddo
         enddo

         ! --------- calc sign(x,y)
         do j = js,je
            dcdy_xmean(j) = 0.
            do i = is,ie
               dcdy_p(i,j) = 0.5*(dcdy(i,j)+dcdy(i,j-1))
               dcdy_xmean(j) = dcdy_xmean(j) + dcdy_p(i,j)
            enddo
            !
            if (dcdy_xmean(j) >= 0.) then
               sign1d(j) = 1.
            elseif (dcdy_xmean(j) < 0.) then
               sign1d(j) = -1.
            endif
         enddo

         ! --------- param EF [cm] ~ delta(c*h)
         dtr_sum = 0.; dtr_abs_sum = 0.
         do j = js,je
            do i = is,ie
               ! ka * del2C [cm] = m2/s * 1/m2 * m*c*s 
               dtr_ka(i,j) = ka(i,j,k) * g%IareaT(i,j) * &
                  ( g%dy_Cu(i,j)*dcdx(i,j) - g%dy_Cu(i-1,j)*dcdx(i-1,j) &
                  + g%dx_Cv(i,j)*dcdy(i,j) - g%dx_Cv(i,j-1)*dcdy(i,j-1) )
               ! -chi*delC [cm] = - chidotn  |h*delC| = [m/s] [c*s]
               hdelc_norm(i,j) = (  (0.5*(dcdx(i,j)+dcdx(i-1,j)))**2. + &
                     (0.5*(dcdy(i,j)+dcdy(i,j-1)))**2.  ) **.5
               dtr_chi(i,j) = - chidotn(i,j,k) * hdelc_norm(i,j) * sign1d(j)
               !
               dtr(i,j) = dtr_ka(i,j) + dtr_chi(i,j)
               
               if(present(forc_ka) .and. present(forc_chi) .and. present(forc_kchi)) then
                  forc_ka(i,j,k,m) = dtr_ka(i,j) / dt
                  forc_chi(i,j,k,m) = dtr_chi(i,j) / dt
                  forc_kchi(i,j,k,m) = dtr(i,j) / dt
               endif

               dtr_sum = dtr_sum + dtr(i,j)
               dtr_abs_sum = dtr_abs_sum + abs(dtr(i,j))
            enddo
         enddo

         ! --------- correct EF --> EF + w*[EF], w = - |EF|/[|EF|]
         do j = js,je
            do i = is,ie
               dtr(i,j) = dtr(i,j) - dtr_sum/dtr_abs_sum * abs(dtr(i,j))
            enddo
         enddo

         ! --------- update tr
         do j = js,je
            do i = is,ie
               q_new = fld(i,j)*h1(i,j,k) + dtr(i,j)
               ! q_new = max( 0.0, q_new)
               q_new = max( fmn(i,j)*h2(i,j,k),       &
                            min(q_new, fmx(i,j)*h2(i,j,k)) )
               reg%Tr(m)%t(i,j,k) = q_new / (h2(i,j,k)+h_neglect) 
               !  Eliminate very small values far from tracer patch
               if (reg%Tr(m)%t(i,j,k) < 1.e-7)  reg%Tr(m)%t(i,j,k) = 0.0
            enddo
         enddo
      enddo !m
   enddo ! k

end subroutine tracer_kachi_chidotn_sign

! ====================================================================
!           Parameterization: feed chi*n (\chi_{\parallel})
! =====================================================================
! update tracer by: 
!   d(ch)/dt = ka*del2C - Chi*  h*delC [m/s*c], 
!            = ka*del2C - chidotn |h*delC|, 
!   where chidotn = |chi*n| is a scalar input
! 
! Note the forcing is calculated using the tracer before adv (Reg0) !
! c2 = ( c1*h1 + dt*EF ) / h2
subroutine tracer_kachi_chidotn_sign2(Reg, Reg0, trsign, dt, ka, chidotn, &
            h0, h1, h2, G, forc_ka, forc_chi, forc_kchi)
   type(tracer_registry_type), intent(inout)    :: reg 
   type(tracer_registry_type), intent(in)       :: reg0 ! 
   type(ocean_grid_type),      intent(in)       :: g      
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke, reg%ntr) &
                         :: trsign ! used to determine sign
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                               intent(in)       :: h1, h2, h0  
   real,                       intent(in)       :: dt   ![s]  
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                               intent(in)       :: ka, chidotn 
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke, reg%ntr), &
                     optional, intent(out)   :: forc_ka, &
                                                forc_chi, &
                                                forc_kchi  

   ! local
   real, parameter :: h_neglect = 1.0e-20 * max(1.0e-10, 1.0e-17)
   integer :: i, j, k, m, is, ie, js, je, nz, ntr &
            , ia, ib, ja, jb
   real, dimension(g%isdb:g%iedb, g%jsd:g%jed)   :: &
               dt_x, coef_x, dcdx
               ! dt times the open face width divided by
               ! the distance between adjacent tracer points [s].
   real, dimension(g%isd:g%ied,   g%jsdb:g%jedb) :: &
               dt_y, coef_y, dcdy, dcdys
   real, dimension(g%isd:g%ied,   g%jsd:g%jed) :: &
               dtr_ka, dtr_chi, dtr         & ! [cm]
              ,fmn, fmx, fld, fld0, fldsign, delc2   &
              , hdelc_norm                  &
              , dcdy_p

   real, dimension(g%jsd:g%jed) :: dcdy_xmean, sign1d

   real :: q_new  ! updated tracer mass at one point [cm]


   is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec ; nz = g%ke
   ntr = reg%ntr

   ! ------------ dt/dx [s/m] 
   do j = js,je
      do i = is-1,ie
         dt_x(i,j) = dt * g%IdxCu(i,j)
      enddo
   enddo
   do j = js-1,je
      do i = is,ie
         dt_y(i,j) = dt * g%IdyCv(i,j)
      enddo
   enddo

   do k = 1,nz

      ! --------- dt/dx*h [s]
      do j = js,je
         do i = is-1,ie
            coef_x(i,j) = dt_x(i,j) *2.0*h0(i,j,k)*h0(i+1,j,k) / &
                          ( h0(i,j,k)+h0(i+1,j,k)+h_neglect ) 
         enddo
      enddo
      do j = js-1,je
         do i = is,ie
            coef_y(i,j) = dt_y(i,j) *2.0*h0(i,j,k)*h0(i,j+1,k) / &
                          ( h0(i,j,k)+h0(i,j+1,k)+h_neglect )
         enddo
      enddo

      ! --------- update tracers
      do m = 1,ntr
         ! assign
         fld(:,:) = reg%Tr(m)%t(:,:,k)
         fld0(:,:) = reg0%Tr(m)%t(:,:,k)
         fldsign(:,:) = trsign(:,:,k,m)

         ! --------- min & max values
         !$OMP parallel do private(i,j,ja,jb,ia,ib)
         do j = js, je
            ja = j - 1
            jb = j + 1
            do i = is,ie
               if (g%mask2dT(i,j) /= 0) then
                  ia = i - 1
                  if (g%mask2dT(ia,j)==0) ia = i
                  ib = i + 1
                  if (g%mask2dT(ib,j)==0) ib = i
                  ja = j - 1
                  if (g%mask2dT(i,ja)==0) ja = j
                  jb = j + 1
                  if (g%mask2dT(i,jb)==0) jb = j
                  fmx(i,j) = max( fld(i,j),     &
                      fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
                  fmn(i,j) = min( fld(i,j),     &
                      fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
               endif
            enddo
          enddo
         !$OMP end parallel do

         ! --------- trac grad: dt/dx*h *dc [s*c], using un-adv TRACER!!! 
         ! loop over all sym-grid points; 0 for boundary points !!
         do j = js,je
            do i = is-1,ie
               dcdx(i,j) = g%mask2dCu(i,j)* coef_x(i,j) * &
                  ( fld0(i+1,j) - fld0(i,j) )
            enddo
         enddo
         do j = js-1,je
            do i = is,ie
               dcdy(i,j) = g%mask2dCv(i,j)* coef_y(i,j) * &
                  ( fld0(i,j+1) - fld0(i,j) )
               dcdys(i,j) = g%mask2dCv(i,j)* coef_y(i,j) * &
                  ( fldsign(i,j+1) - fldsign(i,j) )
            enddo
         enddo

         ! --------- calc sign(x,y)
         do j = js,je
            dcdy_xmean(j) = 0.
            do i = is,ie
               dcdy_p(i,j) = 0.5*(dcdys(i,j)+dcdys(i,j-1))
               dcdy_xmean(j) = dcdy_xmean(j) + dcdy_p(i,j)
            enddo
            !
            if (dcdy_xmean(j) >= 0.) then
               sign1d(j) = 1.
            elseif (dcdy_xmean(j) < 0.) then
               sign1d(j) = -1.
            endif
         enddo

         ! --------- param EF [cm] ~ delta(c*h)
         do j = js,je
            do i = is,ie
               ! ka * del2C [cm] = m2/s * 1/m2 * m*c*s 
               dtr_ka(i,j) = ka(i,j,k) * g%IareaT(i,j) * &
                  ( g%dy_Cu(i,j)*dcdx(i,j) - g%dy_Cu(i-1,j)*dcdx(i-1,j) &
                  + g%dx_Cv(i,j)*dcdy(i,j) - g%dx_Cv(i,j-1)*dcdy(i,j-1) )
               ! -chi*delC [cm] = - chidotn  |h*delC| = [m/s] [c*s]
               hdelc_norm(i,j) = (  (0.5*(dcdx(i,j)+dcdx(i-1,j)))**2. + &
                     (0.5*(dcdy(i,j)+dcdy(i,j-1)))**2.  ) **.5
               dtr_chi(i,j) = - chidotn(i,j,k) * hdelc_norm(i,j) * sign1d(j)
               !
               dtr(i,j) = dtr_ka(i,j) + dtr_chi(i,j)
               
               if(present(forc_ka) .and. present(forc_chi) .and. present(forc_kchi)) then
                  forc_ka(i,j,k,m) = dtr_ka(i,j) / dt
                  forc_chi(i,j,k,m) = dtr_chi(i,j) / dt
                  forc_kchi(i,j,k,m) = dtr(i,j) / dt
               endif
            enddo
         enddo

        ! --------- update tr
         do j = js,je
            do i = is,ie
               q_new = fld(i,j)*h1(i,j,k) + dtr(i,j)
               ! q_new = max( 0.0, q_new)
               q_new = max( fmn(i,j)*h2(i,j,k),       &
                            min(q_new, fmx(i,j)*h2(i,j,k)) )
               reg%Tr(m)%t(i,j,k) = q_new / (h2(i,j,k)+h_neglect) 
               !  Eliminate very small values far from tracer patch
               if (reg%Tr(m)%t(i,j,k) < 1.e-7)  reg%Tr(m)%t(i,j,k) = 0.0
            enddo
         enddo
      enddo !m
   enddo ! k

end subroutine tracer_kachi_chidotn_sign2


! ====================================================================
!           Parameterization: k-Chi model
! =====================================================================
! update tracer ('Reg') by the parameterized forc: 
!   d(ch)/dt = ka*del2C - Chi* h*delC [m/s*c], 
!    where ka [m2/s], Chi [m/s], 
!   Chi = - alpha * u, where gamma ~ [m/s/c]
! -->
!   -Chi* h*delC = alpha* u*hdelC= alpha* (divUc - cdivU)

! Note the forcing is calculated using the tracer before adv (Reg0) !
! c2 = ( c1*h1 + dt*EF ) / h2
! 
subroutine tracer_kachi_alpha(Reg, Reg0, dt, ka, alp, uh, vh, h0, h1, h2 &
            , G, forc_ka, forc_chi, forc_kchi)
   type(tracer_registry_type), intent(inout)    :: reg 
   type(tracer_registry_type), intent(in)       :: reg0 ! 
   type(ocean_grid_type),      intent(in)       :: g      
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                               intent(in)       :: h1, h2, h0  
   real,                       intent(in)       :: dt   ![s]  
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                               intent(in)       :: ka, alp 
   real, dimension(g%isdb:g%iedb, g%jsd:g%jed, g%ke), &
                               intent(in)    :: uh      ! uhL 
   real, dimension(g%isd:g%ied, g%jsdb:g%jedb, g%ke), &
                               intent(in)    :: vh      ! vhL  
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke, reg%ntr), &
                     optional, intent(out)   :: forc_ka, &
                                                forc_chi, &
                                                forc_kchi  

! local
   real, parameter :: h_neglect = 1.0e-20 * max(1.0e-10, 1.0e-17)
   integer :: i, j, k, m, is, ie, js, je, nz, ntr &
            , ia, ib, ja, jb
   real, dimension(g%isdb:g%iedb, g%jsd:g%jed)   :: &
               dt_x, coef_x, dcdx
               ! dt times the open face width divided by
               ! the distance between adjacent tracer points [s].
   real, dimension(g%isd:g%ied,   g%jsdb:g%jedb) :: &
               dt_y, coef_y, dcdy
   real, dimension(g%isd:g%ied,   g%jsd:g%jed) :: &
               dtr_ka, dtr_chi, dtr         & ! [cm]
              ,fmn, fmx, fld, fld0, delc2
   real :: q_new  ! updated tracer mass at one point [cm]
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke, reg%ntr) :: &
                  divuc, cdivu, udelc


   is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec ; nz = g%ke
   ntr = reg%ntr

   ! ------------ u*delC [m*c], uh already timed by dt
    call calc_divforc(Reg0, uh, vh, divuc, cdivu, udelc, g)

   ! ------------ dt/dx [s/m] 
   do j = js,je
      do i = is-1,ie
         dt_x(i,j) = dt * g%IdxCu(i,j)
      enddo
   enddo
   do j = js-1,je
      do i = is,ie
         dt_y(i,j) = dt * g%IdyCv(i,j)
      enddo
   enddo

   do k = 1,nz

      ! --------- dt/dx*h [s]
      do j = js,je
         do i = is-1,ie
            coef_x(i,j) = dt_x(i,j) *2.0*h0(i,j,k)*h0(i+1,j,k) / &
                          ( h0(i,j,k)+h0(i+1,j,k)+h_neglect ) 
         enddo
      enddo
      do j = js-1,je
         do i = is,ie
            coef_y(i,j) = dt_y(i,j) *2.0*h0(i,j,k)*h0(i,j+1,k) / &
                          ( h0(i,j,k)+h0(i,j+1,k)+h_neglect )
         enddo
      enddo

      ! --------- update tracers
      do m = 1,ntr
         ! assign
         fld(:,:) = reg%Tr(m)%t(:,:,k)
         fld0(:,:) = reg0%Tr(m)%t(:,:,k)

         ! --------- min & max values
         !$OMP parallel do private(i,j,ja,jb,ia,ib)
         do j = js, je
            ja = j - 1
            jb = j + 1
            do i = is,ie
               if (g%mask2dT(i,j) /= 0) then
                  ia = i - 1
                  if (g%mask2dT(ia,j)==0) ia = i
                  ib = i + 1
                  if (g%mask2dT(ib,j)==0) ib = i
                  ja = j - 1
                  if (g%mask2dT(i,ja)==0) ja = j
                  jb = j + 1
                  if (g%mask2dT(i,jb)==0) jb = j
                  fmx(i,j) = max( fld(i,j),     &
                      fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
                  fmn(i,j) = min( fld(i,j),     &
                      fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
               endif
            enddo
          enddo
         !$OMP end parallel do

         ! --------- trac grad: dt/dx*h *dc [s*c], using un-adv TRACER!!! 
         ! loop over all sym-grid points; 0 for boundary points !!
         do j = js,je
            do i = is-1,ie
               dcdx(i,j) = g%mask2dCu(i,j)* coef_x(i,j) * &
                  ( fld0(i+1,j) - fld0(i,j) )
            enddo
         enddo
         do j = js-1,je
            do i = is,ie
               dcdy(i,j) = g%mask2dCv(i,j)* coef_y(i,j) * &
                  ( fld0(i,j+1) - fld0(i,j) )
            enddo
         enddo

        ! --------- param EF [cm] ~ delta(c*h)
         do j = js,je
            do i = is,ie
               ! ka * del2C [cm] = m2/s * 1/m2 * m*c*s 
               dtr_ka(i,j) = ka(i,j,k) * g%IareaT(i,j) * &
                  ( g%dy_Cu(i,j)*dcdx(i,j) - g%dy_Cu(i-1,j)*dcdx(i-1,j) &
                  + g%dx_Cv(i,j)*dcdy(i,j) - g%dx_Cv(i,j-1)*dcdy(i,j-1) )
               ! -chi*delC [cm] = alpha* u*delC
               dtr_chi(i,j) = alpha(i,j,k) * udelc(i,j,k,m)
               dtr(i,j) = dtr_ka(i,j) + dtr_chi(i,j)
               
               if(present(forc_ka) .and. present(forc_chi) .and. present(forc_kchi)) then
                  forc_ka(i,j,k,m) = dtr_ka(i,j) / dt
                  forc_chi(i,j,k,m) = dtr_chi(i,j) / dt
                  forc_kchi(i,j,k,m) = dtr(i,j) / dt
               endif
            enddo
         enddo

        ! --------- update tr
         do j = js,je
            do i = is,ie
               q_new = fld(i,j)*h1(i,j,k) + dtr(i,j)
               ! q_new = max( 0.0, q_new)
               q_new = max( fmn(i,j)*h2(i,j,k),       &
                            min(q_new, fmx(i,j)*h2(i,j,k)) )
               reg%Tr(m)%t(i,j,k) = q_new / (h2(i,j,k)+h_neglect) 
               !  Eliminate very small values far from tracer patch
               if (reg%Tr(m)%t(i,j,k) < 1.e-7)  reg%Tr(m)%t(i,j,k) = 0.0
            enddo
         enddo
      enddo !m
   enddo ! k

end subroutine tracer_kachi_alpha

! ====================================================================
!           Parameterization: k-Chi model
! =====================================================================
! update tracer ('Reg') by the parameterized forc: 
!   d(ch)/dt = ka*del2C - Chi* h*delC [m/s*c], 
!    where ka [m2/s], Chi [m/s], 
!   Chi = gamma * h*delC, where gamma ~ [m/s/c]
! -->
!   -Chi* h*delC = -gamma* (h*delC)^2

! Note the forcing is calculated using the tracer before adv (Reg0) !
! c2 = ( c1*h1 + dt*EF ) / h2
! 
subroutine tracer_kachi_gamma(Reg, Reg0, dt, ka, gam, h0, h1, h2, G &
            , forc_ka, forc_chi, forc_kchi)
   type(tracer_registry_type), intent(inout)    :: reg 
   type(tracer_registry_type), intent(in)       :: reg0 ! 
   type(ocean_grid_type),      intent(in)       :: g      
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                               intent(in)       :: h1, h2, h0  
   real,                       intent(in)       :: dt   ![s]  
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                               intent(in)       :: ka, gam   
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke, reg%ntr), &
                     optional, intent(out)   :: forc_ka, &
                                                forc_chi, &
                                                forc_kchi ! forc [m/s*c] 

! local
   real, parameter :: h_neglect = 1.0e-20 * max(1.0e-10, 1.0e-17)
   integer :: i, j, k, m, is, ie, js, je, nz, ntr &
            , ia, ib, ja, jb
   real, dimension(g%isdb:g%iedb, g%jsd:g%jed)   :: &
               dt_x, coef_x, dcdx
               ! dt times the open face width divided by
               ! the distance between adjacent tracer points [s].
   real, dimension(g%isd:g%ied,   g%jsdb:g%jedb) :: &
               dt_y, coef_y, dcdy
   real, dimension(g%isd:g%ied,   g%jsd:g%jed) :: &
               dtr_ka, dtr_chi, dtr         & ! [cm]
              ,fmn, fmx, fld, fld0, delc2
   real :: q_new  ! updated tracer mass at one point [cm]

   is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec ; nz = g%ke
   ntr = reg%ntr

   ! ------------ dt/dx [s/m] 
   do j = js,je
      do i = is-1,ie
         dt_x(i,j) = dt * g%IdxCu(i,j)
      enddo
   enddo
   do j = js-1,je
      do i = is,ie
         dt_y(i,j) = dt * g%IdyCv(i,j)
      enddo
   enddo

   do k = 1,nz

      ! --------- dt/dx*h [s]
      do j = js,je
         do i = is-1,ie
            coef_x(i,j) = dt_x(i,j) *2.0*h0(i,j,k)*h0(i+1,j,k) / &
                          ( h0(i,j,k)+h0(i+1,j,k)+h_neglect ) 
         enddo
      enddo
      do j = js-1,je
         do i = is,ie
            coef_y(i,j) = dt_y(i,j) *2.0*h0(i,j,k)*h0(i,j+1,k) / &
                          ( h0(i,j,k)+h0(i,j+1,k)+h_neglect )
         enddo
      enddo

      ! --------- update tracers
      do m = 1,ntr
         ! assign
         fld(:,:) = reg%Tr(m)%t(:,:,k)
         fld0(:,:) = reg0%Tr(m)%t(:,:,k)

         ! --------- min & max values
         !$OMP parallel do private(i,j,ja,jb,ia,ib)
         do j = js, je
            ja = j - 1
            jb = j + 1
            do i = is,ie
               if (g%mask2dT(i,j) /= 0) then
                  ia = i - 1
                  if (g%mask2dT(ia,j)==0) ia = i
                  ib = i + 1
                  if (g%mask2dT(ib,j)==0) ib = i
                  ja = j - 1
                  if (g%mask2dT(i,ja)==0) ja = j
                  jb = j + 1
                  if (g%mask2dT(i,jb)==0) jb = j
                  fmx(i,j) = max( fld(i,j),     &
                      fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
                  fmn(i,j) = min( fld(i,j),     &
                      fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
               endif
            enddo
          enddo
         !$OMP end parallel do

         ! --------- trac grad: dt/dx*h *dc [s*c], using un-adv TRACER!!! 
         ! loop over all sym-grid points; 0 for boundary points !!
         do j = js,je
            do i = is-1,ie
               dcdx(i,j) = g%mask2dCu(i,j)* coef_x(i,j) * &
                  ( fld0(i+1,j) - fld0(i,j) )
            enddo
         enddo
         do j = js-1,je
            do i = is,ie
               dcdy(i,j) = g%mask2dCv(i,j)* coef_y(i,j) * &
                  ( fld0(i,j+1) - fld0(i,j) )
            enddo
         enddo

        ! --------- param EF [cm] ~ delta(c*h)
         do j = js,je
            do i = is,ie
               ! ka * del2C [cm] = m2/s * 1/m2 * m*c*s 
               dtr_ka(i,j) = ka(i,j,k) * g%IareaT(i,j) * &
                  ( g%dy_Cu(i,j)*dcdx(i,j) - g%dy_Cu(i-1,j)*dcdx(i-1,j) &
                  + g%dx_Cv(i,j)*dcdy(i,j) - g%dx_Cv(i,j-1)*dcdy(i,j-1) )
               ! -chi*delC [cm] = -gamma* (h*delC)^2
               delc2(i,j) = ( 0.5*(dcdx(i,j)+dcdx(i-1,j)) )**2. / dt + &
                  ( 0.5*(dcdy(i,j)+dcdy(i,j-1)) )**2. / dt
               dtr_chi(i,j) = - gam(i,j,k) * delc2(i,j)
               dtr(i,j) = dtr_ka(i,j) + dtr_chi(i,j)
               
               if(present(forc_ka) .and. present(forc_chi) .and. present(forc_kchi)) then
                  forc_ka(i,j,k,m) = dtr_ka(i,j) / dt
                  forc_chi(i,j,k,m) = dtr_chi(i,j) / dt
                  forc_kchi(i,j,k,m) = dtr(i,j) / dt
               endif
            enddo
         enddo

         if (k == 1) then
           !print*, 'Tr:', fld0(22,44:48)
           !print*, 'params:', ka(22,44:48,k), chiu(22,44:48,k), chiv(22,44:48,k)
           !print*, 'EF kchi:', forc_kchi(22,44:48,1,1)
           !print*, 'dcdx:', dcdx(22,44:48) / dt
           !print*, 'dcdy:', dcdy(22,44:48) / dt
         endif
         ! print*, "max dtr_ka(1)", maxval(dtr_ka(is:ie,js:je)), minval(dtr_ka(is:ie,js:je))
         ! print*, "max dtr_chi(1)", maxval(dtr_chi(is:ie,js:je)), minval(dtr_chi(is:ie,js:je))

        ! --------- update tr
         do j = js,je
            do i = is,ie
               q_new = fld(i,j)*h1(i,j,k) + dtr(i,j)
               ! q_new = max( 0.0, q_new)
               q_new = max( fmn(i,j)*h2(i,j,k),       &
                            min(q_new, fmx(i,j)*h2(i,j,k)) )
               reg%Tr(m)%t(i,j,k) = q_new / (h2(i,j,k)+h_neglect) 
               !  Eliminate very small values far from tracer patch
               if (reg%Tr(m)%t(i,j,k) < 1.e-7)  reg%Tr(m)%t(i,j,k) = 0.0
            enddo
         enddo
      enddo !m

   enddo ! k

end subroutine tracer_kachi_gamma

! ====================================================================
!           SUBTOUTINE: calc total amount of c-mass of snapshot
! =====================================================================
subroutine tracer_mass(cmass, reg, h, G)
   type(ocean_grid_type),      intent(in)       :: g    
   type(tracer_registry_type), intent(in)       :: reg     !< registered tracers
   real, dimension(g%isd:g%ied, g%jsd:g%jed, g%ke), &
                              intent(in)       :: h   
   real, dimension(g%ke, reg%ntr),   intent(out):: cmass   ! c*m3, nk-by-ntr

   ! local
   integer :: i, j, k, is, ie, js, je, nz, m 
   is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec ; nz = g%ke

   cmass = 0.
   do m = 1, reg%ntr
      do k = 1, nz
       do j = js,je
         do i = is,ie
            if (g%mask2dT(i,j) /= 0.) cmass(k,m) = cmass(k,m) + &
               reg%tr(m)%t(i,j,k) * g%areaT(i,j) * h(i,j,k) 
         enddo               
       enddo                  
      enddo  !k
   enddo ! ntr

end subroutine tracer_mass

! ====================================================================
!           SUBTOUTINE: update tracer at each point by a ratio !
! =====================================================================
subroutine corr_tracer(ratio,reg,G)
   type(ocean_grid_type),       intent(in)       :: g    
   type(tracer_registry_type),  intent(inout)    :: reg     
   real, dimension(g%ke, reg%ntr),    intent(in) :: ratio  ! r = prev/after 

   ! local
   integer :: i, j, k, is, ie, js, je, nz, m
   is = g%isc ; ie = g%iec ; js = g%jsc ; je = g%jec ; nz = g%ke

   do m = 1,reg%ntr
      do k = 1, nz
         do j = js,je
            do i = is,ie
               if (g%mask2dT(i,j) /= 0.) reg%tr(m)%t(i,j,k) = &
                  reg%tr(m)%t(i,j,k)*ratio(k,m)
            enddo               
         enddo                  
      enddo  !k
   enddo !ntr

end subroutine corr_tracer

! c=====================================================================
! c          SUBTOUTINE: 3-d version of advfct.f
! c=====================================================================
subroutine fct3d(time_interval, tracer_Reg, G, fco, fc, uflx, vflx)
! c
! c --- fully 3-d version of advfct.f
! c     
! c     Main input vars:
! c       iord    - 1: donnor cell; 2: complete with antidiffusive fluxes (default)
! c       fco/fc  - layer thicknesses at previous and new time step
! c       u/vflx  - isopycnal mass fluxes [m3] (times delta T)
! c       fld     - 'tracer', transported mixing ratio, like sal or temp
! c
! c     Some inmportant vars in the adv:
! c       bforej/afterj - (jdm), amount of tracermass in each lat band in 
! c                       current layer at previous/new time steps
! c       bfore/after   - amount of tracermass in the whole domain at
! c                       previous/new time steps
! c
real,                  intent(in)        :: time_interval !< Offline transport time interval
type(ocean_grid_type), intent(in)  :: g  
type(tracer_registry_type), intent(inout)          :: tracer_Reg            !< Control structure for offline module
real, dimension(g%isd:g%ied, g%jsd:g%jed,   g%ke), intent(in)    :: fco !< h at t1 before adv (from onl model)
real, dimension(g%isd:g%ied, g%jsd:g%jed,   g%ke), intent(in)    :: fc !< h at t2 after adv (from onl model)
real, dimension(g%isdb:g%iedb, g%jsd:g%jed, g%ke), intent(in) :: uflx  !< (uhtr) Zonal mass transport [m3]
real, dimension(g%isd:g%ied, g%jsdb:g%jedb, g%ke), intent(in) :: vflx  !< Meridional mass transport

integer :: iord
real, dimension(g%isd:g%ied, g%jsd:g%jed,   g%ke)   :: fld
real, dimension(g%isd:g%ied, g%jsd:g%jed)   ::  flxdiv, fmx, fmn, flp,fln
real, dimension(g%isdb:g%iedb, g%jsd:g%jed) :: flx, uan
real, dimension(g%isd:g%ied, g%jsdb:g%jedb) :: fly, van

! real, dimension(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) :: flp,fln,flx,fly,uan,van,flxdiv,fmx,fmn
! real, dimension(1-nbdy:jj+nbdy) :: clipj, vlumjs

real, dimension(g%jsd:g%jed) :: clipj, vlumj
real, dimension(g%ke)  :: bfore, after, ratio
real a(g%ke),b(g%ke),c(g%ke),athird,dx,fcdx,yl,yr
real onemu,clip,vlume,amount,slab,dslab,thkchg, &
     fluxdv,epsil, q
integer ip1,im1,jp1,jm1,kp,jaa,margin,i1,i2,j1,j2, ia,ib,ja,jb,l
logical wrap,recovr
data (athird = 1./3.)
real,    parameter :: ft14 = 7.0/6.0,  &!4th centered inner coeff
                      ft24 =-1.0/6.0   !4th centered outer coeff

integer :: i, j, k, m, is, ie, js, je, isd, ied, jsd, jed, nz, ntr
integer :: isdb, iedb, jsdb, jedb

! c --- if iord=1, scheme reduces to simple donor cell scheme. usually it
! c     is set to 2 in our model
parameter (epsil = 1.e-11, onemu = 1.e-11)
iord = 2
recovr = .false.

is  = g%isc ; ie  = g%iec ; js  = g%jsc ; je  = g%jec ; nz = g%ke
isd = g%isd ; ied = g%ied ; jsd = g%jsd ; jed = g%jed
isdb = g%IsdB ; iedb = g%IedB ; jsdb = g%JsdB ; jedb = g%JedB
ntr = tracer_Reg%ntr

! c ---------------------------------------------------------------------
! c ---  ADV
! c ---------------------------------------------------------------------
bfore = 0.
after = 0.

m = 1 ! tracers
fld = tracer_Reg%Tr(m)%t

DO k = 1, nz
! c     Tracer mass at each lay
   do j = js,je
      do i = is,ie
         if (g%mask2dT(i,j) /= 0.) bfore(k) = bfore(k) + &
            fld(i,j,k) * fco(i,j,k) * g%areaT(i,j)
      enddo              
   enddo                  
! c
! c ---    Compute antidiffusive tracer fluxes (high- minus low-order)
! c ---    at x, y and z directions [m3*c]
! c   -    flx/fly - low-order c-flux (uflx*c) [m3*c], used for 'flxdiv'
! c   -    q       - 4th order tracer [c]
! c   -    uan/van - high-order flux minus low-order flux
! c   -    fmx/fmn - max/min values among tracers in 5 cells and the 
! c                  adjacent inward vertfx/diaflx (both positive)
! c
   do j = js-1,je+1
      ja = j - 1
      jb = j + 1
      do i = is-1,ie+1
         ! c-flux in x dir
         if (g%mask2dCu(i,j).ne.0.) then
            if (uflx(i,j,k).ge.0.) then
               q = fld(i,j,k)
            else
               q = fld(i+1,j,k)
            endif
            flx(i,j) = uflx(i,j,k) * q
            q = fld(i+1,j,k) + fld(i,j,k) !  2nd order
            if (g%mask2dT(i+2,j)+g%mask2dCu(i-1,j).eq.2.) then !  4th order
               q = 1.125*q - .125*(fld(i+2,j,k) + fld(i-1,j,k)) 
            endif
            uan(i,j) = .5*q*uflx(i,j,k) - flx(i,j)
         endif    
         ! c-flux in y dir
         if (g%mask2dCv(i,j).ne.0) then
            if (vflx(i,j,k).ge.0.) then
               q = fld(i,j,k)
            else
               q = fld(i,j+1,k)
            endif
            fly(i,j) = vflx(i,j,k) * q
            q = fld(i,j+1,k) + fld(i,j,k) !  2nd order
            if (g%mask2dT(i,j+2 )+g%mask2dCv(i,j-1).eq.2.) then !  4th order
               q = 1.125*q - .125*(fld(i,j+2,k) + fld(i,j-1,k)) 
            endif
            van(i,j) = .5*q*vflx(i,j,k) - fly(i,j)
         endif           
         ! c-flux 
         if (g%mask2dT(i,j).ne.0.) then
            ia = i - 1
            if (g%mask2dT(ia,j).eq.0.) ia = i
            ib = i + 1
            if (g%mask2dT(ib,j).eq.0.) ib = i
            ja = j - 1
            if (g%mask2dT(i,ja).eq.0.) ja = j
            jb = j + 1
            if (g%mask2dT(i,jb).eq.0.) jb = j
            fmx(i,j) = max( fld(i,j,k), &
                    fld(ia,j,k),fld(ib,j,k),fld(i,ja,k),fld(i,jb,k) )
            fmn(i,j) = min( fld(i,j,k), &
                    fld(ia,j,k),fld(ib,j,k),fld(i,ja,k),fld(i,jb,k) )
         endif 
      enddo  !i
! c
! c ---       Boundaries for fluxes, u and v
! c

   enddo      !j

! c
! c ---    UPDATE the tracer in current layer with low-order c-flux
! c
   do j = js,je
      if (recovr) vlumj(j) = 0. 
      if (recovr) clipj(j) = 0.
      do i = is,ie
         if (g%mask2dT(i,j).ne.0.) then
            flxdiv(i,j) = (flx(i,j)-flx(i-1,j)+fly(i,j)-fly(i,j-1)) &
                *g%IareaT(i,j)
            ! updated tracer mass
            ! (ch)_2 = (ch)_1 - div_Fc 
            q = fld(i,j,k)*fco(i,j,k) - flxdiv(i,j) 
            ! THIS IS TO KEEP TRACER WITHIN REASONABLE VALUES !
            ! q(q<fmn)=fmn, q(q>fmx)=fmx, q(fmn<q<fmx)=q
            amount = max( fmn(i,j)*fc(i,j,k),  &
                           min(q, fmx(i,j)*fc(i,j,k)) )
            if (recovr) then
               vlumj(j) = vlumj(j) + g%areaT(i,j)*fc(i,j,k)
               clipj(j) = clipj(j) + (q-amount)*g%areaT(i,j)
            endif
            ! amount = q
            ! updated tracer conc, essentially q/h, corrected to 
            ! avoid negative values
            ! fld(i,j,k) = amount / (onemu+fc(i,j,k))
            fld(i,j,k) = (fld(i,j,k)*onemu+amount) / (onemu+fc(i,j,k))
         endif            !ip
      enddo               !i
   enddo                  !j
! c
! c --- at each grid point, determine the ratio of the largest permissible
! c --- pos. (neg.) change in -fld- to the sum of all incoming (outgoing) fluxes
! c
   do j = js,je
      do i = is,ie
         if (g%mask2dT(i,j).ne.0.) then
            flp(i,j) = (fmx(i,j)-fld(i,j,k)) * fc(i,j,k)  &
                 / ( (max(0.,uan(i-1,j))-min(0.,uan(i,j)) &
                     + max(0.,van(i,j-1))-min(0.,van(i,j))+epsil)*g%IareaT(i,j) ) 
! 
            fln(i,j) = (fmn(i,j)-fld(i,j,k))*fc(i,j,k)   &
                 / ( (min(0.,uan(i-1,j))-max(0.,uan(i,j)) &
                     + min(0.,van(i,j-1))-max(0.,van(i,j))-epsil)*g%IareaT(i,j) )
         endif       
      enddo     
   enddo     
! c
! c----    Limit antidiffusive fluxes ("clipping") by multiplying 'clip'
! c
   do j = js,je
      do i = is,ie
         if (g%mask2dCu(i,j).ne.0.) then

            if (uan(i,j).ge.0.) then
               clip = min(1.,flp(i,j),fln(i-1,j))
            else
               clip = min(1.,fln(i,j),flp(i-1,j))
            endif
            flx(i,j) = uan(i,j)*clip
         endif           !iu
         if (g%mask2dCv(i,j).ne.0) then
            if (van(i,j).ge.0.) then
               clip = min(1.,flp(i,j),fln(i,j-1))
            else
               clip = min(1.,fln(i,j),flp(i,j-1))
            endif
            fly(i,j) = van(i,j)*clip
         endif           !iv
      enddo               !i
   enddo                  !j
! c
! c ---    UPDATE (for the 2nd time) with 'clipped' low-order flux
! c
   do j = js,je
      do i = is,ie
         if (g%mask2dT(i,j).ne.0.) then
            flxdiv(i,j) = (flx(i,j)-flx(i-1,j)+fly(i,j)-fly(i,j-1)) &
                *g%IareaT(i,j)
            q = fld(i,j,k)*fc(i,j,k) - flxdiv(i,j)
            amount = max( fmn(i,j)*fc(i,j,k), min(q,fmx(i,j)*fc(i,j,k)) )
            ! fld(i,j,k) = (fld(i,j,k)*onemu+amount) / (onemu+fc(i,j,k))
            if (recovr) clipj(j) = clipj(j)+(q-amount)*g%areaT(i,j)
         endif           !ip
      enddo               !i
   enddo                  !j


! c
! c ---    Recover (if needed) 'clipped' amount and return to field layer by layer
! c
   if (recovr) then
      vlume = 0.
      clip = 0.
      do j = js,je
         vlume = vlume + vlumj(j)
         clip = clip + clipj(j)
      enddo

      if (vlume.ne.0.) then
         clip = clip / vlume
         write(6,'(i2,a,1pe11.3)') k,'  fld drift in fct3d',-clip
         do j = js,je
            do i = is,ie
               if (g%mask2dT(i,j).ne.0) fld(i,j,k) = fld(i,j,k) + clip
            enddo         !i
         enddo            !j
      endif
   endif                 !recover

! c  assign &   Tracer mass at each lay
   do j = js,je
      do i = is,ie
         if (g%mask2dT(i,j) /= 0.) after(k) = after(k) + &
            fld(i,j,k) * fco(i,j,k) * g%areaT(i,j)

         ! assign
         tracer_Reg%Tr(m)%t(i,j,k) = fld(i,j,k)
      enddo              
   enddo                    

ENDDO ! k
! c
! c      if (bfore.ne.0.)
! c     . write (lp,'(a,1p,3e14.6,e11.1)') 'fct3d conservation:',
! c     .  bfore,after,after-bfore,(after-bfore)/bfore
if (after(1).ne.0.) ratio = bfore / after
print*, '-----fct3d: bfore= ',bfore, NEW_LINE('A')// &
         'after= ', after, NEW_LINE('A')// &
         'Err= ',(after-bfore)/bfore
! c
do k = 1,nz
   do j = js,je
      do i = is,ie
         if (g%mask2dT(i,j) /= 0.) tracer_Reg%Tr(m)%t(i,j,k) = &
               tracer_Reg%Tr(m)%t(i,j,k) * ratio(k)
      enddo               !i
   enddo                  !k
enddo                     !j

end subroutine fct3d



end module mom_tracer_aux
