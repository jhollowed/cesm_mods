module held_suarez
  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: Implement idealized Held-Suarez forcings
  !    Held, I. M., and M. J. Suarez, 1994: 'A proposal for the
  !    intercomparison of the dynamical cores of atmospheric general
  !    circulation models.'
  !    Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
  ! 
  !-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8

  ! JH, for writing to log file, not sure if correct way to do this
  use cam_logfile,      only : iulog

  implicit none
  private
  save

  public :: held_suarez_1994_init
  public :: held_suarez_1994

  !!
  !! Forcing parameters
  !!
  real(r8), parameter :: efoldf  =  1._r8  ! efolding time for wind dissipation
  real(r8), parameter :: efolda  = 40._r8  ! efolding time for T dissipation
  real(r8), parameter :: efolds  =  4._r8  ! efolding time for T dissipation
  real(r8), parameter :: sigmab  =  0.7_r8 ! threshold sigma level
  real(r8), parameter :: t00     = 200._r8 ! minimum reference temperature
  real(r8), parameter :: kf      = 1._r8/(86400._r8*efoldf) ! 1./efolding_time for wind dissipation

  real(r8), parameter :: onemsig = 1._r8 - sigmab ! 1. - sigma_reference

  real(r8), parameter :: ka      = 1._r8/(86400._r8 * efolda) ! 1./efolding_time for temperature diss.
  real(r8), parameter :: ks      = 1._r8/(86400._r8 * efolds)

  !!
  !! Model constants, reset in init call
  !!
  real(r8)              :: cappa = 2.0_r8 / 7.0_r8  ! R/Cp
  real(r8)              :: cpair = 1004.0_r8        ! specific heat of dry air (J/K/kg)
! CJ rair, gravit added
  real(r8)              :: rair  = 286.85714_r8     ! gas constant (J/K/kg)
  real(r8)              :: gravit = 9.80665_r8      ! gravity (m/s/s)
  real(r8)              :: psurf_ref = 0.0_r8       ! Surface pressure
  ! pref_mid_norm are layer midpoints normalized by surface pressure ('eta' coordinate)
  real(r8), allocatable :: pref_mid_norm(:)
  integer               :: pver                     ! Num vertical levels



!======================================================================= 
contains
!======================================================================= 

! CJ rair, gravit added
  subroutine held_suarez_1994_init(cappa_in, cpair_in, rair_in, gravit_in, psurf_ref_in, pref_mid_norm_in)
    !! Dummy arguments
    real(r8), intent(in) :: cappa_in
    real(r8), intent(in) :: cpair_in
    real(r8), intent(in) :: rair_in
    real(r8), intent(in) :: gravit_in
    real(r8), intent(in) :: psurf_ref_in
    real(r8), intent(in) :: pref_mid_norm_in(:)

    pver = size(pref_mid_norm_in)
    allocate(pref_mid_norm(pver))
    cappa         = cappa_in
    cpair         = cpair_in
    rair          = rair_in
    gravit        = gravit_in
    psurf_ref     = psurf_ref_in
    pref_mid_norm = pref_mid_norm_in

  end subroutine held_suarez_1994_init

  subroutine held_suarez_1994(pcols, ncol, clat, ztodt, pmid, pint, rpdel,&
       u, v, t, du, dv, s)

    !
    ! Input arguments
    !
    integer,  intent(in)  :: pcols            ! Size of column dimension
    integer,  intent(in)  :: ncol             ! Num active columns
    real(r8), intent(in)  :: clat(pcols)      ! latitudes(radians) for columns
    real(r8), intent(in)  :: ztodt            ! time step
    real(r8), intent(in)  :: pmid(pcols,pver) ! mid-point pressure
    real(r8), intent(in)  :: pint(pcols,pver+1) ! interface pressure
    real(r8), intent(in)  :: rpdel(pcols,pver) ! 1/pdel
    real(r8), intent(in)  :: u(pcols,pver)    ! Zonal wind (m/s)
    real(r8), intent(in)  :: v(pcols,pver)    ! Meridional wind (m/s)
    real(r8), intent(in)  :: t(pcols,pver)    ! Temperature (K)
                                              !
                                              ! Output arguments
                                              !
    real(r8), intent(out) :: du(pcols,pver)   ! Zonal wind tend
    real(r8), intent(out) :: dv(pcols,pver)   ! Meridional wind tend
    real(r8), intent(out) :: s(pcols,pver)    ! Heating rate
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i, k          ! Longitude, level indices

    real(r8) :: kv            ! 1./efolding_time (normalized) for wind
    real(r8) :: kt            ! 1./efolding_time for temperature diss.
    real(r8) :: trefa         ! "radiative equilibrium" T
    real(r8) :: trefc         ! used in calc of "radiative equilibrium" T
    real(r8) :: cossq(ncol)   ! coslat**2
    real(r8) :: cossqsq(ncol) ! coslat**4
    real(r8) :: sinsq(ncol)   ! sinlat**2
    real(r8) :: coslat(ncol)  ! cosine(latitude)

!   CJ add for vertical diffusion
    real(r8), parameter :: p_vdiff_level = 100._r8, &    ! onset of the vertical diffusion 1 hPa
                           nu_v          = 1._r8         ! vertical diffusion coefficient 1 m^2/2     
!                          nu_v          = 0.5_r8        ! vertical diffusion coefficient 0.5 m^2/2     
    real(r8) :: Km(pcols,pver+1)                         ! Eddy diffusivity for vertical diffsuion
    real(r8) :: rho                                      ! density
    real(r8) :: pT                                       ! JH: pressure at top interface
    real(r8) :: CAm(pcols,pver)                          ! Matrix Coefficents for PBL Scheme
    real(r8) :: CCm(pcols,pver)                          ! Matrix Coefficents for PBL Scheme
    real(r8) :: CEm(pcols,pver+1)                        ! Matrix Coefficents for PBL Scheme
    real(r8) :: CFu(pcols,pver+1)                        ! Matrix Coefficents for PBL Scheme

    ! JH
    real(r8) :: pi
    pi = 4._r8*atan(1._r8)

    !
    !-----------------------------------------------------------------------
    !

    do i = 1, ncol
      coslat (i) = cos(clat(i))
      sinsq  (i) = sin(clat(i))*sin(clat(i))
      cossq  (i) = coslat(i)*coslat(i)
      cossqsq(i) = cossq (i)*cossq (i)
    end do

    !
    !-----------------------------------------------------------------------
    !
    ! Held/Suarez IDEALIZED physics algorithm:
    !
    !   Held, I. M., and M. J. Suarez, 1994: A proposal for the
    !   intercomparison of the dynamical cores of atmospheric general
    !   circulation models.
    !   Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
    !
    !-----------------------------------------------------------------------
    !
    ! Compute idealized radiative heating rates (as dry static energy)
    !
    !
    do k = 1, pver
      if (pref_mid_norm(k) > sigmab) then
        do i = 1, ncol
          kt = ka + (ks - ka)*cossqsq(i)*(pref_mid_norm(k) - sigmab)/onemsig
          trefc   = 315._r8 - (60._r8 * sinsq(i))
          trefa = (trefc - 10._r8*cossq(i)*log((pmid(i,k)/psurf_ref)))*(pmid(i,k)/psurf_ref)**cappa
          trefa    = max(t00,trefa)
          s(i,k) = (trefa - t(i,k))*kt*cpair
        end do
      else
        do i = 1, ncol
          trefc   = 315._r8 - 60._r8*sinsq(i)
          trefa = (trefc - 10._r8*cossq(i)*log((pmid(i,k)/psurf_ref)))*(pmid(i,k)/psurf_ref)**cappa
          trefa    = max(t00,trefa)
          s(i,k) = (trefa - t(i,k))*ka*cpair
        end do
      end if
    end do
    !
    ! Add diffusion near the surface for the wind fields
    !
    do k = 1, pver
      do i = 1, pcols
        du(i,k) = 0._r8
        dv(i,k) = 0._r8
      end do
    end do

    !
    do k = 1, pver
      if (pref_mid_norm(k) > sigmab) then
        kv  = kf*(pref_mid_norm(k) - sigmab)/onemsig
        do i = 1, ncol
          du(i,k) = -kv*u(i,k)
          dv(i,k) = -kv*v(i,k)
        end do
      end if
    end do

    !
    ! CJ, no Rayleigh frction in the spone layer
    ! CJ, add explicit vertical diffusion of the zonal wind above a threshold
    !
    
    ! JH: debug
    write(iulog,*) ''
    write(iulog,*) 'JH: Applying vertical diffusion profile'
    
    Km = 0._r8
    pT = pint(1, 1)

    do k=1,pver
       do i=1,ncol  
          ! CJ, set the vertical diffusion coefficient, here a constant
          ! if (pmid(i,k) .le. p_vdiff_level) Km(i,k) = nu_v
          
          ! JH, set the vertical diffusion coefficient, here a sin profile, 
          !     with cutoff pressure p_vdiff_level
          !     diffusion coefficient nu_v
          !     diffusivity Km
          !     midpoint pressure pmid
          !     top interface pressure pT
          if (pmid(i, k) < p_vdiff_level) then
            Km(i,k) = nu_v * sin(0.5*pi * log(p_vdiff_level/pmid(i,k))/log(p_vdiff_level/pT) )**2
            if(i==1) write(iulog,*) Km(i, k)
          end if

       end do
    end do

    ! JH
    if(i==1) write(iulog,*) ''

    !
    ! Calculate Diagonal Variables for Implicit PBL Scheme
    !
    do k=1,pver-1
       do i=1,ncol
          rho        = (pint(i,k+1)/(rair*(t(i,k+1)+t(i,k))/2.0_r8)) ! density
          CAm(i,k)   = rpdel(i,k)*ztodt*gravit*gravit*Km(i,k+1)*rho*rho/(pmid(i,k+1)-pmid(i,k))
          CCm(i,k+1) = rpdel(i,k+1)*ztodt*gravit*gravit*Km(i,k+1)*rho*rho/(pmid(i,k+1)-pmid(i,k))
       end do
    end do
    do i=1,ncol
       CAm(i,pver)   = 0._r8
       CCm(i,1)      = 0._r8
       CEm(i,pver+1) = 0._r8
       CFu(i,pver+1) = 0._r8
    end do
    do i=1,ncol
       do k=pver,1,-1
          CEm(i,k) = CCm(i,k)/(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
          CFu(i,k) = (u(i,k)+CAm(i,k)*CFu(i,k+1))/(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
      end do
    end do
!
! Calculate the updated zonal wind tendencies and add them to other tendencies
!
! First we need to calculate the tendencies at the top model level
!
    do i=1,ncol
       du(i,1)   = du(i,1) + (CFu(i,1)-u(i,1))/ztodt
    end do

    do i=1,ncol
       do k=2,pver
          du(i,k)    = du(i,k) + (CEm(i,k)*u(i,k-1)+CFu(i,k)-u(i,k))/ztodt
       end do
    end do


  end subroutine held_suarez_1994

end module held_suarez
