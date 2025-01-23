abstract type Microphysics end

Base.@kwdef struct SimpleMicrophysics <: Microphysics end


function (::SimpleMicrophysics)(Phys,RhoPos,ThPos,RhoVPos,RhoCPos,RelCloud,Rain)
  @inline function Microphysics(F,U,Thermo)
    FT = eltype(U)
    p = Thermo[1]
    Rho = U[RhoPos]
    RhoTh = U[ThPos]
    RhoV = U[RhoVPos]
    RhoC = U[RhoCPos]
    RhoD = Rho - RhoV - RhoC
    Cpml = Phys.Cpd * RhoD + Phys.Cpv * RhoV + Phys.Cpl * RhoC
    Rm = Phys.Rd * RhoD + Phys.Rv * RhoV
    kappaM = Rm / Cpml
    T = p / Rm
    p_vs = Thermodynamics.fpws(T,Phys.T0)
    a = p_vs / (Phys.Rv * T) - RhoV
    b = RhoC
    FPh = FT(0.5) * RelCloud * (a + b - sqrt(a * a + b * b))
    L = Thermodynamics.LatHeat(T,Phys.L00,Phys.Cpl,Phys.Cpv,Phys.T0)
    FR = -FPh * Rain
    FRhoTh = RhoTh*((-L/(Cpml*T) - log(p / Phys.p0) * (Rm / Cpml) * (Phys.Rv / Rm + 
      (Phys.Cpl - Phys.Cpv) / Cpml)  + Phys.Rv / Rm) * FPh +
      (FT(1) / Rho - log(p/Phys.p0) * (Rm / Cpml) * (Phys.Cpl / Cpml)) * FR)
    F[RhoPos] += -FR
    F[ThPos] += FRhoTh
    F[RhoVPos] += FPh
    F[RhoCPos] += -FPh - FR
  end
  return Microphysics
end
#=
Base.@kwdef struct OneMomentMicrophysics <: Microphysics end

function (::OneMomentMicrophysicsEquil)(Phys,RhoPos,ThPos,RhoTPos,RhoRPos,RelCloud,Rain)
  @inline function Microphysics(F,U,p)
    FT = eltype(U)
    nz = size(U,1)
    @views Rho = U[:,RhoPs]
    @views RhoTh = U[:,ThPos]
    @views RhoV = U[:,RhoVPos]
    @views RhoC = U[:,RhoCPos]
    @views RhoR = U[:,RhoRPos]

    F[RhoPos] += -FR
    F[ThPos] += FRhoTh
    F[RhoVPos] += FPh
    F[RhoCPos] += -FPh - FR
    dt_max = dt
    for k = 1 : nz
      # Liquid water terminal velocity (m/s) following eq. A13
      velqr = FT(36.34) *(RhoR[k] * FT(0.001))^0.1364 * sqrt(rho[1] / rho[k])
      if velqr > FT(0)
        dt_max = min(dt_max, FT(0.8) * (z[k+1]-z[k]) / velqr)
      end
    end
    # Number of subcycles
    rainsplit = ceiling(dt / dt_max)
    dt0 = dt / real(rainsplit,8)
    for i = 1 : rainsplit
      velqrB = FT(36.34) *(RhoR[1] * FT(0.001))^0.1364   
      for k = 1 : nz - 1
        velqrT = FT(36.34) *(RhoR[1] * FT(0.001))^0.1364 * sqrt(rho[1] / rho[k])
        sed = FT(0.001) * dt0*(r(k+1)*qr(k+1)*velqr(k+1)-r(k)*qr(k)*velqr(k))/(r(k)*(z(k+1)-z(k)))
    end    
  end
  return Microphysics
end

!-----------------------------------------------------------------------
!
!  Version:  1.0
!  Modified from https://github.com/ClimateGlobalChange/DCMIP2016/blob/master/interface/kessler.f90
!
!  Date:  September 6th, 2023
!
!  Change log:
!  Added annotations which refer to Hughes and Jablonowski (2023)
!
!  The KESSLER subroutine implements the Kessler (1969) microphysics
!  parameterization as described by Soong and Ogura (1973) and Klemp
!  and Wilhelmson (1978, KW). KESSLER is called at the end of each
!  time step and makes the final adjustments to the potential
!  temperature and moisture variables due to microphysical processes
!  occurring during that time step. KESSLER is called once for each
!  vertical column of grid cells. Increments are computed and added
!  into the respective variables. The Kessler scheme contains three
!  moisture categories: water vapor, cloud water (liquid water that
!  moves with the flow), and rain water (liquid water that falls
!  relative to the surrounding air). There  are no ice categories.
!  Variables in the column are ordered from the surface to the top.

!  Implementation notes:
!     * Vertical levels in this implementation are assumed to decrease in
!           height as index k increases.
!     * Rho passed into the subroutine refers to the dry density.
!     * For legacy reasons, qv, qc, qr in the code correspond to
!           m_v, m_c, m_r in the paper respectively.
!
!  SUBROUTINE KESSLER(theta, qv, qc, qr, rho, pk, dt, z, nz, rainnc)
!
!  Input variables:
!     theta  - potential temperature (K)
!     qv     - water vapor mixing ratio (gm/gm)
!     qc     - cloud water mixing ratio (gm/gm)
!     qr     - rain  water mixing ratio (gm/gm)
!     rho    - dry air density (not mean state as in KW) (kg/m^3)
!     pk     - Exner function  (not mean state as in KW) (p/p0)**(R/cp)
!     dt     - time step (s)
!     z      - heights of thermodynamic levels in the grid column (m)
!     nz     - number of thermodynamic levels in the column
!     precl  - Precipitation rate (m_water/s)
!
! Output variables:
!     Increments are added into t, qv, qc, qr, and rainnc which are
!     returned to the routine from which KESSLER was called. To obtain
!     the total precip qt, after calling the KESSLER routine, compute:
!
!       qt = sum over surface grid cells of (rainnc * cell area)  (kg)
!       [here, the conversion to kg uses (10^3 kg/m^3)*(10^-3 m/mm) = 1]
!
!
!  Authors:
!           Owen Hughes
!           University of Michigan
!           Email: owhughes@umich.edu
!
!           Paul Ullrich
!           University of California, Davis
!           Email: paullrich@ucdavis.edu
!
!           Based on a code by Joseph Klemp
!           (National Center for Atmospheric Research)
!
!  Reference:
!
!    Klemp, J. B., W. C. Skamarock, W. C., and S.-H. Park, 2015:
!    Idealized Global Nonhydrostatic Atmospheric Test Cases on a Reduced
!    Radius Sphere. Journal of Advances in Modeling Earth Systems.
!    doi:10.1002/2015MS000435
!
!=======================================================================

SUBROUTINE KESSLER(theta, qv, qc, qr, rho, pk, dt, z, nz, precl)

  IMPLICIT NONE

  !------------------------------------------------
  !   Input / output parameters
  !------------------------------------------------

  REAL(8), DIMENSION(nz), INTENT(INOUT) :: &
            theta   ,     & ! Potential temperature (K)
            qv      ,     & ! Water vapor (dry) mixing ratio (gm/gm)
            qc      ,     & ! Cloud water (dry) mixing ratio (gm/gm)
            qr              ! Rain  water (dry) mixing ratio (gm/gm)

  REAL(8), DIMENSION(nz), INTENT(IN) :: &
            rho             ! Dry air density (not mean state as in KW) (kg/m^3)

  REAL(8), INTENT(OUT) :: &
            precl          ! Precipitation rate (m_water / s)

  REAL(8), DIMENSION(nz), INTENT(IN) :: &
            z       ,     & ! Heights of thermo. levels in the grid column (m)
            pk              ! Exner function (p/p0)**(R/cp)

  REAL(8), INTENT(IN) :: &
            dt              ! Time step (s)

  INTEGER, INTENT(IN) :: nz ! Number of thermodynamic levels in the column

  !------------------------------------------------
  !   Local variables
  !------------------------------------------------
  REAL, DIMENSION(nz) :: r, rhalf, velqr, sed, pc

  REAL(8) :: f5, f2x, xk, ern, qrprod, prod, qvs, psl, rhoqr, dt_max, dt0

  INTEGER :: k, rainsplit, nt

  !------------------------------------------------
  !   Begin calculation
  !------------------------------------------------
  f2x = 17.27d0
  f5 = 237.3d0 * f2x * 2500000.d0 / 1003.d0
  xk = .2875d0      !  kappa (r/cp)
  psl    = 1000.d0  !  pressure at sea level (mb)
  rhoqr  = 1000.d0  !  density of liquid water (kg/m^3)

  do k=1,nz
    r(k)     = 0.001d0*rho(k)
    rhalf(k) = sqrt(rho(1)/rho(k))
    pc(k)    = 3.8d0/(pk(k)**(1./xk)*psl)

    ! Liquid water terminal velocity (m/s) following eq. A13
    velqr(k)  = 36.34d0*(qr(k)*r(k))**0.1364*rhalf(k)

  end do

  ! Maximum time step size in accordance with CFL condition
  if (dt .le. 0.d0) then
    write(*,*) 'kessler.f90 called with nonpositive dt'
    stop
  end if

  dt_max = dt
  do k=1,nz-1
    if (velqr(k) .ne. 0.d0) then
      dt_max = min(dt_max, 0.8d0*(z(k+1)-z(k))/velqr(k))
    end if
  end do

  ! Number of subcycles
  rainsplit = ceiling(dt / dt_max)
  dt0 = dt / real(rainsplit,8)

  ! Subcycle through rain process
  precl = 0.d0

  do nt=1,rainsplit

    ! Precipitation rate (m/s)
    precl = precl + rho(1) * qr(1) * velqr(1) / rhoqr

    ! Sedimentation term using upstream differencing, following eq. A12
    do k=1,nz-1
      sed(k) = dt0*(r(k+1)*qr(k+1)*velqr(k+1)-r(k)*qr(k)*velqr(k))/(r(k)*(z(k+1)-z(k)))
    end do
    sed(nz)  = -dt0*qr(nz)*velqr(nz)/(.5*(z(nz)-z(nz-1)))

    ! Adjustment terms
    do k=1,nz

      ! Autoconversion and accretion rates following eq. A11
      qrprod = qc(k) - (qc(k)-dt0*amax1(.001*(qc(k)-.001d0),0.))/(1.d0+dt0*2.2d0*qr(k)**.875)
      qc(k) = amax1(qc(k)-qrprod,0.) ! eq. A9
      qr(k) = amax1(qr(k)+qrprod+sed(k),0.) ! eq. A10

      ! Saturation vapor mixing ratio (gm/gm) following eq. B9
      qvs = pc(k)*exp(f2x*(pk(k)*theta(k)-273.d0)   &
             /(pk(k)*theta(k)- 36.d0))
      prod = (qv(k)-qvs)/(1.d0+qvs*f5/(pk(k)*theta(k)-36.d0)**2)

      ! Evaporation rate following eq. A14, A15
      ern = amin1(dt0*(((1.6d0+124.9d0*(r(k)*qr(k))**.2046)  &
            *(r(k)*qr(k))**.525)/(2550000d0*pc(k)            &
            /(3.8d0 *qvs)+540000d0))*(dim(qvs,qv(k))         &
            /(r(k)*qvs)),amax1(-prod-qc(k),0.),qr(k))

      ! Deviation from simple forms given in H&J paper come from
      ! the saturation adjustment following KW eq. 3.10
      theta(k)= theta(k) + 2500000d0/(1003.d0*pk(k))*(amax1( prod,-qc(k))-ern) ! eq. A1
      qv(k) = amax1(qv(k)-max(prod,-qc(k))+ern,0.) ! eq. A2
      qc(k) = qc(k)+max(prod,-qc(k)) ! eq. A3
      qr(k) = qr(k)-ern ! eq. A4
    end do

    ! Recalculate liquid water terminal velocity
    if (nt .ne. rainsplit) then
      do k=1,nz
        velqr(k)  = 36.34d0*(qr(k)*r(k))**0.1364*rhalf(k) ! eq. A13
      end do
    end if
  end do

  precl = precl / dble(rainsplit)

END SUBROUTINE KESSLER

!=======================================================================

=#

