!CTA predicted contstraints likelihood module
!Pat Scott patscott@physics.mcgill.ca
!May 20 2013

module cta_gc

use types
use parameters

implicit none

!Parameters for CTA Array B, Ring Method (backgroundArea is not actually used, only signalArea and alpha)
double precision, public :: signalArea = 9.97e-4, backgroundArea = 4.05e-3, alpha = 0.24612 !areas in sr
double precision, public :: signalJ = 6.569e21, backgroundJ = 7.728e21 !in Gev^2/cm^5

!Various Module-level variables
integer, private :: i, ngam
integer, private, parameter :: num_energy_samples = 300
double precision, private, dimension(num_energy_samples) :: gammas, ekin, y2_gam
double precision, private :: hrs = 3600.d0
double precision, private, parameter :: WYieldMassScaling = 22.449567d0, ZYieldMassScaling = 18.9973d0

!Effective area for CTA, DC-like curve from Fig 2 (right) of arXiv:1211.3181; corresponds to 61 mid-sized 
!telescopes (ie CTA North + South).
integer, private, parameter :: nEAbins = 13
double precision, private, parameter :: CTA_EffArea_logE(nEAbins) = (/(1.3+dble(i)*0.2, i = 1, nEAbins)/) !log10(E/GeV)
double precision, private :: CTA_EffAreas(nEAbins), y2_area(nEAbins) = 0.d0 !effective area is log10(EffA/cm^2)
data CTA_EffAreas / 7.279, 8.509, 9.340, 9.584, 9.704, 9.768, 9.798, 9.772, 9.854, 9.890, 9.954, 9.864, 9.800 /

!Default observing time (in hours)
double precision, private, parameter :: defaultObsTime = 500.d0
!Default systematic error on the background estimate 
double precision, private, parameter :: defaultBGsyserr = 0.d0

!Only the init routine and ctalike function are available from outside the module
public ctalike, init_ctalike
private get_spectrum, get_Ntot, calc_cta_like, int_theory_spectrum, integrand

contains


  subroutine init_ctalike
  !Initialise cubic splines for CTA effective area calculation
  call spline_v(CTA_EffArea_logE,CTA_EffAreas,nEAbins,0.d0,0.d0,y2_area)  
  end subroutine init_ctalike


  double precision function ctalike(BFs_in, sigmav, mass, BGrate, observedHrs, BGsyserr)
  !Returns log-likelihood for specified WIMP model from predicted CTA observations
  !Inputs:  BFs_in      Annihilation branching fractions to different final states
  !         sigmav      Annihilation cross-section (cm^3/s)
  !         mass        WIMP mass (GeV)
  !         BGrate      *Optional* Photons/sr/s expected from background; if not given, calculated automatically assuming E^-3 spectrum
  !         observedHrs *Optional* Number of hours for which the GC has been observed (default is 500)
  Type(BFset), intent(IN) :: BFs_in
  Type(BFset) :: BFs
  double precision, intent(IN) :: sigmav, mass
  double precision, intent(IN), optional :: BGrate, observedHrs, BGsyserr
  double precision :: BGrate_internal, obsTime_internal, BGsyserr_internal, Ntot_Wto4, Ntot_Zto4, Ntot
  Type(Input_params) :: params
  Type(DM) :: dmdata

  if (present(BGrate)) BGrate_internal = BGrate

  if (present(observedHrs)) then
    obsTime_internal = observedHrs
  else
    obsTime_internal = defaultObsTime
  endif

  if (present(BGsyserr)) then
    BGsyserr_internal = BGsyserr
  else
    BGsyserr_internal = defaultBGsyserr
  endif

  call clear(BFs)
  Ntot_Wto4 = 0.d0
  Ntot_Zto4 = 0.d0

  if (BFs_in%Wto4 .ne. 0.d0) then
    BFs%w = BFs_in%Wto4
    !Get the spectrum for annihilation to WW at mass = mW
    call get_spectrum(BFs, sigmav, mw, params, dmdata)
    BFs%w = 0.d0
    !Calculate the total number of photons in the observable energy range, scaled by the effective area
    Ntot_Wto4 = get_Ntot(dmdata)
    !Correct the spectrum to the actual mass, to approximate 4-body final states
    Ntot_Wto4 = Ntot_Wto4 * (1.d0 - WYieldMassScaling * (1.d0 - mass/mw)) 
    if (Ntot_Wto4 .lt. 0.d0) Ntot_Wto4 = 0.d0 !stop('Virtual W yield has gone negative!') 
  endif

  if (BFs_in%Zto4 .ne. 0.d0) then
    BFs%z = BFs_in%Zto4
    !Get the spectrum for annihilation to ZZ at ms = mZ
    call get_spectrum(BFs, sigmav, mz, params, dmdata)
    BFs%z = 0.d0
    !Calculate the total number of photons in the observable energy range, scaled by the effective area
    Ntot_Zto4 = get_Ntot(dmdata)
    !Correct the spectrum to the actual mass, to approximate 4-body final states
    Ntot_Zto4 = Ntot_Zto4 * (1.d0 - ZYieldMassScaling * (1.d0 - mass/mz))
    if (Ntot_Zto4 .lt. 0.d0) Ntot_Zto4 = 0.d0 !stop('Virtual Z yield has gone negative!') 
  endif

  !Get the spectrum for all 2-body final states
  BFs = BFs_in
  BFs%Wto4 = 0.d0
  BFs%Zto4 = 0.d0
  call get_spectrum(BFs, sigmav, mass, params, dmdata)
  Ntot = get_Ntot(dmdata)

  !Combine the total number of photons for the 2-body and 3/4-body final states
  Ntot = Ntot + Ntot_Wto4 + Ntot_Zto4
  
  !Find the expected background event rate for energies equal to the WIMP mass and below, assuming E^-3 power law 
  if (.not. present(BGrate)) BGrate_internal = get_Em3BGrate(dmdata,mass)

  !Get the likeihood
  call calc_cta_like(params, ctalike, Ntot, BGrate_internal, obsTime_internal, BGsyserr_internal)

  end function ctalike


  subroutine get_spectrum(BFs, sigmav, mass, params, dmdata)
  double precision, intent(IN) :: sigmav, mass
  Type(BFset), intent(IN) :: BFs
  Type(Input_params) :: params
  Type(Nuisance_params) :: nuis
  Type(Flags) :: calcflags
  Type(DM) :: dmdata

  !Set the params object up
  params%mx = mass
  params%sigmav = sigmav*1.d27
  params%br_uubar = BFs%u 
  params%br_ddbar = BFs%d 
  params%br_ssbar = BFs%s 
  params%br_ccbar = BFs%c 
  params%br_bbar = BFs%b
  params%br_ttbar = BFs%t
  params%br_ww = BFs%w
  params%br_zz = BFs%z 
  params%br_hh = BFs%h
  params%br_ee = BFs%e 
  params%br_mumu = BFs%mu 
  params%br_tautau = BFs%tau 
  !params%br_Vtoe = BFs%Vtoe 
  !params%br_Vtomu = BFs%Vtomu
  !params%br_Vtotau = BFsh%Vtotau 
  params%br_gg = BFs%g 
  params%br_zgam = BFs%zgam 
  params%br_gamgam = BFs%gamgam
  !params%br_hgam = BFs%hgam 
  params%br_invis1 = 0.d0
  params%br_invis2 = 0.d0

  !Set the nuis object up
  nuis%mtop = mt
  !nuis%mhiggs = mh
  
  !Set the calcflags object up
  calcflags%Debug                    = .true.
  calcflags%use_data                 = 1
  calcflags%ID_predict               = .true.
  calcflags%ID_predict_gamma         = .true.
  calcflags%ID_Flags_gamma%halo_fix  = .true.
  calcflags%ID_Flags_gamma%gadiff    = .false.
  calcflags%ID_Flags_gamma%gac       = .false.
  calcflags%ID_Flags_gamma%GC_region = .false.
  calcflags%ID_Flags_gamma%dwarfs    = .false.
  calcflags%ID_Flags_gamma%cta_gc    = .true.

  !Set the relevant bits of the dm object
  dmdata%ID_in%gammas%nbins = num_energy_samples
  dmdata%ID_in%gammas%efluxes_i = 10.d0**CTA_EffArea_logE(1) !Insert a mass-dependent integration threshold here if you want to optimise limits
  dmdata%ID_in%gammas%efluxes_f = min(10.d0**CTA_EffArea_logE(nEAbins),mass)

  !Get the annihilation spectrum
  call darkmatter(calcflags, params, nuis, dmdata)

  end subroutine get_spectrum

  
  double precision function get_Ntot(DMO)
  !Returns total integrated photon yield over energy range, convolved with effective area (units cm^2)
  Type(DM) :: DMO
  integer :: i
  double precision :: log10Emin, log10Emax

  gammas = 0d0
  ekin= 0d0

  if (feedback > 2) write (*,*) "Calculate CTA GC Limits..."

  ekin = log10(DMO%id%gammas%ekin(1:num_energy_samples))-3.d0 !ekin is log10(E/GeV), DMO%id%gammas%ekin in MeV
  where(DMO%id%gammas%fluxdiff .lt. 1.d-100) DMO%id%gammas%fluxdiff = 1.d-100 !prevent zeros
  gammas = log10(DMO%id%gammas%fluxdiff(1:num_energy_samples))+3.d0 !gammas is log10(N/GeV), DMO%id%gammas%fluxgadiff in MeV^-1

  ngam = 0
  do i=1,num_energy_samples
    if (gammas(i) .gt. -20.d0) ngam=i
  end do

  y2_gam = 0.d0
  call spline_v(ekin(1:ngam),gammas(1:ngam),ngam,0.d0,0.d0,y2_gam)

  ! Total number of photons in the observable energy window, convolved with effective area
  log10Emin = max(CTA_EffArea_logE(1), ekin(1))
  log10Emax = min(CTA_EffArea_logE(nEAbins), ekin(ngam))
  get_Ntot = int_theory_spectrum(log10Emin,log10Emax)
  if (Feedback > 2) write (*,*) "Ntot (",10.d0**log10Emin,"-",10.d0**log10Emax," GeV) = ", get_Ntot

  end function get_Ntot


  double precision function get_Em3BGrate(DMO,mass)
  !Returns total integrated photon yield from background model over energy range, convolved with effective area (units cm^2).
  !Output is in units of per second per steradian
  Type(DM), intent(IN) :: DMO
  double precision, intent(IN) :: mass
  double precision :: log10Emin, log10Emax

  gammas = 0d0
  ekin= 0d0

  ekin = log10(DMO%id%gammas%ekin(1:num_energy_samples))-3.d0 !ekin is log10(E/GeV), DMO%id%gammas%ekin in MeV
  !Default background rate, assuming E^3*(CR e+e- flux) = 1.5e-2 Gev^2 cm^-2 s^-1 sr^-1 spectrum as seen by Fermi
  gammas = log10(1.5d-2) - 3.d0*ekin

  ngam = num_energy_samples
  
  y2_gam = 0.d0
  call spline_v(ekin,gammas,num_energy_samples,0.d0,0.d0,y2_gam)

  !Total CR background rate in the energy window of the analysis, convolved with effective area
  log10Emin = max(CTA_EffArea_logE(1), ekin(1)) 
  log10Emax = min(CTA_EffArea_logE(nEAbins), log10(mass))
  get_Em3BGrate = int_theory_spectrum(log10Emin,log10Emax)
  if (Feedback > 2) write (*,*) "BGrate (",10.d0**log10Emin,"-",10.d0**log10Emax," GeV) = ", get_Em3BGrate

  end function get_Em3BGrate


  subroutine calc_cta_like(HardP, lnlike, Ntot, BG_photons_per_sr_per_second, observingTime, syserr)
  !Do the actual likelihood calculation, returning result in lnlike
  Type(Input_Params):: HardP
  double precision, intent(INOUT) :: lnlike
  double precision, intent(IN) :: Ntot, BG_photons_per_sr_per_second, observingTime, syserr
  double precision :: sigmav, N_galacticBG, N_partial, N_BG, N_signal, dsbessei0
  external dsbessei0

  ! Model sigmav
  sigmav = HardP%sigmav

  !Total number of photons expected from backgrounds in the signal region (same as in scaled BG region)
  N_galacticBG = BG_photons_per_sr_per_second * signalArea * observingTime * hrs

  ! Total expected photon rates from annihilation, divided by J factor
  N_partial = sigmav * 1.d-27 * Ntot / (HardP%mx**2 * 8.d0*pi)
  ! Total expected rate from DM annihilation in the background region, scaled by area
  N_BG = alpha * N_partial * backgroundJ
  ! Total expected rate from DM annihilation in the signal region
  N_signal = N_partial * signalJ
  ! Convert to actual expected numbers of photons, adding background contribution
  N_BG = N_BG * observingTime * hrs + (1.d0 + syserr) * N_galacticBG
  N_signal = N_signal * observingTime * hrs + N_galacticBG
  
  ! Calculate log likelihood with a Skellam distribution, assuming expectation value of k=0 difference has been observed 
  lnlike = log(dsbessei0(2.d0*sqrt(N_signal)*sqrt(N_BG))) + 2.d0*sqrt(N_signal)*sqrt(N_BG) - (N_signal + N_BG)
  ! Form that would be required for dealing with k != 0
  !lnlike = 0.5*k*log(N_signal/N_BG) + log(Bessel_I(abs(k),2.d0*sqrt(N_signal*N_BG))) - (N_signal + N_BG)

  ! Calculate the log likelihood of the background-only model, and use it to construct the log likelihood ratio
  lnlike = lnlike - log(dsbessei0(2.d0*sqrt(1.d0 + syserr)*N_galacticBG)) - 2.d0*sqrt(1.d0 + syserr)*N_galacticBG + (2.d0 + syserr)*N_galacticBG
  ! Form that would be required for dealing with k != 0
  !lnlike = lnlike - log(Bessel_I(abs(k),2.d0*sqrt(1.d0 + syserr)*N_galacticBG)) + (2.d0 + syserr)*N_galacticBG

  end subroutine calc_cta_like


  double precision function int_theory_spectrum(x0, x1)
  !Integrate the theory spectrum from log10(E/GeV) = x0 to x1
  double precision, intent(IN) :: x0, x1
  double precision, external :: dmf_int

  integer i
  double precision :: x, effA, diffspec

  int_theory_spectrum = dmf_int(integrand,x0,x1,1.d-6)

  end function int_theory_spectrum


  double precision function integrand(x) 
  ! Integrand in units of cm^2, x is log10(E/GeV)
  double precision :: x, effA, diffspec
 
  call splint_v(ekin(1:ngam),gammas(1:ngam),y2_gam,ngam,x,diffspec)   
  call splint_v(CTA_EffArea_logE,CTA_EffAreas,y2_area,nEAbins,x,effA)   
  integrand = 10.**(effA + diffspec + x)

  end function integrand


end module cta_gc
