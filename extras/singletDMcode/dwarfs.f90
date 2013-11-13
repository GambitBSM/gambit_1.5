!Fermi-LAT combined dwarfs likelihood module
!Pat Scott patscott@physics.mcgill.ca
!Christoph Weniger (christoph.weniger@gmail.com)
!April 28 2013

module dwarfs_combined

use types
use parameters

implicit none

!Global variables needed by many functions/subroutines
integer, private, parameter :: max_theory_size = 200
real(8), private, dimension(max_theory_size) :: y2, theory_bins, gammas, ekin
integer, private :: ngam
integer, parameter :: original = 1, modified = 2
real(8), private, parameter :: WYieldMassScaling = 0.6724d0, ZYieldMassScaling = 0.68523d0

public dwarflike
private get_spectrum, get_Ntot, calc_dwarf_like, int_theory_spectrum, theory_spectrum2, dwarfLikelihood

contains


  double precision function dwarflike(BFs_in, sigmav, mass, method)

  Type(BFset), intent(IN) :: BFs_in
  Type(BFset) :: BFs
  double precision, intent(IN) :: sigmav, mass
  integer :: method
  double precision :: Ntot, Ntot_Wto4, Ntot_Zto4
  Type(Input_params) :: params
  Type(DM) :: dmdata

  call clear(BFs)
  Ntot_Wto4 = 0.d0
  Ntot_Zto4 = 0.d0

  if (BFs_in%Wto4 .ne. 0.d0) then
    BFs%w = BFs_in%Wto4
    !Get the spectrum for annihilation to WW at ms = mW
    call get_spectrum(BFs, sigmav, mw, params, dmdata)
    BFs%w = 0.d0
    !Calculate the total number of photons between 1 and 100 GeV
    Ntot_Wto4 = get_Ntot(dmdata)
    !Correct the spectrum to the actual mass, to approximate 4-body final states
    Ntot_Wto4 = Ntot_Wto4 * (1.d0 - WYieldMassScaling * (1.d0 - mass/mw)) 
  endif

  if (BFs_in%Zto4 .ne. 0.d0) then
    BFs%z = BFs_in%Zto4
    !Get the spectrum for annihilation to ZZ at ms = mZ
    call get_spectrum(BFs, sigmav, mz, params, dmdata)
    BFs%z = 0.d0
    !Calculate the total number of photons between 1 and 100 GeV
    Ntot_Zto4 = get_Ntot(dmdata)
    !Correct the spectrum to the actual mass, to approximate 4-body final states
    Ntot_Zto4 = Ntot_Zto4 * (1.d0 - ZYieldMassScaling * (1.d0 - mass/mz)) 
  endif

  !Get the spectrum for all 2-body final states
  BFs = BFs_in
  BFs%Wto4 = 0.d0
  BFs%Zto4 = 0.d0
  call get_spectrum(BFs, sigmav, mass, params, dmdata)
  Ntot = get_Ntot(dmdata)

  !Combine the total number of photons for the 2-body and 3/4-body final states
  Ntot = Ntot + Ntot_Wto4 + Ntot_Zto4

  !Get the likeihood
  call calc_dwarf_like(params, dwarflike, Ntot, method)

  dwarflike = -dwarflike

  end function dwarflike


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
  !params%br_Vtotau = BFs%Vtotau 
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
  calcflags%ID_Flags_gamma%dwarfs    = .true.
  calcflags%ID_Flags_gamma%cta_gc    = .false.

  !Set the relevant bits of the dm object
  dmdata%ID_in%gammas%nbins = 100
  dmdata%ID_in%gammas%efluxes_i = 0.1d0 
  dmdata%ID_in%gammas%efluxes_f = 500.d0

  !Get the annihilation spectrum
  call darkmatter(calcflags, params, nuis, dmdata)

  end subroutine get_spectrum

  
  double precision function get_Ntot(DMO)

  Type(DM) :: DMO
  integer :: i, nbins

  gammas = 0d0
  ekin= 0d0
  nbins = DMO%id_in%gammas%nbins

  if (feedback > 2) write (*,*) "Calculate Dwarf Limits..."

  !gammas contains the quantity: spectrum*E^2 in keV (same units as the data)
  !ekin is in keV
  !DMO%id%gammas%ekin(1:nbins) is in MeV
  ekin(1:nbins) = DMO%id%gammas%ekin(1:nbins)*1e-3 !ekin is in GeV !
  !DMO%id%gammas%fluxgadiff is in MeV^-1 !! (if  DMO%ID_in%gammas%delta_gamma  .ne. 0)
  gammas(1:nbins) = DMO%id%gammas%fluxdiff(1:nbins)/1e-3*ekin*ekin
  ! gammas in GeV !
  !ekin, gammas used later on by splines etc... don't remove them!

  ngam = 0
  do i=1,nbins
  if (gammas(i)>0d0) ngam=i
  end do

  y2(:) = 0d0
  call spline_v(log(ekin(1:ngam)),log(gammas(1:ngam)),ngam,0.d0,0.d0,y2(1:ngam))

  ! Total number of photons between 1 and 100 GeV
  get_Ntot = int_theory_spectrum(1.0d0, min(100.0d0, ekin(ngam)))
  if (Feedback > 2) write (*,*) "Ntot (1-100 GeV) = ", get_Ntot

  end function get_Ntot


  subroutine calc_dwarf_like(HardP, lnlike, Ntot, method)

  Type(Input_Params):: HardP
  double precision :: lnlike, Ntot, sigmav, sigmavL2, phi
  integer :: method 

  ! Model sigmav
  sigmav = HardP%sigmav

  !! Koushiappas' limits [1108.2914]

  ! OLD IMPLEMENTATION
  !! To be replaced by Maja´s limits at some point
  !! Koushiappas limit: sigmav < 5e-30 * DMmass**2 / Ntot * 8 * pi      @ 95% CL
  !sigmavL2 = 5.0d-30 * HardP%mx**2 / Ntot * 8 * pi * 1d27    ! 2sigma limit
  !sigmavL3 = 8.0d-30 * HardP%mx**2 / Ntot * 8 * pi * 1d27    ! 3sigma limit (a guess)

  !! Define parameters for cooked-up chi^2 function
  !a = 3*sigmavL2 - 2*sigmavL3
  !b = (sigmavL2 - sigmavL3)**2
  !if (a.lt.0) then
  !  write(*,*) "calc_dwarf_like WARNING: inconsistent chi^2 function"
  !end if

  !if (Feedback > 2) write (*,*) "95% Limit on <sigma v> [10^27 cm^3 s^-1]", sigmavL2

  !! This returns a chi^2/2 function
  !lnlike = 0.5*max(0d0, sign((sigmav-a)**2/b, sigmav-a))

  ! NEW IMPLEMENTATION
  if (Feedback > 2) then
    sigmavL2 = 5.0d-30 * HardP%mx**2 / Ntot * 8 * pi * 1d27    ! 2sigma limit
    write (*,*) "95% Limit on <sigma v> [10^27 cm^3 s^-1]", sigmavL2
  end if

  phi = sigmav/HardP%mx**2*Ntot/8./pi/1d27 * 1d26
  lnlike = 0.5*dwarfLikelihood(phi,method)

  end subroutine calc_dwarf_like


  function int_theory_spectrum(x0, x1)

  real(8) :: int_theory_spectrum, x0, x1, de
  integer :: i

  int_theory_spectrum = 0d0
  de = (x1/x0)**(1d0/100d0)
  do i=1,100
  int_theory_spectrum = int_theory_spectrum + theory_spectrum2(de**i*x0)*(de-1d0)*de**i*x0
  end do

  end function int_theory_spectrum


  function theory_spectrum2(xvalue) 
  ! Theory spectrum is in units of 1/GeV !

  real(8) :: xvalue, theory_spectrum2, yvalue
  !notice: xvalue, theory_spectrum2 are both in linear scale!

  call splint_v(log(ekin(1:ngam)),log(gammas(1:ngam)),y2(1:ngam),ngam,log(xvalue),yvalue)   
  theory_spectrum2 = exp(yvalue)/(xvalue*xvalue)

  end function theory_spectrum2


  function dwarfLikelihood(x,method)
  real(8) :: dwarfLikelihood
  real(8) :: x
  real(8) :: dwarfPhi(101)
  real(8), target :: dwarfL_modified(101)
  real(8), target :: dwarfL_original(101)
  real(8), pointer :: dwarfL(:)
  integer :: method, i

  ! This is the tabulated Phi-Likelihood function from Koushiappas et al.
  ! Above L = 36, we use linear extrapolation up to L = 360000
  data dwarfPhi / 1e-20 , 6.74308086122e-05 , 0.000123192463137 , &
  0.000171713798503 , 0.000215245918518 , 0.000255093268618 , 0.00029207805123 ,&
  0.000326751732695 , 0.000359503469472 , 0.000390620122006 , 0.000420321264006,&
  0.00044878042576 , 0.000476138421008 , 0.000502511975672 , 0.000527999496499,&
  0.000552685056887 , 0.000576641243501 , 0.000599931255273 ,&
  0.000622610497068 , 0.000644727821172 , 0.000666326515638 , 0.000687445105269,&
  0.000708118010141 , 0.000728376093388 , 0.000748247120993 , 0.00076775615078,&
  0.000786925863514 , 0.000805776846231 , 0.000824327835809 , 0.00084259592922,&
  0.000860596765645 , 0.000878344684789 , 0.000895852864914 ,&
  0.000913133443547 , 0.000930197623331 , 0.0009470557651 , 0.000963717469925 ,&
  0.00098019165163 , 0.000996486601006 , 0.00101261004288 , 0.00102856918685 ,&
  0.00104437077256 , 0.00106002111016 , 0.00107552611658 , 0.00109089134805 ,&
  0.00110612202935 , 0.00112122308019 , 0.00113619913897 , 0.00115105458439 ,&
  0.00116579355487 , 0.00118041996631 , 0.00119493752815 , 0.00120934975806 ,&
  0.00122365999528 , 0.00123787141289 , 0.00125198702892 , 0.00126600971667 ,&
  0.00127994221404 , 0.00129378713223 , 0.00130754696367 , 0.00132122408935 ,&
  0.00133482078559 , 0.00134833923028 , 0.00136178150869 , 0.0013751496188 ,&
  0.00138844547626 , 0.00140167091906 , 0.00141482771173 , 0.00142791754942 ,&
  0.00144094206154 , 0.0014539028153 , 0.00146680131887 , 0.00147963902447 ,&
  0.00149241733116 , 0.00150513758749 , 0.00151780109399 , 0.00153040910553 ,&
  0.00154296283341 , 0.00155546344754 , 0.00156791207827 , 0.00158030981824 ,&
  0.00159265772411 , 0.00160495681814 , 0.00161720808976 , 0.00162941249692 ,&
  0.00164157096757 , 0.00165368440081 , 0.00166575366823 , 0.00167777961494 ,&
  0.00168976306076 , 0.00170170480119 , 0.00171360560841 , 0.00172546623219 ,&
  0.00173728740083 , 0.00174906982191 , 0.00176081418314 , 0.00177252115315 ,&
  0.00178419138212 , 0.00179582550256 , 0.00180742412988 , 18.0 /
  data dwarfL_modified / 0.0, & ! Pat's method (normalization w.r.t. p-value of phi=0)
  0.0513551, 0.177438, 0.35228, 0.561353, 0.795726, 1.04953, 1.3187, 1.60032, &
  1.89222, 2.19274, 2.50059, 2.81476, 3.13441, 3.45887, 3.78757, 4.12006,        &
  4.45594, 4.79486, 5.13653, 5.48072, 5.82719, 6.17576, 6.52625, 6.87853,        &
  7.23244, 7.58789, 7.94475, 8.30294, 8.66236, 9.02294, 9.38462, 9.74731,        &
  10.111, 10.4755, 10.841, 11.2072, 11.5742, 11.9419, 12.3104, 12.6795, 13.0492, &
  13.4195, 13.7904, 14.1619, 14.5339, 14.9063, 15.2793, 15.6527, 16.0266,        &
  16.4008, 16.7755, 17.1506, 17.5261, 17.9019, 18.2781, 18.6546, 19.0315,        &
  19.4087, 19.7861, 20.1639, 20.542, 20.9203, 21.2989, 21.6778, 22.0569,         &
  22.4362, 22.8158, 23.1957, 23.5757, 23.956, 24.3365, 24.7171, 25.098, 25.4791, &
  25.8604, 26.2418, 26.6235, 27.0053, 27.3872, 27.7694, 28.1517, 28.5342,        &
  28.9168, 29.2996, 29.6825, 30.0655, 30.4487, 30.8321, 31.2155, 31.5992,        &
  31.9829, 32.3667, 32.7507, 33.1348, 33.519, 33.9034, 34.2878, 34.6724,         &
  35.0571, 350000.0 /
  data dwarfL_original / 0.251157575927 , 0.442789856112 , 0.664135248207 ,&
  0.907852855854 , 1.16899707562 , 1.4440821743 , 1.73056190225 , 2.02652143469,&
  2.33048566777 , 2.64129491971 , 2.95802151481 , 3.27991214086 ,&
  3.60634699937 , 3.93681021383 , 4.27086797558 , 4.60815212568 , 4.94834763107,&
  5.29118289946 , 5.63642219625 , 5.98385964029 , 6.33331440065 ,&
  6.68462681822 , 7.03765524669 , 7.39227345917 , 7.74836850335 , 8.1058389155 ,&
  8.46459322379 , 8.82454868664 , 9.18563022325 , 9.54776950239 , 9.91090416212,&
  10.2749771387 , 10.6399360866 , 11.0057328761 , 11.3723231545 ,&
  11.7396659644 , 12.1077234072 , 12.476460348 , 12.8458441552 , 13.2158444693 ,&
  13.5864329983 , 13.957583335 , 14.3292707946 , 14.7014722687 , 15.0741660946 ,&
  15.4473319381 , 15.8209506877 , 16.1950043588 , 16.5694760077 , 16.9443496528,&
  17.319610204 , 17.6952433975 , 18.0712357373 , 18.4475744412 , 18.8242473915,&
  19.2012430903 , 19.5785506178 , 19.9561595945 , 20.3340601466 ,&
  20.7122428735 , 21.0906988184 , 21.4694194406 , 21.8483965909 , 22.2276224876,&
  22.607089695 , 22.9867911034 , 23.3667199103 , 23.746869603 , 24.1272339425,&
  24.5078069488 , 24.888582886 , 25.2695562501 , 25.6507217564 , 26.0320743279,&
  26.4136090847 , 26.7953213342 , 27.1772065616 , 27.5592604206 ,&
  27.9414787259 , 28.3238574449 , 28.7063926905 , 29.0890807142 , 29.4719178998,&
  29.8549007569 , 30.2380259156 , 30.6212901208 , 31.0046902269 ,&
  31.3882231934 , 31.7718860798 , 32.1556760418 , 32.5395903267 , 32.9236262697,&
  33.3077812904 , 33.6920528891 , 34.0764386433 , 34.4609362049 ,&
  34.8455432972 , 35.2302577118 , 35.6150773059 , 36.0 , 360000.0 /

  select case (method)
  case (original)
    dwarfL => dwarfL_original
  case (modified)
    dwarfL => dwarfL_modified
  case default
    stop 'Unrecognised method in DwarfLikelihood.  Quitting...'
  end select

  if(dwarfPhi(1).ge.x) then
    dwarfLikelihood = dwarfL(1)
    return
  end if
  if(dwarfPhi(101).le.x) then
    dwarfLikelihood = dwarfL(101)
    return
  end if
  do i = 1, 101
    if (dwarfPhi(i).ge.x) exit
  end do
  dwarfLikelihood = (dwarfL(i)-dwarfL(i-1))*(x-dwarfPhi(i-1))/(dwarfPhi(i)-dwarfPhi(i-1)) + dwarfL(i-1)
  return

  end function dwarfLikelihood


end module dwarfs_combined
