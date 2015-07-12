!#######################################################################
MODULE DDCALC0
!#######################################################################

!#######################################################################
! DIRECT DETECTION RATES AND LIKELIHOODS (SIMPLE VERSION)
! Routines to calculate the dark matter direct detection event rates
! and corresponding likelihoods/exclusion levels.
! 
! To compile:
!     gfortran -O3 -fno-range-check -c -o DDCalc0.o DDCalc0.f90
!     ifort -fast -c -o DDCalc0.o DDCalc0.f90
! 
! To compile program 'DDCalc0run':
!     gfortran -O3 -fno-range-check DDCalc0.f90 -o DDCalc0run DDCalc0run.f90
!     ifort -fast DDCalc0.f90 -o DDCalc0run DDCalc0run.f90
! 
! To see usage:
!     ./DDCalc0run --help
! 
! 
!   Created by Chris Savage
!   University of Utah   (2013 - 2014)
!   Nordita              (2014 -     )
! 
!   With contributions from:
!     Andre Scaffidi            University of Adelaide (2014)
!     Lauren Hsu                FermiLab (2015)
! 
! 
! 
!   --------========<<<<<<<< FILE LAYOUT >>>>>>>>========--------
! 
! PUBLIC ROUTINE DECLARATIONS
!   Declares which routines can be accessed externally.  Routines
!   not explicitly listed here are not externally accessible due to
!   globally declared default private scope.
! 
! CONSTANTS/GLOBAL VARIABLES/STRUCTURES
!   Contains various constants, structures definitions, and global
!   variables used throughout the module.
! 
! INTERFACES
!   Define names for accessing routines (useful in some circumstances).
! 
! MODULE ROUTINES
!  * SIMPLE INTERFACE ROUTINES
!    Basic routines that can be used to initialize the module and
!    set WIMP parameters.  Meant to provide a simple, basic interface
!    for external access to the module.
! Experiment specific:
!  * XENON100 2012 ROUTINES
!    Routines for calculating rates/likelihoods for the XENON100 2012
!    analysis.
!  * LUX 2013 ROUTINES
!    Routines for calculating rates/likelihoods for the LUX 2013
!    analysis.
!  * SuperCDMS 2014 ROUTINES
!    Routines for calculating rates/likelihoods for the SuperCDMS 2014
!    low-mass WIMP analysis.
!  * DARWIN ARGON 2015 ROUTINES
!    Routines for calculating rates/likelihoods for the proposed
!    DARWIN argon-based detector, using a 2015 parameters estimate.
!  * DARWIN XENON 2015 ROUTINES
!    Routines for calculating rates/likelihoods for the proposed
!    DARWIN xenon-based detector, using a 2015 parameters estimate.
! Generic:
!  * MAIN ROUTINES
!    Main routines used for each program mode.
!  * INITIALIZATION
!    Various general initialization routines.
!  * OUTPUT ROUTINES
!    Routines for printing out info, parameters, and results.
!  * WIMP ROUTINES
!    Routines for initializing, modifying, or viewing WIMP properties,
!    mainly the mass and couplings.
!  * COUPLING CONVERSION ROUTINES
!    Routines for converting between WIMP-nucleon couplings G, f & a,
!    with different notation corresponding to different normalization.
!  * CROSS-SECTION/COUPLING CONVERSION ROUTINES
!    Routines to convert between WIMP-nucleon couplings G and WIMP-
!    nucleon cross-sections sigma.
!  * HALO ROUTINES
!    Routines for initializing and modifying the dark matter halo
!    distribution, as well as perform the mean inverse speed (eta)
!    calculation needed for recoil rates.
!  * ISOTOPES/NUCLEAR ROUTINES (INCLUDING FORM FACTORS)
!    Routines for determining isotope and nuclear properties, including
!    isotopic compositions, nuclear masses, and form factors (both
!    spin-independent and spin-dependent).
!    [NOTE: LIMITED IMPLEMENTATION FOR SD CASE -- MOSTLY SET TO ZERO.]
!  * DETECTOR SETUP ROUTINES
!    Routines for initializing and modifying the detector setup,
!    notatably by setting isotopes to use, loading efficiencies from
!    file, and tabulating form factors.
!  * RATE ROUTINES
!    Routines to calculate and view rates.  The CalcRates() routine
!    here must be called after changing the WIMP mass and/or couplings.
!  * LIKELIHOOD/P-VALUE ROUTINES
!    Routines to calculate the log-likelihood (Poisson with background)
!    or p-value (Poisson or Maximum Gap method, without background
!    subtraction) of the current result.  Must have called CalcRates()
!    for these routines to return appropriate values for the current
!    WIMP mass and couplings.
!  * USAGE
!    Routines for outputting usage instructions.
!  * COMMAND LINE PROCESSING
!    Utility routines for processing command line arguments, mainly
!    of the form '-x', '--flag', or '--key=value'.
!  * EXPORT/IMPORT DATA
!    Utility routines for import and export of data.
!  * TABULATION ROUTINES
!    Routines for determining tabulation points.
!  * INTERPOLATION
!    Interpolation routines, using several functional forms.
!  * RANDOM NUMBER GENERATOR
!    Internal random number generator routines (safe for use with
!    OpenMP).
!  * MATH FUNCTIONS
!    Various math functions.
!  * PROBABILITY DISTRIBUTIONS
!    Various probability related functions for several distributions.
!  * ARRAY SORT/SEARCH
!    Routines to sort or search arrays.
!  * STRING UTILITY FUNCTIONS
!    Utility functions involving strings.
!  * TIME UTILITY FUNCTIONS
!    Utility functions for examining wall or CPU time.
! 
!#######################################################################


!#######################################################################
! SETUP
!#######################################################################

! All types must be explicit
IMPLICIT NONE

! Set default private.
! Constants/routines must be explicitly set to public for
! external access.
PRIVATE

! Version of this software
CHARACTER*(*), PUBLIC, PARAMETER :: VERSION_STRING = '0.26'



!#######################################################################
! PUBLIC ROUTINE DECLARATIONS
!#######################################################################

! Here, we declare which of the routines defined below should have
! public access.  Due to the global private default scope, only
! routines explicitly listed here are externally accessible.
! 
! Interfaces are used to provide an exernally accessible name that
! differs from the actual function name (the latter will remain only
! privately accessible).  The externally accessible names all include
! a 'DDCalc0_' prefix or an experiment-specific prefix (e.g.
! 'LUX_2013_').
! 
! Public constants are declared in the constants section below.


! C++ INTERFACE --------------------------------------------------------

! For ease of use in linking to these routines from C++ code, a second
! (wrapper) version of each public routine is defined that uses
! C-compatible types and an explicitly specified object file symbolic
! name.  These routines are declared here with an extra 'C_' prefix.
! An interface file 'DDCalc0.hpp' is provided that declares functions,
! with the same name and signature as the Fortran ones below, that wrap
! the actual function calls.
! 
! TYPE COMPATIBILITY
! The C++ interface uses the following types in place of the given
! Fortran type:
!   LOGICAL -> bool
!   INTEGER -> int
!   REAL*8  -> double
! As the related quantities may in some cases have different machine
! representations (e.g. number of bits), which can be platform/
! compiler-dependent, BIND() statements are used to declare Fortran
! arguments that are type-compatible with the C++ types.  These
! Fortran types may have different KIND specifications to the ones
! shown above, so conversion to the appropriate types are performed
! in the routines below.
! 
! NAME MANGLING
! The external symbol names in compiled objects/libraries is compiler-
! dependent in some cases.  The gfortran compiler module procedure
! symbol names are of the form:
!   __<modulename>_MOD_<routinename>
! with <modulename> and <routinename> the name of the module and the
! subroutine, respectively, converted to lower case.  The ifort
! compiler instead uses:
!   <modulename>_mp_<routinename>_
! again with the module and subroutine names in lower case.  To allow
! for easier interfacing with C/C++, which must call routines by this
! compiler-dependent external symbol name, we use the Fortran 2003
! BIND() feature to explicitly specify the symbol name for some of the
! public subroutines; the BIND() directive appears in the function
! definitions later in this file.  We choose as our naming convention
! the following:
!   C_DDCALC0_<routinename>
! 
! Name mangling is mainly an issue for modules.  For non-module
! routines, both gfortran and ifort use:
!   <routinename>_
! where the routine name is in lower case.  Use of BIND() for such
! routines requires an interface to be declared in any Fortran routine
! that calls the BIND() routine.


! SIMPLE INTERFACE ROUTINES --------------------------------------------

! These are basic routines for externally accessing the module that
! can be used for initialization and setting WIMP parameters.  These
! are meant to allow for a simpler interface to this module; other
! routines are more robust and provide more capabilities.

! Module initialization:
!     SUBROUTINE DDCalc0_Init()
PUBLIC :: DDCalc0_Init

! Standard Halo Model (SHM) settings:
!     SUBROUTINE DDCalc0_SetSHM(rho,vrot,v0,vesc)
! where the arguments are the local dark matter density [GeV/cm^3],
! the local disk rotation speed [km/s], the most probable dark matter
! speed in the galactic rest frame [km/s], and the galactic escape
! speed [km/s].  This routine need only be called if non-default
! values are desired (0.4, 235, 235, and 550, respectively).
PUBLIC :: DDCalc0_SetSHM

! WIMP parameter settings:
!     SUBROUTINE DDCalc0_SetWIMP_mfa(m,fp,fn,ap,an)
!     SUBROUTINE DDCalc0_SetWIMP_mG(m,GpSI,GnSI,GpSD,GnSD)
!     SUBROUTINE DDCalc0_SetWIMP_msigma(m,sigmapSI,sigmanSI,sigmapSD,sigmanSD)
! where m is the WIMP mass [GeV], f is spin-independent (SI) WIMP-
! nucleon coupling [GeV^-2], a is the spin-dependent (SD) WIMP-nucleon
! coupling [unitless], G are the effective 4-fermion vertex couplings,
! related by:
!     GpSI = 2 fp        GpSD = 2\sqrt{2} G_F ap
!     GnSI = 2 fn        GnSD = 2\sqrt{2} G_F an
! and sigma is the WIMP-nucleon scattering cross-section [pb].
! Negative cross-sections indicate the corresponding coupling should
! be negative.  In all cases, 'p' refers to proton and 'n' to neutron.
PUBLIC :: DDCalc0_SetWIMP_mfa,DDCalc0_SetWIMP_mG,DDCalc0_SetWIMP_msigma

! Same as above, but retrieves current WIMP parameters.
! The only difference: cross-sections are always positive, regardless
! of the sign of the corresponding coupling.
PUBLIC :: DDCalc0_GetWIMP_mfa,DDCalc0_GetWIMP_mG,DDCalc0_GetWIMP_msigma


! EXPERIMENT-SPECIFIC ROUTINES -----------------------------------------

! To use any of the routines for a particular experiment, the
! corresponding init() routine must be called once.  The corresponding
! calc() routine must be called after each change in the WIMP
! parameters.

! Detector initialization:
!     SUBROUTINE Init(intervals)
! where intervals is logical indicating if calculations should be
! performed for analysis sub-intervals (i.e. intervals between
! observed events).  This is only necessary for maximum gap
! calculations and can be set to .FALSE. for likelihood analyses.
PUBLIC :: XENON100_2012_Init,LUX_2013_Init,                             &
          SuperCDMS_2014_Init,SIMPLE_2014_Init,                         &
          DARWIN_Ar_2015_Init,DARWIN_Xe_2015_Init

! Set minimum recoil energy Emin to consider [keV] (initially set to
! 0 keV):
!     SUBROUTINE SetEmin(Emin)
! Note the efficiency curves already account for detector and analysis
! thresholds regardless of this setting, so setting this to 0 keV (the
! default behavior when initialization is performed) does not imply
! that very low energy recoils actually contribute to the signal.
PUBLIC :: XENON100_2012_SetEmin,LUX_2013_SetEmin,                       &
          SuperCDMS_2014_SetEMin,SIMPLE_2014_SetEMin,                   &
          DARWIN_Ar_2015_SetEmin,DARWIN_Xe_2015_SetEmin

! Rate calculation:
!     SUBROUTINE CalcRates()
! Performs the rate calculations used for likelihoods/confidence
! intervals.  Must be called after any changes to the WIMP parameters.
! Actual calculation values are accessed through other routines.
PUBLIC :: XENON100_2012_CalcRates,LUX_2013_CalcRates,                   &
          SuperCDMS_2014_CalcRates,SIMPLE_2014_CalcRates,               &
          DARWIN_Ar_2015_CalcRates,DARWIN_Xe_2015_CalcRates

! Number of observed events in the analysis:
!     INTEGER FUNCTION Events()
PUBLIC :: XENON100_2012_Events,LUX_2013_Events,                         &
          SuperCDMS_2014_Events,SIMPLE_2014_Events,                     &
          DARWIN_Ar_2015_Events,DARWIN_Xe_2015_Events

! Average expected number of background events in the analysis:
!     REAL*8 FUNCTION Background()
PUBLIC :: XENON100_2012_Background,LUX_2013_Background,                 &
          SuperCDMS_2014_Background,SIMPLE_2014_Background,             &
          DARWIN_Ar_2015_Background,DARWIN_Xe_2015_Background

! Average expected number of signal events in the analysis:
!     REAL*8 FUNCTION Signal()
! Or the separate spin-independent and spin-dependent contributions:
!     REAL*8 FUNCTION SignalSI()
!     REAL*8 FUNCTION SignalSD()
PUBLIC :: XENON100_2012_Signal,   LUX_2013_Signal,                      &
          SuperCDMS_2014_Signal,  SIMPLE_2014_Signal,                   &
          DARWIN_Ar_2015_Signal,  DARWIN_Xe_2015_Signal
PUBLIC :: XENON100_2012_SignalSI, LUX_2013_SignalSI,                    &
          SuperCDMS_2014_SignalSI,SIMPLE_2014_SignalSI,                 &
          DARWIN_Ar_2015_SignalSI,DARWIN_Xe_2015_SignalSI
PUBLIC :: XENON100_2012_SignalSD, LUX_2013_SignalSD,                    &
          SuperCDMS_2014_SignalSD,SIMPLE_2014_SignalSD,                 &
          DARWIN_Ar_2015_SignalSD,DARWIN_Xe_2015_SignalSD

! Log-likelihood for current WIMP:
!     REAL*8 FUNCTION LogLikelihood()
! Based upon a Poisson distribution in the number of observed events
! given the expected background+signal.  Calc() must be called first.
PUBLIC :: XENON100_2012_LogLikelihood,LUX_2013_LogLikelihood,           &
          SuperCDMS_2014_LogLikelihood,SIMPLE_2014_LogLikelihood,       &
          DARWIN_Xe_2015_LogLikelihood,DARWIN_Ar_2015_LogLikelihood

! Logarithm of the p-value for current WIMP:
!     REAL*8 FUNCTION LogPValue()
! Based upon the maximum gap method if Init() called with intervals =
! .TRUE., a signal-only (no-background) Poisson distribution otherwise.
! Calc() must be called first.
PUBLIC :: XENON100_2012_LogPValue,LUX_2013_LogPValue,                   &
          SuperCDMS_2014_LogPValue,SIMPLE_2014_LogPValue,               &
          DARWIN_Xe_2015_LogPValue,DARWIN_Ar_2015_LogPValue

! Scale by which the current WIMP cross-sections must be multiplied to
! achieve the given p-value:
!     REAL*8 FUNCTION ScaleToPValue(lnp)
! where lnp is the logarithm of the desired p-value (p=1-CL).
! Calc() must be called first.
PUBLIC :: XENON100_2012_ScaleToPValue,LUX_2013_ScaleToPValue,           &
          SuperCDMS_2014_ScaleToPValue,SIMPLE_2014_ScaleToPValue,       &
          DARWIN_Xe_2015_ScaleToPValue,DARWIN_Ar_2015_ScaleToPValue

! Detector structure initialization (ADVANCED USAGE ONLY):
!     SUBROUTINE InitTo(D,intervals)
! where D is a DetectorRateStruct structure and intervals is as
! described for the initialization routines above.  The D structure
! can be used in the generic routines below.  These routines are
! not intended for standard usage, so ignore them unless you are
! familiar with the internals of this code.
PUBLIC :: XENON100_2012_InitTo,LUX_2013_InitTo,                         &
          SuperCDMS_2014_InitTo,SIMPLE_2014_InitTo,                     &
          DARWIN_Ar_2015_InitTo,DARWIN_Xe_2015_InitTo


! MAIN ROUTINES --------------------------------------------------------

! Main routine
PUBLIC :: DDCalc0_Main
INTERFACE DDCalc0_Main
  MODULE PROCEDURE main
END INTERFACE

! Main routines for specific running modes
PUBLIC :: DDCalc0_MainEventsAndLikelihoods,                             &
          DDCalc0_MainEventsAndLikelihoodsInteractive,                  &
          DDCalc0_MainLogLikelihood,DDCalc0_MainLogLikelihoodInteractive,&
          DDCalc0_MainLogPValue,DDCalc0_MainLogPValueInteractive,       &
          DDCalc0_MainSpectrum,DDCalc0_MainEventsByMass,                &
          DDCalc0_MainConstraintsSI,DDCalc0_MainConstraintsSD,          &
          DDCalc0_MainLimitsSI,DDCalc0_MainLimitsSD
INTERFACE DDCalc0_MainEventsAndLikelihoods
  MODULE PROCEDURE MainEventsAndLikelihoods
END INTERFACE
INTERFACE DDCalc0_MainEventsAndLikelihoodsInteractive
  MODULE PROCEDURE MainEventsAndLikelihoodsInteractive
END INTERFACE
INTERFACE DDCalc0_MainLogLikelihood
  MODULE PROCEDURE MainLogLikelihood
END INTERFACE
INTERFACE DDCalc0_MainLogLikelihoodInteractive
  MODULE PROCEDURE MainLogLikelihoodInteractive
END INTERFACE
INTERFACE DDCalc0_MainLogPValue
  MODULE PROCEDURE MainLogPValue
END INTERFACE
INTERFACE DDCalc0_MainLogPValueInteractive
  MODULE PROCEDURE MainLogPValueInteractive
END INTERFACE
INTERFACE DDCalc0_MainSpectrum
  MODULE PROCEDURE MainSpectrum
END INTERFACE
INTERFACE DDCalc0_MainEventsByMass
  MODULE PROCEDURE MainEventsByMass
END INTERFACE
INTERFACE DDCalc0_MainConstraintsSI
  MODULE PROCEDURE MainConstraintsSI
END INTERFACE
INTERFACE DDCalc0_MainConstraintsSD
  MODULE PROCEDURE MainConstraintsSD
END INTERFACE
INTERFACE DDCalc0_MainLimitsSI
  MODULE PROCEDURE MainLimitsSI
END INTERFACE
INTERFACE DDCalc0_MainLimitsSD
  MODULE PROCEDURE MainLimitsSD
END INTERFACE


! GENERIC ROUTINES -----------------------------------------------------

! Initialization routines
PUBLIC :: DDCalc0_Initialize,DDCalc0_InitializeCL
INTERFACE DDCalc0_Initialize
  MODULE PROCEDURE Initialize
END INTERFACE
INTERFACE DDCalc0_InitializeCL
  MODULE PROCEDURE InitializeCL
END INTERFACE

! WIMP mass & couplings routines
PUBLIC :: DDCalc0_GetWIMP,DDCalc0_SetWIMP,                              &
          DDCalc0_InitWIMP,DDCalc0_InitWIMPCL
INTERFACE DDCalc0_GetWIMP
  MODULE PROCEDURE GetWIMP
END INTERFACE
INTERFACE DDCalc0_SetWIMP
  MODULE PROCEDURE SetWIMP
END INTERFACE
INTERFACE DDCalc0_InitWIMP
  MODULE PROCEDURE InitWIMP
END INTERFACE
INTERFACE DDCalc0_InitWIMPCL
  MODULE PROCEDURE InitWIMPCL
END INTERFACE

! Dark matter halo distribution routines
PUBLIC :: DDCalc0_GetHalo,DDCalc0_SetHalo,                              &
          DDCalc0_InitHalo,DDCalc0_InitHaloCL
INTERFACE DDCalc0_GetHalo
  MODULE PROCEDURE GetHalo
END INTERFACE
INTERFACE DDCalc0_SetHalo
  MODULE PROCEDURE SetHalo
END INTERFACE
INTERFACE DDCalc0_InitHalo
  MODULE PROCEDURE InitHalo
END INTERFACE
INTERFACE DDCalc0_InitHaloCL
  MODULE PROCEDURE InitHaloCL
END INTERFACE

! Detector properties routines
PUBLIC :: DDCalc0_GetDetector,DDCalc0_SetDetector,                      &
          DDCalc0_InitDetector,DDCalc0_InitDetectorCL
INTERFACE DDCalc0_GetDetector
  MODULE PROCEDURE GetDetector
END INTERFACE
INTERFACE DDCalc0_SetDetector
  MODULE PROCEDURE SetDetector
END INTERFACE
INTERFACE DDCalc0_InitDetector
  MODULE PROCEDURE InitDetector
END INTERFACE
INTERFACE DDCalc0_InitDetectorCL
  MODULE PROCEDURE InitDetectorCL
END INTERFACE

! Rate calculation routines
PUBLIC :: DDCalc0_GetRates,DDCalc0_CalcRates
INTERFACE DDCalc0_GetRates
  MODULE PROCEDURE GetRates
END INTERFACE
INTERFACE DDCalc0_CalcRates
  MODULE PROCEDURE CalcRates
END INTERFACE

! Event routines (utility)
PUBLIC :: DDCalc0_Events,DDCalc0_Background
PUBLIC :: DDCalc0_Signal,DDCalc0_SignalSI,DDCalc0_SignalSD
INTERFACE DDCalc0_Events
  MODULE PROCEDURE GetEvents
END INTERFACE
INTERFACE DDCalc0_Background
  MODULE PROCEDURE GetBackground
END INTERFACE
INTERFACE DDCalc0_Signal
  MODULE PROCEDURE GetSignal
END INTERFACE
INTERFACE DDCalc0_SignalSI
  MODULE PROCEDURE GetSignalSI
END INTERFACE
INTERFACE DDCalc0_SignalSD
  MODULE PROCEDURE GetSignalSD
END INTERFACE

! Likelihood/p-value routines
PUBLIC :: DDCalc0_LogLikelihood,DDCalc0_LogPValue,DDCalc0_ScaleToPValue
INTERFACE DDCalc0_LogLikelihood
  MODULE PROCEDURE LogLikelihood
END INTERFACE
INTERFACE DDCalc0_LogPValue
  MODULE PROCEDURE LogPValue
END INTERFACE
INTERFACE DDCalc0_ScaleToPValue
  MODULE PROCEDURE ScaleToPValue
END INTERFACE



!#######################################################################
! CONSTANTS/GLOBAL VARIABLES/STRUCTURES
!#######################################################################

! MATH CONSTANTS -------------------------------------------------------

REAL*8, PUBLIC, PARAMETER :: PI         = 3.1415926535897932d0   ! Pi
REAL*8, PUBLIC, PARAMETER :: TWOPI      = 6.2831853071795865d0   ! 2*Pi
REAL*8, PUBLIC, PARAMETER :: HALFPI     = 1.5707963267948966d0   ! Pi/2
REAL*8, PUBLIC, PARAMETER :: FOURPI     =12.5663706143591730d0   ! 4*Pi
REAL*8, PUBLIC, PARAMETER :: SQRTPI     = 1.7724538509055160d0   ! Sqrt(Pi)
REAL*8, PUBLIC, PARAMETER :: SQRT2PI    = 2.5066282746310005d0   ! Sqrt(2*Pi)
REAL*8, PUBLIC, PARAMETER :: INVPI      = 0.31830988618379067d0  ! 1/Pi
REAL*8, PUBLIC, PARAMETER :: INV2PI     = 0.15915494309189534d0  ! 1/2Pi
REAL*8, PUBLIC, PARAMETER :: INVSQRTPI  = 0.56418958354775629d0  ! 1/Sqrt(Pi)
REAL*8, PUBLIC, PARAMETER :: INVSQRT2PI = 0.39894228040143268d0  ! 1/Sqrt(2*Pi)

REAL*8, PUBLIC, PARAMETER :: SQRT2  = 1.4142135623730950d0   ! Sqrt(2)
REAL*8, PUBLIC, PARAMETER :: SQRT3  = 1.7320508075688773d0   ! Sqrt(3)


! PHYSICS CONSTANTS ----------------------------------------------------

! Speed of light [m/s]
REAL*8, PUBLIC, PARAMETER :: SPEED_OF_LIGHT = 2.99792458d8

! Planck constant times speed of light [GeV fm]
REAL*8, PUBLIC, PARAMETER :: HBARC = 0.1973269718d0

! Fermi coupling constant, in units of /(hbar c)^3 [GeV^-2]
REAL*8, PUBLIC, PARAMETER :: FERMI_COUPLING_CONSTANT = 1.1663787d-5

! Proton and neutron masses [GeV]
REAL*8, PUBLIC, PARAMETER :: PROTON_MASS    = 0.9382720d0
REAL*8, PUBLIC, PARAMETER :: NEUTRON_MASS   = 0.9395654d0


! ARGUMENTS ------------------------------------------------------------

! Will store some command line arguments globally to make parsing
! easier.
!   options     Arguments that begin with '--', usually of the form
!               --<flag> or --<flag>=<value>
!   parameters  Arguments that are not options
!   values      Conversion of parameters to floating point (if possible)
! 
! Currently, only the parameters are used (options are parsed directly).
! 
TYPE, PRIVATE :: ArgumentStructure
  INTEGER :: Noptions    = -1
  INTEGER :: Nparameters = -1
  CHARACTER*256, ALLOCATABLE :: options(:)
  CHARACTER*256, ALLOCATABLE :: parameters(:)
  REAL*8, ALLOCATABLE :: values(:)
END TYPE
  
! For easier parsing, will store command line arguments here
TYPE (ArgumentStructure), PRIVATE :: Arguments


! VERBOSITY ------------------------------------------------------------

! Various parameters that are not accounted for elsewhere

! Verbosity
! Affects the level of output.  Generally, higher magnitude
! means more output; positive will include headers, while
! negative will not.
! 
INTEGER, PRIVATE :: VerbosityLevel = 1


! WIMP -----------------------------------------------------------------

! Structure containing WIMP mass and couplings.
! 
TYPE, PRIVATE :: WIMPStruct
  ! WIMP mass [GeV]
  REAL*8 :: m
  
  ! Effective couplings to the proton and neutron in units of [GeV^-2].
  ! In terms of more commonly used notation:
  !   SI (scalar):        G = 2f
  !   SD (axial-vector):  G = 2\sqrt{2} G_F a
  ! where G, f, and a have 'p' and 'n' subscripts.
  REAL*8 :: GpSI,GnSI,GpSD,GnSD
  
  ! Couplings that yield \sigma = 1 pb in each case.
  REAL*8 :: GpSI0,GnSI0,GpSD0,GnSD0
END TYPE

! Default (internal) WIMP structure
TYPE(WIMPStruct), PRIVATE :: WIMP


! HALO -----------------------------------------------------------------

! Parameters describing the dark matter halo.  Only the Standard Halo
! Model (SHM) can be used for the velocity distribution (i.e. Maxwell-
! Boltzmann distribution with a finite cutoff).

! Local dark matter halo density [GeV/cm^3]:
!   0.3 [standard (old)]
! * Catena & Ullio, JCAP 1008, 004 (2010) [arXiv:0907.0018]
!   For spherical halos, not including structure
!     Einasto profile: 0.385 +/- 0.027
!     NFW profile:     0.389 +/- 0.025
! * Weber & de Boer, arXiv:0910.4272
!     0.2 - 0.4 (depending on model)
! * Salucci et al., arXiv:1003.3101
!   Model independent technique?
!     0.430 +/- 0.113 (alpha) +/- 0.096 (r)
! * Pato et al., PRD 82, 023531 (2010) [arXiv:1006.1322]
!   Density at disk may be 1.01 - 1.41 times larger than shell
!   averaged quantity, so above measurements might be underestimates
!   of local density.
! DEFAULT: 0.4 GeV/cm^3

! Sun's peculiar velocity [km/s]:
! motion relative to local standard of rest (LSR)
! * Mignard, Astron. Astrophys. 354, 522 (2000)
! * Schoenrich, Binney & Dehnen, arXiv:0912.3693
! DEFAULT: (11,12,7) km/s

! Disk rotation velocity [km/s]:
! * Kerr & Lynden-Bell, MNRAS 221, 1023 (1986)
!     220 [standard]
!     222 +/- 20 (average over multiple measurements, could be biased
!                 by systematics)
! * Reid et al., Astrophys. J. 700, 137 (2009) [arXiv:0902.3913]
!   Estimate based on masers.
!     254 +/- 16
! * McMillan & Binney, MNRAS 402, 934 (2010) [arXiv:0907.4685]
!   Reanalysis of Reid et al. masers.  Range of estimates based on
!   various models; suggest Sun's velocity with respect to LSR should
!   be modified.
!     200 +/- 20 to 279 +/- 33
! * Bovy, Hogg & Rix, ApJ 704, 1704 (2009) [arXiv:0907.5423]
!     244 +/- 13 (masers only)
!     236 +/- 11 (combined estimate)
! DEFAULT: 235 km/s

! The Local Standard of Rest (LSR) [km/s], which we take to be
! (0,vrot,0).
! DEFAULT: (0,235,0) km/s

! Sun's velocity vector relative to the galactic rest frame [km/s],
! sum of LSR and peculiar velocity, where LSR = (0,vrot,0):
! DEFAULT: (0,235,0) + (11,12,7) km/s

! Sun's speed relative to the galactic rest frame [km/s].
! Equal to magnitude of Sun's velocity.
! DEFAULT: sqrt{11^2 + (235+12)^2 + 7^2} km/s

! Most probable speed (characterizing velocity dispersion) [km/s]:
! Typically similar to rotation velocity.
!     vrms = sqrt(3/2) v0    [rms velocity]
!     vmp  = v0              [most probably velocity]
!     vave = sqrt(4/pi) v0   [mean velocity]
! DEFAULT: 235 km/s

! Local escape velocity [km/s]:
!   650 [standard (old)]
! * Smith et al., MNRAS 379, 755 (2007) [astro-ph/0611671]
!   Note from Fig 7 that the distribution is asymmetric.  The following
!   results assume vrot = 220.
!     544 (mean), 498 - 608 (90% CL)
!     462 - 640 (90% CL when also fitting parameter k)
! DEFAULT: 550 km/s

! Structure containing halo parameters
TYPE, PRIVATE :: HaloStruct
  ! Galactic motions ---------------------------
  ! Local galactic disk rotation speed [km/s]
  REAL*8 :: vrot = 235d0
  ! Local standard of rest velocity vector [km/s], defined relative to
  ! galactic rest frame.
  REAL*8 :: vlsr(3) = (/ 0d0, 235d0, 0d0 /)
  ! Sun's peculiar velocity vector [km/s], defined relative to local
  ! standard of rest.
  REAL*8 :: vpec(3) = (/ 11d0, 12d0, 7d0 /)
  ! Sun's velocity vector [km/s], defined relative to galactic rest
  ! frame.
  REAL*8 :: vsun(3) = (/ 0d0, 235d0, 0d0 /) + (/ 11d0, 12d0, 7d0 /)
  ! Sun's speed (or observer's speed) [km/s], defined relative to
  ! galactic rest frame.
  REAL*8 :: vobs = SQRT(11d0**2 + (235d0+12d0)**2 + 7d0**2)
  
  ! Local DM density ---------------------------
  ! Local dark matter density [GeV/cm^3]
  REAL*8 :: rho = 0.4d0
  
  ! DM distribution (SHM) ----------------------
  ! Truncated Maxwell-Boltzmann ("MB") distribution.
  
  ! Bulk velocity of the dark matter [km/s] (i.e. the velocity of
  ! the MB rest frame), defined relative to the galactic rest frame.
  REAL*8 :: vbulk(3) = (/ 0d0, 0d0, 0d0 /)
  ! Most probable speed [km/s] in the MB rest frame.
  REAL*8 :: v0 = 235d0
  ! Escape speed [km/s] in the MB rest frame.
  REAL*8 :: vesc = 550d0
  
  ! DM distribution (tabulated) ----------------
  ! Instead of being calculated for SHM above, mean inverse speed
  ! can be given explicitly as a tabulation.
  LOGICAL :: tabulated = .FALSE.
  CHARACTER(LEN=1024) :: eta_file = ''
  INTEGER :: Nvmin = -1
  REAL*8, ALLOCATABLE :: vmin(:)
  REAL*8, ALLOCATABLE :: eta(:)
END TYPE

! Default (internal) WIMP structure
TYPE(HaloStruct), PRIVATE :: Halo


! DETECTOR EFFICIENCY --------------------------------------------------

! <<<<FOR FUTURE USE>>>>
! Structure to contain tabulated detection efficiencies as a
! function of energy, for the overall analysis range and possibly
! for subintervals/bins.
! 
TYPE, PUBLIC :: DetectorEfficiencyStruct
  
  ! File containing efficiencies
  CHARACTER(LEN=1024) :: file = ''
  
  ! Number of tabulation points (energies).
  INTEGER :: NE = -1
  
  ! Tabulated energies [keV].  Array of size [1:NE].
  REAL*8, ALLOCATABLE :: E(:)
  
  ! Number of S1 bins/sub-intervals with efficiencies (does not
  ! include total interval). Will calculate rates for each bin/interval
  ! plus total.
  INTEGER :: Neff = -1
  
  ! Array of size [1:NE,0:Neff] with the second index for the S1
  ! bin/interval (zero for full range)
  REAL*8, ALLOCATABLE :: eff(:,:)
  
END TYPE


! DETECTOR PARAMETERS --------------------------------------------------

! <<<<FOR FUTURE USE>>>>
! Structure to contain various detector parameters.
! 
TYPE, PUBLIC :: DetectorParametersStruct
  
  ! Exposure -----------------------------------
  ! Detector fiducial mass [kg]
  REAL*8 :: mass = 118d0
  
  ! Detector exposure time [day]
  REAL*8 :: time = 85.3d0
  
  ! Total detector exposure [kg*day]
  REAL*8 :: exposure = 118d0*85.3d0
  
  ! Isotopes -----------------------------------
  ! Number of isotopes
  INTEGER :: Niso = -1
  
  ! Detector isotopes, their mass fractions, and nuclear masses [GeV]
  INTEGER, ALLOCATABLE :: Ziso(:)
  INTEGER, ALLOCATABLE :: Aiso(:)
  REAL*8, ALLOCATABLE  :: fiso(:)
  REAL*8, ALLOCATABLE  :: Miso(:)  ! Calculated internally
  
END TYPE


! DETECTOR SPECTRA -----------------------------------------------------

! <<<<FOR FUTURE USE>>>>
! Structure to contain tabulated differential rates dR/dE as a function
! of energy.
! 
TYPE, PUBLIC :: DetectorSpectraStruct
  
  ! Tabulation ---------------------------------
  ! Number of tabulation points (energies).
  ! NOTE: This tabulation is fixed to that used by the efficiency data.
  INTEGER :: NE = -1
  
  ! Tabulated energies [keV].  Array of size [1:NE].
  REAL*8, ALLOCATABLE :: E(:)
  
  ! Efficiencies -------------------------------
  ! Tabulated detection efficiencies.  Here tabulated at desired
  ! E for dR/dE calculations.
  
  ! Number of S1 bins/intervals with efficiencies.
  ! Will calculate rates for each bin/interval plus total.
  INTEGER :: Neff = -1
  
  ! Array of size [1:NE,0:Neff] with the second index for the S1
  ! bin/interval (zero for full range)
  REAL*8, ALLOCATABLE :: eff(:,:)
  
  ! Form factors -------------------------------
  ! Tabulated spin-independent or spin-dependent form factors combined
  ! with prefactors.  Arrays of size [-1:1,1:NE,1:Niso].  Defined as
  ! [unitless]:
  !   Wsi(+1,:,:) = (1/pi) Z^2 F^2(q)        ! SI proton
  !   Wsi( 0,:,:) = (1/pi) 2*Z*(A-Z) F^2(q)  ! SI crossterm
  !   Wsi(-1,:,:) = (1/pi) (A-Z)^2 F^2(q)    ! SI neutron
  !   Wsd(+1,:,:) = 4/(2J+1) Spp(q)          ! SD proton
  !   Wsd( 0,:,:) = 4/(2J+1) Spn(q)          ! SD crossterm
  !   Wsd(-1,:,:) = 4/(2J+1) Snn(q)          ! SD neutron
  ! The above definitions give for the conventional SI and SD
  ! cross-sections:
  !   \sigma(q) = \mu^2 (hbar c)^2 [W(1)*Gp^2 + W(0)*Gp*Gn + W(-1)*Gn^2]
  ! where Gp and Gn are the effective proton and neutron couplings
  ! in units of [GeV^-2] and \mu is the reduced mass.  In terms of
  ! more commonly used notation:
  !   SI (scalar):        G = 2f
  !   SD (axial-vector):  G = 2\sqrt{2} G_F a
  ! where G, f, and a have 'p' and 'n' subscripts.  While form factors
  ! are often a function of the momentum transfer, we tabulate them
  ! here as a function of recoil energy E = q^2/2M.
  ! NOTE: Need only be calculated once.
  REAL*8, ALLOCATABLE :: Wsi(:,:,:),Wsd(:,:,:)
  
  ! Halo ---------------------------------------
  ! The minimum velocity for producing a recoil of energy E, given
  ! by vmin = sqrt{M E/(2\mu^2)} [km/s].
  ! Array of size [1:NE,1:Niso] that needs to be recalculated when the
  ! WIMP mass changes.
  REAL*8, ALLOCATABLE :: vmin(:,:)
  
  ! Tabulated mean inverse speed (eta) [s/km] at the above vmin.
  REAL*8, ALLOCATABLE :: eta(:,:)
  
  ! Differential rates (reference couplings) ---
  ! Reference differential rates calculated at \sigma = 1 pb for each
  ! coupling and tabulated by energy.  Differential arrays are of size
  ! [-1:1,1:NE] and are given separately for SI and SD couplings.  The
  ! first index is for the proton & neutron components at sigma = 1 pb
  ! in each case and the second index is that of the energy (given by
  ! the E array at the same index).  These represent raw rates prior to
  ! any efficiency cuts. [cpd/kg/keV]
  REAL*8, ALLOCATABLE :: dRdEsi0(:,:),dRdEsd0(:,:)
  
  ! Differential rates (actual couplings) ------
  ! Differential rates at given couplings, tabulated by energy.
  ! Arrays are of size [1:NE], given separately for SI and SD couplings
  ! as well as the SI+SD total. [cpd/kg/keV]
  REAL*8, ALLOCATABLE :: dRdEsi(:),dRdEsd(:),dRdE(:)
  
END TYPE


! DETECTOR RATES -------------------------------------------------------

! <<<<FOR FUTURE USE>>>>
! Structure to contain tabulated rates as a function of energy.
! 
TYPE, PUBLIC :: DetectorRateStruct
  
  ! Integrated rates (reference couplings) -----
  ! Reference integrated rates calculated at \sigma = 1 pb for each
  ! coupling, with the efficiency-weighted integral done separately
  ! for each available efficiency curve.  Arrays are of size
  ! [-1:1,0:Neff].  The first index is for the proton & neutron
  ! components at sigma = 1 pb in each case and the second index is for
  ! the S1 bin/interval (0 for full range). [cpd/kg]
  REAL*8, ALLOCATABLE :: Rsi0(:,:),Rsd0(:,:)
  
  ! Integrated rates (actual couplings) --------
  ! Efficiency-corrected rates at given couplings.  Arrays are of size
  ! [0:Neff] with the index being that of the S1 bin/interval
  ! efficiency curve used in the integral (0 for full range).  Given
  ! separately for SI and SD components as well as the SI+SD total.
  ! [cpd/kg]
  REAL*8, ALLOCATABLE :: Rsi(:),Rsd(:),R(:)
  
  ! Events -------------------------------------
  ! Expected number of signal events at reference couplings.
  ! Arrays of size [-1:1,0:Neff].
  REAL*8, ALLOCATABLE :: MuSignalSI0(:,:),MuSignalSD0(:,:)
  ! Expected number of signal events.  Arrays of size [0:Neff].
  REAL*8, ALLOCATABLE :: MuSignalSI(:),MuSignalSD(:),MuSignal(:)
  
  ! Average expected background events
  REAL*8 :: MuBackground = 0d0
  
  ! Observed number of events
  INTEGER :: Nevents = -1
  
END TYPE


! DETECTORS ------------------------------------------------------------

! Structure to contain detector characteristics and tabulated rates as
! a function of energy.
! 
TYPE, PUBLIC :: DetectorStruct
  
  ! Label --------------------------------------
  ! Label for the experimental result contained in this structure.
  ! Must be at most 12 characters as it will be used in column headers.
  CHARACTER(LEN=12) :: label = ''
  ! More detailed description.
  CHARACTER(LEN=1024) :: description = ''
  
  ! FUTURE IMPLEMENTATION ----------------------
  ! Flag to indicate if array sizes and values are outdated
  ! and need to be reinitialized.
  LOGICAL :: stale = .TRUE.
  
  ! Detector parameters
  TYPE(DetectorParametersStruct) :: parameters
  
  ! Detector efficiencies
  TYPE(DetectorEfficiencyStruct) :: efficiency
  
  ! Detector rates
  TYPE(DetectorRateStruct) :: rates
  
  ! Exposure -----------------------------------
  ! Detector fiducial mass [kg]
  REAL*8 :: mass = 118d0
  
  ! Detector exposure time [day]
  REAL*8 :: time = 85.3d0
  
  ! Total detector exposure [kg*day]
  REAL*8 :: exposure = 118d0*85.3d0
  
  ! Events -------------------------------------
  ! Observed number of events
  INTEGER :: Nevents = -1
  
  ! Average expected background events
  REAL*8 :: MuBackground = 0d0
  
  ! Isotopes -----------------------------------
  ! Number of isotopes
  INTEGER :: Niso = -1
  
  ! Detector isotopes, their mass fractions, and nuclear masses [GeV]
  INTEGER, ALLOCATABLE :: Ziso(:)
  INTEGER, ALLOCATABLE :: Aiso(:)
  REAL*8, ALLOCATABLE  :: fiso(:)
  REAL*8, ALLOCATABLE  :: Miso(:)  ! Calculated internally
  
  ! Tabulation ---------------------------------
  ! Number of tabulation points (energies).
  INTEGER :: NE = -1
  
  ! Tabulated energies [keV].  Array of size [1:NE].
  REAL*8, ALLOCATABLE :: E(:)
  
  ! Cached energy tabulation as changes to E array will be made
  ! for internal purposes.
  REAL*8, ALLOCATABLE :: E_cache(:)
  
  ! Efficiencies -------------------------------
  ! Tabulated detection efficiencies.
  ! File containing efficiencies
  CHARACTER(LEN=1024) :: eff_file = ''
  
  ! Number of energies for efficiency tabulation.
  INTEGER :: NEeff = -1
  
  ! Energies for efficiency tabulation [keV].  Array of size [1:NEeff].
  REAL*8, ALLOCATABLE :: Eeff(:)
  
  ! Number of S1 bins/intervals with efficiencies.
  ! Will calculate rates for each bin/interval plus total.
  INTEGER :: Neff = -1
  
  ! Array of size [1:NEeff,0:Neff] with the second index for the S1
  ! bin/interval (zero for full range)
  REAL*8, ALLOCATABLE :: eff(:,:)
  
  ! Indicates if rates for intervals/bins are to also be calculated
  ! in addition to the total rate.  Needed for maximum gap analysis,
  ! but unnecessary for likelihood.  This flag is ignored if the
  ! efficiencies for intervals/bins are not provided.
  LOGICAL :: intervals = .TRUE.
  
  ! Array of size [1:NE,0:Neff] or [1:NE,0:0] containing efficiencies
  ! retabulated to that of the E array.  Sub-interval efficiencies are
  ! dropped if intervals=.FALSE.
  REAL*8, ALLOCATABLE :: eff0(:,:)
  
  ! Form factors -------------------------------
  ! Tabulated spin-independent or spin-dependent form factors combined
  ! with prefactors.  Arrays of size [-1:1,1:NE,1:Niso].  Defined as
  ! [unitless]:
  !   Wsi(+1,:,:) = (1/pi) Z^2 F^2(q)        ! SI proton
  !   Wsi( 0,:,:) = (1/pi) 2*Z*(A-Z) F^2(q)  ! SI crossterm
  !   Wsi(-1,:,:) = (1/pi) (A-Z)^2 F^2(q)    ! SI neutron
  !   Wsd(+1,:,:) = 4/(2J+1) Spp(q)          ! SD proton
  !   Wsd( 0,:,:) = 4/(2J+1) Spn(q)          ! SD crossterm
  !   Wsd(-1,:,:) = 4/(2J+1) Snn(q)          ! SD neutron
  ! The above definitions give for the conventional SI and SD
  ! cross-sections:
  !   \sigma(q) = \mu^2 (hbar c)^2 [W(1)*Gp^2 + W(0)*Gp*Gn + W(-1)*Gn^2]
  ! where Gp and Gn are the effective proton and neutron couplings
  ! in units of [GeV^-2] and \mu is the reduced mass.  In terms of
  ! more commonly used notation:
  !   SI (scalar):        G = 2f
  !   SD (axial-vector):  G = 2\sqrt{2} G_F a
  ! where G, f, and a have 'p' and 'n' subscripts.  While form factors
  ! are often a function of the momentum transfer, we tabulate them
  ! here as a function of recoil energy E = q^2/2M.
  ! NOTE: Need only be calculated once.
  REAL*8, ALLOCATABLE :: Wsi(:,:,:),Wsd(:,:,:)
  
  ! Halo ---------------------------------------
  ! The minimum velocity for producing a recoil of energy E, given
  ! by vmin = sqrt{M E/(2\mu^2)} [km/s].
  ! Array of size [1:NE,1:Niso] that needs to be recalculated when the
  ! WIMP mass changes.
  REAL*8, ALLOCATABLE :: vmin(:,:)
  
  ! Tabulated mean inverse speed (eta) [s/km] at the above vmin.
  REAL*8, ALLOCATABLE :: eta(:,:)
  
  ! Differential rates (reference couplings) ---
  ! Reference differential rates calculated at \sigma = 1 pb for each
  ! coupling and tabulated by energy.  Differential arrays are of size
  ! [-1:1,1:NE] and are given separately for SI and SD couplings.  The
  ! first index is for the proton & neutron components at sigma = 1 pb
  ! in each case and the second index is that of the energy (given by
  ! the E array at the same index).  These represent raw rates prior to
  ! any efficiency cuts. [cpd/kg/keV]
  REAL*8, ALLOCATABLE :: dRdEsi0(:,:),dRdEsd0(:,:)
  
  ! Integrated rates (reference couplings) -----
  ! Reference integrated rates calculated at \sigma = 1 pb for each
  ! coupling, with the efficiency-weighted integral done separately
  ! for each available efficiency curve.  Arrays are of size
  ! [-1:1,0:Neff].  The first index is for the proton & neutron
  ! components at sigma = 1 pb in each case and the second index is for
  ! the S1 bin/interval (0 for full range). [cpd/kg]
  REAL*8, ALLOCATABLE :: Rsi0(:,:),Rsd0(:,:)
  
  ! Differential rates (actual couplings) ------
  ! Differential rates at given couplings, tabulated by energy.
  ! Arrays are of size [1:NE], given separately for SI and SD couplings
  ! as well as the SI+SD total. [cpd/kg/keV]
  REAL*8, ALLOCATABLE :: dRdEsi(:),dRdEsd(:),dRdE(:)
  
  ! Integrated rates (actual couplings) --------
  ! Efficiency-corrected rates at given couplings.  Arrays are of size
  ! [0:Neff] with the index being that of the S1 bin/interval
  ! efficiency curve used in the integral (0 for full range).  Given
  ! separately for SI and SD components as well as the SI+SD total.
  ! [cpd/kg]
  REAL*8, ALLOCATABLE :: Rsi(:),Rsd(:),R(:)
  
  ! Events -------------------------------------
  ! Expected number of signal events at reference couplings.
  ! Arrays of size [-1:1,0:Neff].
  REAL*8, ALLOCATABLE :: MuSignalSI0(:,:),MuSignalSD0(:,:)
  ! Expected number of signal events.  Arrays of size [0:Neff].
  REAL*8, ALLOCATABLE :: MuSignalSI(:),MuSignalSD(:),MuSignal(:)
  
END TYPE

! Default (internal) detector structure.
TYPE(DetectorStruct), TARGET, PRIVATE :: DefaultDetector

! Structures for specific experiments
! NOTE: These will be initialized to internally defined states and
! are not externally modifiable.
TYPE(DetectorStruct), PRIVATE :: XENON100_2012
TYPE(DetectorStruct), PRIVATE :: LUX_2013
TYPE(DetectorStruct), PRIVATE :: SuperCDMS_2014
TYPE(DetectorStruct), PRIVATE :: SIMPLE_2014
TYPE(DetectorStruct), PRIVATE :: DARWIN_Ar_2015,DARWIN_Xe_2015


! TABULATION -----------------------------------------------------------

! Structure to contain fixed spacing linear or logarithmic tabulation
! parametrization.
! 
TYPE, PRIVATE :: TabulationStruct
  ! Linear or logarithmic spacing?
  LOGICAL :: logarithmic = .FALSE.
  ! Number of intervals between tabulation points
  INTEGER :: N
  ! Range of tabulation
  REAL*8 :: xmin,xmax
  ! Logarithm of tabulation range values
  REAL*8 :: lnxmin,lnxmax
  ! Spacing between tabulation points:
  !   linear:       delta = x_{k+1} - x_k
  !   logarithmic:  delta = ln(x_{k+1}) - ln(x_k) = ln(x_{k+1}/x_k)
  REAL*8 :: delta
END TYPE


! FORMATTING -----------------------------------------------------------

! Prefix to place at beginning of comment and data lines
CHARACTER*(*), PARAMETER :: COMMENT_PREFIX = '# '
CHARACTER*(*), PARAMETER :: DATA_PREFIX    = '  '
CHARACTER*(*), PARAMETER :: COMMENT_LINE   = &
    '# ----------------------------------'   &
    // '------------------------------------'


! RANDOM NUMBER GENERATOR STATE ----------------------------------------

! The structure here stores state data for the uniform random number
! generator based on the algorithm of Marsaglia & Zaman (and later
! Tsang).  See the 'RANDOM NUMBER GENERATOR' section below for
! references.
! 
TYPE, PRIVATE :: RandState
  LOGICAL :: initialized = .FALSE.
  INTEGER :: I,J      ! I97,J97
  REAL*8 ::  U(1:97)
  REAL*8 ::  C,CD,CM
END TYPE

! Instantiation for default
! OpenMP directive ensures each thread maintains own copy
TYPE(RandState), PRIVATE :: DEFAULT_RandState
!$OMP THREADPRIVATE(DEFAULT_RandState)



!#######################################################################
! INTERFACES
!#######################################################################

! Here define interfaces, which allow for the following:
!  *) Provide a new name for accessing a routine.
!  *) Allow private routines to be externally accessed
!     under the new name (if interface is public).
!  *) Allow for multiple routines to be accessed in the same
!     name, as long as signatures are not ambiguous.

! Interpolation functions: scalar and array-valued
INTERFACE LinearInterpolate
  MODULE PROCEDURE LinearInterpolate_S,LinearInterpolate_A
END INTERFACE

! Convert E to vmin given the scalar or vector arguments.
! Access all versions under function name 'EToVmin'.
INTERFACE EToVmin
  MODULE PROCEDURE EToVmin0,EToVmin1,EToVmin2
END INTERFACE

! Calculate mean inverse speed from vmin given the scalar or vector
! arguments. Access all versions under function name 'MeanInverseSpeed'.
INTERFACE MeanInverseSpeed
  MODULE PROCEDURE MeanInverseSpeed0,MeanInverseSpeed1,MeanInverseSpeed2
END INTERFACE



!#######################################################################
! ROUTINES
!#######################################################################

CONTAINS


!=======================================================================
! SIMPLE INTERFACE ROUTINES
! Basic versions of routines that are required for using this module
! externally.  These are meant to allow for a simpler interface to
! this module; other routines are more robust and provide more
! capabilities.
! 
! BIND() is used to specify compiler-independent object file symbol
! names to allow for easier interfacing with C/C++.  
!=======================================================================

!-----------------------------------------------------------------------
! Initializes the module.  Must be run once before use of any
! module routines.
! 
! For more detailed initialization, see:
!   Initialize() [interface name: DDCalc0_Initialize]
! 
SUBROUTINE DDCalc0_Init()
  IMPLICIT NONE
  CALL Initialize(.TRUE.)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_DDCalc0_Init() &
           BIND(C,NAME='C_DDCALC0_ddcalc0_init')
  IMPLICIT NONE
  CALL DDCalc0_Init()
END SUBROUTINE


!-----------------------------------------------------------------------
! Sets the dark matter halo to the Standard Halo Model (SHM) with the
! given parameters.  Need only be run if non-default parameters are
! to be used.
! 
! For more detailed halo settings, see:
!   SetHalo() [interface name: DDCalc0_SetHalo]
! 
! Input arguments:
!   rho         Local dark matter density [GeV/cm^3].  Default is
!               0.4 GeV/cm^3.
!   vrot        Local galactic disk rotation speed [km/s].  Default is
!               235 km/s.
!   v0          Most probable speed [km/s] in the galactic rest frame.
!               For the conventional isothermal sphere, this should be
!               the same as vrot.  Default is 235 km/s.
!   vesc        Galactic escape speed [km/s] in the galactic rest
!               frame.  Default is 550 km/s.
! 
SUBROUTINE DDCalc0_SetSHM(rho,vrot,v0,vesc)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: rho,vrot,v0,vesc
  CALL SetHalo(rho=rho,vrot=vrot,v0=v0,vesc=vesc)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_DDCalc0_SetSHM(rho,vrot,v0,vesc) &
           BIND(C,NAME='C_DDCALC0_ddcalc0_setshm')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: rho,vrot,v0,vesc
  CALL DDCalc0_SetSHM(rho=REAL(rho,KIND=8),vrot=REAL(vrot,KIND=8),      &
                      v0=REAL(v0,KIND=8),vesc=REAL(vesc,KIND=8))
END SUBROUTINE


!-----------------------------------------------------------------------
! Sets/gets the WIMP mass and couplings.  Couplings are specified via
! the commonly used fp/fn (spin-independent) and ap/an (spin-dependent)
! normalizations.
! 
! For more detailed WIMP settings, see:
!   SetWIMP() [interface name: DDCalc0_SetWIMP]
!   GetWIMP() [interface name: DDCalc0_GetWIMP]
! 
! Input/output arguments:
!   m           WIMP mass [GeV].
!   fp          Spin-independent WIMP-proton coupling [GeV^-2].
!               Related by GpSI = 2 fp.
!   fn          Spin-independent WIMP-neutron coupling [GeV^-2].
!               Related by GnSI = 2 fn.
!   ap          Spin-dependent WIMP-proton coupling [unitless].
!               Related by GpSD = 2\sqrt{2} G_F ap.
!   an          Spin-dependent WIMP-neutron coupling [unitless].
!               Related by GnSD = 2\sqrt{2} G_F an.
! 
SUBROUTINE DDCalc0_SetWIMP_mfa(m,fp,fn,ap,an)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: m,fp,fn,ap,an
  CALL SetWIMP(m=m,fp=fp,fn=fn,ap=ap,an=an)
END SUBROUTINE

SUBROUTINE DDCalc0_GetWIMP_mfa(m,fp,fn,ap,an)
  IMPLICIT NONE
  REAL*8, INTENT(OUT) :: m,fp,fn,ap,an
  CALL GetWIMP(m=m,fp=fp,fn=fn,ap=ap,an=an)
END SUBROUTINE


! C++ interface wrappers
SUBROUTINE C_DDCalc0_SetWIMP_mfa(m,fp,fn,ap,an) &
           BIND(C,NAME='C_DDCALC0_ddcalc0_setwimp_mfa')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: m,fp,fn,ap,an
  CALL DDCalc0_SetWIMP_mfa(m=REAL(m,KIND=8),                            &
               fp=REAL(fp,KIND=8),fn=REAL(fn,KIND=8),                   &
               ap=REAL(ap,KIND=8),an=REAL(an,KIND=8))
END SUBROUTINE

SUBROUTINE C_DDCalc0_GetWIMP_mfa(m,fp,fn,ap,an) &
           BIND(C,NAME='C_DDCALC0_ddcalc0_getwimp_mfa')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: m,fp,fn,ap,an
  REAL*8 :: m0,fp0,fn0,ap0,an0
  CALL DDCalc0_GetWIMP_mfa(m=m0,fp=fp0,fn=fn0,ap=ap0,an=an0)
  ! Automatic type conversions here
  m  = m0
  fp = fp0
  fn = fn0
  ap = ap0
  an = an0
END SUBROUTINE


!-----------------------------------------------------------------------
! Sets/gets the WIMP mass and couplings.  Couplings are specified via
! their effective 4-fermion vertex couplings 'G'.
! 
! For more detailed WIMP settings, see:
!   SetWIMP() [interface name: DDCalc0_SetWIMP]
!   GetWIMP() [interface name: DDCalc0_GetWIMP]
! 
! Input/output arguments:
!   m           WIMP mass [GeV].
!   GpSI        Spin-independent WIMP-proton coupling [GeV^-2].
!               Related by GpSI = 2 fp.
!   GnSI        Spin-independent WIMP-neutron coupling [GeV^-2].
!               Related by GnSI = 2 fn.
!   GpSD        Spin-dependent WIMP-proton coupling [GeV^-2].
!               Related by GpSD = 2\sqrt{2} G_F ap.
!   GnSD        Spin-dependent WIMP-neutron coupling [GeV^-2].
!               Related by GnSD = 2\sqrt{2} G_F an.
! 
SUBROUTINE DDCalc0_SetWIMP_mG(m,GpSI,GnSI,GpSD,GnSD)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: m,GpSI,GnSI,GpSD,GnSD
  CALL SetWIMP(m=m,GpSI=GpSI,GnSI=GnSI,GpSD=GpSD,GnSD=GnSD)
END SUBROUTINE

SUBROUTINE DDCalc0_GetWIMP_mG(m,GpSI,GnSI,GpSD,GnSD)
  IMPLICIT NONE
  REAL*8, INTENT(OUT) :: m,GpSI,GnSI,GpSD,GnSD
  CALL GetWIMP(m=m,GpSI=GpSI,GnSI=GnSI,GpSD=GpSD,GnSD=GnSD)
END SUBROUTINE


! C++ interface wrappers
SUBROUTINE C_DDCalc0_SetWIMP_mG(m,GpSI,GnSI,GpSD,GnSD) &
           BIND(C,NAME='C_DDCALC0_ddcalc0_setwimp_mg')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: m,GpSI,GnSI,GpSD,GnSD
  CALL DDCalc0_SetWIMP_mG(m=REAL(m,KIND=8),                             &
               GpSI=REAL(GpSI,KIND=8),GnSI=REAL(GnSI,KIND=8),           &
               GpSD=REAL(GpSD,KIND=8),GnSD=REAL(GnSD,KIND=8))
END SUBROUTINE

SUBROUTINE C_DDCalc0_GetWIMP_mG(m,GpSI,GnSI,GpSD,GnSD) &
           BIND(C,NAME='C_DDCALC0_ddcalc0_getwimp_mg')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: m,GpSI,GnSI,GpSD,GnSD
  REAL*8 :: m0,GpSI0,GnSI0,GpSD0,GnSD0
  CALL DDCalc0_GetWIMP_mG(m=m0,GpSI=GpSI0,GnSI=GnSI0,GpSD=GpSD0,GnSD=GnSD0)
  ! Automatic type conversions here
  m    = m0
  GpSI = GpSI0
  GnSI = GnSI0
  GpSD = GpSD0
  GnSD = GnSD0
END SUBROUTINE


!-----------------------------------------------------------------------
! Sets/gets the WIMP mass and couplings.  Couplings are specified via
! their corresponding WIMP-nucleon cross-sections.
! 
! For more detailed WIMP settings, see:
!   SetWIMP() [interface name: DDCalc0_SetWIMP].
!   GetWIMP() [interface name: DDCalc0_GetWIMP]
! 
! Input/output arguments.  As input, give negative value to cross-
! sections to set the corresponding coupling negative:
!   m           WIMP mass [GeV].
!   sigmapSI    Spin-independent WIMP-proton cross-section [pb].
!   sigmanSI    Spin-independent WIMP-neutron cross-section [pb].
!   sigmapSD    Spin-dependent WIMP-proton cross-section [pb].
!   sigmanSD    Spin-dependent WIMP-neutron cross-section [pb].
! 
SUBROUTINE DDCalc0_SetWIMP_msigma(m,sigmapSI,sigmanSI,sigmapSD,sigmanSD)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: m,sigmapSI,sigmanSI,sigmapSD,sigmanSD
  CALL SetWIMP(m=m,sigmapSI=sigmapSI,sigmanSI=sigmanSI,                 &
               sigmapSD=sigmapSD,sigmanSD=sigmanSD)
END SUBROUTINE

SUBROUTINE DDCalc0_GetWIMP_msigma(m,sigmapSI,sigmanSI,sigmapSD,sigmanSD)
  IMPLICIT NONE
  REAL*8, INTENT(OUT) :: m,sigmapSI,sigmanSI,sigmapSD,sigmanSD
  CALL GetWIMP(m=m,sigmapSI=sigmapSI,sigmanSI=sigmanSI,                 &
               sigmapSD=sigmapSD,sigmanSD=sigmanSD)
END SUBROUTINE


! C++ interface wrappers
SUBROUTINE C_DDCalc0_SetWIMP_msigma(m,sigmapSI,sigmanSI,sigmapSD,sigmanSD) &
           BIND(C,NAME='C_DDCALC0_ddcalc0_setwimp_msigma')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: m,sigmapSI,sigmanSI,sigmapSD,sigmanSD
  CALL DDCalc0_SetWIMP_msigma(m=REAL(m,KIND=8),                         &
               sigmapSI=REAL(sigmapSI,KIND=8),sigmanSI=REAL(sigmanSI,KIND=8),&
               sigmapSD=REAL(sigmapSD,KIND=8),sigmanSD=REAL(sigmanSD,KIND=8))
END SUBROUTINE

SUBROUTINE C_DDCalc0_GetWIMP_msigma(m,sigmapSI,sigmanSI,sigmapSD,sigmanSD) &
           BIND(C,NAME='C_DDCALC0_ddcalc0_getwimp_msigma')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(OUT) :: m,sigmapSI,sigmanSI,sigmapSD,sigmanSD
  REAL*8 :: m0,sigmapSI0,sigmanSI0,sigmapSD0,sigmanSD0
  CALL DDCalc0_GetWIMP_msigma(m=m0,sigmapSI=sigmapSI0,sigmanSI=sigmanSI0,&
                              sigmapSD=sigmapSD0,sigmanSD=sigmanSD0)
  ! Automatic type conversions here
  m        = m0
  sigmapSI = sigmapSI0
  sigmanSI = sigmanSI0
  sigmapSD = sigmapSD0
  sigmanSD = sigmanSD0
END SUBROUTINE



!=======================================================================
! XENON100 2012 ANALYSIS ROUTINES
! Based upon the XENON100 2012 analysis [1207.5988].  Two events seen in
! the analysis region.
! 
! BIND() is used to specify compiler-independent object file symbol
! names to allow for easier interfacing with C/C++.  
!=======================================================================

!-----------------------------------------------------------------------
! Initializes the module to perform calculations for the XENON100 2012
! analysis.  This must be called if any of the following XENON100 2012
! routines are to be used.
! 
! Required input arguments:
!     intervals   Indicates if sub-intervals should be included.
!                 Only necessary if confidence intervals using the
!                 maximum gap method are desired.
! 
SUBROUTINE XENON100_2012_Init(intervals)
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: intervals
  CALL XENON100_2012_InitTo(XENON100_2012,intervals)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_XENON100_2012_Init(intervals) &
           BIND(C,NAME='C_DDCALC0_xenon100_2012_init')
  USE ISO_C_BINDING, only: C_BOOL
  IMPLICIT NONE
  LOGICAL(KIND=C_BOOL), INTENT(IN) :: intervals
  CALL XENON100_2012_Init(LOGICAL(intervals))
END SUBROUTINE


! ----------------------------------------------------------------------
! Sets the minimum recoil energy to be included in the calculations.
! Note the efficiency curves already account for detector and analysis
! thresholds regardless of this setting, so setting this to 0 keV (the
! default behavior when initialization is performed) does not imply
! that very low energy recoils actually contribute to the signal.
! 
! Required input arguments:
!     Emin        The minimum recoil energy to consider [keV]
! 
SUBROUTINE XENON100_2012_SetEmin(Emin)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: Emin
  CALL SetDetector(XENON100_2012,Emin=Emin)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_XENON100_2012_SetEmin(Emin) &
           BIND(C,NAME='C_DDCALC0_xenon100_2012_setemin')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: Emin
  CALL XENON100_2012_SetEmin(REAL(Emin,KIND=8))
END SUBROUTINE


! ----------------------------------------------------------------------
! Calculates various rate quantities using the current WIMP.
! Must be called each time the WIMP parameters are modified.
! 
SUBROUTINE XENON100_2012_CalcRates()
  IMPLICIT NONE
  CALL CalcRates(XENON100_2012)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_XENON100_2012_CalcRates() &
           BIND(C,NAME='C_DDCALC0_xenon100_2012_calcrates')
  IMPLICIT NONE
  CALL XENON100_2012_CalcRates()
END SUBROUTINE


! ----------------------------------------------------------------------
! Returns the observed number of events.
! 
FUNCTION XENON100_2012_Events() RESULT(N)
  IMPLICIT NONE
  INTEGER :: N
  CALL GetRates(XENON100_2012,Nevents=N)
END FUNCTION

! C++ interface wrapper
FUNCTION C_XENON100_2012_Events() RESULT(N) &
         BIND(C,NAME='C_DDCALC0_xenon100_2012_events')
  USE ISO_C_BINDING, only: C_INT
  IMPLICIT NONE
  INTEGER(KIND=C_INT) :: N
  ! Automatic type conversions here
  N = XENON100_2012_Events()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of background events.
! 
FUNCTION XENON100_2012_Background() RESULT(b)
  IMPLICIT NONE
  REAL*8 :: b
  CALL GetRates(XENON100_2012,background=b)
END FUNCTION

! C++ interface wrapper
FUNCTION C_XENON100_2012_Background() RESULT(b) &
         BIND(C,NAME='C_DDCALC0_xenon100_2012_background')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: b
  ! Automatic type conversions here
  b = XENON100_2012_Background()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of signal events for the
! current WIMP.
! 
FUNCTION XENON100_2012_Signal() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(XENON100_2012,signal=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_XENON100_2012_Signal() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_xenon100_2012_signal')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = XENON100_2012_Signal()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of spin-independent signal events
! for the current WIMP.
! 
FUNCTION XENON100_2012_SignalSI() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(XENON100_2012,signal_si=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_XENON100_2012_SignalSI() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_xenon100_2012_signalsi')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = XENON100_2012_SignalSI()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of spin-dependent signal events
! for the current WIMP.
! 
FUNCTION XENON100_2012_SignalSD() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(XENON100_2012,signal_sd=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_XENON100_2012_SignalSD() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_xenon100_2012_signalsd')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = XENON100_2012_SignalSD()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the log-likelihood for the current WIMP mass and couplings.
! Uses a Poisson distribution in the number of observed events N:
!    P(N|s+b)
! where s is the average expected signal and b is the average expected
! background.
! 
FUNCTION XENON100_2012_LogLikelihood() RESULT(lnlike)
  IMPLICIT NONE
  REAL*8 :: lnlike
  lnlike = LogLikelihood(XENON100_2012)
END FUNCTION

! C++ interface wrapper
FUNCTION C_XENON100_2012_LogLikelihood() RESULT(lnlike) &
         BIND(C,NAME='C_DDCALC0_xenon100_2012_loglikelihood')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: lnlike
  ! Automatic type conversions here
  lnlike = XENON100_2012_LogLikelihood()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the log of the p-value for the current WIMP mass and
! couplings (NO BACKGROUND SUBTRACTION).  Uses the maximum gap method
! if XENON100_2012_Init was called with argument intervals=.TRUE.,
! otherwise uses a Poisson distribution in the number of observed
! events N:
!    P(N|s)
! where s is the average expected signal (background contributions are
! ignored).
! 
FUNCTION XENON100_2012_LogPValue() RESULT(lnp)
  IMPLICIT NONE
  REAL*8 :: lnp
  lnp = LogPValue(XENON100_2012)
END FUNCTION

! C++ interface wrapper
FUNCTION C_XENON100_2012_LogPValue() RESULT(lnp) &
         BIND(C,NAME='C_DDCALC0_xenon100_2012_logpvalue')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: lnp
  ! Automatic type conversions here
  lnp = XENON100_2012_LogPValue()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the factor x by which the cross-sections must be scaled
! (sigma -> x*sigma) to achieve the desired p-value (given as log(p)).
! See LogPValue() above for a description of the statistics.
! 
! Required input argument:
!   lnp         The logarithm of the desired p-value (p = 1-CL).
! 
FUNCTION XENON100_2012_ScaleToPValue(lnp) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN) :: lnp
  x = ScaleToPValue(XENON100_2012,lnp)
END FUNCTION

! C++ interface wrapper
FUNCTION C_XENON100_2012_ScaleToPValue(lnp) RESULT(x) &
         BIND(C,NAME='C_DDCALC0_xenon100_2012_scaletopvalue')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: x
  REAL(KIND=C_DOUBLE), INTENT(IN) :: lnp
  ! Automatic type conversions here
  x = XENON100_2012_ScaleToPValue(REAL(lnp,KIND=8))
END FUNCTION


!-----------------------------------------------------------------------
! INTERNAL ROUTINE.
! Initializes the given DetectorStruct to the XENON100 2012 analysis.
! This is meant as an internal routine; external access should be
! through XENON100_2012_Init instead.
! 
! The efficiencies used here were generated using TPCMC.
! 
! Required input arguments:
!     D           The DetectorStruct to initialize
!     intervals   Indicates if sub-intervals should be included
! 
SUBROUTINE XENON100_2012_InitTo(D,intervals)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(OUT) :: D
  LOGICAL, INTENT(IN) :: intervals
  INTEGER, PARAMETER :: NE = 151
  INTEGER, PARAMETER :: NEFF = 3
  ! Efficiency curves energy tabulation points
  REAL*8, PARAMETER :: E(NE)                                            &
      =       (/ 0.10000d0, 0.10471d0, 0.10965d0, 0.11482d0, 0.12023d0, &
      0.12589d0, 0.13183d0, 0.13804d0, 0.14454d0, 0.15136d0, 0.15849d0, &
      0.16596d0, 0.17378d0, 0.18197d0, 0.19055d0, 0.19953d0, 0.20893d0, &
      0.21878d0, 0.22909d0, 0.23988d0, 0.25119d0, 0.26303d0, 0.27542d0, &
      0.28840d0, 0.30200d0, 0.31623d0, 0.33113d0, 0.34674d0, 0.36308d0, &
      0.38019d0, 0.39811d0, 0.41687d0, 0.43652d0, 0.45709d0, 0.47863d0, &
      0.50119d0, 0.52481d0, 0.54954d0, 0.57544d0, 0.60256d0, 0.63096d0, &
      0.66069d0, 0.69183d0, 0.72444d0, 0.75858d0, 0.79433d0, 0.83176d0, &
      0.87096d0, 0.91201d0, 0.95499d0, 1.0000d0,  1.0471d0,  1.0965d0,  &
      1.1482d0,  1.2023d0,  1.2589d0,  1.3183d0,  1.3804d0,  1.4454d0,  &
      1.5136d0,  1.5849d0,  1.6596d0,  1.7378d0,  1.8197d0,  1.9055d0,  &
      1.9953d0,  2.0893d0,  2.1878d0,  2.2909d0,  2.3988d0,  2.5119d0,  &
      2.6303d0,  2.7542d0,  2.8840d0,  3.0200d0,  3.1623d0,  3.3113d0,  &
      3.4674d0,  3.6308d0,  3.8019d0,  3.9811d0,  4.1687d0,  4.3652d0,  &
      4.5709d0,  4.7863d0,  5.0119d0,  5.2481d0,  5.4954d0,  5.7544d0,  &
      6.0256d0,  6.3096d0,  6.6069d0,  6.9183d0,  7.2444d0,  7.5858d0,  &
      7.9433d0,  8.3176d0,  8.7096d0,  9.1201d0,  9.5499d0, 10.000d0,   &
     10.471d0,  10.965d0,  11.482d0,  12.023d0,  12.589d0,  13.183d0,   &
     13.804d0,  14.454d0,  15.136d0,  15.849d0,  16.596d0,  17.378d0,   &
     18.197d0,  19.055d0,  19.953d0,  20.893d0,  21.878d0,  22.909d0,   &
     23.988d0,  25.119d0,  26.303d0,  27.542d0,  28.840d0,  30.200d0,   &
     31.623d0,  33.113d0,  34.674d0,  36.308d0,  38.019d0,  39.811d0,   &
     41.687d0,  43.652d0,  45.709d0,  47.863d0,  50.119d0,  52.481d0,   &
     54.954d0,  57.544d0,  60.256d0,  63.096d0,  66.069d0,  69.183d0,   &
     72.444d0,  75.858d0,  79.433d0,  83.176d0,  87.096d0,  91.201d0,   &
     95.499d0, 100.00d0 /)
  ! Efficiency (total)
  REAL*8, PARAMETER :: EFF0(NE)                                         &
      =       (/ 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      4.00000d-6,1.00000d-6,4.00000d-6,2.00000d-6,5.00000d-6,1.20000d-5,&
      1.00000d-5,1.40000d-5,2.60000d-5,3.60000d-5,3.50000d-5,7.60000d-5,&
      8.70000d-5,1.16000d-4,1.76000d-4,2.20000d-4,3.49000d-4,4.83000d-4,&
      7.35000d-4,9.10000d-4,1.33200d-3,1.65400d-3,2.34300d-3,2.92200d-3,&
      4.00500d-3,5.06900d-3,6.38100d-3,7.86100d-3,9.81400d-3,1.27250d-2,&
      1.56120d-2,1.92370d-2,2.36900d-2,2.82510d-2,3.36450d-2,4.02370d-2,&
      4.66420d-2,5.47130d-2,6.41220d-2,7.29530d-2,8.46030d-2,9.60320d-2,&
      1.09300d-1,1.23320d-1,1.39660d-1,1.55620d-1,1.74690d-1,1.94520d-1,&
      2.16200d-1,2.38030d-1,2.61850d-1,2.85090d-1,3.09390d-1,3.29750d-1,&
      3.52810d-1,3.70130d-1,3.84230d-1,3.95210d-1,4.00710d-1,4.03610d-1,&
      4.01040d-1,3.94630d-1,3.84650d-1,3.72750d-1,3.59150d-1,3.43050d-1,&
      3.27380d-1,3.13490d-1,2.98190d-1,2.83900d-1,2.72670d-1,2.60230d-1,&
      2.50270d-1,2.39970d-1,2.30500d-1,2.21070d-1,2.10540d-1,1.98930d-1,&
      1.86480d-1,1.70950d-1,1.53000d-1,1.32510d-1,1.09770d-1,8.64570d-2,&
      6.48630d-2,4.59830d-2,3.03230d-2,1.85540d-2,1.06490d-2,5.61100d-3,&
      2.66100d-3,1.22700d-3,5.37000d-4,2.15000d-4,7.20000d-5,2.70000d-5,&
      6.00000d-6,5.00000d-6,2.00000d-6,0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0 /)
  ! Efficiency (first interval)
  REAL*8, PARAMETER :: EFF1(NE)                                         &
      =       (/ 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      4.00000d-6,1.00000d-6,2.00000d-6,2.00000d-6,2.00000d-6,7.00000d-6,&
      7.00000d-6,1.20000d-5,1.80000d-5,1.90000d-5,1.80000d-5,4.50000d-5,&
      5.00000d-5,6.60000d-5,9.90000d-5,1.31000d-4,2.04000d-4,2.62000d-4,&
      3.86000d-4,4.72000d-4,6.91000d-4,8.25000d-4,1.12100d-3,1.32800d-3,&
      1.81300d-3,2.17300d-3,2.63700d-3,3.16700d-3,4.01900d-3,5.00000d-3,&
      5.84600d-3,7.08400d-3,8.07700d-3,9.52700d-3,1.07220d-2,1.22500d-2,&
      1.38680d-2,1.54030d-2,1.73800d-2,1.91490d-2,2.11010d-2,2.34970d-2,&
      2.55000d-2,2.75090d-2,2.98440d-2,3.21020d-2,3.46350d-2,3.63250d-2,&
      3.87330d-2,3.98790d-2,4.15350d-2,4.20320d-2,4.20200d-2,4.11880d-2,&
      3.96450d-2,3.72390d-2,3.35650d-2,2.97950d-2,2.55980d-2,2.09740d-2,&
      1.67220d-2,1.26450d-2,9.26300d-3,6.46100d-3,4.18800d-3,2.65400d-3,&
      1.49700d-3,9.08000d-4,4.94000d-4,2.52000d-4,1.15000d-4,5.30000d-5,&
      1.80000d-5,1.00000d-5,3.00000d-6,2.00000d-6,0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0 /)
  ! Efficiency (second interval)
  REAL*8, PARAMETER :: EFF2(NE)                                         &
     =        (/ 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 2.00000d-6,0.00000d0, 3.00000d-6,5.00000d-6,&
      2.00000d-6,2.00000d-6,7.00000d-6,1.60000d-5,1.30000d-5,2.40000d-5,&
      2.70000d-5,4.40000d-5,6.60000d-5,7.10000d-5,1.23000d-4,1.81000d-4,&
      2.80000d-4,3.18000d-4,4.81000d-4,6.01000d-4,8.92000d-4,1.12500d-3,&
      1.57100d-3,1.97900d-3,2.53400d-3,3.11900d-3,3.76600d-3,4.85500d-3,&
      5.99400d-3,7.26200d-3,8.97600d-3,1.06020d-2,1.24910d-2,1.47470d-2,&
      1.67300d-2,1.92640d-2,2.21710d-2,2.49040d-2,2.81090d-2,3.10880d-2,&
      3.45400d-2,3.84720d-2,4.16550d-2,4.56960d-2,4.96820d-2,5.38700d-2,&
      5.85350d-2,6.20200d-2,6.60030d-2,6.81040d-2,7.06770d-2,7.07810d-2,&
      7.09020d-2,6.90490d-2,6.44180d-2,5.90480d-2,5.28870d-2,4.62860d-2,&
      3.85470d-2,3.06790d-2,2.35860d-2,1.74240d-2,1.25220d-2,8.14600d-3,&
      5.31700d-3,3.32900d-3,1.80300d-3,1.01000d-3,5.14000d-4,2.38000d-4,&
      1.05000d-4,5.40000d-5,1.80000d-5,7.00000d-6,3.00000d-6,0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0 /)
  ! Efficiency (third interval)
  REAL*8, PARAMETER :: EFF3(NE)                                         &
      =       (/ 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      1.00000d-6,0.00000d0, 1.00000d-6,1.00000d-6,4.00000d-6,7.00000d-6,&
      1.00000d-5,6.00000d-6,1.10000d-5,1.80000d-5,2.20000d-5,4.00000d-5,&
      6.90000d-5,1.20000d-4,1.60000d-4,2.28000d-4,3.30000d-4,4.69000d-4,&
      6.21000d-4,9.17000d-4,1.21000d-3,1.57500d-3,2.02900d-3,2.87000d-3,&
      3.77200d-3,4.89100d-3,6.63700d-3,8.12200d-3,1.04320d-2,1.32400d-2,&
      1.60440d-2,2.00460d-2,2.45710d-2,2.89000d-2,3.53930d-2,4.14470d-2,&
      4.92590d-2,5.73430d-2,6.81650d-2,7.78210d-2,9.03720d-2,1.04320d-1,&
      1.18930d-1,1.36130d-1,1.54310d-1,1.74960d-1,1.96690d-1,2.17780d-1,&
      2.42270d-1,2.63840d-1,2.86240d-1,3.06370d-1,3.22220d-1,3.36350d-1,&
      3.45770d-1,3.51300d-1,3.51800d-1,3.48860d-1,3.42440d-1,3.32250d-1,&
      3.20570d-1,3.09260d-1,2.95900d-1,2.82640d-1,2.72040d-1,2.59930d-1,&
      2.50150d-1,2.39910d-1,2.30480d-1,2.21060d-1,2.10540d-1,1.98930d-1,&
      1.86480d-1,1.70950d-1,1.53000d-1,1.32510d-1,1.09770d-1,8.64570d-2,&
      6.48630d-2,4.59830d-2,3.03230d-2,1.85540d-2,1.06490d-2,5.61100d-3,&
      2.66100d-3,1.22700d-3,5.37000d-4,2.15000d-4,7.20000d-5,2.70000d-5,&
      6.00000d-6,5.00000d-6,2.00000d-6,0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0 /)
  ! Efficiencies array (2D)
  REAL*8, PARAMETER :: EFF(NE,0:NEFF)                                   &
      = RESHAPE( (/ EFF0(:), EFF1(:), EFF2(:), EFF3(:) /) ,SHAPE(EFF))
  
  ! One call for all settings.
  ! Most of these _must_ be there to ensure everything get initialized.
  CALL SetDetector(D,mass=34d0,time=224.6d0,Nevents=2,                  &
                   background=1.0d0,Nelem=1,Zelem=(/54/),               &
                   NEeff=NE,Eeff=E,Neff=NEFF,eff=EFF,                   &
                   intervals=intervals)
  D%eff_file = '[XENON100 2012]'
  
END SUBROUTINE



!=======================================================================
! LUX 2013 ANALYSIS ROUTINES
! Based upon the LUX 2013 analysis [1310.8214].  One event seen in
! the analysis region.
! 
! BIND() is used to specify compiler-independent object file symbol
! names to allow for easier interfacing with C/C++.  
!=======================================================================

!-----------------------------------------------------------------------
! Initializes the module to perform calculations for the LUX 2013
! analysis.  This must be called if any of the following LUX 2013
! routines are to be used.
! 
! Required input arguments:
!     intervals   Indicates if sub-intervals should be included.
!                 Only necessary if confidence intervals using the
!                 maximum gap method are desired.
! 
SUBROUTINE LUX_2013_Init(intervals)
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: intervals
  CALL LUX_2013_InitTo(LUX_2013,intervals)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_LUX_2013_Init(intervals) &
           BIND(C,NAME='C_DDCALC0_lux_2013_init')
  USE ISO_C_BINDING, only: C_BOOL
  IMPLICIT NONE
  LOGICAL(KIND=C_BOOL), INTENT(IN) :: intervals
  CALL LUX_2013_Init(LOGICAL(intervals))
END SUBROUTINE


! ----------------------------------------------------------------------
! Sets the minimum recoil energy to be included in the calculations.
! Note the efficiency curves already account for detector and analysis
! thresholds regardless of this setting, so setting this to 0 keV (the
! default behavior when initialization is performed) does not imply
! that very low energy recoils actually contribute to the signal.
! 
! Required input arguments:
!     Emin        The minimum recoil energy to consider [keV]
! 
SUBROUTINE LUX_2013_SetEmin(Emin)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: Emin
  CALL SetDetector(LUX_2013,Emin=Emin)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_LUX_2013_SetEmin(Emin) &
           BIND(C,NAME='C_DDCALC0_lux_2013_setemin')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: Emin
  CALL LUX_2013_SetEmin(REAL(Emin,KIND=8))
END SUBROUTINE


! ----------------------------------------------------------------------
! Calculates various rate quantities using the current WIMP.
! Must be called each time the WIMP parameters are modified.
! 
SUBROUTINE LUX_2013_CalcRates()
  IMPLICIT NONE
  CALL CalcRates(LUX_2013)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_LUX_2013_CalcRates() &
           BIND(C,NAME='C_DDCALC0_lux_2013_calcrates')
  IMPLICIT NONE
  CALL LUX_2013_CalcRates()
END SUBROUTINE


! ----------------------------------------------------------------------
! Returns the observed number of events.
! 
FUNCTION LUX_2013_Events() RESULT(N)
  IMPLICIT NONE
  INTEGER :: N
  CALL GetRates(LUX_2013,Nevents=N)
END FUNCTION

! C++ interface wrapper
FUNCTION C_LUX_2013_Events() RESULT(N) &
         BIND(C,NAME='C_DDCALC0_lux_2013_events')
  USE ISO_C_BINDING, only: C_INT
  IMPLICIT NONE
  INTEGER(KIND=C_INT) :: N
  ! Automatic type conversions here
  N = LUX_2013_Events()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of background events.
! 
FUNCTION LUX_2013_Background() RESULT(b)
  IMPLICIT NONE
  REAL*8 :: b
  CALL GetRates(LUX_2013,background=b)
END FUNCTION

! C++ interface wrapper
FUNCTION C_LUX_2013_Background() RESULT(b) &
         BIND(C,NAME='C_DDCALC0_lux_2013_background')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: b
  ! Automatic type conversions here
  b = LUX_2013_Background()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of signal events for the
! current WIMP.
! 
FUNCTION LUX_2013_Signal() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(LUX_2013,signal=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_LUX_2013_Signal() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_lux_2013_signal')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = LUX_2013_Signal()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of spin-independent signal events
! for the current WIMP.
! 
FUNCTION LUX_2013_SignalSI() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(LUX_2013,signal_si=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_LUX_2013_SignalSI() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_lux_2013_signalsi')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = LUX_2013_SignalSI()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of spin-dependent signal events
! for the current WIMP.
! 
FUNCTION LUX_2013_SignalSD() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(LUX_2013,signal_sd=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_LUX_2013_SignalSD() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_lux_2013_signalsd')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = LUX_2013_SignalSD()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the log-likelihood for the current WIMP mass and couplings.
! Uses a Poisson distribution in the number of observed events N:
!    P(N|s+b)
! where s is the average expected signal and b is the average expected
! background.
! 
FUNCTION LUX_2013_LogLikelihood() RESULT(lnlike)
  IMPLICIT NONE
  REAL*8 :: lnlike
  lnlike = LogLikelihood(LUX_2013)
END FUNCTION

! C++ interface wrapper
FUNCTION C_LUX_2013_LogLikelihood() RESULT(lnlike) &
         BIND(C,NAME='C_DDCALC0_lux_2013_loglikelihood')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: lnlike
  ! Automatic type conversions here
  lnlike = LUX_2013_LogLikelihood()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the log of the p-value for the current WIMP mass and
! couplings (NO BACKGROUND SUBTRACTION).  Uses the maximum gap method
! if LUX_2013_Init was called with argument intervals=.TRUE.,
! otherwise uses a Poisson distribution in the number of observed
! events N:
!    P(N|s)
! where s is the average expected signal (background contributions are
! ignored).
! 
FUNCTION LUX_2013_LogPValue() RESULT(lnp)
  IMPLICIT NONE
  REAL*8 :: lnp
  lnp = LogPValue(LUX_2013)
END FUNCTION

! C++ interface wrapper
FUNCTION C_LUX_2013_LogPValue() RESULT(lnp) &
         BIND(C,NAME='C_DDCALC0_lux_2013_logpvalue')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: lnp
  ! Automatic type conversions here
  lnp = LUX_2013_LogPValue()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the factor x by which the cross-sections must be scaled
! (sigma -> x*sigma) to achieve the desired p-value (given as log(p)).
! See LogPValue() above for a description of the statistics.
! 
! Required input argument:
!   lnp         The logarithm of the desired p-value (p = 1-CL).
! 
FUNCTION LUX_2013_ScaleToPValue(lnp) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN) :: lnp
  x = ScaleToPValue(LUX_2013,lnp)
END FUNCTION

! C++ interface wrapper
FUNCTION C_LUX_2013_ScaleToPValue(lnp) RESULT(x) &
         BIND(C,NAME='C_DDCALC0_lux_2013_scaletopvalue')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: x
  REAL(KIND=C_DOUBLE), INTENT(IN) :: lnp
  ! Automatic type conversions here
  x = LUX_2013_ScaleToPValue(REAL(lnp,KIND=8))
END FUNCTION


!-----------------------------------------------------------------------
! INTERNAL ROUTINE.
! Initializes the given DetectorStruct to the LUX 2013 analysis.
! This is meant as an internal routine; external access should be
! through LUX_2013_Init instead.
! 
! The efficiencies used here were generated using TPCMC.
! 
! Required input arguments:
!     D           The DetectorStruct to initialize
!     intervals   Indicates if sub-intervals should be included
! 
SUBROUTINE LUX_2013_InitTo(D,intervals)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(OUT) :: D
  LOGICAL, INTENT(IN) :: intervals
  INTEGER, PARAMETER :: NE = 151
  INTEGER, PARAMETER :: NEFF = 2
  ! Efficiency curves energy tabulation points
  REAL*8, PARAMETER :: E(NE)                                            &
      =       (/ 0.10000d0, 0.10471d0, 0.10965d0, 0.11482d0, 0.12023d0, &
      0.12589d0, 0.13183d0, 0.13804d0, 0.14454d0, 0.15136d0, 0.15849d0, &
      0.16596d0, 0.17378d0, 0.18197d0, 0.19055d0, 0.19953d0, 0.20893d0, &
      0.21878d0, 0.22909d0, 0.23988d0, 0.25119d0, 0.26303d0, 0.27542d0, &
      0.28840d0, 0.30200d0, 0.31623d0, 0.33113d0, 0.34674d0, 0.36308d0, &
      0.38019d0, 0.39811d0, 0.41687d0, 0.43652d0, 0.45709d0, 0.47863d0, &
      0.50119d0, 0.52481d0, 0.54954d0, 0.57544d0, 0.60256d0, 0.63096d0, &
      0.66069d0, 0.69183d0, 0.72444d0, 0.75858d0, 0.79433d0, 0.83176d0, &
      0.87096d0, 0.91201d0, 0.95499d0, 1.0000d0,  1.0471d0,  1.0965d0,  &
      1.1482d0,  1.2023d0,  1.2589d0,  1.3183d0,  1.3804d0,  1.4454d0,  &
      1.5136d0,  1.5849d0,  1.6596d0,  1.7378d0,  1.8197d0,  1.9055d0,  &
      1.9953d0,  2.0893d0,  2.1878d0,  2.2909d0,  2.3988d0,  2.5119d0,  &
      2.6303d0,  2.7542d0,  2.8840d0,  3.0200d0,  3.1623d0,  3.3113d0,  &
      3.4674d0,  3.6308d0,  3.8019d0,  3.9811d0,  4.1687d0,  4.3652d0,  &
      4.5709d0,  4.7863d0,  5.0119d0,  5.2481d0,  5.4954d0,  5.7544d0,  &
      6.0256d0,  6.3096d0,  6.6069d0,  6.9183d0,  7.2444d0,  7.5858d0,  &
      7.9433d0,  8.3176d0,  8.7096d0,  9.1201d0,  9.5499d0, 10.000d0,   &
     10.471d0,  10.965d0,  11.482d0,  12.023d0,  12.589d0,  13.183d0,   &
     13.804d0,  14.454d0,  15.136d0,  15.849d0,  16.596d0,  17.378d0,   &
     18.197d0,  19.055d0,  19.953d0,  20.893d0,  21.878d0,  22.909d0,   &
     23.988d0,  25.119d0,  26.303d0,  27.542d0,  28.840d0,  30.200d0,   &
     31.623d0,  33.113d0,  34.674d0,  36.308d0,  38.019d0,  39.811d0,   &
     41.687d0,  43.652d0,  45.709d0,  47.863d0,  50.119d0,  52.481d0,   &
     54.954d0,  57.544d0,  60.256d0,  63.096d0,  66.069d0,  69.183d0,   &
     72.444d0,  75.858d0,  79.433d0,  83.176d0,  87.096d0,  91.201d0,   &
     95.499d0, 100.00d0 /)
  ! LOWER 50% NR BAND >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Efficiency (total)
  REAL*8, PARAMETER :: EFF0(NE)                                         &
      =       (/ 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      2.00000d-7,2.00000d-7,2.00000d-7,3.00000d-7,6.00000d-7,1.20000d-6,&
      1.40000d-6,1.90000d-6,4.20000d-6,4.80000d-6,6.50000d-6,1.08000d-5,&
      1.74000d-5,2.38000d-5,3.42000d-5,4.58000d-5,6.88000d-5,9.74000d-5,&
      1.37700d-4,1.96200d-4,2.84900d-4,3.97300d-4,5.54700d-4,7.81600d-4,&
      1.08660d-3,1.54220d-3,2.12830d-3,2.95810d-3,4.36510d-3,6.05560d-3,&
      8.07280d-3,1.08272d-2,1.40706d-2,1.84873d-2,2.39837d-2,3.08424d-2,&
      3.95926d-2,4.99352d-2,6.27487d-2,7.66224d-2,9.42524d-2,1.15720d-1,&
      1.40320d-1,1.67930d-1,1.98780d-1,2.32650d-1,2.68550d-1,3.06370d-1,&
      3.45110d-1,3.83720d-1,4.20080d-1,4.55080d-1,4.86810d-1,5.14040d-1,&
      5.37080d-1,5.54520d-1,5.66770d-1,5.73020d-1,5.74600d-1,5.71510d-1,&
      5.64900d-1,5.54790d-1,5.43390d-1,5.31050d-1,5.17910d-1,5.05330d-1,&
      4.94050d-1,4.83350d-1,4.74620d-1,4.67740d-1,4.61250d-1,4.56760d-1,&
      4.52870d-1,4.49950d-1,4.47830d-1,4.46160d-1,4.44270d-1,4.42740d-1,&
      4.41050d-1,4.39550d-1,4.38810d-1,4.37010d-1,4.33900d-1,4.28450d-1,&
      4.18620d-1,4.02840d-1,3.78120d-1,3.42460d-1,2.94850d-1,2.37970d-1,&
      1.77820d-1,1.21260d-1,7.52971d-2,4.19808d-2,2.08447d-2,9.20250d-3,&
      3.55500d-3,1.23880d-3,3.67200d-4,9.76000d-5,2.22000d-5,5.50000d-6,&
      9.00000d-7,0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0 /)
  ! Efficiency (first interval)
  REAL*8, PARAMETER :: EFF1(NE)                                         &
      =       (/ 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      1.00000d-7,2.00000d-7,2.00000d-7,3.00000d-7,6.00000d-7,1.10000d-6,&
      1.10000d-6,1.70000d-6,3.60000d-6,4.50000d-6,5.70000d-6,9.30000d-6,&
      1.54000d-5,2.07000d-5,2.90000d-5,3.92000d-5,5.91000d-5,8.30000d-5,&
      1.17700d-4,1.62500d-4,2.37000d-4,3.26000d-4,4.56700d-4,6.23200d-4,&
      8.61500d-4,1.21160d-3,1.64560d-3,2.25640d-3,3.25450d-3,4.43340d-3,&
      5.82990d-3,7.69370d-3,9.85520d-3,1.27067d-2,1.61984d-2,2.03405d-2,&
      2.56110d-2,3.15482d-2,3.86713d-2,4.64675d-2,5.56741d-2,6.63380d-2,&
      7.77153d-2,8.96590d-2,1.02090d-1,1.14360d-1,1.25880d-1,1.36440d-1,&
      1.45460d-1,1.51910d-1,1.55090d-1,1.55910d-1,1.53090d-1,1.46780d-1,&
      1.37830d-1,1.26400d-1,1.12860d-1,9.78608d-2,8.25226d-2,6.74473d-2,&
      5.35181d-2,4.07813d-2,3.01037d-2,2.13387d-2,1.44853d-2,9.50610d-3,&
      5.88720d-3,3.47590d-3,1.97510d-3,1.07430d-3,5.52700d-4,2.55800d-4,&
      1.12100d-4,5.05000d-5,2.12000d-5,6.90000d-6,2.10000d-6,8.00000d-7,&
      0.00000d0, 2.00000d-7,0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0 /)
  ! Efficiency (second interval)
  REAL*8, PARAMETER :: EFF2(NE)                                         &
      =       (/ 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      1.00000d-7,0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 1.00000d-7,&
      3.00000d-7,2.00000d-7,6.00000d-7,3.00000d-7,8.00000d-7,1.50000d-6,&
      2.00000d-6,3.10000d-6,5.20000d-6,6.60000d-6,9.70000d-6,1.44000d-5,&
      2.00000d-5,3.37000d-5,4.79000d-5,7.13000d-5,9.80000d-5,1.58400d-4,&
      2.25100d-4,3.30600d-4,4.82700d-4,7.01700d-4,1.11060d-3,1.62220d-3,&
      2.24290d-3,3.13350d-3,4.21540d-3,5.78060d-3,7.78530d-3,1.05019d-2,&
      1.39816d-2,1.83870d-2,2.40774d-2,3.01549d-2,3.85783d-2,4.93773d-2,&
      6.26065d-2,7.82662d-2,9.66959d-2,1.18290d-1,1.42670d-1,1.69920d-1,&
      1.99650d-1,2.31820d-1,2.64980d-1,2.99170d-1,3.33710d-1,3.67250d-1,&
      3.99250d-1,4.28120d-1,4.53910d-1,4.75160d-1,4.92080d-1,5.04070d-1,&
      5.11390d-1,5.14010d-1,5.13290d-1,5.09720d-1,5.03430d-1,4.95830d-1,&
      4.88170d-1,4.79880d-1,4.72650d-1,4.66660d-1,4.60700d-1,4.56500d-1,&
      4.52750d-1,4.49900d-1,4.47800d-1,4.46160d-1,4.44260d-1,4.42740d-1,&
      4.41050d-1,4.39550d-1,4.38810d-1,4.37010d-1,4.33900d-1,4.28450d-1,&
      4.18620d-1,4.02840d-1,3.78120d-1,3.42460d-1,2.94850d-1,2.37970d-1,&
      1.77820d-1,1.21260d-1,7.52971d-2,4.19808d-2,2.08447d-2,9.20250d-3,&
      3.55500d-3,1.23880d-3,3.67200d-4,9.76000d-5,2.22000d-5,5.50000d-6,&
      9.00000d-7,0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0 /)
  ! END: LOWER 50% NR BAND <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ! BEGIN 10-50% NR BAND >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! Efficiency (total)
  !REAL*8, PARAMETER :: EFF0(NE)                                         &
  !    =       (/ 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 2.00000d-7,1.00000d-7,2.00000d-7,5.00000d-7,5.00000d-7,&
  !    7.00000d-7,7.00000d-7,1.20000d-6,1.80000d-6,2.30000d-6,3.90000d-6,&
  !    5.40000d-6,8.50000d-6,1.23000d-5,1.70000d-5,2.38000d-5,3.31000d-5,&
  !    4.61000d-5,5.94000d-5,9.26000d-5,1.21600d-4,1.72500d-4,2.34500d-4,&
  !    3.33300d-4,4.56700d-4,6.38400d-4,8.84300d-4,1.33540d-3,1.84720d-3,&
  !    2.51600d-3,3.31170d-3,4.43180d-3,5.93570d-3,7.78490d-3,1.02271d-2,&
  !    1.33162d-2,1.72063d-2,2.22009d-2,2.73286d-2,3.45888d-2,4.36821d-2,&
  !    5.45799d-2,6.77302d-2,8.32638d-2,1.00680d-1,1.21170d-1,1.43630d-1,&
  !    1.68140d-1,1.94250d-1,2.21300d-1,2.48970d-1,2.76390d-1,3.02480d-1,&
  !    3.26300d-1,3.48160d-1,3.66100d-1,3.79720d-1,3.89880d-1,3.95620d-1,&
  !    3.97780d-1,3.96640d-1,3.92580d-1,3.87090d-1,3.79880d-1,3.73050d-1,&
  !    3.65850d-1,3.58740d-1,3.53390d-1,3.48230d-1,3.43970d-1,3.41220d-1,&
  !    3.38810d-1,3.37650d-1,3.36940d-1,3.36130d-1,3.35760d-1,3.34690d-1,&
  !    3.33280d-1,3.32160d-1,3.31870d-1,3.32150d-1,3.32200d-1,3.31440d-1,&
  !    3.28070d-1,3.20960d-1,3.07340d-1,2.84620d-1,2.50940d-1,2.07520d-1,&
  !    1.58600d-1,1.10370d-1,6.96978d-2,3.92375d-2,1.98434d-2,8.79110d-3,&
  !    3.44320d-3,1.18140d-3,3.65300d-4,9.83000d-5,2.27000d-5,5.00000d-6,&
  !    3.00000d-7,0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0 /)
  ! Efficiency (first interval)
  !REAL*8, PARAMETER :: EFF1(NE)                                         &
  !    =       (/ 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 2.00000d-7,1.00000d-7,2.00000d-7,5.00000d-7,5.00000d-7,&
  !    6.00000d-7,7.00000d-7,1.20000d-6,1.80000d-6,2.30000d-6,3.80000d-6,&
  !    5.20000d-6,7.90000d-6,1.18000d-5,1.66000d-5,2.28000d-5,3.16000d-5,&
  !    4.33000d-5,5.57000d-5,8.52000d-5,1.12700d-4,1.59900d-4,2.15100d-4,&
  !    3.04800d-4,4.11400d-4,5.74500d-4,7.88300d-4,1.16420d-3,1.59150d-3,&
  !    2.16230d-3,2.78910d-3,3.71100d-3,4.89860d-3,6.29290d-3,8.15860d-3,&
  !    1.03830d-2,1.31709d-2,1.66408d-2,2.02664d-2,2.50647d-2,3.07883d-2,&
  !    3.72727d-2,4.46930d-2,5.28245d-2,6.11301d-2,7.03714d-2,7.93184d-2,&
  !    8.75027d-2,9.49859d-2,1.01120d-1,1.05340d-1,1.07160d-1,1.06570d-1,&
  !    1.03200d-1,9.76218d-2,8.97488d-2,7.98012d-2,6.90608d-2,5.76738d-2,&
  !    4.66902d-2,3.63217d-2,2.71772d-2,1.96357d-2,1.34951d-2,8.86060d-3,&
  !    5.57630d-3,3.35120d-3,1.90010d-3,1.02690d-3,5.29500d-4,2.54600d-4,&
  !    1.14600d-4,4.86000d-5,1.93000d-5,6.70000d-6,2.90000d-6,1.40000d-6,&
  !    0.00000d0, 1.00000d-7,0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0 /)
  ! Efficiency (second interval)
  !REAL*8, PARAMETER :: EFF2(NE)                                         &
  !    =       (/ 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    1.00000d-7,0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 1.00000d-7,&
  !    2.00000d-7,6.00000d-7,5.00000d-7,4.00000d-7,1.00000d-6,1.50000d-6,&
  !    2.80000d-6,3.70000d-6,7.40000d-6,8.90000d-6,1.26000d-5,1.94000d-5,&
  !    2.85000d-5,4.53000d-5,6.39000d-5,9.60000d-5,1.71200d-4,2.55700d-4,&
  !    3.53700d-4,5.22600d-4,7.20800d-4,1.03710d-3,1.49200d-3,2.06850d-3,&
  !    2.93320d-3,4.03540d-3,5.56010d-3,7.06220d-3,9.52410d-3,1.28938d-2,&
  !    1.73072d-2,2.30372d-2,3.04393d-2,3.95450d-2,5.07938d-2,6.43164d-2,&
  !    8.06401d-2,9.92666d-2,1.20180d-1,1.43630d-1,1.69230d-1,1.95910d-1,&
  !    2.23100d-1,2.50540d-1,2.76350d-1,2.99920d-1,3.20820d-1,3.37950d-1,&
  !    3.51090d-1,3.60320d-1,3.65400d-1,3.67460d-1,3.66390d-1,3.64190d-1,&
  !    3.60270d-1,3.55390d-1,3.51490d-1,3.47210d-1,3.43440d-1,3.40970d-1,&
  !    3.38690d-1,3.37600d-1,3.36920d-1,3.36130d-1,3.35750d-1,3.34690d-1,&
  !    3.33280d-1,3.32160d-1,3.31870d-1,3.32150d-1,3.32200d-1,3.31440d-1,&
  !    3.28070d-1,3.20960d-1,3.07340d-1,2.84620d-1,2.50940d-1,2.07520d-1,&
  !    1.58600d-1,1.10370d-1,6.96978d-2,3.92375d-2,1.98434d-2,8.79110d-3,&
  !    3.44320d-3,1.18140d-3,3.65300d-4,9.83000d-5,2.27000d-5,5.00000d-6,&
  !    3.00000d-7,0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
  !    0.00000d0, 0.00000d0 /)
  ! END: 10-50% NR BAND <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ! Efficiencies array (2D)
  REAL*8, PARAMETER :: EFF(NE,0:NEFF)                                   &
      = RESHAPE( (/ EFF0(:), EFF1(:), EFF2(:) /) ,SHAPE(EFF))
  
  ! One call for all settings.
  ! Most of these _must_ be there to ensure everything get initialized.
  CALL SetDetector(D,mass=118d0,time=85.3d0,Nevents=1,                  &
                   background=0.64d0,Nelem=1,Zelem=(/54/),              &
                   NEeff=NE,Eeff=E,Neff=NEFF,eff=EFF,                   &
                   intervals=intervals)
  D%eff_file = '[LUX 2013]'
  
END SUBROUTINE



!=======================================================================
! SuperCDMS 2014 ANALYSIS ROUTINES
! added by LLH, Jan. 2015
! Based upon the SuperCDMS low-mass WIMP search: PRL 112, 241302 (2014)
! [1402.7137].  11 candidate events seen in signal region (consistent
! with background).
! 
! BIND() is used to specify compiler-independent object file symbol
! names to allow for easier interfacing with C/C++.  
!=======================================================================

!-----------------------------------------------------------------------
! Initializes the module to perform calculations for the SuperCDMS 2014
! analysis.  This must be called if any of the following SuperCDMS 2014
! routines are to be used.
! 
! Required input arguments:
!     intervals   Indicates if sub-intervals should be included.
!                 Only necessary if confidence intervals using the
!                 maximum gap method are desired.
! 
SUBROUTINE SuperCDMS_2014_Init(intervals)
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: intervals
  CALL SuperCDMS_2014_InitTo(SuperCDMS_2014,intervals)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_SuperCDMS_2014_Init(intervals) &
           BIND(C,NAME='C_DDCALC0_supercdms_2014_init')
  USE ISO_C_BINDING, only: C_BOOL
  IMPLICIT NONE
  LOGICAL(KIND=C_BOOL), INTENT(IN) :: intervals
  CALL SuperCDMS_2014_Init(LOGICAL(intervals))
END SUBROUTINE


! ----------------------------------------------------------------------
! Sets the minimum recoil energy to be included in the calculations.
! Note the efficiency curves already account for detector and analysis
! thresholds regardless of this setting, so setting this to 0 keV (the
! default behavior when initialization is performed) does not imply
! that very low energy recoils actually contribute to the signal.
! 
! Required input arguments:
!     Emin        The minimum recoil energy to consider [keV]
! 
SUBROUTINE SuperCDMS_2014_SetEmin(Emin)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: Emin
  CALL SetDetector(SuperCDMS_2014,Emin=Emin)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_SuperCDMS_2014_SetEmin(Emin) &
           BIND(C,NAME='C_DDCALC0_supercdms_2014_setemin')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: Emin
  CALL SuperCDMS_2014_SetEmin(REAL(Emin,KIND=8))
END SUBROUTINE


! ----------------------------------------------------------------------
! Calculates various rate quantities using the current WIMP.
! Must be called each time the WIMP parameters are modified.
! 
SUBROUTINE SuperCDMS_2014_CalcRates()
  IMPLICIT NONE
  CALL CalcRates(SuperCDMS_2014)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_SuperCDMS_2014_CalcRates() &
           BIND(C,NAME='C_DDCALC0_supercdms_2014_calcrates')
  IMPLICIT NONE
  CALL SuperCDMS_2014_CalcRates()
END SUBROUTINE


! ----------------------------------------------------------------------
! Returns the observed number of events.
! 
FUNCTION SuperCDMS_2014_Events() RESULT(N)
  IMPLICIT NONE
  INTEGER :: N
  CALL GetRates(SuperCDMS_2014,Nevents=N)
END FUNCTION

! C++ interface wrapper
FUNCTION C_SuperCDMS_2014_Events() RESULT(N) &
         BIND(C,NAME='C_DDCALC0_supercdms_2014_events')
  USE ISO_C_BINDING, only: C_INT
  IMPLICIT NONE
  INTEGER(KIND=C_INT) :: N
  ! Automatic type conversions here
  N = SuperCDMS_2014_Events()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of background events.
! 
FUNCTION SuperCDMS_2014_Background() RESULT(b)
  IMPLICIT NONE
  REAL*8 :: b
  CALL GetRates(SuperCDMS_2014,background=b)
END FUNCTION

! C++ interface wrapper
FUNCTION C_SuperCDMS_2014_Background() RESULT(b) &
         BIND(C,NAME='C_DDCALC0_supercdms_2014_background')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: b
  ! Automatic type conversions here
  b = SuperCDMS_2014_Background()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of signal events for the
! current WIMP.
! 
FUNCTION SuperCDMS_2014_Signal() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(SuperCDMS_2014,signal=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_SuperCDMS_2014_Signal() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_supercdms_2014_signal')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = SuperCDMS_2014_Signal()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of spin-independent signal events
! for the current WIMP.
! 
FUNCTION SuperCDMS_2014_SignalSI() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(SuperCDMS_2014,signal_si=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_SuperCDMS_2014_SignalSI() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_supercdms_2014_signalsi')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = SuperCDMS_2014_SignalSI()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of spin-dependent signal events
! for the current WIMP.
! 
FUNCTION SuperCDMS_2014_SignalSD() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(SuperCDMS_2014,signal_sd=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_SuperCDMS_2014_SignalSD() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_supercdms_2014_signalsd')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = SuperCDMS_2014_SignalSD()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the log-likelihood for the current WIMP mass and couplings.
! Uses a Poisson distribution in the number of observed events N:
!    P(N|s+b)
! where s is the average expected signal and b is the average expected
! background.
! 
FUNCTION SuperCDMS_2014_LogLikelihood() RESULT(lnlike)
  IMPLICIT NONE
  REAL*8 :: lnlike
  lnlike = LogLikelihood(SuperCDMS_2014)
END FUNCTION

! C++ interface wrapper
FUNCTION C_SuperCDMS_2014_LogLikelihood() RESULT(lnlike) &
         BIND(C,NAME='C_DDCALC0_supercdms_2014_loglikelihood')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: lnlike
  ! Automatic type conversions here
  lnlike = SuperCDMS_2014_LogLikelihood()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the log of the p-value for the current WIMP mass and
! couplings (NO BACKGROUND SUBTRACTION).  Uses the maximum gap method
! if SuperCDMS_2014_Init was called with argument intervals=.TRUE.,
! otherwise uses a Poisson distribution in the number of observed
! events N:
!    P(N|s)
! where s is the average expected signal (background contributions are
! ignored).
! 
FUNCTION SuperCDMS_2014_LogPValue() RESULT(lnp)
  IMPLICIT NONE
  REAL*8 :: lnp
  lnp = LogPValue(SuperCDMS_2014)
END FUNCTION

! C++ interface wrapper
FUNCTION C_SuperCDMS_2014_LogPValue() RESULT(lnp) &
         BIND(C,NAME='C_DDCALC0_supercdms_2014_logpvalue')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: lnp
  ! Automatic type conversions here
  lnp = SuperCDMS_2014_LogPValue()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the factor x by which the cross-sections must be scaled
! (sigma -> x*sigma) to achieve the desired p-value (given as log(p)).
! See LogPValue() above for a description of the statistics.
! 
! Required input argument:
!   lnp         The logarithm of the desired p-value (p = 1-CL).
! 
FUNCTION SuperCDMS_2014_ScaleToPValue(lnp) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN) :: lnp
  x = ScaleToPValue(SuperCDMS_2014,lnp)
END FUNCTION

! C++ interface wrapper
FUNCTION C_SuperCDMS_2014_ScaleToPValue(lnp) RESULT(x) &
         BIND(C,NAME='C_DDCALC0_supercdms_2014_scaletopvalue')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: x
  REAL(KIND=C_DOUBLE), INTENT(IN) :: lnp
  ! Automatic type conversions here
  x = SuperCDMS_2014_ScaleToPValue(REAL(lnp,KIND=8))
END FUNCTION


!-----------------------------------------------------------------------
! INTERNAL ROUTINE.
! Initializes the given DetectorStruct to the SuperCDMS 2014 analysis.
! This is meant as an internal routine; external access should be
! through SuperCDMS_2014_Init instead.
! 
! Required input arguments:
!     D           The DetectorStruct to initialize
!     intervals   Indicates if sub-intervals should be included
! 
SUBROUTINE SuperCDMS_2014_InitTo(D,intervals)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(OUT) :: D
  LOGICAL, INTENT(IN) :: intervals
  INTEGER :: K
  INTEGER, PARAMETER :: NE = 1112
  ! Efficiency curves energy tabulation points
  ! NOTE: Converted from phonon energies using Lindhard
  REAL*8, PARAMETER :: E(NE)                                            &
      =       (/ 1.60867d0, 1.61651d0, 1.62434d0, 1.63218d0, 1.64002d0, &
      1.64785d0, 1.65568d0, 1.66351d0, 1.67134d0, 1.67917d0, 1.68700d0, &
      1.69483d0, 1.70265d0, 1.71048d0, 1.71830d0, 1.72613d0, 1.73395d0, &
      1.74177d0, 1.74959d0, 1.75740d0, 1.76522d0, 1.77304d0, 1.78085d0, &
      1.78867d0, 1.79648d0, 1.80429d0, 1.81210d0, 1.81991d0, 1.82772d0, &
      1.83553d0, 1.84333d0, 1.85114d0, 1.85895d0, 1.86675d0, 1.87455d0, &
      1.88235d0, 1.89015d0, 1.89795d0, 1.90575d0, 1.91355d0, 1.92135d0, &
      1.92914d0, 1.93694d0, 1.94473d0, 1.95253d0, 1.96032d0, 1.96811d0, &
      1.97590d0, 1.98369d0, 1.99148d0, 1.99927d0, 2.00705d0, 2.01484d0, &
      2.02263d0, 2.03041d0, 2.03819d0, 2.04598d0, 2.05376d0, 2.06154d0, &
      2.06932d0, 2.07710d0, 2.08488d0, 2.09265d0, 2.10043d0, 2.10821d0, &
      2.11598d0, 2.12376d0, 2.13153d0, 2.13930d0, 2.14707d0, 2.15484d0, &
      2.16261d0, 2.17038d0, 2.17815d0, 2.18592d0, 2.19368d0, 2.20145d0, &
      2.20921d0, 2.21698d0, 2.22474d0, 2.23250d0, 2.24027d0, 2.24803d0, &
      2.25579d0, 2.26355d0, 2.27130d0, 2.27906d0, 2.28682d0, 2.29457d0, &
      2.30233d0, 2.31008d0, 2.31784d0, 2.32559d0, 2.33334d0, 2.34109d0, &
      2.34884d0, 2.35659d0, 2.36434d0, 2.37209d0, 2.37984d0, 2.38758d0, &
      2.39533d0, 2.40307d0, 2.41082d0, 2.41856d0, 2.42630d0, 2.43404d0, &
      2.44179d0, 2.44953d0, 2.45727d0, 2.46500d0, 2.47274d0, 2.48048d0, &
      2.48822d0, 2.49595d0, 2.50369d0, 2.51142d0, 2.51915d0, 2.52689d0, &
      2.53462d0, 2.54235d0, 2.55008d0, 2.55781d0, 2.56554d0, 2.57327d0, &
      2.58100d0, 2.58872d0, 2.59645d0, 2.60417d0, 2.61190d0, 2.61962d0, &
      2.62735d0, 2.63507d0, 2.64279d0, 2.65051d0, 2.65823d0, 2.66595d0, &
      2.67367d0, 2.68139d0, 2.68911d0, 2.69682d0, 2.70454d0, 2.71225d0, &
      2.71997d0, 2.72768d0, 2.73540d0, 2.74311d0, 2.75082d0, 2.75853d0, &
      2.76624d0, 2.77395d0, 2.78166d0, 2.78937d0, 2.79708d0, 2.80479d0, &
      2.81249d0, 2.82020d0, 2.82790d0, 2.83561d0, 2.84331d0, 2.85101d0, &
      2.85872d0, 2.86642d0, 2.87412d0, 2.88182d0, 2.88952d0, 2.89722d0, &
      2.90492d0, 2.91262d0, 2.92031d0, 2.92801d0, 2.93571d0, 2.94340d0, &
      2.95110d0, 2.95879d0, 2.96648d0, 2.97418d0, 2.98187d0, 2.98956d0, &
      2.99725d0, 3.00494d0, 3.01263d0, 3.02032d0, 3.02801d0, 3.03570d0, &
      3.04338d0, 3.05107d0, 3.05876d0, 3.06644d0, 3.07413d0, 3.08181d0, &
      3.08949d0, 3.09718d0, 3.10486d0, 3.11254d0, 3.12022d0, 3.12790d0, &
      3.13558d0, 3.14326d0, 3.15094d0, 3.15862d0, 3.16629d0, 3.17397d0, &
      3.18165d0, 3.18932d0, 3.19700d0, 3.20467d0, 3.21235d0, 3.22002d0, &
      3.22769d0, 3.23536d0, 3.24303d0, 3.25071d0, 3.25838d0, 3.26605d0, &
      3.27371d0, 3.28138d0, 3.28905d0, 3.29672d0, 3.30438d0, 3.31205d0, &
      3.31972d0, 3.32738d0, 3.33505d0, 3.34271d0, 3.35037d0, 3.35804d0, &
      3.36570d0, 3.37336d0, 3.38102d0, 3.38868d0, 3.39634d0, 3.40400d0, &
      3.41166d0, 3.41932d0, 3.42697d0, 3.43463d0, 3.44229d0, 3.44994d0, &
      3.45760d0, 3.46525d0, 3.47291d0, 3.48056d0, 3.48822d0, 3.49587d0, &
      3.50352d0, 3.51117d0, 3.51882d0, 3.52647d0, 3.53412d0, 3.54177d0, &
      3.54942d0, 3.55707d0, 3.56472d0, 3.57236d0, 3.58001d0, 3.58766d0, &
      3.59530d0, 3.60295d0, 3.61059d0, 3.61824d0, 3.62588d0, 3.63352d0, &
      3.64117d0, 3.64881d0, 3.65645d0, 3.66409d0, 3.67173d0, 3.67937d0, &
      3.68701d0, 3.69465d0, 3.70228d0, 3.70992d0, 3.71756d0, 3.72520d0, &
      3.73283d0, 3.74047d0, 3.74810d0, 3.75574d0, 3.76337d0, 3.77100d0, &
      3.77864d0, 3.78627d0, 3.79390d0, 3.80153d0, 3.80916d0, 3.81679d0, &
      3.82442d0, 3.83205d0, 3.83968d0, 3.84731d0, 3.85494d0, 3.86257d0, &
      3.87019d0, 3.87782d0, 3.88545d0, 3.89307d0, 3.90070d0, 3.90832d0, &
      3.91594d0, 3.92357d0, 3.93119d0, 3.93881d0, 3.94643d0, 3.95406d0, &
      3.96168d0, 3.96930d0, 3.97692d0, 3.98454d0, 3.99216d0, 3.99977d0, &
      4.00739d0, 4.01501d0, 4.02263d0, 4.03024d0, 4.03786d0, 4.04547d0, &
      4.05309d0, 4.06070d0, 4.06832d0, 4.07593d0, 4.08355d0, 4.09116d0, &
      4.09877d0, 4.10638d0, 4.11399d0, 4.12160d0, 4.12921d0, 4.13682d0, &
      4.14443d0, 4.15204d0, 4.15965d0, 4.16726d0, 4.17487d0, 4.18247d0, &
      4.19008d0, 4.19769d0, 4.20529d0, 4.21290d0, 4.22050d0, 4.22811d0, &
      4.23571d0, 4.24331d0, 4.25091d0, 4.25852d0, 4.26612d0, 4.27372d0, &
      4.28132d0, 4.28892d0, 4.29652d0, 4.30412d0, 4.31172d0, 4.31932d0, &
      4.32692d0, 4.33452d0, 4.34211d0, 4.34971d0, 4.35731d0, 4.36490d0, &
      4.37250d0, 4.38009d0, 4.38769d0, 4.39528d0, 4.40287d0, 4.41047d0, &
      4.41806d0, 4.42565d0, 4.43325d0, 4.44084d0, 4.44843d0, 4.45602d0, &
      4.46361d0, 4.47120d0, 4.47879d0, 4.48638d0, 4.49396d0, 4.50155d0, &
      4.50914d0, 4.51673d0, 4.52431d0, 4.53190d0, 4.53949d0, 4.54707d0, &
      4.55466d0, 4.56224d0, 4.56982d0, 4.57741d0, 4.58499d0, 4.59257d0, &
      4.60016d0, 4.60774d0, 4.61532d0, 4.62290d0, 4.63048d0, 4.63806d0, &
      4.64564d0, 4.65322d0, 4.66080d0, 4.66838d0, 4.67596d0, 4.68353d0, &
      4.69111d0, 4.69869d0, 4.70626d0, 4.71384d0, 4.72142d0, 4.72899d0, &
      4.73657d0, 4.74414d0, 4.75171d0, 4.75929d0, 4.76686d0, 4.77443d0, &
      4.78200d0, 4.78958d0, 4.79715d0, 4.80472d0, 4.81229d0, 4.81986d0, &
      4.82743d0, 4.83500d0, 4.84257d0, 4.85014d0, 4.85770d0, 4.86527d0, &
      4.87284d0, 4.88041d0, 4.88797d0, 4.89554d0, 4.90310d0, 4.91067d0, &
      4.91823d0, 4.92580d0, 4.93336d0, 4.94093d0, 4.94849d0, 4.95605d0, &
      4.96361d0, 4.97118d0, 4.97874d0, 4.98630d0, 4.99386d0, 5.00142d0, &
      5.00898d0, 5.01654d0, 5.02410d0, 5.03166d0, 5.03922d0, 5.04677d0, &
      5.05433d0, 5.06189d0, 5.06944d0, 5.07700d0, 5.08456d0, 5.09211d0, &
      5.09967d0, 5.10722d0, 5.11478d0, 5.12233d0, 5.12988d0, 5.13744d0, &
      5.14499d0, 5.15254d0, 5.16010d0, 5.16765d0, 5.17520d0, 5.18275d0, &
      5.19030d0, 5.19785d0, 5.20540d0, 5.21295d0, 5.22050d0, 5.22805d0, &
      5.23560d0, 5.24314d0, 5.25069d0, 5.25824d0, 5.26578d0, 5.27333d0, &
      5.28088d0, 5.28842d0, 5.29597d0, 5.30351d0, 5.31106d0, 5.31860d0, &
      5.32614d0, 5.33369d0, 5.34123d0, 5.34877d0, 5.35632d0, 5.36386d0, &
      5.37140d0, 5.37894d0, 5.38648d0, 5.39402d0, 5.40156d0, 5.40910d0, &
      5.41664d0, 5.42418d0, 5.43172d0, 5.43925d0, 5.44679d0, 5.45433d0, &
      5.46187d0, 5.46940d0, 5.47694d0, 5.48447d0, 5.49201d0, 5.49955d0, &
      5.50708d0, 5.51461d0, 5.52215d0, 5.52968d0, 5.53722d0, 5.54475d0, &
      5.55228d0, 5.55981d0, 5.56735d0, 5.57488d0, 5.58241d0, 5.58994d0, &
      5.59747d0, 5.60500d0, 5.61253d0, 5.62006d0, 5.62759d0, 5.63512d0, &
      5.64264d0, 5.65017d0, 5.65770d0, 5.66523d0, 5.67275d0, 5.68028d0, &
      5.68781d0, 5.69533d0, 5.70286d0, 5.71038d0, 5.71791d0, 5.72543d0, &
      5.73295d0, 5.74048d0, 5.74800d0, 5.75552d0, 5.76305d0, 5.77057d0, &
      5.77809d0, 5.78561d0, 5.79313d0, 5.80065d0, 5.80817d0, 5.81569d0, &
      5.82321d0, 5.83073d0, 5.83825d0, 5.84577d0, 5.85329d0, 5.86081d0, &
      5.86833d0, 5.87584d0, 5.88336d0, 5.89088d0, 5.89839d0, 5.90591d0, &
      5.91342d0, 5.92094d0, 5.92845d0, 5.93597d0, 5.94348d0, 5.95100d0, &
      5.95851d0, 5.96602d0, 5.97354d0, 5.98105d0, 5.98856d0, 5.99607d0, &
      6.00358d0, 6.01110d0, 6.01861d0, 6.02612d0, 6.03363d0, 6.04114d0, &
      6.04865d0, 6.05615d0, 6.06366d0, 6.07117d0, 6.07868d0, 6.08619d0, &
      6.09369d0, 6.10120d0, 6.10871d0, 6.11622d0, 6.12372d0, 6.13123d0, &
      6.13873d0, 6.14624d0, 6.15374d0, 6.16125d0, 6.16875d0, 6.17625d0, &
      6.18376d0, 6.19126d0, 6.19876d0, 6.20627d0, 6.21377d0, 6.22127d0, &
      6.22877d0, 6.23627d0, 6.24377d0, 6.25127d0, 6.25877d0, 6.26627d0, &
      6.27377d0, 6.28127d0, 6.28877d0, 6.29627d0, 6.30377d0, 6.31127d0, &
      6.31876d0, 6.32626d0, 6.33376d0, 6.34125d0, 6.34875d0, 6.35625d0, &
      6.36374d0, 6.37124d0, 6.37873d0, 6.38623d0, 6.39372d0, 6.40121d0, &
      6.40871d0, 6.41620d0, 6.42370d0, 6.43119d0, 6.43868d0, 6.44617d0, &
      6.45366d0, 6.46115d0, 6.46865d0, 6.47614d0, 6.48363d0, 6.49112d0, &
      6.49861d0, 6.50610d0, 6.51359d0, 6.52107d0, 6.52856d0, 6.53605d0, &
      6.54354d0, 6.55103d0, 6.55851d0, 6.56600d0, 6.57349d0, 6.58098d0, &
      6.58846d0, 6.59595d0, 6.60343d0, 6.61092d0, 6.61840d0, 6.62589d0, &
      6.63337d0, 6.64085d0, 6.64834d0, 6.65582d0, 6.66330d0, 6.67079d0, &
      6.67827d0, 6.68575d0, 6.69323d0, 6.70071d0, 6.70819d0, 6.71568d0, &
      6.72316d0, 6.73064d0, 6.73812d0, 6.74560d0, 6.75308d0, 6.76055d0, &
      6.76803d0, 6.77551d0, 6.78299d0, 6.79047d0, 6.79794d0, 6.80542d0, &
      6.81290d0, 6.82037d0, 6.82785d0, 6.83533d0, 6.84280d0, 6.85028d0, &
      6.85775d0, 6.86523d0, 6.87270d0, 6.88018d0, 6.88765d0, 6.89512d0, &
      6.90260d0, 6.91007d0, 6.91754d0, 6.92501d0, 6.93249d0, 6.93996d0, &
      6.94743d0, 6.95490d0, 6.96237d0, 6.96984d0, 6.97731d0, 6.98478d0, &
      6.99225d0, 6.99972d0, 7.00719d0, 7.01466d0, 7.02213d0, 7.02959d0, &
      7.03706d0, 7.04453d0, 7.05200d0, 7.05946d0, 7.06693d0, 7.07440d0, &
      7.08186d0, 7.08933d0, 7.09679d0, 7.10426d0, 7.11172d0, 7.11919d0, &
      7.12665d0, 7.13412d0, 7.14158d0, 7.14904d0, 7.15651d0, 7.16397d0, &
      7.17143d0, 7.17889d0, 7.18636d0, 7.19382d0, 7.20128d0, 7.20874d0, &
      7.21620d0, 7.22366d0, 7.23112d0, 7.23858d0, 7.24604d0, 7.25350d0, &
      7.26096d0, 7.26842d0, 7.27588d0, 7.28334d0, 7.29079d0, 7.29825d0, &
      7.30571d0, 7.31317d0, 7.32062d0, 7.32808d0, 7.33554d0, 7.34299d0, &
      7.35045d0, 7.35790d0, 7.36536d0, 7.37281d0, 7.38027d0, 7.38772d0, &
      7.39518d0, 7.40263d0, 7.41008d0, 7.41754d0, 7.42499d0, 7.43244d0, &
      7.43989d0, 7.44735d0, 7.45480d0, 7.46225d0, 7.46970d0, 7.47715d0, &
      7.48460d0, 7.49205d0, 7.49950d0, 7.50695d0, 7.51440d0, 7.52185d0, &
      7.52930d0, 7.53675d0, 7.54420d0, 7.55164d0, 7.55909d0, 7.56654d0, &
      7.57399d0, 7.58143d0, 7.58888d0, 7.59633d0, 7.60377d0, 7.61122d0, &
      7.61866d0, 7.62611d0, 7.63355d0, 7.64100d0, 7.64844d0, 7.65589d0, &
      7.66333d0, 7.67077d0, 7.67822d0, 7.68566d0, 7.69310d0, 7.70055d0, &
      7.70799d0, 7.71543d0, 7.72287d0, 7.73032d0, 7.73776d0, 7.74520d0, &
      7.75264d0, 7.76008d0, 7.76752d0, 7.77496d0, 7.78240d0, 7.78984d0, &
      7.79728d0, 7.80471d0, 7.81215d0, 7.81959d0, 7.82703d0, 7.83447d0, &
      7.84190d0, 7.84934d0, 7.85678d0, 7.86421d0, 7.87165d0, 7.87909d0, &
      7.88652d0, 7.89396d0, 7.90139d0, 7.90883d0, 7.91626d0, 7.92370d0, &
      7.93113d0, 7.93857d0, 7.94600d0, 7.95343d0, 7.96087d0, 7.96830d0, &
      7.97573d0, 7.98316d0, 7.99060d0, 7.99803d0, 8.00546d0, 8.01289d0, &
      8.02032d0, 8.02775d0, 8.03518d0, 8.04261d0, 8.05004d0, 8.05747d0, &
      8.06490d0, 8.07233d0, 8.07976d0, 8.08719d0, 8.09462d0, 8.10205d0, &
      8.10947d0, 8.11690d0, 8.12433d0, 8.13176d0, 8.13918d0, 8.14661d0, &
      8.15404d0, 8.16146d0, 8.16889d0, 8.17631d0, 8.18374d0, 8.19116d0, &
      8.19859d0, 8.20601d0, 8.21344d0, 8.22086d0, 8.22828d0, 8.23571d0, &
      8.24313d0, 8.25056d0, 8.25798d0, 8.26540d0, 8.27282d0, 8.28024d0, &
      8.28767d0, 8.29509d0, 8.30251d0, 8.30993d0, 8.31735d0, 8.32477d0, &
      8.33219d0, 8.33961d0, 8.34703d0, 8.35445d0, 8.36187d0, 8.36929d0, &
      8.37671d0, 8.38412d0, 8.39154d0, 8.39896d0, 8.40638d0, 8.41380d0, &
      8.42121d0, 8.42863d0, 8.43605d0, 8.44346d0, 8.45088d0, 8.45829d0, &
      8.46571d0, 8.47313d0, 8.48054d0, 8.48796d0, 8.49537d0, 8.50278d0, &
      8.51020d0, 8.51761d0, 8.52503d0, 8.53244d0, 8.53985d0, 8.54727d0, &
      8.55468d0, 8.56209d0, 8.56950d0, 8.57691d0, 8.58433d0, 8.59174d0, &
      8.59915d0, 8.60656d0, 8.61397d0, 8.62138d0, 8.62879d0, 8.63620d0, &
      8.64361d0, 8.65102d0, 8.65843d0, 8.66584d0, 8.67325d0, 8.68065d0, &
      8.68806d0, 8.69547d0, 8.70288d0, 8.71029d0, 8.71769d0, 8.72510d0, &
      8.73251d0, 8.73991d0, 8.74732d0, 8.75472d0, 8.76213d0, 8.76954d0, &
      8.77694d0, 8.78435d0, 8.79175d0, 8.79916d0, 8.80656d0, 8.81396d0, &
      8.82137d0, 8.82877d0, 8.83617d0, 8.84358d0, 8.85098d0, 8.85838d0, &
      8.86578d0, 8.87319d0, 8.88059d0, 8.88799d0, 8.89539d0, 8.90279d0, &
      8.91019d0, 8.91759d0, 8.92499d0, 8.93239d0, 8.93979d0, 8.94719d0, &
      8.95459d0, 8.96199d0, 8.96939d0, 8.97679d0, 8.98419d0, 8.99159d0, &
      8.99899d0, 9.00638d0, 9.01378d0, 9.02118d0, 9.02857d0, 9.03597d0, &
      9.04337d0, 9.05076d0, 9.05816d0, 9.06556d0, 9.07295d0, 9.08035d0, &
      9.08774d0, 9.09514d0, 9.10253d0, 9.10993d0, 9.11732d0, 9.12471d0, &
      9.13211d0, 9.13950d0, 9.14689d0, 9.15429d0, 9.16168d0, 9.16907d0, &
      9.17647d0, 9.18386d0, 9.19125d0, 9.19864d0, 9.20603d0, 9.21342d0, &
      9.22081d0, 9.22820d0, 9.23559d0, 9.24299d0, 9.25037d0, 9.25776d0, &
      9.26515d0, 9.27254d0, 9.27993d0, 9.28732d0, 9.29471d0, 9.30210d0, &
      9.30949d0, 9.31687d0, 9.32426d0, 9.33165d0, 9.33904d0, 9.34642d0, &
      9.35381d0, 9.36120d0, 9.36858d0, 9.37597d0, 9.38336d0, 9.39074d0, &
      9.39813d0, 9.40551d0, 9.41290d0, 9.42028d0, 9.42766d0, 9.43505d0, &
      9.44243d0, 9.44982d0, 9.45720d0, 9.46458d0, 9.47197d0, 9.47935d0, &
      9.48673d0, 9.49411d0, 9.50150d0, 9.50888d0, 9.51626d0, 9.52364d0, &
      9.53102d0, 9.53840d0, 9.54578d0, 9.55317d0, 9.56055d0, 9.56793d0, &
      9.57531d0, 9.58269d0, 9.59007d0, 9.59744d0, 9.60482d0, 9.61220d0, &
      9.61958d0, 9.62696d0, 9.63434d0, 9.64172d0, 9.64909d0, 9.65647d0, &
      9.66385d0, 9.67122d0, 9.67860d0, 9.68598d0, 9.69335d0, 9.70073d0, &
      9.70811d0, 9.71548d0, 9.72286d0, 9.73023d0, 9.73761d0, 9.74498d0, &
      9.75236d0, 9.75973d0, 9.76710d0, 9.77448d0, 9.78185d0, 9.78922d0, &
      9.79660d0, 9.80397d0, 9.81134d0, 9.81872d0, 9.82609d0, 9.83346d0, &
      9.84083d0, 9.84820d0, 9.85558d0, 9.86295d0, 9.87032d0, 9.87769d0, &
      9.88506d0, 9.89243d0, 9.89980d0, 9.90717d0, 9.91454d0, 9.92191d0, &
      9.92928d0, 9.93665d0, 9.94401d0, 9.95138d0, 9.95875d0, 9.96612d0, &
      9.97349d0, 9.98086d0, 9.98086d0 /)
  ! Efficiency (total), including T5Z3
  REAL*8, PARAMETER :: EFF0_ALL_TZ(NE)                                  &
      =       (/ 4.21628d-2,4.56586d-2,4.64347d-2,4.72082d-2,4.79792d-2,&
      4.87478d-2,4.95142d-2,5.02783d-2,5.10401d-2,5.17996d-2,5.25568d-2,&
      5.33114d-2,5.41480d-2,5.48988d-2,5.56465d-2,5.63907d-2,5.71493d-2,&
      5.78875d-2,5.86215d-2,5.93508d-2,6.00752d-2,6.07943d-2,6.15076d-2,&
      6.22149d-2,6.29158d-2,6.38929d-2,6.45838d-2,6.52673d-2,6.59429d-2,&
      6.66103d-2,6.72693d-2,6.79195d-2,6.85606d-2,6.91924d-2,6.98146d-2,&
      7.04270d-2,7.10294d-2,7.16216d-2,7.22034d-2,7.27747d-2,7.33354d-2,&
      7.38854d-2,7.44246d-2,7.49529d-2,7.54704d-2,7.59770d-2,7.64727d-2,&
      7.69576d-2,7.74317d-2,7.78950d-2,7.85230d-2,7.84758d-2,7.84249d-2,&
      7.83712d-2,7.83155d-2,7.82586d-2,7.82014d-2,7.81445d-2,7.80887d-2,&
      7.80347d-2,7.79832d-2,7.79347d-2,7.78899d-2,7.78494d-2,7.78136d-2,&
      7.77830d-2,7.77582d-2,7.77395d-2,7.77274d-2,7.77223d-2,7.77245d-2,&
      7.77344d-2,7.77523d-2,7.77785d-2,7.78132d-2,7.80296d-2,7.80833d-2,&
      7.81463d-2,7.82187d-2,7.83007d-2,7.83924d-2,7.84941d-2,7.86057d-2,&
      7.87276d-2,7.88597d-2,7.90021d-2,7.91550d-2,7.93183d-2,7.94922d-2,&
      7.96767d-2,7.98718d-2,8.00776d-2,8.02940d-2,8.05212d-2,8.07590d-2,&
      8.10076d-2,8.12668d-2,8.15368d-2,8.18174d-2,8.21087d-2,8.25698d-2,&
      8.28828d-2,8.32064d-2,8.35404d-2,8.38849d-2,8.42399d-2,8.46051d-2,&
      8.49808d-2,8.53666d-2,8.57627d-2,8.61689d-2,8.65853d-2,8.70116d-2,&
      8.74479d-2,8.78941d-2,8.83502d-2,8.88160d-2,8.92915d-2,8.97766d-2,&
      9.02712d-2,9.07753d-2,9.12888d-2,9.18117d-2,9.23437d-2,9.28849d-2,&
      9.37831d-2,9.43442d-2,9.49143d-2,9.54932d-2,9.60809d-2,9.66772d-2,&
      9.72821d-2,9.78955d-2,9.85173d-2,9.91473d-2,9.97855d-2,1.00432d-1,&
      1.01086d-1,1.01748d-1,1.02418d-1,1.03096d-1,1.03781d-1,1.04473d-1,&
      1.05173d-1,1.05881d-1,1.06595d-1,1.07316d-1,1.08045d-1,1.08780d-1,&
      1.09521d-1,1.10459d-1,1.11214d-1,1.11976d-1,1.12744d-1,1.13518d-1,&
      1.14298d-1,1.15084d-1,1.15876d-1,1.16673d-1,1.17476d-1,1.18284d-1,&
      1.19098d-1,1.19916d-1,1.20740d-1,1.21568d-1,1.22401d-1,1.23239d-1,&
      1.24082d-1,1.24928d-1,1.25779d-1,1.26634d-1,1.27493d-1,1.28356d-1,&
      1.29222d-1,1.30092d-1,1.31117d-1,1.31994d-1,1.32873d-1,1.33756d-1,&
      1.34642d-1,1.35530d-1,1.36421d-1,1.37314d-1,1.38209d-1,1.39106d-1,&
      1.40005d-1,1.40906d-1,1.41808d-1,1.42712d-1,1.43617d-1,1.44523d-1,&
      1.45430d-1,1.46338d-1,1.47246d-1,1.48155d-1,1.49065d-1,1.49974d-1,&
      1.50884d-1,1.51794d-1,1.52703d-1,1.53593d-1,1.54503d-1,1.55413d-1,&
      1.56322d-1,1.57231d-1,1.58139d-1,1.59045d-1,1.59951d-1,1.60856d-1,&
      1.61759d-1,1.62662d-1,1.63562d-1,1.64462d-1,1.65359d-1,1.66255d-1,&
      1.67149d-1,1.68042d-1,1.68932d-1,1.69820d-1,1.70706d-1,1.71590d-1,&
      1.72471d-1,1.73350d-1,1.74226d-1,1.75100d-1,1.76170d-1,1.77038d-1,&
      1.77904d-1,1.78766d-1,1.79626d-1,1.80481d-1,1.81334d-1,1.82183d-1,&
      1.83028d-1,1.83870d-1,1.84707d-1,1.85541d-1,1.86371d-1,1.87196d-1,&
      1.88017d-1,1.88834d-1,1.89646d-1,1.90454d-1,1.91257d-1,1.92055d-1,&
      1.92847d-1,1.93635d-1,1.94418d-1,1.95195d-1,1.95967d-1,1.96544d-1,&
      1.97305d-1,1.98059d-1,1.98808d-1,1.99552d-1,2.00291d-1,2.01024d-1,&
      2.01751d-1,2.02474d-1,2.03192d-1,2.03904d-1,2.04612d-1,2.05314d-1,&
      2.06012d-1,2.06705d-1,2.07394d-1,2.08077d-1,2.08757d-1,2.09431d-1,&
      2.10102d-1,2.10768d-1,2.11430d-1,2.12087d-1,2.12741d-1,2.13390d-1,&
      2.14779d-1,2.15422d-1,2.16061d-1,2.16696d-1,2.17328d-1,2.17955d-1,&
      2.18580d-1,2.19200d-1,2.19817d-1,2.20431d-1,2.21041d-1,2.21649d-1,&
      2.22252d-1,2.22853d-1,2.23451d-1,2.24140d-1,2.24738d-1,2.25333d-1,&
      2.25926d-1,2.26516d-1,2.27103d-1,2.27689d-1,2.28272d-1,2.28853d-1,&
      2.29432d-1,2.30646d-1,2.31222d-1,2.31797d-1,2.32369d-1,2.32940d-1,&
      2.33509d-1,2.34077d-1,2.34644d-1,2.35209d-1,2.35773d-1,2.36336d-1,&
      2.36897d-1,2.37458d-1,2.38018d-1,2.38576d-1,2.39135d-1,2.39692d-1,&
      2.40249d-1,2.40805d-1,2.41361d-1,2.41917d-1,2.42472d-1,2.43028d-1,&
      2.43583d-1,2.44138d-1,2.44654d-1,2.45208d-1,2.45762d-1,2.46317d-1,&
      2.46873d-1,2.47429d-1,2.47985d-1,2.48542d-1,2.49101d-1,2.49660d-1,&
      2.50219d-1,2.50780d-1,2.51343d-1,2.51906d-1,2.52471d-1,2.53037d-1,&
      2.53605d-1,2.54174d-1,2.54745d-1,2.55318d-1,2.55892d-1,2.56469d-1,&
      2.57048d-1,2.57629d-1,2.58212d-1,2.59283d-1,2.59874d-1,2.60466d-1,&
      2.61062d-1,2.68432d-1,2.69150d-1,2.69871d-1,2.70596d-1,2.71325d-1,&
      2.72057d-1,2.72792d-1,2.73531d-1,2.74273d-1,2.75019d-1,2.75768d-1,&
      2.76521d-1,2.77277d-1,2.78036d-1,2.78799d-1,2.79566d-1,2.80335d-1,&
      2.81108d-1,2.81885d-1,2.82665d-1,2.83448d-1,2.84380d-1,2.85170d-1,&
      2.85963d-1,2.86759d-1,2.87559d-1,2.88361d-1,2.89167d-1,2.89977d-1,&
      2.90789d-1,2.91605d-1,2.92423d-1,2.93245d-1,2.94070d-1,2.94898d-1,&
      2.95729d-1,2.96563d-1,2.97400d-1,2.98240d-1,2.99083d-1,2.99929d-1,&
      3.00777d-1,3.01629d-1,3.02483d-1,3.03340d-1,3.04200d-1,3.05048d-1,&
      3.05914d-1,3.06783d-1,3.07655d-1,3.08530d-1,3.09407d-1,3.10286d-1,&
      3.11168d-1,3.12052d-1,3.12939d-1,3.13828d-1,3.14720d-1,3.15614d-1,&
      3.16510d-1,3.17409d-1,3.18310d-1,3.19213d-1,3.20118d-1,3.21025d-1,&
      3.21934d-1,3.22846d-1,3.23759d-1,3.24675d-1,3.25592d-1,3.26512d-1,&
      3.27892d-1,3.28817d-1,3.29744d-1,3.30672d-1,3.31603d-1,3.32535d-1,&
      3.33469d-1,3.34404d-1,3.35341d-1,3.36280d-1,3.37220d-1,3.38162d-1,&
      3.39106d-1,3.40051d-1,3.40997d-1,3.42894d-1,3.43763d-1,3.43845d-1,&
      3.44797d-1,3.45751d-1,3.46706d-1,3.47662d-1,3.48620d-1,3.49578d-1,&
      3.50539d-1,3.51547d-1,3.52510d-1,3.53475d-1,3.54440d-1,3.55407d-1,&
      3.56374d-1,3.57343d-1,3.58312d-1,3.59283d-1,3.60254d-1,3.61226d-1,&
      3.62198d-1,3.63171d-1,3.64144d-1,3.65118d-1,3.66093d-1,3.67068d-1,&
      3.68043d-1,3.69018d-1,3.69994d-1,3.70969d-1,3.71945d-1,3.72920d-1,&
      3.73896d-1,3.74871d-1,3.75274d-1,3.77636d-1,3.78639d-1,3.79643d-1,&
      3.80647d-1,3.81651d-1,3.82654d-1,3.83658d-1,3.84662d-1,3.85665d-1,&
      3.86668d-1,3.87670d-1,3.88672d-1,3.89674d-1,3.90675d-1,3.91675d-1,&
      3.92675d-1,3.93674d-1,3.94672d-1,3.95670d-1,3.96666d-1,3.97661d-1,&
      3.98655d-1,3.99648d-1,4.00640d-1,4.01866d-1,4.02858d-1,4.03849d-1,&
      4.04838d-1,4.05825d-1,4.06811d-1,4.07795d-1,4.08777d-1,4.09758d-1,&
      4.10736d-1,4.11712d-1,4.12687d-1,4.13659d-1,4.14629d-1,4.15597d-1,&
      4.16562d-1,4.17525d-1,4.18485d-1,4.19443d-1,4.20398d-1,4.21350d-1,&
      4.22300d-1,4.23247d-1,4.24190d-1,4.25131d-1,4.27238d-1,4.28175d-1,&
      4.29108d-1,4.30038d-1,4.30965d-1,4.31888d-1,4.32808d-1,4.33724d-1,&
      4.34637d-1,4.35545d-1,4.36450d-1,4.37351d-1,4.38248d-1,4.39140d-1,&
      4.40029d-1,4.40913d-1,4.41793d-1,4.42669d-1,4.43540d-1,4.44407d-1,&
      4.45269d-1,4.46127d-1,4.46979d-1,4.47827d-1,4.48670d-1,4.48208d-1,&
      4.49038d-1,4.49863d-1,4.50684d-1,4.51499d-1,4.52309d-1,4.53114d-1,&
      4.53914d-1,4.54710d-1,4.55500d-1,4.56285d-1,4.57065d-1,4.57840d-1,&
      4.58610d-1,4.59375d-1,4.60135d-1,4.60891d-1,4.61641d-1,4.62386d-1,&
      4.63127d-1,4.63862d-1,4.64593d-1,4.65319d-1,4.66040d-1,4.66756d-1,&
      4.68439d-1,4.69146d-1,4.69849d-1,4.70547d-1,4.71240d-1,4.71928d-1,&
      4.72612d-1,4.73291d-1,4.73965d-1,4.74634d-1,4.75299d-1,4.75959d-1,&
      4.76615d-1,4.77266d-1,4.77912d-1,4.78554d-1,4.79191d-1,4.79824d-1,&
      4.80453d-1,4.81077d-1,4.81696d-1,4.82311d-1,4.82922d-1,4.83529d-1,&
      4.84131d-1,4.84724d-1,4.85320d-1,4.85911d-1,4.86498d-1,4.87081d-1,&
      4.87660d-1,4.88234d-1,4.88805d-1,4.89372d-1,4.89934d-1,4.90493d-1,&
      4.91047d-1,4.91598d-1,4.92145d-1,4.92688d-1,4.93227d-1,4.93762d-1,&
      4.94294d-1,4.94821d-1,4.95345d-1,4.95866d-1,4.96383d-1,4.96896d-1,&
      4.97406d-1,4.97912d-1,4.99263d-1,4.99763d-1,5.00259d-1,5.00752d-1,&
      5.01241d-1,5.01728d-1,5.02211d-1,5.02690d-1,5.03167d-1,5.03641d-1,&
      5.04111d-1,5.04578d-1,5.05043d-1,5.05504d-1,5.05963d-1,5.06418d-1,&
      5.06871d-1,5.07321d-1,5.07768d-1,5.08212d-1,5.08654d-1,5.09093d-1,&
      5.09530d-1,5.09963d-1,5.10395d-1,5.09351d-1,5.09777d-1,5.10200d-1,&
      5.10621d-1,5.11039d-1,5.11455d-1,5.11868d-1,5.12279d-1,5.12687d-1,&
      5.13093d-1,5.13496d-1,5.13897d-1,5.14295d-1,5.14690d-1,5.15083d-1,&
      5.15473d-1,5.15861d-1,5.16246d-1,5.16628d-1,5.17008d-1,5.17385d-1,&
      5.17760d-1,5.18131d-1,5.18501d-1,5.18867d-1,5.19634d-1,5.19993d-1,&
      5.20350d-1,5.20704d-1,5.21056d-1,5.21405d-1,5.21750d-1,5.22094d-1,&
      5.22434d-1,5.22772d-1,5.23107d-1,5.23439d-1,5.23768d-1,5.24094d-1,&
      5.24418d-1,5.24739d-1,5.25057d-1,5.25372d-1,5.25684d-1,5.25993d-1,&
      5.26300d-1,5.26603d-1,5.26904d-1,5.27202d-1,5.27497d-1,5.28945d-1,&
      5.29235d-1,5.29522d-1,5.29806d-1,5.30086d-1,5.30364d-1,5.30639d-1,&
      5.30910d-1,5.31179d-1,5.31445d-1,5.31707d-1,5.31967d-1,5.32223d-1,&
      5.32476d-1,5.32726d-1,5.32973d-1,5.33217d-1,5.33458d-1,5.33695d-1,&
      5.33929d-1,5.34160d-1,5.34388d-1,5.34613d-1,5.34834d-1,5.35053d-1,&
      5.34267d-1,5.34478d-1,5.34685d-1,5.34889d-1,5.35090d-1,5.35287d-1,&
      5.35481d-1,5.35672d-1,5.35859d-1,5.36043d-1,5.36223d-1,5.36401d-1,&
      5.36574d-1,5.36744d-1,5.36911d-1,5.37075d-1,5.37234d-1,5.37391d-1,&
      5.37544d-1,5.37693d-1,5.37839d-1,5.37981d-1,5.38119d-1,5.38254d-1,&
      5.38386d-1,5.36868d-1,5.36992d-1,5.37113d-1,5.37229d-1,5.37343d-1,&
      5.37453d-1,5.37560d-1,5.37663d-1,5.37763d-1,5.37860d-1,5.37954d-1,&
      5.38044d-1,5.38132d-1,5.38216d-1,5.38298d-1,5.38377d-1,5.38452d-1,&
      5.38525d-1,5.38596d-1,5.38663d-1,5.38728d-1,5.38790d-1,5.38850d-1,&
      5.38907d-1,5.38962d-1,5.39884d-1,5.39935d-1,5.39983d-1,5.40030d-1,&
      5.40074d-1,5.40115d-1,5.40155d-1,5.40192d-1,5.40228d-1,5.40261d-1,&
      5.40293d-1,5.40322d-1,5.40350d-1,5.40376d-1,5.40400d-1,5.40422d-1,&
      5.40443d-1,5.40462d-1,5.40480d-1,5.40496d-1,5.40510d-1,5.40523d-1,&
      5.40535d-1,5.40546d-1,5.40555d-1,5.39900d-1,5.39907d-1,5.39912d-1,&
      5.39916d-1,5.39919d-1,5.39922d-1,5.39923d-1,5.39923d-1,5.39923d-1,&
      5.39922d-1,5.39920d-1,5.39917d-1,5.39914d-1,5.39910d-1,5.39905d-1,&
      5.39900d-1,5.39895d-1,5.39889d-1,5.39883d-1,5.39876d-1,5.39870d-1,&
      5.39862d-1,5.39855d-1,5.39848d-1,5.39840d-1,5.40458d-1,5.40451d-1,&
      5.40444d-1,5.40437d-1,5.40430d-1,5.40423d-1,5.40416d-1,5.40410d-1,&
      5.40404d-1,5.40399d-1,5.40393d-1,5.40389d-1,5.40385d-1,5.40381d-1,&
      5.40378d-1,5.40376d-1,5.40374d-1,5.40374d-1,5.40373d-1,5.40374d-1,&
      5.40376d-1,5.40379d-1,5.40382d-1,5.40387d-1,5.40393d-1,5.39703d-1,&
      5.39710d-1,5.39719d-1,5.39729d-1,5.39740d-1,5.39752d-1,5.39765d-1,&
      5.39779d-1,5.39794d-1,5.39811d-1,5.39828d-1,5.39846d-1,5.39864d-1,&
      5.39884d-1,5.39905d-1,5.39926d-1,5.39948d-1,5.39971d-1,5.39995d-1,&
      5.40019d-1,5.40044d-1,5.40070d-1,5.40096d-1,5.40123d-1,5.40151d-1,&
      5.40620d-1,5.40649d-1,5.40679d-1,5.40709d-1,5.40740d-1,5.40771d-1,&
      5.40802d-1,5.40834d-1,5.40866d-1,5.40898d-1,5.40931d-1,5.40964d-1,&
      5.40997d-1,5.41031d-1,5.41064d-1,5.41098d-1,5.41132d-1,5.41166d-1,&
      5.41200d-1,5.41234d-1,5.41268d-1,5.41302d-1,5.41336d-1,5.41371d-1,&
      5.41404d-1,5.40845d-1,5.40879d-1,5.40914d-1,5.40947d-1,5.40981d-1,&
      5.41015d-1,5.41048d-1,5.41081d-1,5.41113d-1,5.41146d-1,5.41178d-1,&
      5.41209d-1,5.41240d-1,5.41271d-1,5.41301d-1,5.41331d-1,5.41360d-1,&
      5.41389d-1,5.41417d-1,5.41445d-1,5.41472d-1,5.41498d-1,5.41524d-1,&
      5.41549d-1,5.41573d-1,5.41584d-1,5.41606d-1,5.41628d-1,5.41649d-1,&
      5.41669d-1,5.41688d-1,5.41706d-1,5.41723d-1,5.41740d-1,5.41755d-1,&
      5.41769d-1,5.41783d-1,5.41795d-1,5.41806d-1,5.41816d-1,5.41825d-1,&
      5.41833d-1,5.41840d-1,5.41845d-1,5.41850d-1,5.41853d-1,5.41854d-1,&
      5.41855d-1,5.41854d-1,5.41851d-1,5.41989d-1,5.41984d-1,5.41977d-1,&
      5.41969d-1,5.41960d-1,5.41948d-1,5.41936d-1,5.41921d-1,5.41905d-1,&
      5.41888d-1,5.41869d-1,5.41848d-1,5.41825d-1,5.41801d-1,5.41775d-1,&
      5.41747d-1,5.41718d-1,5.41686d-1,5.41653d-1,5.41618d-1,5.41581d-1,&
      5.41542d-1,5.41501d-1,5.41458d-1,5.41414d-1,5.41711d-1,5.41663d-1,&
      5.41613d-1,5.41561d-1,5.41506d-1,5.41450d-1,5.41391d-1,5.41331d-1,&
      5.41268d-1,5.41202d-1,5.41135d-1,5.41065d-1,5.40993d-1,5.40918d-1,&
      5.40841d-1,5.40762d-1,5.40680d-1,5.40596d-1,5.40509d-1,5.40419d-1,&
      5.40327d-1,5.40233d-1,5.40136d-1,5.40036d-1,5.39934d-1,5.37776d-1,&
      5.37668d-1,5.37558d-1,5.37445d-1,5.37329d-1,5.37210d-1,5.37088d-1,&
      5.36963d-1,5.36836d-1,5.36706d-1,5.36573d-1,5.36436d-1,5.36297d-1,&
      5.36155d-1,5.36010d-1,5.35862d-1,5.35710d-1,5.35556d-1,5.35398d-1,&
      5.35238d-1,5.35074d-1,5.34907d-1,5.34737d-1,5.34563d-1,5.34386d-1,&
      5.36406d-1,5.36223d-1,5.36036d-1,5.35845d-1,5.35652d-1,5.35454d-1,&
      5.35254d-1,5.35049d-1,5.34842d-1,5.34630d-1,5.34415d-1,5.34197d-1,&
      5.33975d-1,5.33749d-1,5.33520d-1,5.33287d-1,5.33050d-1,5.32809d-1,&
      5.32565d-1,5.32317d-1,5.32065d-1,5.31809d-1,5.31550d-1,5.31286d-1,&
      5.31018d-1,5.30784d-1,5.30784d-1,5.30784d-1,5.30784d-1,5.30784d-1,&
      5.30784d-1,5.30784d-1,5.30784d-1,5.30784d-1,5.30784d-1,5.30784d-1,&
      5.30784d-1,5.30784d-1,5.30784d-1,5.30784d-1,5.30784d-1,5.30784d-1,&
      5.30784d-1,5.30784d-1,5.30784d-1,5.30784d-1,5.30784d-1,5.30784d-1,&
      5.30784d-1,5.30784d-1,5.29415d-1,5.29415d-1,5.29415d-1,5.29415d-1,&
      5.29415d-1,5.29415d-1,5.29415d-1,5.29415d-1,5.29415d-1,5.29415d-1,&
      5.29415d-1,5.29415d-1,5.29415d-1,5.29415d-1,5.29415d-1,5.29415d-1,&
      5.29415d-1,5.29415d-1,5.29415d-1,5.29415d-1,5.29415d-1,5.29415d-1,&
      5.29415d-1,5.29415d-1,5.29415d-1,5.29802d-1,5.29802d-1,5.29802d-1,&
      5.29802d-1,5.29802d-1,5.29802d-1,5.29802d-1,5.29802d-1,5.29802d-1,&
      5.29802d-1,5.29802d-1,5.29802d-1 /)
  ! Efficiency (total), excluding T5Z3
  REAL*8, PARAMETER :: EFF0_NO_T5Z3(NE)                                 &
      =       (/ 4.03332d-2,4.38198d-2,4.45862d-2,4.53504d-2,4.61127d-2,&
      4.68732d-2,4.76319d-2,4.83888d-2,4.91439d-2,4.98972d-2,5.06485d-2,&
      5.13976d-2,5.21445d-2,5.28888d-2,5.36302d-2,5.43685d-2,5.51216d-2,&
      5.58545d-2,5.65833d-2,5.73079d-2,5.80276d-2,5.87423d-2,5.94514d-2,&
      6.01547d-2,6.08517d-2,6.18145d-2,6.25019d-2,6.31818d-2,6.38541d-2,&
      6.45183d-2,6.51741d-2,6.58212d-2,6.64594d-2,6.70882d-2,6.77076d-2,&
      6.83172d-2,6.89168d-2,6.95063d-2,7.00855d-2,7.06542d-2,7.12123d-2,&
      7.17597d-2,7.22963d-2,7.28221d-2,7.33370d-2,7.38411d-2,7.43343d-2,&
      7.48167d-2,7.52882d-2,7.57491d-2,7.63636d-2,7.62879d-2,7.62091d-2,&
      7.61278d-2,7.60449d-2,7.59614d-2,7.58778d-2,7.57951d-2,7.57138d-2,&
      7.56348d-2,7.55586d-2,7.54859d-2,7.54173d-2,7.53533d-2,7.52944d-2,&
      7.52412d-2,7.51941d-2,7.51535d-2,7.51199d-2,7.50936d-2,7.50751d-2,&
      7.50646d-2,7.50624d-2,7.50689d-2,7.50842d-2,7.52904d-2,7.53255d-2,&
      7.53703d-2,7.54248d-2,7.54893d-2,7.55639d-2,7.56487d-2,7.57439d-2,&
      7.58495d-2,7.59657d-2,7.60927d-2,7.62303d-2,7.63787d-2,7.65381d-2,&
      7.67083d-2,7.68894d-2,7.70815d-2,7.72846d-2,7.74987d-2,7.77238d-2,&
      7.79599d-2,7.82069d-2,7.84649d-2,7.87339d-2,7.90138d-2,7.94429d-2,&
      7.97451d-2,8.00580d-2,8.03818d-2,8.07162d-2,8.10613d-2,8.14170d-2,&
      8.17833d-2,8.21601d-2,8.25473d-2,8.29450d-2,8.33529d-2,8.37711d-2,&
      8.41995d-2,8.46380d-2,8.50865d-2,8.55451d-2,8.60135d-2,8.64917d-2,&
      8.69797d-2,8.74773d-2,8.79845d-2,8.85012d-2,8.90274d-2,8.95628d-2,&
      9.04395d-2,9.09952d-2,9.15600d-2,9.21337d-2,9.27164d-2,9.33079d-2,&
      9.39080d-2,9.45168d-2,9.51341d-2,9.57598d-2,9.63937d-2,9.70359d-2,&
      9.76861d-2,9.83443d-2,9.90104d-2,9.96842d-2,1.00366d-1,1.01055d-1,&
      1.01751d-1,1.02454d-1,1.03165d-1,1.03883d-1,1.04608d-1,1.05339d-1,&
      1.06077d-1,1.07014d-1,1.07766d-1,1.08524d-1,1.09288d-1,1.10058d-1,&
      1.10834d-1,1.11616d-1,1.12403d-1,1.13195d-1,1.13993d-1,1.14796d-1,&
      1.15604d-1,1.16417d-1,1.17235d-1,1.18057d-1,1.18883d-1,1.19714d-1,&
      1.20549d-1,1.21389d-1,1.22232d-1,1.23078d-1,1.23929d-1,1.24782d-1,&
      1.25639d-1,1.26500d-1,1.27496d-1,1.28363d-1,1.29232d-1,1.30104d-1,&
      1.30978d-1,1.31855d-1,1.32734d-1,1.33614d-1,1.34497d-1,1.35382d-1,&
      1.36268d-1,1.37156d-1,1.38045d-1,1.38935d-1,1.39827d-1,1.40719d-1,&
      1.41613d-1,1.42507d-1,1.43402d-1,1.44297d-1,1.45192d-1,1.46088d-1,&
      1.46984d-1,1.47880d-1,1.48776d-1,1.49650d-1,1.50546d-1,1.51442d-1,&
      1.52338d-1,1.53232d-1,1.54126d-1,1.55019d-1,1.55911d-1,1.56802d-1,&
      1.57691d-1,1.58579d-1,1.59466d-1,1.60350d-1,1.61233d-1,1.62114d-1,&
      1.62994d-1,1.63871d-1,1.64745d-1,1.65618d-1,1.66488d-1,1.67355d-1,&
      1.68220d-1,1.69082d-1,1.69941d-1,1.70797d-1,1.71826d-1,1.72676d-1,&
      1.73522d-1,1.74366d-1,1.75205d-1,1.76041d-1,1.76873d-1,1.77701d-1,&
      1.78524d-1,1.79344d-1,1.80159d-1,1.80970d-1,1.81776d-1,1.82578d-1,&
      1.83374d-1,1.84166d-1,1.84953d-1,1.85734d-1,1.86511d-1,1.87281d-1,&
      1.88046d-1,1.88806d-1,1.89560d-1,1.90308d-1,1.91049d-1,1.91608d-1,&
      1.92337d-1,1.93060d-1,1.93776d-1,1.94487d-1,1.95192d-1,1.95890d-1,&
      1.96583d-1,1.97270d-1,1.97951d-1,1.98626d-1,1.99295d-1,1.99959d-1,&
      2.00618d-1,2.01270d-1,2.01918d-1,2.02560d-1,2.03196d-1,2.03827d-1,&
      2.04453d-1,2.05074d-1,2.05690d-1,2.06300d-1,2.06906d-1,2.07507d-1,&
      2.08857d-1,2.09449d-1,2.10037d-1,2.10620d-1,2.11198d-1,2.11771d-1,&
      2.12340d-1,2.12904d-1,2.13464d-1,2.14019d-1,2.14570d-1,2.15117d-1,&
      2.15659d-1,2.16197d-1,2.16731d-1,2.17356d-1,2.17888d-1,2.18416d-1,&
      2.18940d-1,2.19461d-1,2.19978d-1,2.20492d-1,2.21002d-1,2.21509d-1,&
      2.22013d-1,2.23162d-1,2.23661d-1,2.24156d-1,2.24648d-1,2.25138d-1,&
      2.25624d-1,2.26108d-1,2.26589d-1,2.27068d-1,2.27543d-1,2.28017d-1,&
      2.28488d-1,2.28956d-1,2.29422d-1,2.29886d-1,2.30348d-1,2.30808d-1,&
      2.31265d-1,2.31721d-1,2.32175d-1,2.32627d-1,2.33077d-1,2.33526d-1,&
      2.33973d-1,2.34418d-1,2.34752d-1,2.35193d-1,2.35632d-1,2.36071d-1,&
      2.36508d-1,2.36943d-1,2.37378d-1,2.37812d-1,2.38246d-1,2.38678d-1,&
      2.39110d-1,2.39541d-1,2.39971d-1,2.40401d-1,2.40831d-1,2.41260d-1,&
      2.41689d-1,2.42118d-1,2.42547d-1,2.42976d-1,2.43405d-1,2.43834d-1,&
      2.44263d-1,2.44692d-1,2.45122d-1,2.46053d-1,2.46487d-1,2.46921d-1,&
      2.47356d-1,2.54565d-1,2.55119d-1,2.55675d-1,2.56232d-1,2.56792d-1,&
      2.57353d-1,2.57916d-1,2.58482d-1,2.59049d-1,2.59618d-1,2.60189d-1,&
      2.60762d-1,2.61337d-1,2.61914d-1,2.62493d-1,2.63075d-1,2.63658d-1,&
      2.64244d-1,2.64832d-1,2.65422d-1,2.66014d-1,2.66774d-1,2.67371d-1,&
      2.67970d-1,2.68571d-1,2.69175d-1,2.69781d-1,2.70390d-1,2.71000d-1,&
      2.71613d-1,2.72229d-1,2.72847d-1,2.73467d-1,2.74090d-1,2.74715d-1,&
      2.75342d-1,2.75972d-1,2.76604d-1,2.77239d-1,2.77876d-1,2.78515d-1,&
      2.79157d-1,2.79802d-1,2.80449d-1,2.81098d-1,2.81750d-1,2.82328d-1,&
      2.82986d-1,2.83645d-1,2.84307d-1,2.84972d-1,2.85639d-1,2.86309d-1,&
      2.86981d-1,2.87655d-1,2.88332d-1,2.89012d-1,2.89694d-1,2.90378d-1,&
      2.91065d-1,2.91754d-1,2.92446d-1,2.93140d-1,2.93837d-1,2.94536d-1,&
      2.95237d-1,2.95941d-1,2.96648d-1,2.97357d-1,2.98068d-1,2.98782d-1,&
      2.99834d-1,3.00553d-1,3.01275d-1,3.02000d-1,3.02727d-1,3.03456d-1,&
      3.04188d-1,3.04922d-1,3.05659d-1,3.06398d-1,3.07140d-1,3.07884d-1,&
      3.08631d-1,3.09381d-1,3.10132d-1,3.10887d-1,3.11643d-1,3.12403d-1,&
      3.13165d-1,3.13929d-1,3.14696d-1,3.15466d-1,3.16238d-1,3.17013d-1,&
      3.17790d-1,3.18636d-1,3.19419d-1,3.20205d-1,3.20994d-1,3.21785d-1,&
      3.22578d-1,3.23373d-1,3.24170d-1,3.24970d-1,3.25771d-1,3.26574d-1,&
      3.27378d-1,3.28184d-1,3.28992d-1,3.29801d-1,3.30611d-1,3.31422d-1,&
      3.32235d-1,3.33048d-1,3.33863d-1,3.34678d-1,3.35493d-1,3.36310d-1,&
      3.37127d-1,3.37944d-1,3.38338d-1,3.40542d-1,3.41389d-1,3.42237d-1,&
      3.43085d-1,3.43934d-1,3.44782d-1,3.45631d-1,3.46479d-1,3.47328d-1,&
      3.48176d-1,3.49024d-1,3.49871d-1,3.50718d-1,3.51564d-1,3.52410d-1,&
      3.53254d-1,3.54098d-1,3.54941d-1,3.55782d-1,3.56622d-1,3.57461d-1,&
      3.58298d-1,3.59133d-1,3.59967d-1,3.61131d-1,3.61964d-1,3.62795d-1,&
      3.63624d-1,3.64451d-1,3.65276d-1,3.66098d-1,3.66917d-1,3.67734d-1,&
      3.68548d-1,3.69359d-1,3.70166d-1,3.70971d-1,3.71773d-1,3.72571d-1,&
      3.73366d-1,3.74157d-1,3.74944d-1,3.75727d-1,3.76507d-1,3.77283d-1,&
      3.78054d-1,3.78821d-1,3.79584d-1,3.80342d-1,3.81987d-1,3.82737d-1,&
      3.83483d-1,3.84223d-1,3.84958d-1,3.85688d-1,3.86413d-1,3.87132d-1,&
      3.87846d-1,3.88554d-1,3.89256d-1,3.89953d-1,3.90643d-1,3.91328d-1,&
      3.92006d-1,3.92677d-1,3.93343d-1,3.94002d-1,3.94654d-1,3.95299d-1,&
      3.95937d-1,3.96569d-1,3.97193d-1,3.97810d-1,3.98419d-1,3.97986d-1,&
      3.98579d-1,3.99165d-1,3.99743d-1,4.00315d-1,4.00878d-1,4.01435d-1,&
      4.01985d-1,4.02528d-1,4.03063d-1,4.03593d-1,4.04115d-1,4.04631d-1,&
      4.05141d-1,4.05644d-1,4.06141d-1,4.06632d-1,4.07116d-1,4.07595d-1,&
      4.08068d-1,4.08535d-1,4.08997d-1,4.09453d-1,4.09903d-1,4.10348d-1,&
      4.11741d-1,4.12176d-1,4.12607d-1,4.13033d-1,4.13454d-1,4.13870d-1,&
      4.14281d-1,4.14688d-1,4.15091d-1,4.15489d-1,4.15883d-1,4.16273d-1,&
      4.16660d-1,4.17042d-1,4.17420d-1,4.17795d-1,4.18167d-1,4.18535d-1,&
      4.18899d-1,4.19261d-1,4.19619d-1,4.19975d-1,4.20328d-1,4.20678d-1,&
      4.21025d-1,4.21213d-1,4.21557d-1,4.21898d-1,4.22238d-1,4.22575d-1,&
      4.22911d-1,4.23245d-1,4.23577d-1,4.23907d-1,4.24237d-1,4.24565d-1,&
      4.24891d-1,4.25217d-1,4.25542d-1,4.25866d-1,4.26189d-1,4.26512d-1,&
      4.26834d-1,4.27156d-1,4.27477d-1,4.27799d-1,4.28121d-1,4.28443d-1,&
      4.28765d-1,4.29087d-1,4.30424d-1,4.30749d-1,4.31074d-1,4.31400d-1,&
      4.31727d-1,4.32056d-1,4.32385d-1,4.32717d-1,4.33049d-1,4.33384d-1,&
      4.33720d-1,4.34058d-1,4.34398d-1,4.34741d-1,4.35085d-1,4.35433d-1,&
      4.35782d-1,4.36135d-1,4.36490d-1,4.36848d-1,4.37210d-1,4.37574d-1,&
      4.37942d-1,4.38313d-1,4.38688d-1,4.37505d-1,4.37887d-1,4.38273d-1,&
      4.38662d-1,4.39055d-1,4.39451d-1,4.39850d-1,4.40252d-1,4.40657d-1,&
      4.41066d-1,4.41477d-1,4.41890d-1,4.42306d-1,4.42724d-1,4.43145d-1,&
      4.43568d-1,4.43992d-1,4.44419d-1,4.44847d-1,4.45277d-1,4.45708d-1,&
      4.46141d-1,4.46575d-1,4.47010d-1,4.47446d-1,4.48392d-1,4.48829d-1,&
      4.49265d-1,4.49702d-1,4.50140d-1,4.50578d-1,4.51016d-1,4.51453d-1,&
      4.51891d-1,4.52329d-1,4.52766d-1,4.53202d-1,4.53638d-1,4.54073d-1,&
      4.54508d-1,4.54941d-1,4.55373d-1,4.55804d-1,4.56233d-1,4.56661d-1,&
      4.57088d-1,4.57512d-1,4.57935d-1,4.58356d-1,4.58774d-1,4.60348d-1,&
      4.60762d-1,4.61175d-1,4.61584d-1,4.61991d-1,4.62395d-1,4.62795d-1,&
      4.63193d-1,4.63587d-1,4.63978d-1,4.64365d-1,4.64749d-1,4.65129d-1,&
      4.65505d-1,4.65877d-1,4.66244d-1,4.66608d-1,4.66967d-1,4.67321d-1,&
      4.67671d-1,4.68015d-1,4.68355d-1,4.68690d-1,4.69020d-1,4.69344d-1,&
      4.68767d-1,4.69079d-1,4.69386d-1,4.69686d-1,4.69980d-1,4.70269d-1,&
      4.70551d-1,4.70827d-1,4.71096d-1,4.71359d-1,4.71614d-1,4.71863d-1,&
      4.72106d-1,4.72340d-1,4.72568d-1,4.72788d-1,4.73001d-1,4.73206d-1,&
      4.73404d-1,4.73593d-1,4.73774d-1,4.73948d-1,4.74113d-1,4.74269d-1,&
      4.74417d-1,4.73115d-1,4.73245d-1,4.73366d-1,4.73480d-1,4.73584d-1,&
      4.73681d-1,4.73769d-1,4.73849d-1,4.73921d-1,4.73986d-1,4.74043d-1,&
      4.74092d-1,4.74134d-1,4.74169d-1,4.74196d-1,4.74217d-1,4.74231d-1,&
      4.74237d-1,4.74238d-1,4.74231d-1,4.74219d-1,4.74200d-1,4.74174d-1,&
      4.74143d-1,4.74106d-1,4.74698d-1,4.74650d-1,4.74596d-1,4.74537d-1,&
      4.74473d-1,4.74404d-1,4.74329d-1,4.74250d-1,4.74166d-1,4.74077d-1,&
      4.73983d-1,4.73885d-1,4.73783d-1,4.73676d-1,4.73566d-1,4.73451d-1,&
      4.73333d-1,4.73211d-1,4.73086d-1,4.72957d-1,4.72824d-1,4.72689d-1,&
      4.72550d-1,4.72409d-1,4.72264d-1,4.71711d-1,4.71561d-1,4.71409d-1,&
      4.71255d-1,4.71099d-1,4.70940d-1,4.70780d-1,4.70618d-1,4.70454d-1,&
      4.70288d-1,4.70122d-1,4.69953d-1,4.69784d-1,4.69613d-1,4.69441d-1,&
      4.69269d-1,4.69096d-1,4.68922d-1,4.68748d-1,4.68573d-1,4.68398d-1,&
      4.68223d-1,4.68049d-1,4.67874d-1,4.67699d-1,4.67890d-1,4.67716d-1,&
      4.67542d-1,4.67370d-1,4.67198d-1,4.67027d-1,4.66857d-1,4.66689d-1,&
      4.66522d-1,4.66356d-1,4.66192d-1,4.66030d-1,4.65869d-1,4.65711d-1,&
      4.65554d-1,4.65400d-1,4.65248d-1,4.65099d-1,4.64952d-1,4.64808d-1,&
      4.64666d-1,4.64528d-1,4.64393d-1,4.64260d-1,4.64132d-1,4.63422d-1,&
      4.63300d-1,4.63181d-1,4.63066d-1,4.62955d-1,4.62847d-1,4.62742d-1,&
      4.62640d-1,4.62542d-1,4.62447d-1,4.62355d-1,4.62266d-1,4.62181d-1,&
      4.62098d-1,4.62019d-1,4.61942d-1,4.61869d-1,4.61798d-1,4.61731d-1,&
      4.61666d-1,4.61604d-1,4.61545d-1,4.61488d-1,4.61434d-1,4.61383d-1,&
      4.61748d-1,4.61702d-1,4.61659d-1,4.61619d-1,4.61580d-1,4.61545d-1,&
      4.61511d-1,4.61480d-1,4.61452d-1,4.61425d-1,4.61401d-1,4.61378d-1,&
      4.61358d-1,4.61340d-1,4.61324d-1,4.61310d-1,4.61298d-1,4.61288d-1,&
      4.61280d-1,4.61273d-1,4.61269d-1,4.61266d-1,4.61265d-1,4.61265d-1,&
      4.61267d-1,4.60858d-1,4.60864d-1,4.60871d-1,4.60880d-1,4.60890d-1,&
      4.60902d-1,4.60915d-1,4.60929d-1,4.60945d-1,4.60961d-1,4.60979d-1,&
      4.60998d-1,4.61019d-1,4.61040d-1,4.61062d-1,4.61086d-1,4.61110d-1,&
      4.61135d-1,4.61161d-1,4.61188d-1,4.61216d-1,4.61244d-1,4.61273d-1,&
      4.61303d-1,4.61334d-1,4.61426d-1,4.61457d-1,4.61489d-1,4.61521d-1,&
      4.61554d-1,4.61587d-1,4.61621d-1,4.61655d-1,4.61689d-1,4.61723d-1,&
      4.61758d-1,4.61792d-1,4.61827d-1,4.61862d-1,4.61897d-1,4.61932d-1,&
      4.61967d-1,4.62002d-1,4.62037d-1,4.62071d-1,4.62106d-1,4.62140d-1,&
      4.62174d-1,4.62207d-1,4.62241d-1,4.62189d-1,4.62221d-1,4.62253d-1,&
      4.62284d-1,4.62315d-1,4.62345d-1,4.62374d-1,4.62403d-1,4.62431d-1,&
      4.62459d-1,4.62485d-1,4.62511d-1,4.62536d-1,4.62560d-1,4.62583d-1,&
      4.62605d-1,4.62626d-1,4.62646d-1,4.62665d-1,4.62682d-1,4.62699d-1,&
      4.62714d-1,4.62728d-1,4.62741d-1,4.62752d-1,4.63147d-1,4.63157d-1,&
      4.63165d-1,4.63171d-1,4.63176d-1,4.63180d-1,4.63181d-1,4.63181d-1,&
      4.63180d-1,4.63177d-1,4.63171d-1,4.63164d-1,4.63156d-1,4.63145d-1,&
      4.63132d-1,4.63117d-1,4.63101d-1,4.63082d-1,4.63061d-1,4.63038d-1,&
      4.63013d-1,4.62985d-1,4.62955d-1,4.62923d-1,4.62889d-1,4.60928d-1,&
      4.60889d-1,4.60847d-1,4.60802d-1,4.60755d-1,4.60706d-1,4.60654d-1,&
      4.60599d-1,4.60541d-1,4.60481d-1,4.60418d-1,4.60352d-1,4.60284d-1,&
      4.60212d-1,4.60138d-1,4.60060d-1,4.59980d-1,4.59896d-1,4.59809d-1,&
      4.59720d-1,4.59627d-1,4.59531d-1,4.59431d-1,4.59329d-1,4.59223d-1,&
      4.61360d-1,4.61247d-1,4.61131d-1,4.61011d-1,4.60887d-1,4.60760d-1,&
      4.60629d-1,4.60495d-1,4.60357d-1,4.60215d-1,4.60070d-1,4.59920d-1,&
      4.59767d-1,4.59610d-1,4.59449d-1,4.59284d-1,4.59115d-1,4.58942d-1,&
      4.58765d-1,4.58584d-1,4.58398d-1,4.58209d-1,4.58015d-1,4.57817d-1,&
      4.57614d-1,4.57470d-1,4.57470d-1,4.57470d-1,4.57470d-1,4.57470d-1,&
      4.57470d-1,4.57470d-1,4.57470d-1,4.57470d-1,4.57470d-1,4.57470d-1,&
      4.57470d-1,4.57470d-1,4.57470d-1,4.57470d-1,4.57470d-1,4.57470d-1,&
      4.57470d-1,4.57470d-1,4.57470d-1,4.57470d-1,4.57470d-1,4.57470d-1,&
      4.57470d-1,4.57470d-1,4.56133d-1,4.56133d-1,4.56133d-1,4.56133d-1,&
      4.56133d-1,4.56133d-1,4.56133d-1,4.56133d-1,4.56133d-1,4.56133d-1,&
      4.56133d-1,4.56133d-1,4.56133d-1,4.56133d-1,4.56133d-1,4.56133d-1,&
      4.56133d-1,4.56133d-1,4.56133d-1,4.56133d-1,4.56133d-1,4.56133d-1,&
      4.56133d-1,4.56133d-1,4.56133d-1,4.56458d-1,4.56458d-1,4.56458d-1,&
      4.56458d-1,4.56458d-1,4.56458d-1,4.56458d-1,4.56458d-1,4.56458d-1,&
      4.56458d-1,4.56458d-1,4.56458d-1/)
  
  ! Number of events for case of all detectors
  INTEGER, PARAMETER :: NINTERVALS_ALL_TZ  = 12
  ! Number of events for case without T5Z3
  INTEGER, PARAMETER :: NINTERVALS_NO_T5Z3 = 9
  
  ! Array containing interval bounds [E_k,E_{k+1}]
  ! i.e. list of event energies plus bounds of overall energy range
  ! recoil energies converted from total phonon energy using Lindhard
  ! rather than SCDMS charge model to be consistent with efficiencies.
  ! For case of all detectors
  REAL*8, PARAMETER :: E_INTERVALS_ALL_TZ(0:NINTERVALS_ALL_TZ)          &
      = (/ 1.60867d0, 1.67648d0, 1.72268d0, 1.84366d0, 1.9301d0,        &
           2.20435d0, 2.57638d0, 2.82073d0, 5.57222d0, 6.74554d0,       &
           7.53282d0, 9.16426d0, 9.98086d0 /)
  ! for case without T5Z3
  REAL*8, PARAMETER :: E_INTERVALS_NO_T5Z3(0:NINTERVALS_NO_T5Z3)        &
      = (/ 1.60867d0, 1.67648d0, 1.72268d0, 1.84366d0, 1.9301d0,        &
           2.20435d0, 2.57638d0, 2.82073d0, 5.57222d0, 9.98086d0 /)
  
  ! Include T5Z3?
  LOGICAL, PARAMETER :: INCLUDE_T5Z3 = .FALSE.
  
  ! Will build array of efficiencies below
  INTEGER :: Nintervals
  REAL*8, ALLOCATABLE :: eff(:,:)
  
  ! Fill in efficiencies
  ! Interval efficiencies are just total efficiency with interval
  ! energy range and zero elsewhere (ignoring energy resolution).
  IF (INCLUDE_T5Z3) THEN
    Nintervals = NINTERVALS_ALL_TZ
    ALLOCATE(eff(NE,0:Nintervals))
    eff(:,0) = EFF0_ALL_TZ
    DO K=1,Nintervals
      WHERE ((E .GE. E_INTERVALS_ALL_TZ(K-1)) .AND. (E .LE. E_INTERVALS_ALL_TZ(K)))
        eff(:,K) = EFF0_ALL_TZ
      ELSE WHERE
        eff(:,K) = 0d0
      END WHERE
    END DO
  ELSE
    Nintervals = NINTERVALS_NO_T5Z3
    ALLOCATE(eff(NE,0:Nintervals))
    eff(:,0) = EFF0_ALL_TZ
    DO K=1,Nintervals
      WHERE ((E .GE. E_INTERVALS_NO_T5Z3(K-1)) .AND. (E .LE. E_INTERVALS_NO_T5Z3(K)))
        eff(:,K) = EFF0_NO_T5Z3
      ELSE WHERE
        eff(:,K) = 0d0
      END WHERE
    END DO
  END IF
  
  ! One call for all settings.
  ! Most of these _must_ be there to ensure everything get initialized.
  IF (INCLUDE_T5Z3) THEN
    ! These settings are for the analysis with all detectors
    CALL SetDetector(D,mass=4.2d0,time=137.4d0,Nevents=11,              &
                     background=6.1d0,Nelem=1,Zelem=(/32/),             &
                     NEeff=NE,Eeff=E,Neff=Nintervals,eff=eff,           &
                     intervals=intervals)
  ELSE
    ! These settings are for the analysis without t5z3
    CALL SetDetector(D,mass=3.6d0,time=137.4d0,Nevents=8,               &
                     background=6.07d0,Nelem=1,Zelem=(/32/),            &
                     NEeff=NE,Eeff=E,Neff=Nintervals,eff=eff,           &
                     intervals=intervals)
  END IF
  
  D%eff_file = '[SuperCDMS 2014]'
  
END SUBROUTINE



!=======================================================================
! SIMPLE 2014 ANALYSIS ROUTINES
! Based upon the SIMPLE Phase II WIMP search: PRD 89, 072013 (2014)
! [1404.4309].  7+1 candidate events seen with 10.8+1.9 expected
! background events.
! 
! BIND() is used to specify compiler-independent object file symbol
! names to allow for easier interfacing with C/C++.  
!=======================================================================

!-----------------------------------------------------------------------
! Initializes the module to perform calculations for the SIMPLE 2014
! analysis.  This must be called if any of the following SIMPLE 2014
! routines are to be used.
! 
! Required input arguments:
!     intervals   Indicates if sub-intervals should be included.
!                 Only necessary if confidence intervals using the
!                 maximum gap method are desired.
! 
SUBROUTINE SIMPLE_2014_Init(intervals)
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: intervals
  CALL SIMPLE_2014_InitTo(SIMPLE_2014,intervals)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_SIMPLE_2014_Init(intervals) &
           BIND(C,NAME='C_DDCALC0_simple_2014_init')
  USE ISO_C_BINDING, only: C_BOOL
  IMPLICIT NONE
  LOGICAL(KIND=C_BOOL), INTENT(IN) :: intervals
  CALL SIMPLE_2014_Init(LOGICAL(intervals))
END SUBROUTINE


! ----------------------------------------------------------------------
! Sets the minimum recoil energy to be included in the calculations.
! Note the efficiency curves already account for detector and analysis
! thresholds regardless of this setting, so setting this to 0 keV (the
! default behavior when initialization is performed) does not imply
! that very low energy recoils actually contribute to the signal.
! 
! Required input arguments:
!     Emin        The minimum recoil energy to consider [keV]
! 
SUBROUTINE SIMPLE_2014_SetEmin(Emin)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: Emin
  CALL SetDetector(SIMPLE_2014,Emin=Emin)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_SIMPLE_2014_SetEmin(Emin) &
           BIND(C,NAME='C_DDCALC0_simple_2014_setemin')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: Emin
  CALL SIMPLE_2014_SetEmin(REAL(Emin,KIND=8))
END SUBROUTINE


! ----------------------------------------------------------------------
! Calculates various rate quantities using the current WIMP.
! Must be called each time the WIMP parameters are modified.
! 
SUBROUTINE SIMPLE_2014_CalcRates()
  IMPLICIT NONE
  CALL CalcRates(SIMPLE_2014)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_SIMPLE_2014_CalcRates() &
           BIND(C,NAME='C_DDCALC0_simple_2014_calcrates')
  IMPLICIT NONE
  CALL SIMPLE_2014_CalcRates()
END SUBROUTINE


! ----------------------------------------------------------------------
! Returns the observed number of events.
! 
FUNCTION SIMPLE_2014_Events() RESULT(N)
  IMPLICIT NONE
  INTEGER :: N
  CALL GetRates(SIMPLE_2014,Nevents=N)
END FUNCTION

! C++ interface wrapper
FUNCTION C_SIMPLE_2014_Events() RESULT(N) &
         BIND(C,NAME='C_DDCALC0_simple_2014_events')
  USE ISO_C_BINDING, only: C_INT
  IMPLICIT NONE
  INTEGER(KIND=C_INT) :: N
  ! Automatic type conversions here
  N = SIMPLE_2014_Events()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of background events.
! 
FUNCTION SIMPLE_2014_Background() RESULT(b)
  IMPLICIT NONE
  REAL*8 :: b
  CALL GetRates(SIMPLE_2014,background=b)
END FUNCTION

! C++ interface wrapper
FUNCTION C_SIMPLE_2014_Background() RESULT(b) &
         BIND(C,NAME='C_DDCALC0_simple_2014_background')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: b
  ! Automatic type conversions here
  b = SIMPLE_2014_Background()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of signal events for the
! current WIMP.
! 
FUNCTION SIMPLE_2014_Signal() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(SIMPLE_2014,signal=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_SIMPLE_2014_Signal() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_simple_2014_signal')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = SIMPLE_2014_Signal()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of spin-independent signal events
! for the current WIMP.
! 
FUNCTION SIMPLE_2014_SignalSI() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(SIMPLE_2014,signal_si=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_SIMPLE_2014_SignalSI() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_simple_2014_signalsi')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = SIMPLE_2014_SignalSI()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of spin-dependent signal events
! for the current WIMP.
! 
FUNCTION SIMPLE_2014_SignalSD() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(SIMPLE_2014,signal_sd=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_SIMPLE_2014_SignalSD() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_simple_2014_signalsd')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = SIMPLE_2014_SignalSD()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the log-likelihood for the current WIMP mass and couplings.
! Uses a Poisson distribution in the number of observed events N:
!    P(N|s+b)
! where s is the average expected signal and b is the average expected
! background.
! 
FUNCTION SIMPLE_2014_LogLikelihood() RESULT(lnlike)
  IMPLICIT NONE
  REAL*8 :: lnlike
  lnlike = LogLikelihood(SIMPLE_2014)
END FUNCTION

! C++ interface wrapper
FUNCTION C_SIMPLE_2014_LogLikelihood() RESULT(lnlike) &
         BIND(C,NAME='C_DDCALC0_simple_2014_loglikelihood')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: lnlike
  ! Automatic type conversions here
  lnlike = SIMPLE_2014_LogLikelihood()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the log of the p-value for the current WIMP mass and
! couplings (NO BACKGROUND SUBTRACTION).  Uses the maximum gap method
! if SIMPLE_2014_Init was called with argument intervals=.TRUE.,
! otherwise uses a Poisson distribution in the number of observed
! events N:
!    P(N|s)
! where s is the average expected signal (background contributions are
! ignored).
! 
FUNCTION SIMPLE_2014_LogPValue() RESULT(lnp)
  IMPLICIT NONE
  REAL*8 :: lnp
  lnp = LogPValue(SIMPLE_2014)
END FUNCTION

! C++ interface wrapper
FUNCTION C_SIMPLE_2014_LogPValue() RESULT(lnp) &
         BIND(C,NAME='C_DDCALC0_simple_2014_logpvalue')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: lnp
  ! Automatic type conversions here
  lnp = SIMPLE_2014_LogPValue()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the factor x by which the cross-sections must be scaled
! (sigma -> x*sigma) to achieve the desired p-value (given as log(p)).
! See LogPValue() above for a description of the statistics.
! 
! Required input argument:
!   lnp         The logarithm of the desired p-value (p = 1-CL).
! 
FUNCTION SIMPLE_2014_ScaleToPValue(lnp) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN) :: lnp
  x = ScaleToPValue(SIMPLE_2014,lnp)
END FUNCTION

! C++ interface wrapper
FUNCTION C_SIMPLE_2014_ScaleToPValue(lnp) RESULT(x) &
         BIND(C,NAME='C_DDCALC0_simple_2014_scaletopvalue')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: x
  REAL(KIND=C_DOUBLE), INTENT(IN) :: lnp
  ! Automatic type conversions here
  x = SIMPLE_2014_ScaleToPValue(REAL(lnp,KIND=8))
END FUNCTION


!-----------------------------------------------------------------------
! INTERNAL ROUTINE.
! Initializes the given DetectorStruct to the SIMPLE 2014 analysis.
! This is meant as an internal routine; external access should be
! through SIMPLE_2014_Init instead.
! 
! Required input arguments:
!     D           The DetectorStruct to initialize
!     intervals   Indicates if sub-intervals should be included.
!                 Currently ignored as there is no event energy
!                 information.
! 
SUBROUTINE SIMPLE_2014_InitTo(D,intervals)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(OUT) :: D
  LOGICAL, INTENT(IN) :: intervals
  INTEGER :: K,Kmin,Kmax,Niso,Niso0
  INTEGER, ALLOCATABLE :: Ziso(:),Ziso0(:),Aiso(:),Aiso0(:)
  REAL*8, ALLOCATABLE :: fiso(:),fiso0(:),Miso(:),Miso0(:)
  
  ! Will build array of efficiencies below
  INTEGER :: NE
  REAL*8 :: Ethresh,Gamma
  REAL*8, ALLOCATABLE :: E(:),eff(:,:)
  
  ! Uses C2ClF5, but C threshold is ~ 100 keV, so we will drop
  ! the C components.  First get full set of isotopes to get mass
  ! fractions correct.  C is placed last in the stoichiometry to
  ! allow easier extraction of F & Cl parts.
  CALL CompoundIsotopeList(3,(/17,9,6/),(/1,5,2/),                      &
                           Niso0,Ziso0,Aiso0,fiso0,Miso0)
  Niso = COUNT(Ziso0 .NE. 6)
  ALLOCATE(Ziso(Niso),Aiso(Niso),fiso(Niso),Miso(Niso))
  Ziso = Ziso0(1:Niso)
  Aiso = Aiso0(1:Niso)
  fiso = fiso0(1:Niso)
  Miso = Miso0(1:Niso)
  
  ! Set efficiencies.
  ! No event energies, so no intervals for max gap calculations.
  ! Tabulation set at 100 per decade, up to 1000 keV.
  Ethresh = 8d0      ! F & Cl only, C is ~ 100 keV
  Gamma   = 4.2d0    ! 4.2 +/- 0.3
  Kmin = INT(100*LOG10(Ethresh))
  Kmax = NINT(100*LOG10(1d3))
  NE = (Kmax-Kmin)+1
  ALLOCATE(E(1:NE),eff(1:NE,0:0))
  DO K=Kmin,Kmax
    E(K-Kmin+1) = 10**(K/100d0)
  END DO
  eff(:,0) = 1d0 - EXP(-Gamma*(1d0-Ethresh/MAX(E,Ethresh)))
  
  ! One call for all settings.
  ! Most of these _must_ be there to ensure everything get initialized.
  CALL SetDetector(D,exposure=18.24d0,Nevents=8,background=12.7d0,      &
                   Niso=Niso,Ziso=Ziso,Aiso=Aiso,fiso=fiso,             &
                   NEeff=NE,Eeff=E,Neff=0,eff=eff,                      &
                   intervals=.FALSE.)
  D%eff_file = '[SIMPLE 2014]'
  
END SUBROUTINE



!=======================================================================
! DARWIN ARGON ANALYSIS ROUTINES
! Based upon a DARWIN argon-based analysis (2015 estimated parameters)
! [15MM.NNNNN].  Zero events assumed in the analysis region.
! 
! BIND() is used to specify compiler-independent object file symbol
! names to allow for easier interfacing with C/C++.  
!=======================================================================

!-----------------------------------------------------------------------
! Initializes the module to perform calculations for a DARWIN
! argon-based analysis (2015 estimates).  This must be called if any
! of the following DARWIN routines are to be used.
! 
! Required input arguments:
!     intervals   Indicates if sub-intervals should be included.
!                 Only necessary if confidence intervals using the
!                 maximum gap method are desired.
! 
SUBROUTINE DARWIN_Ar_2015_Init(intervals)
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: intervals
  CALL DARWIN_Ar_2015_InitTo(DARWIN_Ar_2015,intervals)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_DARWIN_Ar_2015_Init(intervals) &
           BIND(C,NAME='C_DDCALC0_darwin_ar_2015_init')
  USE ISO_C_BINDING, only: C_BOOL
  IMPLICIT NONE
  LOGICAL(KIND=C_BOOL), INTENT(IN) :: intervals
  CALL DARWIN_Ar_2015_Init(LOGICAL(intervals))
END SUBROUTINE


! ----------------------------------------------------------------------
! Sets the minimum recoil energy to be included in the calculations.
! Note the efficiency curves already account for detector and analysis
! thresholds regardless of this setting, so setting this to 0 keV (the
! default behavior when initialization is performed) does not imply
! that very low energy recoils actually contribute to the signal.
! 
! Required input arguments:
!     Emin        The minimum recoil energy to consider [keV]
! 
SUBROUTINE DARWIN_Ar_2015_SetEmin(Emin)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: Emin
  CALL SetDetector(DARWIN_Ar_2015,Emin=Emin)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_DARWIN_Ar_2015_SetEmin(Emin) &
           BIND(C,NAME='C_DDCALC0_darwin_ar_2015_setemin')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: Emin
  CALL DARWIN_Ar_2015_SetEmin(REAL(Emin,KIND=8))
END SUBROUTINE


! ----------------------------------------------------------------------
! Calculates various rate quantities using the current WIMP.
! Must be called each time the WIMP parameters are modified.
! 
SUBROUTINE DARWIN_Ar_2015_CalcRates()
  IMPLICIT NONE
  CALL CalcRates(DARWIN_Ar_2015)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_DARWIN_Ar_2015_CalcRates() &
           BIND(C,NAME='C_DDCALC0_darwin_ar_2015_calcrates')
  IMPLICIT NONE
  CALL DARWIN_Ar_2015_CalcRates()
END SUBROUTINE


! ----------------------------------------------------------------------
! Returns the observed number of events.
! 
FUNCTION DARWIN_Ar_2015_Events() RESULT(N)
  IMPLICIT NONE
  INTEGER :: N
  CALL GetRates(DARWIN_Ar_2015,Nevents=N)
END FUNCTION

! C++ interface wrapper
FUNCTION C_DARWIN_Ar_2015_Events() RESULT(N) &
         BIND(C,NAME='C_DDCALC0_darwin_ar_2015_events')
  USE ISO_C_BINDING, only: C_INT
  IMPLICIT NONE
  INTEGER(KIND=C_INT) :: N
  ! Automatic type conversions here
  N = DARWIN_Ar_2015_Events()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of background events.
! 
FUNCTION DARWIN_Ar_2015_Background() RESULT(b)
  IMPLICIT NONE
  REAL*8 :: b
  CALL GetRates(DARWIN_Ar_2015,background=b)
END FUNCTION

! C++ interface wrapper
FUNCTION C_DARWIN_Ar_2015_Background() RESULT(b) &
         BIND(C,NAME='C_DDCALC0_darwin_ar_2015_background')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: b
  ! Automatic type conversions here
  b = DARWIN_Ar_2015_Background()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of signal events for the
! current WIMP.
! 
FUNCTION DARWIN_Ar_2015_Signal() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(DARWIN_Ar_2015,signal=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_DARWIN_Ar_2015_Signal() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_darwin_ar_2015_signal')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = DARWIN_Ar_2015_Signal()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of spin-independent signal events
! for the current WIMP.
! 
FUNCTION DARWIN_Ar_2015_SignalSI() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(DARWIN_Ar_2015,signal_si=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_DARWIN_Ar_2015_SignalSI() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_darwin_ar_2015_signalsi')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = DARWIN_Ar_2015_SignalSI()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of spin-dependent signal events
! for the current WIMP.
! 
FUNCTION DARWIN_Ar_2015_SignalSD() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(DARWIN_Ar_2015,signal_sd=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_DARWIN_Ar_2015_SignalSD() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_darwin_ar_2015_signalsd')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = DARWIN_Ar_2015_SignalSD()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the log-likelihood for the current WIMP mass and couplings.
! Uses a Poisson distribution in the number of observed events N:
!    P(N|s+b)
! where s is the average expected signal and b is the average expected
! background.
! 
FUNCTION DARWIN_Ar_2015_LogLikelihood() RESULT(lnlike)
  IMPLICIT NONE
  REAL*8 :: lnlike
  lnlike = LogLikelihood(DARWIN_Ar_2015)
END FUNCTION

! C++ interface wrapper
FUNCTION C_DARWIN_Ar_2015_LogLikelihood() RESULT(lnlike) &
         BIND(C,NAME='C_DDCALC0_darwin_ar_2015_loglikelihood')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: lnlike
  ! Automatic type conversions here
  lnlike = DARWIN_Ar_2015_LogLikelihood()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the log of the p-value for the current WIMP mass and
! couplings (NO BACKGROUND SUBTRACTION).  Uses the maximum gap method
! if DARWIN_Ar_2015_Init was called with argument intervals=.TRUE.,
! otherwise uses a Poisson distribution in the number of observed
! events N:
!    P(N|s)
! where s is the average expected signal (background contributions are
! ignored).
! 
FUNCTION DARWIN_Ar_2015_LogPValue() RESULT(lnp)
  IMPLICIT NONE
  REAL*8 :: lnp
  lnp = LogPValue(DARWIN_Ar_2015)
END FUNCTION

! C++ interface wrapper
FUNCTION C_DARWIN_Ar_2015_LogPValue() RESULT(lnp) &
         BIND(C,NAME='C_DDCALC0_darwin_ar_2015_logpvalue')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: lnp
  ! Automatic type conversions here
  lnp = DARWIN_Ar_2015_LogPValue()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the factor x by which the cross-sections must be scaled
! (sigma -> x*sigma) to achieve the desired p-value (given as log(p)).
! See LogPValue() above for a description of the statistics.
! 
! Required input argument:
!   lnp         The logarithm of the desired p-value (p = 1-CL).
! 
FUNCTION DARWIN_Ar_2015_ScaleToPValue(lnp) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN) :: lnp
  x = ScaleToPValue(DARWIN_Ar_2015,lnp)
END FUNCTION

! C++ interface wrapper
FUNCTION C_DARWIN_Ar_2015_ScaleToPValue(lnp) RESULT(x) &
         BIND(C,NAME='C_DDCALC0_darwin_ar_2015_scaletopvalue')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: x
  REAL(KIND=C_DOUBLE), INTENT(IN) :: lnp
  ! Automatic type conversions here
  x = DARWIN_Ar_2015_ScaleToPValue(REAL(lnp,KIND=8))
END FUNCTION


!-----------------------------------------------------------------------
! INTERNAL ROUTINE.
! Initializes the given DetectorStruct to the proposed DARWIN
! argon-based detector (2015 estimate).  This is meant as an internal
! routine; external access should be through
! DARWIN_Ar_2015_Init instead.
! 
! The efficiencies used here were generated using TPCMC.
! 
! Required input arguments:
!     D           The DetectorStruct to initialize
!     intervals   Indicates if sub-intervals should be included
! 
SUBROUTINE DARWIN_Ar_2015_InitTo(D,intervals)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(OUT) :: D
  LOGICAL, INTENT(IN) :: intervals
  INTEGER, PARAMETER :: NE = 151
  INTEGER, PARAMETER :: NEFF = 1
  ! Efficiency curves energy tabulation points
  REAL*8, PARAMETER :: E(NE)                                            &
      =       (/ 0.10000d0, 0.10471d0, 0.10965d0, 0.11482d0, 0.12023d0, &
      0.12589d0, 0.13183d0, 0.13804d0, 0.14454d0, 0.15136d0, 0.15849d0, &
      0.16596d0, 0.17378d0, 0.18197d0, 0.19055d0, 0.19953d0, 0.20893d0, &
      0.21878d0, 0.22909d0, 0.23988d0, 0.25119d0, 0.26303d0, 0.27542d0, &
      0.28840d0, 0.30200d0, 0.31623d0, 0.33113d0, 0.34674d0, 0.36308d0, &
      0.38019d0, 0.39811d0, 0.41687d0, 0.43652d0, 0.45709d0, 0.47863d0, &
      0.50119d0, 0.52481d0, 0.54954d0, 0.57544d0, 0.60256d0, 0.63096d0, &
      0.66069d0, 0.69183d0, 0.72444d0, 0.75858d0, 0.79433d0, 0.83176d0, &
      0.87096d0, 0.91201d0, 0.95499d0, 1.0000d0,  1.0471d0,  1.0965d0,  &
      1.1482d0,  1.2023d0,  1.2589d0,  1.3183d0,  1.3804d0,  1.4454d0,  &
      1.5136d0,  1.5849d0,  1.6596d0,  1.7378d0,  1.8197d0,  1.9055d0,  &
      1.9953d0,  2.0893d0,  2.1878d0,  2.2909d0,  2.3988d0,  2.5119d0,  &
      2.6303d0,  2.7542d0,  2.8840d0,  3.0200d0,  3.1623d0,  3.3113d0,  &
      3.4674d0,  3.6308d0,  3.8019d0,  3.9811d0,  4.1687d0,  4.3652d0,  &
      4.5709d0,  4.7863d0,  5.0119d0,  5.2481d0,  5.4954d0,  5.7544d0,  &
      6.0256d0,  6.3096d0,  6.6069d0,  6.9183d0,  7.2444d0,  7.5858d0,  &
      7.9433d0,  8.3176d0,  8.7096d0,  9.1201d0,  9.5499d0, 10.000d0,   &
     10.471d0,  10.965d0,  11.482d0,  12.023d0,  12.589d0,  13.183d0,   &
     13.804d0,  14.454d0,  15.136d0,  15.849d0,  16.596d0,  17.378d0,   &
     18.197d0,  19.055d0,  19.953d0,  20.893d0,  21.878d0,  22.909d0,   &
     23.988d0,  25.119d0,  26.303d0,  27.542d0,  28.840d0,  30.200d0,   &
     31.623d0,  33.113d0,  34.674d0,  36.308d0,  38.019d0,  39.811d0,   &
     41.687d0,  43.652d0,  45.709d0,  47.863d0,  50.119d0,  52.481d0,   &
     54.954d0,  57.544d0,  60.256d0,  63.096d0,  66.069d0,  69.183d0,   &
     72.444d0,  75.858d0,  79.433d0,  83.176d0,  87.096d0,  91.201d0,   &
     95.499d0, 100.00d0 /)
  ! Efficiency (total)
  REAL*8, PARAMETER :: EFF0(NE)                                         &
      =       (/ 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      1.00000d-6,1.00000d-6,3.00000d-6,5.00000d-6,4.00000d-6,4.00000d-6,&
      6.00000d-6,9.00000d-6,1.20000d-5,1.70000d-5,3.60000d-5,5.40000d-5,&
      7.00000d-5,1.06000d-4,1.48000d-4,2.48000d-4,3.74000d-4,5.11000d-4,&
      6.61000d-4,9.53000d-4,1.32400d-3,1.91900d-3,2.57200d-3,3.42800d-3,&
      4.78200d-3,6.31100d-3,8.37200d-3,1.11320d-2,1.43600d-2,1.88560d-2,&
      2.39690d-2,3.04640d-2,3.85180d-2,4.79960d-2,5.92610d-2,7.20370d-2,&
      8.75540d-2,1.04670d-1,1.24360d-1,1.44260d-1,1.67600d-1,1.92550d-1,&
      2.19830d-1,2.45220d-1,2.71260d-1,2.95940d-1,3.20590d-1,3.44590d-1,&
      3.66380d-1,3.85290d-1,4.02490d-1,4.17410d-1,4.30070d-1,4.40120d-1,&
      4.47510d-1,4.55340d-1,4.59120d-1,4.63940d-1,4.66210d-1,4.67680d-1,&
      4.70120d-1,4.70570d-1,4.72180d-1,4.73500d-1,4.75110d-1,4.75610d-1,&
      4.76660d-1,4.77150d-1,4.77350d-1,4.77890d-1,4.78970d-1,4.79310d-1,&
      4.80410d-1,4.81470d-1,4.81580d-1,4.83040d-1,4.82920d-1,4.82570d-1,&
      4.83820d-1,4.84220d-1,4.85050d-1,4.84690d-1,4.85190d-1,4.84650d-1,&
      4.84290d-1,4.83720d-1,4.84160d-1,4.83160d-1,4.81950d-1,4.79380d-1,&
      4.76860d-1,4.73220d-1,4.62790d-1,4.43970d-1,4.14520d-1,3.73500d-1,&
      3.19880d-1,2.59970d-1,1.96280d-1,1.36860d-1,8.72900d-2,5.09020d-2,&
      2.65530d-2,1.23380d-2,5.26700d-3,1.86600d-3,6.02000d-4,1.50000d-4,&
      4.60000d-5,3.00000d-6,3.00000d-6,0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0 /)
  ! Efficiency (first and only interval)
  REAL*8, PARAMETER :: EFF1(NE) = EFF0
  ! Efficiencies array (2D)
  REAL*8, PARAMETER :: EFF(NE,0:NEFF)                                   &
      = RESHAPE( (/ EFF0(:), EFF1(:) /) ,SHAPE(EFF))
  
  ! One call for all settings.
  ! Most of these _must_ be there to ensure everything get initialized.
  CALL SetDetector(D,mass=20d3,time=2d0*365d0,Nevents=0,                &
                   background=0.5d0,Nelem=1,Zelem=(/18/),               &
                   NEeff=NE,Eeff=E,Neff=NEFF,eff=EFF,                   &
                   intervals=intervals)
  D%eff_file = '[DARWIN Ar 2015]'
  
END SUBROUTINE



!=======================================================================
! DARWIN XENON ANALYSIS ROUTINES
! Based upon a DARWIN xenon-based analysis (2015 estimated parameters)
! [15MM.NNNNN].  Zero events assumed in the analysis region.
!=======================================================================

!-----------------------------------------------------------------------
! Initializes the module to perform calculations for a DARWIN
! xenon-based analysis (2015 estimates).  This must be called if any
! of the following DARWIN routines are to be used.
! 
! Required input arguments:
!     intervals   Indicates if sub-intervals should be included.
!                 Only necessary if confidence intervals using the
!                 maximum gap method are desired.
! 
SUBROUTINE DARWIN_Xe_2015_Init(intervals)
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: intervals
  CALL DARWIN_Xe_2015_InitTo(DARWIN_Xe_2015,intervals)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_DARWIN_Xe_2015_Init(intervals) &
           BIND(C,NAME='C_DDCALC0_darwin_xe_2015_init')
  USE ISO_C_BINDING, only: C_BOOL
  IMPLICIT NONE
  LOGICAL(KIND=C_BOOL), INTENT(IN) :: intervals
  CALL DARWIN_Xe_2015_Init(LOGICAL(intervals))
END SUBROUTINE


! ----------------------------------------------------------------------
! Sets the minimum recoil energy to be included in the calculations.
! Note the efficiency curves already account for detector and analysis
! thresholds regardless of this setting, so setting this to 0 keV (the
! default behavior when initialization is performed) does not imply
! that very low energy recoils actually contribute to the signal.
! 
! Required input arguments:
!     Emin        The minimum recoil energy to consider [keV]
! 
SUBROUTINE DARWIN_Xe_2015_SetEmin(Emin)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: Emin
  CALL SetDetector(DARWIN_Xe_2015,Emin=Emin)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_DARWIN_Xe_2015_SetEmin(Emin) &
           BIND(C,NAME='C_DDCALC0_darwin_xe_2015_setemin')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE), INTENT(IN) :: Emin
  CALL DARWIN_Xe_2015_SetEmin(REAL(Emin,KIND=8))
END SUBROUTINE


! ----------------------------------------------------------------------
! Calculates various rate quantities using the current WIMP.
! Must be called each time the WIMP parameters are modified.
! 
SUBROUTINE DARWIN_Xe_2015_CalcRates()
  IMPLICIT NONE
  CALL CalcRates(DARWIN_Xe_2015)
END SUBROUTINE

! C++ interface wrapper
SUBROUTINE C_DARWIN_Xe_2015_CalcRates() &
           BIND(C,NAME='C_DDCALC0_darwin_xe_2015_calcrates')
  IMPLICIT NONE
  CALL DARWIN_Xe_2015_CalcRates()
END SUBROUTINE


! ----------------------------------------------------------------------
! Returns the observed number of events.
! 
FUNCTION DARWIN_Xe_2015_Events() RESULT(N)
  IMPLICIT NONE
  INTEGER :: N
  CALL GetRates(DARWIN_Xe_2015,Nevents=N)
END FUNCTION

! C++ interface wrapper
FUNCTION C_DARWIN_Xe_2015_Events() RESULT(N) &
         BIND(C,NAME='C_DDCALC0_darwin_xe_2015_events')
  USE ISO_C_BINDING, only: C_INT
  IMPLICIT NONE
  INTEGER(KIND=C_INT) :: N
  ! Automatic type conversions here
  N = DARWIN_Xe_2015_Events()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of background events.
! 
FUNCTION DARWIN_Xe_2015_Background() RESULT(b)
  IMPLICIT NONE
  REAL*8 :: b
  CALL GetRates(DARWIN_Xe_2015,background=b)
END FUNCTION

! C++ interface wrapper
FUNCTION C_DARWIN_Xe_2015_Background() RESULT(b) &
         BIND(C,NAME='C_DDCALC0_darwin_xe_2015_background')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: b
  ! Automatic type conversions here
  b = DARWIN_Xe_2015_Background()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of signal events for the
! current WIMP.
! 
FUNCTION DARWIN_Xe_2015_Signal() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(DARWIN_Xe_2015,signal=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_DARWIN_Xe_2015_Signal() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_darwin_xe_2015_signal')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = DARWIN_Xe_2015_Signal()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of spin-independent signal events
! for the current WIMP.
! 
FUNCTION DARWIN_Xe_2015_SignalSI() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(DARWIN_Xe_2015,signal_si=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_DARWIN_Xe_2015_SignalSI() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_darwin_xe_2015_signalsi')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = DARWIN_Xe_2015_SignalSI()
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of spin-dependent signal events
! for the current WIMP.
! 
FUNCTION DARWIN_Xe_2015_SignalSD() RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  CALL GetRates(DARWIN_Xe_2015,signal_sd=s)
END FUNCTION

! C++ interface wrapper
FUNCTION C_DARWIN_Xe_2015_SignalSD() RESULT(s) &
         BIND(C,NAME='C_DDCALC0_darwin_xe_2015_signalsd')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: s
  ! Automatic type conversions here
  s = DARWIN_Xe_2015_SignalSD()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the log-likelihood for the current WIMP mass and couplings.
! Uses a Poisson distribution in the number of observed events N:
!    P(N|s+b)
! where s is the average expected signal and b is the average expected
! background.
! 
FUNCTION DARWIN_Xe_2015_LogLikelihood() RESULT(lnlike)
  IMPLICIT NONE
  REAL*8 :: lnlike
  lnlike = LogLikelihood(DARWIN_Xe_2015)
END FUNCTION

! C++ interface wrapper
FUNCTION C_DARWIN_Xe_2015_LogLikelihood() RESULT(lnlike) &
         BIND(C,NAME='C_DDCALC0_darwin_xe_2015_loglikelihood')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: lnlike
  ! Automatic type conversions here
  lnlike = DARWIN_Xe_2015_LogLikelihood()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the log of the p-value for the current WIMP mass and
! couplings (NO BACKGROUND SUBTRACTION).  Uses the maximum gap method
! if DARWIN_Xe_2015_Init was called with argument intervals=.TRUE.,
! otherwise uses a Poisson distribution in the number of observed
! events N:
!    P(N|s)
! where s is the average expected signal (background contributions are
! ignored).
! 
FUNCTION DARWIN_Xe_2015_LogPValue() RESULT(lnp)
  IMPLICIT NONE
  REAL*8 :: lnp
  lnp = LogPValue(DARWIN_Xe_2015)
END FUNCTION

! C++ interface wrapper
FUNCTION C_DARWIN_Xe_2015_LogPValue() RESULT(lnp) &
         BIND(C,NAME='C_DDCALC0_darwin_xe_2015_logpvalue')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: lnp
  ! Automatic type conversions here
  lnp = DARWIN_Xe_2015_LogPValue()
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the factor x by which the cross-sections must be scaled
! (sigma -> x*sigma) to achieve the desired p-value (given as log(p)).
! See LogPValue() above for a description of the statistics.
! 
! Required input argument:
!   lnp         The logarithm of the desired p-value (p = 1-CL).
! 
FUNCTION DARWIN_Xe_2015_ScaleToPValue(lnp) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN) :: lnp
  x = ScaleToPValue(DARWIN_Xe_2015,lnp)
END FUNCTION

! C++ interface wrapper
FUNCTION C_DARWIN_Xe_2015_ScaleToPValue(lnp) RESULT(x) &
         BIND(C,NAME='C_DDCALC0_darwin_xe_2015_scaletopvalue')
  USE ISO_C_BINDING, only: C_DOUBLE
  IMPLICIT NONE
  REAL(KIND=C_DOUBLE) :: x
  REAL(KIND=C_DOUBLE), INTENT(IN) :: lnp
  ! Automatic type conversions here
  x = DARWIN_Xe_2015_ScaleToPValue(REAL(lnp,KIND=8))
END FUNCTION


!-----------------------------------------------------------------------
! INTERNAL ROUTINE.
! Initializes the given DetectorStruct to the proposed DARWIN
! xenon-based detector (2015 estimate).  This is meant as an internal
! routine; external access should be through
! DARWIN_Xe_2015_Init instead.
! 
! The efficiencies used here were generated using TPCMC.
! 
! Required input arguments:
!     D           The DetectorStruct to initialize
!     intervals   Indicates if sub-intervals should be included
! 
SUBROUTINE DARWIN_Xe_2015_InitTo(D,intervals)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(OUT) :: D
  LOGICAL, INTENT(IN) :: intervals
  INTEGER, PARAMETER :: NE = 151
  INTEGER, PARAMETER :: NEFF = 1
  ! Efficiency curves energy tabulation points
  REAL*8, PARAMETER :: E(NE)                                            &
      =       (/ 0.10000d0, 0.10471d0, 0.10965d0, 0.11482d0, 0.12023d0, &
      0.12589d0, 0.13183d0, 0.13804d0, 0.14454d0, 0.15136d0, 0.15849d0, &
      0.16596d0, 0.17378d0, 0.18197d0, 0.19055d0, 0.19953d0, 0.20893d0, &
      0.21878d0, 0.22909d0, 0.23988d0, 0.25119d0, 0.26303d0, 0.27542d0, &
      0.28840d0, 0.30200d0, 0.31623d0, 0.33113d0, 0.34674d0, 0.36308d0, &
      0.38019d0, 0.39811d0, 0.41687d0, 0.43652d0, 0.45709d0, 0.47863d0, &
      0.50119d0, 0.52481d0, 0.54954d0, 0.57544d0, 0.60256d0, 0.63096d0, &
      0.66069d0, 0.69183d0, 0.72444d0, 0.75858d0, 0.79433d0, 0.83176d0, &
      0.87096d0, 0.91201d0, 0.95499d0, 1.0000d0,  1.0471d0,  1.0965d0,  &
      1.1482d0,  1.2023d0,  1.2589d0,  1.3183d0,  1.3804d0,  1.4454d0,  &
      1.5136d0,  1.5849d0,  1.6596d0,  1.7378d0,  1.8197d0,  1.9055d0,  &
      1.9953d0,  2.0893d0,  2.1878d0,  2.2909d0,  2.3988d0,  2.5119d0,  &
      2.6303d0,  2.7542d0,  2.8840d0,  3.0200d0,  3.1623d0,  3.3113d0,  &
      3.4674d0,  3.6308d0,  3.8019d0,  3.9811d0,  4.1687d0,  4.3652d0,  &
      4.5709d0,  4.7863d0,  5.0119d0,  5.2481d0,  5.4954d0,  5.7544d0,  &
      6.0256d0,  6.3096d0,  6.6069d0,  6.9183d0,  7.2444d0,  7.5858d0,  &
      7.9433d0,  8.3176d0,  8.7096d0,  9.1201d0,  9.5499d0, 10.000d0,   &
     10.471d0,  10.965d0,  11.482d0,  12.023d0,  12.589d0,  13.183d0,   &
     13.804d0,  14.454d0,  15.136d0,  15.849d0,  16.596d0,  17.378d0,   &
     18.197d0,  19.055d0,  19.953d0,  20.893d0,  21.878d0,  22.909d0,   &
     23.988d0,  25.119d0,  26.303d0,  27.542d0,  28.840d0,  30.200d0,   &
     31.623d0,  33.113d0,  34.674d0,  36.308d0,  38.019d0,  39.811d0,   &
     41.687d0,  43.652d0,  45.709d0,  47.863d0,  50.119d0,  52.481d0,   &
     54.954d0,  57.544d0,  60.256d0,  63.096d0,  66.069d0,  69.183d0,   &
     72.444d0,  75.858d0,  79.433d0,  83.176d0,  87.096d0,  91.201d0,   &
     95.499d0, 100.00d0 /)
  ! Efficiency (total)
  REAL*8, PARAMETER :: EFF0(NE)                                         &
      =       (/ 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 1.00000d-5,2.00000d-5,0.00000d0, &
      1.00000d-5,2.00000d-5,1.00000d-5,3.00000d-5,6.00000d-5,4.00000d-5,&
      1.60000d-4,2.50000d-4,3.30000d-4,4.30000d-4,5.80000d-4,8.30000d-4,&
      1.29000d-3,1.73000d-3,2.61000d-3,3.38000d-3,4.71000d-3,6.67000d-3,&
      8.93000d-3,1.18500d-2,1.66700d-2,2.11800d-2,2.76600d-2,3.70100d-2,&
      4.62800d-2,5.57800d-2,7.03000d-2,8.48500d-2,1.00920d-1,1.19770d-1,&
      1.40370d-1,1.59900d-1,1.83670d-1,2.13730d-1,2.34620d-1,2.60470d-1,&
      2.90490d-1,3.10700d-1,3.37130d-1,3.59870d-1,3.81220d-1,3.97590d-1,&
      4.13050d-1,4.29350d-1,4.33960d-1,4.43460d-1,4.50340d-1,4.54710d-1,&
      4.58700d-1,4.64210d-1,4.63840d-1,4.65210d-1,4.67910d-1,4.69700d-1,&
      4.71550d-1,4.70310d-1,4.72880d-1,4.73510d-1,4.77040d-1,4.76320d-1,&
      4.78180d-1,4.74830d-1,4.77390d-1,4.79300d-1,4.80260d-1,4.81980d-1,&
      4.79980d-1,4.83400d-1,4.85070d-1,4.83660d-1,4.86320d-1,4.83460d-1,&
      4.84910d-1,4.85100d-1,4.81940d-1,4.83840d-1,4.83400d-1,4.78140d-1,&
      4.77990d-1,4.70180d-1,4.48760d-1,4.18240d-1,3.67580d-1,3.03320d-1,&
      2.31100d-1,1.54890d-1,9.42000d-2,5.02400d-2,2.18100d-2,8.18000d-3,&
      2.44000d-3,7.40000d-4,1.70000d-4,1.00000d-5,1.00000d-5,0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, 0.00000d0, &
      0.00000d0, 0.00000d0 /)
  ! Efficiency (first and only interval)
  REAL*8, PARAMETER :: EFF1(NE) = EFF0
  ! Efficiencies array (2D)
  REAL*8, PARAMETER :: EFF(NE,0:NEFF)                                   &
      = RESHAPE( (/ EFF0(:), EFF1(:) /) ,SHAPE(EFF))
  
  ! One call for all settings.
  ! Most of these _must_ be there to ensure everything get initialized.
  CALL SetDetector(D,mass=12d3,time=2d0*365d0,Nevents=0,                &
                   background=0.5d0,Nelem=1,Zelem=(/54/),               &
                   NEeff=NE,Eeff=E,Neff=NEFF,eff=EFF,                   &
                   intervals=intervals)
  D%eff_file = '[DARWIN Xe 2015]'
  
END SUBROUTINE



!=======================================================================
! MAIN ROUTINES
!=======================================================================

! ----------------------------------------------------------------------
! Generic main routine.
! 
SUBROUTINE main()
  IMPLICIT NONE
  LOGICAL :: interactive
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Load arguments
  CALL ReadArguments()
  
  ! Interactive mode?
  ! True if --interactive given or if the WIMP mass is not specified
  ! through the command-line (either with --m=<val> option or non-option
  ! argument).
  interactive = GetLongArg('interactive') .OR.                          &
                ((Arguments%Nparameters .EQ. 0) .AND. .NOT. GetLongArg('m'))
  
  ! Determine program mode =============
  
  ! ------------------------------------
  ! Calculates the log-likelihood (w/ background)
  IF (GetLongArg('log-likelihood')) THEN
    IF (interactive) THEN
      CALL MainLogLikelihoodInteractive()
    ELSE
      CALL MainLogLikelihood()
    END IF
  ! ------------------------------------
  ! Calculates the log of the p-value (no background subtraction)
  ELSE IF (GetLongArg('log-pvalue')) THEN
    IF (interactive) THEN
      CALL MainLogPValueInteractive()
    ELSE
      CALL MainLogPValue()
    END IF
  ! ------------------------------------
  ! Prints the raw recoil spectrum dR/dE
  ELSE IF (GetLongArg('spectrum')) THEN
    CALL MainSpectrum()
  ! ------------------------------------
  ! Calculates expected events, tabulated by WIMP mass
  ELSE IF (GetLongArg('events-by-mass')) THEN
    CALL MainEventsByMass()
  ! ------------------------------------
  ! Calculates spin-independent likelihood contraints,
  ! tabulated by WIMP mass
  ELSE IF (GetLongArg('constraints-SI')) THEN
    CALL MainConstraintsSI()
  ! ------------------------------------
  ! Calculates spin-dependent likelihood contraints,
  ! tabulated by WIMP mass
  ELSE IF (GetLongArg('constraints-SD')) THEN
    CALL MainConstraintsSD()
  ! ------------------------------------
  ! Calculates spin-independent no-background-subtraction limits,
  ! tabulated by WIMP mass
  ELSE IF (GetLongArg('limits-SI')) THEN
    CALL MainLimitsSI()
  ! ------------------------------------
  ! Calculates spin-dependent no-background-subtraction limits,
  ! tabulated by WIMP mass
  ELSE IF (GetLongArg('limits-SD')) THEN
    CALL MainLimitsSD()
  ! ------------------------------------
  ! Default case:
  ! Calculate both expected events and likelihoods
  ELSE
    IF (interactive) THEN
      CALL MainEventsAndLikelihoodsInteractive()
    ELSE IF (Arguments%Nparameters .GE. 1) THEN
      CALL MainEventsAndLikelihoods()
    ELSE
      CALL MainEventsAndLikelihoodsInteractive()
    END IF
  END IF
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate log-likelihood.
! 
SUBROUTINE MainLogLikelihood()
  IMPLICIT NONE
  REAL*8 :: lnLike
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes various global structures.
  ! Will calculate only total rates (not sub-intervals).
  CALL InitializeCL(intervals=.FALSE.)
  
  ! Print header
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteWIMPHeader()
    CALL WriteHaloHeader()
    CALL WriteDetectorHeader()
    CALL WriteLogLikelihoodHeader()
  END IF
  
  ! Do rate calculations
  CALL CalcRates()
  
  ! Get log-likelihood
  lnLike = LogLikelihood()
  
  ! Print results
  WRITE(*,*)  lnLike
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate log-likelihood (interactive mode).
! Reads one line of input containing WIMP parameters and writes out the
! corresponding log-likelihood, terminating when the input stream ends
! (EOF) or a blank line is given.
! 
SUBROUTINE MainLogLikelihoodInteractive()
  IMPLICIT NONE
  REAL*8 :: lnLike
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes various global structures.
  ! Will calculate only total rates (not sub-intervals).
  CALL InitializeCL(intervals=.FALSE.)
  
  ! Print header
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteHaloHeader()
    CALL WriteDetectorHeader()
    CALL WriteLogLikelihoodHeader()
  END IF
  
  ! Print instructions
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteInteractiveHeader()
  END IF
  
  ! Cycle over input.
  ! ParseWIMPInput() reads line containing WIMP parameters from
  ! standard input and parses it, returning false if a blank line
  ! or EOF is found.  See ParseWIMPInput() for description of input
  ! line.
  DO WHILE (ParseWIMPInput())
    ! Do rate calculations
    CALL CalcRates()
    ! Get log-likelihood
    lnLike = LogLikelihood()
    ! Print results to standard output
    WRITE(*,*)  lnLike
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate log of p-value using the maximum gap method.
! See S. Yellin, PRD 66, 032005 (2002) [physics/0203002].  Uses Poisson
! distribution if efficiencies for sub-intervals are not available (but
! no background subtraction!).
! 
SUBROUTINE MainLogPValue()
  IMPLICIT NONE
  REAL*8 :: lnP
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes various global structures.
  ! Will calculate rates for sub-intervals as these are used by
  ! maximum gap method.
  CALL InitializeCL(intervals=.TRUE.)
  
  ! Print header
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteWIMPHeader()
    CALL WriteHaloHeader()
    CALL WriteDetectorHeader()
    CALL WriteLogPValueHeader()
  END IF
  
  ! Do rate calculations
  CALL CalcRates()
  
  ! Get log of p-value
  lnP = LogPValue()
  
  ! Print results
  WRITE(*,*)  lnP
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate log of p-value using the maximum gap method
! (interactive mode).  See S. Yellin, PRD 66, 032005 (2002)
! [physics/0203002].  Uses Poisson distribution if efficiencies for
! sub-intervals are not available (but no background subtraction!).
! Reads one line of input containing WIMP parameters and writes out the
! corresponding log of the p-value, terminating when the input stream
! ends (EOF) or a blank line is given.
! 
SUBROUTINE MainLogPValueInteractive()
  IMPLICIT NONE
  REAL*8 :: lnP
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes various global structures.
  ! Will calculate rates for sub-intervals as these are used by
  ! maximum gap method.
  CALL InitializeCL(intervals=.TRUE.)
  
  ! Print header
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteHaloHeader()
    CALL WriteDetectorHeader()
    CALL WriteLogPValueHeader()
  END IF
  
  ! Print instructions
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteInteractiveHeader()
  END IF
  
  ! Cycle over input.
  ! ParseWIMPInput() reads line containing WIMP parameters from
  ! standard input and parses it, returning false if a blank line
  ! or EOF is found.  See ParseWIMPInput() for description of input
  ! line.
  DO WHILE (ParseWIMPInput())
    ! Do rate calculations
    CALL CalcRates()
    ! Get log of p-value
    lnP = LogPValue()
    ! Print results to standard output
    WRITE(*,*)  lnP
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to print expected events and likelihoods.
! 
SUBROUTINE MainEventsAndLikelihoods()
  IMPLICIT NONE
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes various global structures.
  CALL InitializeCL()
  
  ! Print header
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    IF (VerbosityLevel .EQ. 2) CALL WriteWIMPHeader()
    CALL WriteHaloHeader()
    CALL WriteDetectorHeader()
  END IF
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteEventsAndLikelihoodsHeader()
    CALL WriteEventsAndLikelihoodsColumnHeader()
  END IF
  
  ! Do rate calculations
  CALL CalcRates()
  
  ! Print results
  CALL WriteEventsAndLikelihoodsData()
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to print expected events and likelihoods.
! 
SUBROUTINE MainEventsAndLikelihoodsInteractive()
  IMPLICIT NONE
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes various global structures.
  CALL InitializeCL()
  
  ! Print header
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteHaloHeader()
    CALL WriteDetectorHeader()
    CALL WriteEventsAndLikelihoodsHeader()
  END IF
  
  ! Print instructions
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteInteractiveHeader(1)
  END IF
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteEventsAndLikelihoodsColumnHeader()
  END IF
  
  ! Cycle over input.
  ! ParseWIMPInput() reads line containing WIMP parameters from
  ! standard input and parses it, returning false if a blank line
  ! or EOF is found.  See ParseWIMPInput() for description of input
  ! line.
  DO WHILE (ParseWIMPInput())
    ! Do rate calculations
    CALL CalcRates()
    ! Print results to standard output
    CALL WriteEventsAndLikelihoodsData()
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate raw differential rates as a function of
! energy, i.e. the raw recoil energy spectrum (not including
! efficiencies or energy resolution).
! 
SUBROUTINE MainSpectrum()
  IMPLICIT NONE
  LOGICAL :: use_log
  INTEGER :: NE,K
  REAL*8 :: Emin,Emax
  REAL*8, ALLOCATABLE :: E(:),Eeff(:),eff(:,:)
  TYPE(TabulationStruct) :: TS
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes various global structures.
  CALL InitializeCL()
  
  ! Set tabulation energies; set efficiencies to 1.
  ! Note same results can be achieved by simply giving NE=-1
  ! and NEeff=-1 arguments to SetDetector (apart from any
  ! --E-tabulation command line specification).
  Emin = 0.1d0
  Emax = 1000d0
  NE   = -50
  use_log = .TRUE.
  CALL GetTabulationArgs('E-tabulation',Emin,Emax,NE,use_log)
  CALL InitTabulation(TS,Emin,Emax,NE,.TRUE.)
  NE = TS%N+2
  ALLOCATE(E(1:NE),Eeff(1:2),eff(1:2,0:0))
  E(1) = 0d0
  DO K = 2,NE
    E(K) = TabulationValue(TS,K-2)
  END DO
  Eeff = (/ 0d0, HUGE(Eeff) /)
  eff = 1d0
  CALL SetDetector(NE=NE,E=E,NEeff=2,Eeff=Eeff,Neff=0,eff=eff)
  
  ! For high verbosity level, we are printing reference rates, so
  ! set "actual" rates to same.
  IF (VerbosityLevel .GE. 4) CALL SetWIMP(sigmaSI=1d0,sigmaSD=1d0)
  
  ! Do rate calculations
  CALL CalcRates()
  
  ! Write out header.
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteWIMPHeader()
    CALL WriteHaloHeader()
    CALL WriteDetectorHeader()
    CALL WriteSpectrumHeader()
  END IF
  
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteSpectrumColumnHeader()
  END IF
  
  ! Write out table.
  CALL WriteSpectrumData()
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate events as a function of mass.
! Given for fixed WIMP-nucleon cross-sections.
! 
SUBROUTINE MainEventsByMass()
  IMPLICIT NONE
  LOGICAL :: use_log
  INTEGER :: Nm,I
  REAL*8 :: m,mmin,mmax,sigmapSI,sigmanSI,sigmapSD,sigmanSD
  TYPE(TabulationStruct) :: TS
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes various global structures.
  CALL InitializeCL()
  
  ! Determine WIMP mass tabulation
  mmin    =    1d0
  mmax    = 1000d0
  Nm      =  -20
  use_log = .TRUE.
  CALL GetTabulationArgs('m-tabulation',mmin,mmax,Nm,use_log)
  CALL InitTabulation(TS,mmin,mmax,Nm,use_log)
  
  ! Get fixed cross-sections; will need to be reset at each mass.
  ! Set to negative if coupling is negative (intended meaning for
  ! negative cross-sections as input).
  CALL GetWIMP(sigmapSI=sigmapSI,sigmanSI=sigmanSI,sigmapSD=sigmapSD,sigmanSD=sigmanSD)
  IF (WIMP%GpSI .LT. 0d0) sigmapSI = -ABS(sigmapSI)
  IF (WIMP%GnSI .LT. 0d0) sigmanSI = -ABS(sigmanSI)
  IF (WIMP%GpSD .LT. 0d0) sigmapSD = -ABS(sigmapSD)
  IF (WIMP%GnSD .LT. 0d0) sigmanSD = -ABS(sigmanSD)
  
  ! Write out header.
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteHaloHeader()
    CALL WriteDetectorHeader()
    CALL WriteEventsByMassHeader()
  END IF
  
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteEventsByMassColumnHeader()
  END IF
  
  ! Cycle over masses
  DO I = 0,TS%N
    ! Set WIMP mass and cross-sections
    m = TabulationValue(TS,I)
    CALL SetWIMP(m=m,sigmapSI=sigmapSI,sigmanSI=sigmanSI,sigmapSD=sigmapSD,sigmanSD=sigmanSD)
    ! Do rate calculations
    CALL CalcRates()
    ! Write out table data line.
    CALL WriteEventsByMassData()
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate likelihood constraints as a function of
! mass for SI couplings.
! 
SUBROUTINE MainConstraintsSI()
  IMPLICIT NONE
  LOGICAL :: use_log
  INTEGER :: Nm,I
  REAL*8 :: m,mmin,mmax,x,lnp,thetaG,Gp,Gn,s1,s2
  TYPE(TabulationStruct) :: TS
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes various global structures.
  ! Do not need sub-intervals as we only use the total events.
  CALL InitializeCL(intervals=.FALSE.)
  
  ! Determine WIMP mass tabulation
  mmin    =    1d0
  mmax    = 1000d0
  Nm      =  -20
  use_log = .TRUE.
  CALL GetTabulationArgs('m-tabulation',mmin,mmax,Nm,use_log)
  CALL InitTabulation(TS,mmin,mmax,Nm,use_log)
  
  ! Get p-value for exclusion limit
  CALL ReadLogPValue(lnp)
  
  ! Get angle of (Gp,Gn), which will be kept fixed.
  IF (GetLongArgReal('theta-SI',x)) THEN
    thetaG = x
  ELSE IF (GetLongArgReal('theta-SI-pi',x)) THEN
    thetaG = PI*x
  ELSE
    thetaG = 0.25d0*PI
  END IF
  
  ! Will initially use these couplings at every mass;
  ! constraint calculations will perform appropriate rescaling.
  Gp = 1d0*COS(thetaG)
  Gn = 1d0*SIN(thetaG)
  IF (ABS(Gp) .LT. 1d-8) Gp = 0d0
  IF (ABS(Gn) .LT. 1d-8) Gn = 0d0
  
  ! Calculate allowed signal rates
  CALL FeldmanCousinsPoissonCI(lnp,DefaultDetector%Nevents,             &
                               DefaultDetector%MuBackground,s1,s2)
  
  ! Write out header.
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteHaloHeader()
    CALL WriteDetectorHeader()
    CALL WriteConstraintsSIHeader(lnp,thetaG,s1,s2)
  END IF
  
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteConstraintsSIColumnHeader()
  END IF
  
  ! Cycle over masses
  DO I = 0,TS%N
    ! Set WIMP mass and cross-sections
    m = TabulationValue(TS,I)
    CALL SetWIMP(m=m,GpSI=Gp,GnSI=Gn,GpSD=0d0,GnSD=0d0)
    ! Do rate calculations
    CALL CalcRates()
    ! Write out table data line.
    CALL WriteConstraintsSIData(s1,s2)
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate likelihood constraints as a function of
! mass for SD couplings.
! 
SUBROUTINE MainConstraintsSD()
  IMPLICIT NONE
  LOGICAL :: use_log
  INTEGER :: Nm,I
  REAL*8 :: m,mmin,mmax,x,lnp,thetaG,Gp,Gn,s1,s2
  TYPE(TabulationStruct) :: TS
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes various global structures.
  ! Do not need sub-intervals as we only use the total events.
  CALL InitializeCL(intervals=.FALSE.)
  
  ! Determine WIMP mass tabulation
  mmin    =    1d0
  mmax    = 1000d0
  Nm      =  -20
  use_log = .TRUE.
  CALL GetTabulationArgs('m-tabulation',mmin,mmax,Nm,use_log)
  CALL InitTabulation(TS,mmin,mmax,Nm,use_log)
  
  ! Get p-value for exclusion limit
  CALL ReadLogPValue(lnp)
  
  ! Get angle of (Gp,Gn), which will be kept fixed.
  IF (GetLongArgReal('theta-SD',x)) THEN
    thetaG = x
  ELSE IF (GetLongArgReal('theta-SD-pi',x)) THEN
    thetaG = PI*x
  ELSE
    thetaG = 0.25d0*PI
  END IF
  
  ! Will initially use these couplings at every mass;
  ! constraint calculations will perform appropriate rescaling.
  Gp = 1d0*COS(thetaG)
  Gn = 1d0*SIN(thetaG)
  IF (ABS(Gp) .LT. 1d-8) Gp = 0d0
  IF (ABS(Gn) .LT. 1d-8) Gn = 0d0
  
  ! Calculate allowed signal rates
  CALL FeldmanCousinsPoissonCI(lnp,DefaultDetector%Nevents,             &
                               DefaultDetector%MuBackground,s1,s2)
  
  ! Write out header.
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteHaloHeader()
    CALL WriteDetectorHeader()
    CALL WriteConstraintsSDHeader(lnp,thetaG,s1,s2)
  END IF
  
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteConstraintsSDColumnHeader()
  END IF
  
  ! Cycle over masses
  DO I = 0,TS%N
    ! Set WIMP mass and cross-sections
    m = TabulationValue(TS,I)
    CALL SetWIMP(m=m,GpSI=0d0,GnSI=0d0,GpSD=Gp,GnSD=Gn)
    ! Do rate calculations
    CALL CalcRates()
    ! Write out table data line.
    CALL WriteConstraintsSDData(s1,s2)
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate exclusion limits as a function of mass for
! SI couplings.
! 
SUBROUTINE MainLimitsSI()
  IMPLICIT NONE
  LOGICAL :: use_log
  INTEGER :: Nm,I
  REAL*8 :: m,mmin,mmax,x,lnp,thetaG,Gp,Gn
  TYPE(TabulationStruct) :: TS
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes various global structures.
  ! Will calculate rates for sub-intervals as these are used by
  ! maximum gap method.
  CALL InitializeCL(intervals=.TRUE.)
  
  ! Determine WIMP mass tabulation
  mmin    =    1d0
  mmax    = 1000d0
  Nm      =  -20
  use_log = .TRUE.
  CALL GetTabulationArgs('m-tabulation',mmin,mmax,Nm,use_log)
  CALL InitTabulation(TS,mmin,mmax,Nm,use_log)
  
  ! Get p-value for exclusion limit
  CALL ReadLogPValue(lnp)
  
  ! Get angle of (Gp,Gn), which will be kept fixed.
  IF (GetLongArgReal('theta-SI',x)) THEN
    thetaG = x
  ELSE IF (GetLongArgReal('theta-SI-pi',x)) THEN
    thetaG = PI*x
  ELSE
    thetaG = 0.25d0*PI
  END IF
  
  ! Will initially use these couplings at every mass;
  ! limit calculations will perform appropriate rescaling.
  Gp = 1d0*COS(thetaG)
  Gn = 1d0*SIN(thetaG)
  IF (ABS(Gp) .LT. 1d-8) Gp = 0d0
  IF (ABS(Gn) .LT. 1d-8) Gn = 0d0
  
  ! Write out header.
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteHaloHeader()
    CALL WriteDetectorHeader()
    CALL WriteLimitsSIHeader(lnp,thetaG)
  END IF
  
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteLimitsSIColumnHeader()
  END IF
  
  ! Cycle over masses
  DO I = 0,TS%N
    ! Set WIMP mass and cross-sections
    m = TabulationValue(TS,I)
    CALL SetWIMP(m=m,GpSI=Gp,GnSI=Gn,GpSD=0d0,GnSD=0d0)
    ! Do rate calculations
    CALL CalcRates()
    ! Write out table data line.
    CALL WriteLimitsSIData(lnp)
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Main routine to calculate exclusion limits as a function of mass for
! SD couplings.
! 
SUBROUTINE MainLimitsSD()
  IMPLICIT NONE
  LOGICAL :: use_log
  INTEGER :: Nm,I
  REAL*8 :: m,mmin,mmax,x,lnp,thetaG,Gp,Gn
  TYPE(TabulationStruct) :: TS
  
  ! Show usage and exit
  IF (ShowUsageRequested()) THEN
    CALL ShowUsage()
    STOP
  END IF
  
  ! Initializes various global structures.
  ! Will calculate rates for sub-intervals as these are used by
  ! maximum gap method.
  CALL InitializeCL(intervals=.TRUE.)
  
  ! Determine WIMP mass tabulation
  mmin    =    1d0
  mmax    = 1000d0
  Nm      =  -20
  use_log = .TRUE.
  CALL GetTabulationArgs('m-tabulation',mmin,mmax,Nm,use_log)
  CALL InitTabulation(TS,mmin,mmax,Nm,use_log)
  
  ! Get p-value for exclusion limit
  CALL ReadLogPValue(lnp)
  
  ! Get angle of (Gp,Gn), which will be kept fixed.
  IF (GetLongArgReal('theta-SD',x)) THEN
    thetaG = x
  ELSE IF (GetLongArgReal('theta-SD-pi',x)) THEN
    thetaG = PI*x
  ELSE
    thetaG = 0.25d0*PI
  END IF
  
  ! Will initially use these couplings at every mass;
  ! limit calculations will perform appropriate rescaling.
  Gp = 1d0*COS(thetaG)
  Gn = 1d0*SIN(thetaG)
  IF (ABS(Gp) .LT. 1d-8) Gp = 0d0
  IF (ABS(Gn) .LT. 1d-8) Gn = 0d0
  
  ! Write out header.
  IF (VerbosityLevel .GE. 2) THEN
    CALL WriteCommandHeader()
    CALL WriteHaloHeader()
    CALL WriteDetectorHeader()
    CALL WriteLimitsSDHeader(lnp,thetaG)
  END IF
  
  IF (VerbosityLevel .GE. 1) THEN
    CALL WriteLimitsSDColumnHeader()
  END IF
  
  ! Cycle over masses
  DO I = 0,TS%N
    ! Set WIMP mass and cross-sections
    m = TabulationValue(TS,I)
    CALL SetWIMP(m=m,GpSI=0d0,GnSI=0d0,GpSD=Gp,GnSD=Gn)
    ! Do rate calculations
    CALL CalcRates()
    ! Write out table data line.
    CALL WriteLimitsSDData(lnp)
  END DO
  
END SUBROUTINE



!=======================================================================
! INITIALIZATION
!=======================================================================

! ----------------------------------------------------------------------
! Initialization routine that sets various parameters to default values.
! This must be called before using the module (unless one of the other
! initialization routines is called).
! 
! The current default initialization is to that of the LUX 2013
! analysis.  If a different setup is desired, the various Set() methods
! must be used.
! 
! Optional input arguments:
!   intervals   Indicates if rates for sub-intervals should be
!               calculated as well as the total rate (default: true).
!               These are unnecessary for some calculations.
! 
SUBROUTINE Initialize(intervals)
  IMPLICIT NONE
  LOGICAL, INTENT(IN), OPTIONAL :: intervals
  LOGICAL :: intervals0
  
  intervals0 = .TRUE.
  IF (PRESENT(intervals)) intervals0 = intervals
  
  ! Initialize WIMP mass and couplings.
  CALL InitWIMP()
  
  ! Initialize halo.
  CALL InitHalo()
  
  ! Initialize detector isotopes, efficiencies, array sizing, etc.
  CALL InitDetector(intervals=intervals0)
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Initialization routine that sets various parameters to default values
! or values specified on the command line.  This must be called before
! using the module (unless one of the other initialization routines is
! called).
! 
! This version of the routine uses command-line options to set up the
! various structures.  It should _only_ be used in programs that intend
! to use such command-line options for setup; otherwise, the
! Initialize() routine above should be used instead.
! 
! Optional input arguments:
!   cmdline     Specifies if the command-line should be checked for
!               values (default: false).  Should be set false if this
!               module is used by external programs.
!   intervals   Indicates if rates for any sub-intervals available from
!               the efficiencies file should be calculated as well as
!               the total rate (default: true).  These are unnecessary
!               for some calculations.
! 
SUBROUTINE InitializeCL(intervals)
  IMPLICIT NONE
  LOGICAL, INTENT(IN), OPTIONAL :: intervals
  LOGICAL :: intervals0
  
  intervals0 = .TRUE.
  IF (PRESENT(intervals)) intervals0 = intervals
  
  ! Extract parameters and options from command-line.
  CALL ReadArguments()
  
  ! Set the output verbosity level.
  CALL InitVerbosity()
  
  ! Initialize WIMP mass and couplings.
  ! Requires ReadArguments to have been called.
  CALL InitWIMPCL()
  
  ! Initialize halo.
  CALL InitHaloCL()
  
  ! Initialize detector isotopes, efficiencies, array sizing, etc.
  CALL InitDetectorCL(intervals=intervals0)
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Reads the command line arguments and separates them into options
! and parameters.  See the definition of the ArgumentStructure
! structure above for details.
! 
SUBROUTINE ReadArguments()
  IMPLICIT NONE
  LOGICAL :: L
  CHARACTER*2 :: firstchars
  INTEGER :: Narg,Nopt,Nparam,I,Iopt,Iparam,ios
  REAL*8 :: x
  
  Narg   = IARGC()
  Nopt   = 0
  Nparam = 0
  
  ! NOTE: cannot use flags of form -<flag> as negative numbers
  ! will be improperly parsed
  
  ! Count argument types
  DO I=1,Narg
    CALL GETARG(I,firstchars)
    IF (firstchars .EQ. '--') THEN
      Nopt = Nopt + 1
    ELSE
      Nparam = Nparam + 1
    END IF
  END DO
  
  Arguments%Noptions = Nopt
  IF (ALLOCATED(Arguments%options)) DEALLOCATE(Arguments%options)
  ALLOCATE(Arguments%options(1:Nopt))
  
  Arguments%Nparameters = Nparam
  IF (ALLOCATED(Arguments%parameters)) DEALLOCATE(Arguments%parameters)
  ALLOCATE(Arguments%parameters(1:Nparam))
  IF (ALLOCATED(Arguments%values)) DEALLOCATE(Arguments%values)
  ALLOCATE(Arguments%values(1:Nparam))
  
  ! Divide arguments up
  Iopt   = 0
  Iparam = 0
  DO I=1,Narg
    CALL GETARG(I,firstchars)
    IF (firstchars .EQ. '--') THEN
      Iopt = Iopt + 1
      CALL GETARG(I,Arguments%options(Iopt))
    ELSE
      Iparam = Iparam + 1
      CALL GETARG(I,Arguments%parameters(Iparam))
    END IF
  END DO
  
  ! Try to parse parameters to floating point or logical -> floating
  ! point (T=1,F=0).  Set to -HUGE(1d0) if cannot be parsed.
  DO I=1,Nparam
    READ(UNIT=Arguments%parameters(I),FMT=*,IOSTAT=ios) x
    IF (ios .EQ. 0) THEN
      Arguments%values(I) = x
    ELSE
      READ(UNIT=Arguments%parameters(I),FMT=*,IOSTAT=ios) L
      IF (ios .EQ. 0) THEN
        IF (L) THEN
          Arguments%values(I) = 1d0
        ELSE
          Arguments%values(I) = 0d0
        END IF
      ELSE
        Arguments%values(I) = -HUGE(1d0)
      END IF
    END IF
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Initializes the verbosity level using command line arguments or
! default values.
! 
! Possible options:
!   --verbosity=<value>  ! Level of verbosity (default is 1)
!   --verbose            ! Equivalent to --verbosity=2
!   --quiet              ! Equivalent to --verbosity=0
! 
SUBROUTINE InitVerbosity()
  IMPLICIT NONE
  INTEGER :: verb
  
  IF (GetLongArgInteger('verbosity',verb)) THEN
    VerbosityLevel = verb
  ELSE IF (GetLongArg('verbose')) THEN
    VerbosityLevel = 2
  ELSE IF (GetLongArg('quiet')) THEN
    VerbosityLevel = 0
  ELSE
    VerbosityLevel = 1
  END IF
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Reads the p-value or CL from command line arguments or sets them to
! a default value [p=0.1 or 90% CL].  We use p+CL = 1 and return only
! the p-value (as the logarithm of its value).
! 
! Possible options:
!   --p-value=<val>         ! The p-value
!   --log-p-value=<val>     ! The logarithm of the p-value
!   --p-value-sigma=<val>   ! The p-value corresponding to the given number of s.d.'s
!                           ! in the normal distribution
!  --confidence-level=<val> ! The confidence level (1-p)
!  --confidence-level-sigma=<val>
!                           ! Equivalent to --p-value-sigma
! 
! Output argument:
!   lnp             The logarithm of the p-value (CL = 1-p)
! 
SUBROUTINE ReadLogPValue(lnp)
  IMPLICIT NONE
  REAL*8, INTENT(OUT) :: lnp
  LOGICAL :: calc_nsigma
  REAL*8 :: x,y,p,nsigma
  
  calc_nsigma = .FALSE.
  
  ! Default
  lnp = LOG(0.1d0)
  
  ! Process arguments
  IF (GetLongArgReal('confidence-level',x)) THEN
    lnp = LOG(MAX(1-x,TINY(1d0)))
  ELSE IF (GetLongArgReal('p-value',x)) THEN
    lnp = LOG(MAX(x,TINY(1d0)))
  ELSE IF (GetLongArgReal('log-p-value',x)) THEN
    lnp = x
  ELSE IF (GetLongArgReal('p-value-sigma',x)) THEN
    calc_nsigma = .TRUE.
    nsigma = x
  ELSE IF (GetLongArgReal('confidence-level-sigma',x)) THEN
    calc_nsigma = .TRUE.
    nsigma = x
  END IF
  
  ! calculate p-value in terms of normal distribution at nsigma s.d.'s.
  ! In that case:
  !    p = erfc(nsigma/\sqrt{2})
  IF (calc_nsigma) THEN
    IF (x .LT. 25d0) THEN
      lnp = LOG(ERFC(x/SQRT2))
    ELSE
      ! For large nsigma, use asymptotic expansion of erfc.
      ! Unnecessarily calculated to near full double precision....
      y = 1d0 / x**2
      lnp = -0.5d0*x**2 - LOG(SQRTPI*x/SQRT2) + LOGp1(-y*(1-3*y*(1-5*y*(1-7*y*(1-9*y*(1-11*y))))))
    END IF
  END IF
  
END SUBROUTINE



!=======================================================================
! OUTPUT ROUTINES
!=======================================================================

!-----------------------------------------------------------------------
! Prints the given number of empty comment lines (just the comment
! prefix).  Utility function.
! 
SUBROUTINE WriteEmptyCommentLines(N)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  INTEGER :: I
  DO I = 1,N
    WRITE(*,'(A)') COMMENT_PREFIX
  END DO
END SUBROUTINE


!-----------------------------------------------------------------------
! Write command header.
! Outputs the command used and date.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteCommandHeader(extra_lines)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  CHARACTER*1024 :: cmd
  CHARACTER*10 :: date,time
  CHARACTER*100 :: datetime
  
  ! Command
  CALL GetFullCommand(cmd)
  WRITE(*,'(A)') COMMENT_PREFIX &
      //  'Command: ' // TRIM(cmd)
  
  ! Date & time
  CALL DATE_AND_TIME(date,time)
  datetime = 'Generated on '                                            &
      // date(1:4) // '-' // date(5:6) // '-' // date(7:8)  // ' at '   &
      // time(1:2) // ':' // time(3:4) // ':' // time(5:10)
  WRITE(*,'(A)') COMMENT_PREFIX &
      // TRIM(datetime)
  
  ! Version
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // 'DDCalc0 version: ' // TRIM(VERSION_STRING)
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // 'See/cite C. Savage et al., arxiv:15MM.NNNNN.'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'DDCalc0 version ' // TRIM(VERSION_STRING) // '.  See/cite C. Savage et al., arxiv:15MM.NNNNN.'
  
  WRITE(*,'(A)') COMMENT_PREFIX
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output information regarding the WIMP,
! notably the mass and couplings.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteWIMPHeader(extra_lines)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  ! Mass
  WRITE(*,'(A,F11.4)') COMMENT_PREFIX &
      // 'WIMP mass [GeV]                       =',WIMP%m
  WRITE(*,'(A)') COMMENT_PREFIX
  
  ! Couplings
  WRITE(*,'(A,2(2X,A12))') COMMENT_PREFIX &
      // 'WIMP-nucleon couplings:  ','  G [GeV^-2]','  sigma [pb]'
  WRITE(*,'(A,2(2X,1PG12.4))') COMMENT_PREFIX &
      // '  proton SI              ',WIMP%GpSI,GpToSigmapSI(WIMP%m,WIMP%GpSI)
  WRITE(*,'(A,2(2X,1PG12.4))') COMMENT_PREFIX &
      // '  neutron SI             ',WIMP%GnSI,GnToSigmanSI(WIMP%m,WIMP%GnSI)
  WRITE(*,'(A,2(2X,1PG12.4))') COMMENT_PREFIX &
      // '  proton SD              ',WIMP%GpSD,GpToSigmapSD(WIMP%m,WIMP%GpSD)
  WRITE(*,'(A,2(2X,1PG12.4))') COMMENT_PREFIX &
      // '  neutron SD             ',WIMP%GnSD,GnToSigmanSD(WIMP%m,WIMP%GnSD)
  WRITE(*,'(A)') COMMENT_PREFIX
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output information regarding the dark matter halo,
! notably the density and velocity distribution.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteHaloHeader(extra_lines)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  ! Solar motion
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'Sun''s velocity in the Galactic rest frame [km/s] in Galactic'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'coordinates (U,V,W), where U is anti-radial (towards the'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'Galactic center), V is in the direction of the disk rotation,'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'and W is vertical (out of the disk plane):'
  WRITE(*,'(A,3(1X,F8.2),3X,A)') COMMENT_PREFIX &
      // '  ',Halo%vsun
  WRITE(*,'(A)') COMMENT_PREFIX
  
  ! Density
  WRITE(*,'(A,F11.4)') COMMENT_PREFIX &
      // 'Local dark matter density [GeV/cm^3]  =',Halo%rho
  WRITE(*,'(A)') COMMENT_PREFIX
  
  ! Halo velocity distribution
  IF (.NOT. Halo%tabulated) THEN
    ! Maxwell-Boltzmann parameters
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'SHM-like velocity distribution (Maxwell-Boltzmann with finite cutoff):'
    WRITE(*,'(A,3(1X,F8.2),3X,A)') COMMENT_PREFIX &
        // '  Bulk motion (galactic frame) [km/s] =',Halo%vbulk
    WRITE(*,'(A,F9.2)') COMMENT_PREFIX &
        // '  Most probable speed (v0) [km/s]     =',Halo%v0
    IF (Halo%vesc .GE. 1e6) THEN
      WRITE(*,'(A)') COMMENT_PREFIX &
        // '  Escape speed (vesc) [km/s]          = (infinite)'
    ELSE
      WRITE(*,'(A,F9.2)') COMMENT_PREFIX &
        // '  Escape speed (vesc) [km/s]          =',Halo%vesc
    END IF
  ELSE IF (Halo%eta_file .NE. '') THEN
    ! Tabulated from file
    !WRITE(*,'(A)') COMMENT_PREFIX &
    !    // 'The mean inverse speed was provided in tabulated form in the following'
    !WRITE(*,'(A)') COMMENT_PREFIX &
    !    // 'file:'
    !WRITE(*,'(A)') COMMENT_PREFIX &
    !    // '    ' // TRIM(Halo%eta_file)
    WRITE(*,'(A,A)') COMMENT_PREFIX &
        // 'Mean inverse speed tabulation file    = ',TRIM(Halo%eta_file)
  ELSE
    ! User provided tabulation
    ! This really shouldn't occur here because there is no way
    ! to provide the tabulation this way in program modes!
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The mean inverse speed was provided in tabulated form by the user.'
  END IF
  WRITE(*,'(A)') COMMENT_PREFIX
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output information regarding the detector.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteDetectorHeader(extra_lines)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  ! Exposure, observed events, expected background events
  WRITE(*,'(A,1PG12.4)') COMMENT_PREFIX &
      // 'Detector exposure [kg day]            =',DefaultDetector%exposure
  IF (DefaultDetector%Nevents .GE. 0) THEN
    WRITE(*,'(A,I6)') COMMENT_PREFIX &
      // 'Observed events                       =',DefaultDetector%Nevents
  END IF
  WRITE(*,'(A,F11.4)') COMMENT_PREFIX &
      // 'Average expected background events    =',DefaultDetector%MuBackground
  WRITE(*,'(A)') COMMENT_PREFIX
  
  ! Detector isotopes
  WRITE(*,'(A,I6)') COMMENT_PREFIX &
      // 'Isotopes                              =',DefaultDetector%Niso
  WRITE(*,'(A,99(3X,I4,2X))') COMMENT_PREFIX &
      // '  Atomic number Z       ',DefaultDetector%Ziso
  WRITE(*,'(A,99(3X,I4,2X))') COMMENT_PREFIX &
      // '  Atomic mass number A  ',DefaultDetector%Aiso
  WRITE(*,'(A,99(1X,F8.5))') COMMENT_PREFIX &
      // '  Mass fraction         ',DefaultDetector%fiso
  WRITE(*,'(A)') COMMENT_PREFIX
  
  ! Efficiencies and intervals/bins
  WRITE(*,'(A,A)') COMMENT_PREFIX &
      // 'Efficiency file                       = ',TRIM(DefaultDetector%eff_file)
  !IF (DefaultDetector%Neff .GT. 0) THEN
  IF (DefaultDetector%intervals .AND. (DefaultDetector%Neff .GT. 0)) THEN
    WRITE(*,'(A,I6)') COMMENT_PREFIX &
      // 'Number of bins/sub-intervals          =',DefaultDetector%Neff
  END IF
  WRITE(*,'(A)') COMMENT_PREFIX
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Write the interactive-mode instruction header.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteInteractiveHeader(extra_lines)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'Enter WIMP parameters in one of the following formats:'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // '  $>  m'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // '  $>  m sigmaSI'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // '  $>  m sigmaSI sigmaSD'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // '  $>  m sigmaSI sigmapSD sigmanSD'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // '  $>  m sigmapSI sigmanSI sigmapSD sigmanSD'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'where m is the WIMP mass [GeV] and sigma is a WIMP-nucleon cross-section'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // '[pb].  In the first case, all cross-sections are set to 1 pb.  In the'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'second case, SD couplings are set to zero.  Negative cross-sections may'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'be given to indicate the corresponding coupling is negative.  A blank'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'line will terminate the program (as will an invalid format).'
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // '  $>  m GpSI GnSI GpSD GnSD'
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // 'where m is the WIMP mass [GeV], sigma is a WIMP-nucleon cross-section'
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // '[pb], and G is a WIMP-nucleon coupling [GeV^-2].  In the first case,'
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // 'all cross-sections are set to 1 pb.  In the second case, SD couplings'
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // 'are set to zero.  Negative cross-sections may be given to indicate'
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // 'the corresponding coupling is negative.  A blank line will terminate'
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // 'the program (as will an invalid format).'
  !WRITE(*,'(A)') COMMENT_PREFIX &
  !    // ''
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Write log-likelihood header.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteLogLikelihoodHeader(extra_lines)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'The log-likelihood for the given parameters is given below,'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'defined as L = P(N|s+b) where P is the Poisson distribution of'
  WRITE(*,'(A)') COMMENT_PREFIX &
      // 'N observed events given an expected signal s and background b.'
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Write log p-value header.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteLogPValueHeader(extra_lines)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  IF (DefaultDetector%intervals .AND. (DefaultDetector%Neff .EQ. DefaultDetector%Nevents+1)) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The log of the p-value for the given parameters is given below,'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'calculated using the maximum gap method.  See:'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  S. Yellin, PRD 66, 032005 (2002) [physics/0203002]'
  ELSE
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The log of the p-value for the given parameters are given below,'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'calculated using a Poisson distribution on the number of events.'
  END IF
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a header containing information regarding
! the events and likelihoods data to follow.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteEventsAndLikelihoodsHeader(extra_lines)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The columns below contain the following data:'
  END IF
  
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  m            WIMP mass [GeV].'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  GpSI,GnSI,GpSD,GnSD'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               WIMP-nucleon spin-independent and spin-dependent couplings [GeV^2].'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmapSI,sigmanSI,sigmapSD,sigmanSD'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               WIMP-nucleon SI and SD scattering cross-sections [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  observed     The observed number of events.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  background   Average expected background events.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  signal(SI)   Average expected spin-independent signal events.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  signal(SD)   Average expected spin-dependent signal events.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  log(L)       Log-likelihood using the Poisson distribution (signal+background).'
    IF (DefaultDetector%intervals .AND. (DefaultDetector%Neff .EQ. DefaultDetector%Nevents+1)) THEN
      WRITE(*,'(A)') COMMENT_PREFIX &
        // '  log(p)       Log of the p-value determined using the maximum gap method;'
      WRITE(*,'(A)') COMMENT_PREFIX &
        // '               see Yellin, Phys. Rev. D 66, 032005 (2002) [physics/0203002].'
    ELSE
      WRITE(*,'(A)') COMMENT_PREFIX &
        // '  log(p)       Log of the p-value determined using the Poisson distribution'
      WRITE(*,'(A)') COMMENT_PREFIX &
        // '               (signal only: no background subtraction).'
    END IF
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a column header for the events and
! likelihoods data to follow.
! 
SUBROUTINE WriteEventsAndLikelihoodsColumnHeader()
  IMPLICIT NONE
  
  SELECT CASE (VerbosityLevel)
  CASE (:0)
    CONTINUE
  CASE (1:2)
    WRITE(*,'(A,1(1X,A8),3(1X,A11),2(1X,A11))') COMMENT_PREFIX,         &
        'observed','background ','signal(SI) ','signal(SD) ',           &
        '  log(L)   ','  log(p)   '
  CASE (3:)
    WRITE(*,'(A,1(1X,A10),4(1X,A10),4(1X,A10),1(1X,A8),3(1X,A11),2(1X,A11))') &
        COMMENT_PREFIX,'    mass  ',                                    &
        '   GpSI   ','   GnSI   ','   GpSD   ','   GnSD   ',            &
        ' sigmapSI ',' sigmanSI ',' sigmapSD ',' sigmanSD ',            &
        'observed','background ','signal(SI) ','signal(SD) ',           &
        '  log(L)   ','  log(p)   '
  END SELECT
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output the events and likelihoods data for the
! current WIMP mass and couplings.
! 
SUBROUTINE WriteEventsAndLikelihoodsData()
  IMPLICIT NONE
  REAL*8 :: lnLike
  REAL*8 :: lnp
  
  ! Get log-likelihood and p-value
  lnLike = LogLikelihood()
  lnp    = LogPValue()
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  SELECT CASE (ABS(VerbosityLevel))
  CASE (:2)
    WRITE(*,'(A,1(2X,I5,2X),3(1X,1PG11.4),2(1X,1PG11.4))')              &
        DATA_PREFIX,                                                    &
        DefaultDetector%Nevents,DefaultDetector%MuBackground,           &
        DefaultDetector%MuSignalSI(0),DefaultDetector%MuSignalSD(0),    &
        lnLike,lnp
  CASE (3:)
    WRITE(*,'(A,1(1X,F10.3),4(1X,1PG10.3),4(1X,1PG10.3),1(2X,I5,2X),3(1X,1PG11.4),2(1X,1PG11.4))') &
        DATA_PREFIX,                                                    &
        WIMP%m,WIMP%GpSI,WIMP%GnSI,WIMP%GpSD,WIMP%GnSD,                 &
        GpToSigmapSI(WIMP%m,WIMP%GpSI),GnToSigmanSI(WIMP%m,WIMP%GnSI),  &
        GpToSigmapSD(WIMP%m,WIMP%GpSD),GnToSigmanSD(WIMP%m,WIMP%GnSD),  &
        DefaultDetector%Nevents,DefaultDetector%MuBackground,           &
        DefaultDetector%MuSignalSI(0),DefaultDetector%MuSignalSD(0),    &
        lnLike,lnp
  END SELECT
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a header containing information regarding
! the events and likelihoods data to follow.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteSpectrumHeader(extra_lines)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  IF (VerbosityLevel .GE. 4) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The table below gives differential rate components at reference'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'couplings that yield WIMP-nucleon scattering cross-sections of 1 pb.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'Separate columns are given for the spectrum contribution arising from'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'the proton-proton [dRdEpp0], neutron-neutron [dRdEnn0], and proton-'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'neutron [dRdEpn0] components of the cross-section/coupling formula.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'This allows rates for arbitrary couplings to be constructed as follows:'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '    dRdE = (sigmap/[pb])*dRdEpp0 + (sigman/[pb])*dRdEnn0'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '              +/- (sqrt{sigmap*sigman}/[pb])*dRdEpn0'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'where sigmap and sigman are the WIMP-nucleon scattering cross-sections'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'in pb.  The sign of the cross-term should be the same as the sign of'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'the product of couplings Gp*Gn.'
    !WRITE(*,'(A)') COMMENT_PREFIX &
    !    // ''
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The columns below contain the following data:'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  E            Recoil energy [keV].'
  END IF
  
  SELECT CASE (VerbosityLevel)
  CASE (2)
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  dRdE         Differential recoil spectrum [cpd/kg/keV].'
  CASE (3)
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  dRdE         Differential recoil spectrum [cpd/kg/keV].'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               One column for each of the total, SI, and SD spectra.'
  CASE (4:)
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  dRdEpp0,dRdEpn0,dRdEnn0'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               Differential recoil spectrum components [cpd/kg/keV],'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               given separately for SI and SD interactions.'
  END SELECT
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a column header for the tabulated
! differential rate spectrum to follow.
! 
SUBROUTINE WriteSpectrumColumnHeader()
  IMPLICIT NONE
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  SELECT CASE (ABS(VerbosityLevel))
  CASE (1)
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'  E [keV] '
    WRITE(*,'(1(1X,A12))',ADVANCE='NO')                                 &
        ' dR/dE [dru]'
  CASE (2)
    ! Differential rate
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'    E     '
    WRITE(*,'(1(1X,A12))',ADVANCE='NO')                                 &
        '   dR/dE    '
  CASE (3)
    ! Combined, si, sd
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'    E     '
    WRITE(*,'(3(1X,A12))',ADVANCE='NO')                                 &
        '   dR/dE    ',' dR/dE(SI)  ',' dR/dE(SD)  '
  CASE (4:)
    ! Reference rate components
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'    E     '
    WRITE(*,'(6(1X,A12))',ADVANCE='NO')                                 &
        'dRdEpp0(SI) ','dRdEpn0(SI) ','dRdEnn0(SI) ',                   &
        'dRdEpp0(SD) ','dRdEpn0(SD) ','dRdEnn0(SD) '
  END SELECT
  IF (VerbosityLevel .GE. 1) WRITE(*,'(A)') ''
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output the tabulated differential rate spectrum
! for the current WIMP mass and couplings.
! 
SUBROUTINE WriteSpectrumData()
  IMPLICIT NONE
  INTEGER :: KE
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  DO KE = 1,DefaultDetector%NE
    WRITE(*,'(A,1(1X,F10.4))',ADVANCE='NO') DATA_PREFIX,                &
        DefaultDetector%E(KE)
    SELECT CASE (ABS(VerbosityLevel))
    CASE (:2)
      ! Differential rate
      WRITE(*,'(1(1X,1PG12.5))',ADVANCE='NO')                           &
          CoerceExponent(DefaultDetector%dRdE(KE),2,5)
    CASE (3)
      ! Combined, si, sd
      WRITE(*,'(3(1X,1PG12.5))',ADVANCE='NO')                           &
          CoerceExponent(DefaultDetector%dRdE(KE),2,5),                 &
          CoerceExponent(DefaultDetector%dRdEsi(KE),2,5),               &
          CoerceExponent(DefaultDetector%dRdEsd(KE),2,5)
    CASE (4:)
      ! Reference rate components
      WRITE(*,'(6(1X,1PG12.5))',ADVANCE='NO')                           &
          CoerceExponent(DefaultDetector%dRdEsi0(+1,KE),2,5),           &
          CoerceExponent(DefaultDetector%dRdEsi0( 0,KE),2,5),           &
          CoerceExponent(DefaultDetector%dRdEsi0(-1,KE),2,5),           &
          CoerceExponent(DefaultDetector%dRdEsd0(+1,KE),2,5),           &
          CoerceExponent(DefaultDetector%dRdEsd0( 0,KE),2,5),           &
          CoerceExponent(DefaultDetector%dRdEsd0(-1,KE),2,5)
    END SELECT
    WRITE(*,'(A)') ''
  END DO
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a header containing information regarding
! the events table (tabulated by mass) to follow.
! 
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteEventsByMassHeader(extra_lines)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  REAL*8 :: sigmapSI,sigmanSI,sigmapSD,sigmanSD
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  ! Will write out the fixed WIMP-nucleon cross-sections.
  ! Negative for negative couplings.
  CALL GetWIMP(sigmapSI=sigmapSI,sigmanSI=sigmanSI,sigmapSD=sigmapSD,sigmanSD=sigmanSD)
  IF (WIMP%GpSI .LT. 0d0) sigmapSI = -ABS(sigmapSI)
  IF (WIMP%GnSI .LT. 0d0) sigmanSI = -ABS(sigmanSI)
  IF (WIMP%GpSD .LT. 0d0) sigmapSD = -ABS(sigmapSD)
  IF (WIMP%GnSD .LT. 0d0) sigmanSD = -ABS(sigmanSD)
  
  ! Description and fixed cross-sections
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The table below gives the expected spin-independent (SI) and spin-'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'dependent (SD) interaction events, tabulated by WIMP mass, for fixed'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'WIMP-nucleon cross-sections.  The fixed cross-sections are [pb]:'
    WRITE(*,'(A,1(2X,1PG12.4))') COMMENT_PREFIX &
        // '  sigmapSI (proton SI)  =',CoerceExponent(sigmapSI,2,4)
    WRITE(*,'(A,1(2X,1PG12.4))') COMMENT_PREFIX &
        // '  sigmanSI (neutron SI) =',CoerceExponent(sigmanSI,2,4)
    WRITE(*,'(A,1(2X,1PG12.4))') COMMENT_PREFIX &
        // '  sigmapSD (proton SD)  =',CoerceExponent(sigmapSD,2,4)
    WRITE(*,'(A,1(2X,1PG12.4))') COMMENT_PREFIX &
        // '  sigmanSD (neutron SD) =',CoerceExponent(sigmanSD,2,4)
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'A negative cross-section means the corresponding WIMP-nucleon coupling'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'will be taken to be negative.'
    !WRITE(*,'(A)') COMMENT_PREFIX &
    !    // ''
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The columns below contain the following data:'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  mass         WIMP mass [GeV].'
  END IF
  
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  G            WIMP-nucleon couplings for spin-independent (SI) and'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               spin-dependent (SD) interactions [GeV^-2].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  events(SI)   Average expected spin-independent events.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  events(SD)   Average expected spin-dependent events.'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a column header for the table of expected
! events tabulated by mass.
! 
SUBROUTINE WriteEventsByMassColumnHeader()
  IMPLICIT NONE
  INTEGER :: Keff
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  ! Multiple lines in some cases
  IF (VerbosityLevel .GE. 3) THEN
    ! Mass column
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,''
    ! Coupling columns
    IF (VerbosityLevel .GE. 4) THEN
      WRITE(*,'(4(1X,A10))',ADVANCE='NO') '','','',''
    END IF
    ! Events for full range
    WRITE(*,'(1X,A1)',ADVANCE='NO') ' '
    WRITE(*,'(1(1X,A23))',ADVANCE='NO')                                 &
        '----- full range ------'
    ! Events for sub-intervals
    IF (DefaultDetector%intervals) THEN
      DO Keff = 1,DefaultDetector%Neff
        WRITE(*,'(1X,A1)',ADVANCE='NO') '|'
        WRITE(*,'(1(1X,A14,I3,A6))',ADVANCE='NO')                       &
            '----- interval',Keff,' -----'
      END DO
    END IF
    WRITE(*,'(A)') ''
  END IF
  
  ! Main column header line below
  
  ! Mass column
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'mass [GeV]'
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'    mass  '
  END IF
  
  ! Coupling columns
  IF (VerbosityLevel .GE. 4) THEN
    WRITE(*,'(4(1X,A10))',ADVANCE='NO')                                 &
        '   GpSI   ','   GnSI   ','   GpSD   ','   GnSD   '
  END IF
  
  ! Events for full range
  IF (VerbosityLevel .GE. 3) WRITE(*,'(1X,A1)',ADVANCE='NO') ' '
  IF (VerbosityLevel .GE. 1) THEN
    WRITE(*,'(2(1X,A11))',ADVANCE='NO') ' events(SI)',' events(SD)'
  END IF
  
  ! Events for sub-intervals
  IF (VerbosityLevel .GE. 3) THEN
    IF (DefaultDetector%intervals) THEN
      DO Keff = 1,DefaultDetector%Neff
        WRITE(*,'(1X,A1)',ADVANCE='NO') ' '
        WRITE(*,'(2(1X,A11))',ADVANCE='NO') ' events(SI)',' events(SD)'
      END DO
    END IF
  END IF
  
  IF (VerbosityLevel .GE. 1) WRITE(*,'(A)') ''
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output the mass and expected events for the
! current WIMP.  Used for tabulation of events by WIMP mass.
! 
SUBROUTINE WriteEventsByMassData()
  IMPLICIT NONE
  INTEGER :: Keff
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  ! Mass
  WRITE(*,'(A,1(1X,F10.4))',ADVANCE='NO') DATA_PREFIX,                  &
      CoerceNumber(WIMP%m,10,4)
  
  ! Couplings
  IF (ABS(VerbosityLevel) .GE. 4) THEN
    WRITE(*,'(4(1X,1PG10.3))',ADVANCE='NO')                             &
        CoerceExponent(WIMP%GpSI,2,3),CoerceExponent(WIMP%GnSI,2,3),    &
        CoerceExponent(WIMP%GpSD,2,3),CoerceExponent(WIMP%GnSD,2,3)
  END IF
  
  ! Events (full range)
  IF (ABS(VerbosityLevel) .GE. 3) WRITE(*,'(1X,A1)',ADVANCE='NO') ' '
  WRITE(*,'(2(1X,1PG11.4))',ADVANCE='NO')                               &
      CoerceExponent(DefaultDetector%MuSignalSI(0),2,4),                &
      CoerceExponent(DefaultDetector%MuSignalSD(0),2,4)
  
  ! Events for sub-intervals
  IF (ABS(VerbosityLevel) .GE. 3) THEN
    IF (DefaultDetector%intervals) THEN
      DO Keff = 1,DefaultDetector%Neff
        WRITE(*,'(1X,A1)',ADVANCE='NO') ' '
        WRITE(*,'(2(1X,1PG11.4))',ADVANCE='NO')                         &
            CoerceExponent(DefaultDetector%MuSignalSI(Keff),2,4),       &
            CoerceExponent(DefaultDetector%MuSignalSD(Keff),2,4)
      END DO
    END IF
  END IF
  
  WRITE(*,'(A)') ''
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a header containing information regarding
! the SI cross-section constraints table (tabulated by mass) to follow.
! 
! Required input arguments:
!   lnp             The logarithm of the p-value to use for the
!                   constraint CL.
!   thetaG          The angle of (Gp,Gn) that is fixed for each
!                   constraint determination.
!   s1,s2           The range of allowed signal expectation values.
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteConstraintsSIHeader(lnp,thetaG,s1,s2,extra_lines)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: lnp,thetaG,s1,s2
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  ! Description
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The table below gives the range of allowed spin-independent (SI) cross-'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'sections.  The allowed cross-sections are those that predict a mean'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'signal compatible with the observed number of events and estimated'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'background.  The confidence interval (CI) for the mean signal at the'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'given confidence level CL is determined using a Poisson distribution'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'with Feldman-Cousins ordering; see Feldman & Cousins, Phys. Rev. D 57,'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '3873 (1998) [physics/9711021].'
    WRITE(*,'(A,1(2X,1PG12.4))') COMMENT_PREFIX &
        // '  1-CL                  =',CoerceExponent(EXP(lnp),2,4)
    WRITE(*,'(A,1(2X,1PG12.4),2X,A,1(2X,1PG12.4))') COMMENT_PREFIX &
        // '  signal events CI      =',CoerceExponent(s1,2,4),'-',CoerceExponent(s2,2,4)
    WRITE(*,'(A)') COMMENT_PREFIX
    
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The constraints are determined using a fixed ratio of the two WIMP-nucleon'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'couplings Gp & Gn, described by the angle theta s.t. tan(theta) = Gn/Gp.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The conventional isospin-invariant case Gn=Gp corresponds to theta=PI/4.'
    WRITE(*,'(A,1(2X,F9.5))') COMMENT_PREFIX &
        // '  theta/PI              =',CoerceNumber(thetaG/PI,10,5)
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The columns below contain the following data:'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  mass         WIMP mass [GeV].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmapSI     The WIMP-proton spin-independent cross-section lower and'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               upper limits [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmanSI     The WIMP-neutron spin-independent cross-section lower and'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               upper limits [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a column header for the table of SI cross-
! section constraints tabulated by mass.
! 
SUBROUTINE WriteConstraintsSIColumnHeader()
  IMPLICIT NONE
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  ! Mass column
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'mass [GeV]'
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'    mass  '
  END IF
  
  ! WIMP-proton cross-section columns
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(1(1X,A25))',ADVANCE='NO') ' --sigmapSI range [pb]-- '
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(1(1X,A25))',ADVANCE='NO') ' --- sigmapSI range ---- '
  END IF
  
  ! WIMP-neutron cross-section column
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(1(1X,A25))',ADVANCE='NO') ' --- sigmanSI range ---- '
  END IF
  
  IF (VerbosityLevel .GE. 1) WRITE(*,'(A)') ''
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output the mass and SI cross-section limit for the
! current WIMP.  Used for tabulation of limits by WIMP mass.
! 
! Input argument:
!   s1,s2       The range of allowed signal expectation values.
! 
SUBROUTINE WriteConstraintsSIData(s1,s2)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s1,s2
  REAL*8 :: mu,x1,x2,sigmapSI,sigmanSI
  
  ! Need scale factors x s.t. sigma -> x*sigma gives desired
  ! number of events.
  CALL GetRates(signal_si=mu)
  ! Empty set case
  IF (s2 .EQ. 0d0) THEN
    x1 = 0d0
    x2 = 0d0
  ! General case
  ELSE IF (mu .GT. 0d0) THEN
    x1 = s1/mu
    x2 = s2/mu
  ! No events case (at any scale)
  ELSE
    x1 = 0d0
    x2 = HUGE(1d0)
  END IF
  
  ! Cross-sections (multiply by x for limit)
  CALL GetWIMP(sigmapSI=sigmapSI,sigmanSI=sigmanSI)
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  ! Mass
  WRITE(*,'(A,1(1X,F10.4))',ADVANCE='NO') DATA_PREFIX,                  &
      CoerceNumber(WIMP%m,10,4)
  
  ! WIMP-proton cross-section
  WRITE(*,'(2(1X,1PG12.5))',ADVANCE='NO')                               &
      CoerceExponent(x1*sigmapSI,2,5),CoerceExponent(x2*sigmapSI,2,5)
  
  ! WIMP-neutron cross-section
  IF (ABS(VerbosityLevel) .GE. 3) THEN
    WRITE(*,'(2(1X,1PG12.5))',ADVANCE='NO')                             &
        CoerceExponent(x1*sigmanSI,2,5),CoerceExponent(x2*sigmanSI,2,5)
  END IF
  
  WRITE(*,'(A)') ''
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a header containing information regarding
! the SD cross-section constraints table (tabulated by mass) to follow.
! 
! Required input arguments:
!   lnp             The logarithm of the p-value to use for the
!                   constraint CL.
!   thetaG          The angle of (Gp,Gn) that is fixed for each
!                   constraint determination.
!   s1,s2           The range of allowed signal expectation values.
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteConstraintsSDHeader(lnp,thetaG,s1,s2,extra_lines)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: lnp,thetaG,s1,s2
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  ! Description
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The table below gives the range of allowed spin-dependent (SD) cross-'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'sections.  The allowed cross-sections are those that predict a mean'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'signal compatible with the observed number of events and estimated'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'background.  The confidence interval (CI) for the mean signal at the'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'given confidence level CL is determined using a Poisson distribution'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'with Feldman-Cousins ordering; see Feldman & Cousins, Phys. Rev. D 57,'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '3873 (1998) [physics/9711021].'
    WRITE(*,'(A,1(2X,1PG12.4))') COMMENT_PREFIX &
        // '  1-CL                  =',CoerceExponent(EXP(lnp),2,4)
    WRITE(*,'(A,1(2X,1PG12.4),2X,A,1(2X,1PG12.4))') COMMENT_PREFIX &
        // '  signal events CI      =',CoerceExponent(s1,2,4),'-',CoerceExponent(s2,2,4)
    WRITE(*,'(A)') COMMENT_PREFIX
    
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The constraints are determined using a fixed ratio of the two WIMP-nucleon'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'couplings Gp & Gn, described by the angle theta s.t. tan(theta) = Gn/Gp.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The isospin-invariant case Gn=Gp corresponds to theta=PI/4, proton-only'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'coupling is theta=0, and neutron-only coupling is theta=PI/2.'
    WRITE(*,'(A,1(2X,F9.5))') COMMENT_PREFIX &
        // '  theta/PI              =',CoerceNumber(thetaG/PI,10,5)
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The columns below contain the following data:'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  mass         WIMP mass [GeV].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmapSD     The WIMP-proton spin-dependent cross-section lower and'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               upper limits [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmanSD     The WIMP-neutron spin-dependent cross-section lower and'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '               upper limits [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a column header for the table of SD cross-
! section constraints tabulated by mass.
! 
SUBROUTINE WriteConstraintsSDColumnHeader()
  IMPLICIT NONE
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  ! Mass column
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'mass [GeV]'
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'    mass  '
  END IF
  
  ! WIMP-proton cross-section columns
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(1(1X,A25))',ADVANCE='NO') ' --sigmapSD range [pb]-- '
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(1(1X,A25))',ADVANCE='NO') ' --- sigmapSD range ---- '
  END IF
  
  ! WIMP-neutron cross-section column
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(1(1X,A25))',ADVANCE='NO') ' --- sigmanSD range ---- '
  END IF
  
  IF (VerbosityLevel .GE. 1) WRITE(*,'(A)') ''
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output the mass and SD cross-section limit for the
! current WIMP.  Used for tabulation of limits by WIMP mass.
! 
! Input argument:
!   s1,s2       The range of allowed signal expectation values.
! 
SUBROUTINE WriteConstraintsSDData(s1,s2)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s1,s2
  REAL*8 :: mu,x1,x2,sigmapSD,sigmanSD
  
  ! Need scale factors x s.t. sigma -> x*sigma gives desired
  ! number of events.
  CALL GetRates(signal_sd=mu)
  ! Empty set case
  IF (s2 .EQ. 0d0) THEN
    x1 = 0d0
    x2 = 0d0
  ! General case
  ELSE IF (mu .GT. 0d0) THEN
    x1 = s1/mu
    x2 = s2/mu
  ! No events case (at any scale)
  ELSE
    x1 = 0d0
    x2 = HUGE(1d0)
  END IF
  
  ! Cross-sections (multiply by x for limit)
  CALL GetWIMP(sigmapSD=sigmapSD,sigmanSD=sigmanSD)
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  ! Mass
  WRITE(*,'(A,1(1X,F10.4))',ADVANCE='NO') DATA_PREFIX,                  &
      CoerceNumber(WIMP%m,10,4)
  
  ! WIMP-proton cross-section
  WRITE(*,'(2(1X,1PG12.5))',ADVANCE='NO')                               &
      CoerceExponent(x1*sigmapSD,2,5),CoerceExponent(x2*sigmapSD,2,5)
  
  ! WIMP-neutron cross-section
  IF (ABS(VerbosityLevel) .GE. 3) THEN
    WRITE(*,'(2(1X,1PG12.5))',ADVANCE='NO')                             &
        CoerceExponent(x1*sigmanSD,2,5),CoerceExponent(x2*sigmanSD,2,5)
  END IF
  
  WRITE(*,'(A)') ''
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a header containing information regarding
! the SI cross-section limits table (tabulated by mass) to follow.
! 
! Required input arguments:
!   lnp             The logarithm of the p-value to use for the limit
!                   CL.
!   thetaG          The angle of (Gp,Gn) that is fixed for each limit
!                   determination.
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteLimitsSIHeader(lnp,thetaG,extra_lines)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: lnp,thetaG
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  ! Description
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The table below gives the upper limit on spin-independent (SI) cross-'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'section(s) that are not excluded.  Cross-sections are excluded if their'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'p-value is smaller than the given p-value, where the p-value is'
    IF (DefaultDetector%intervals .AND. (DefaultDetector%Neff .EQ. DefaultDetector%Nevents+1)) THEN
      WRITE(*,'(A)') COMMENT_PREFIX &
        // 'determined using the maximum gap method; see Yellin, Phys. Rev. D 66,'
      WRITE(*,'(A)') COMMENT_PREFIX &
        // '032005 (2002) [physics/0203002].'
    ELSE
      WRITE(*,'(A)') COMMENT_PREFIX &
        // 'determined using the Poisson distribution (signal only: no background).'
    END IF
    WRITE(*,'(A,1(2X,1PG12.4))') COMMENT_PREFIX &
        // '  p-value               =',CoerceExponent(EXP(lnp),2,4)
    WRITE(*,'(A)') COMMENT_PREFIX
    
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The limit is determined using a fixed ratio of the two WIMP-nucleon'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'couplings Gp & Gn, described by the angle theta s.t. tan(theta) = Gn/Gp.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The conventional isospin-invariant case Gn=Gp corresponds to theta=PI/4.'
    WRITE(*,'(A,1(2X,F9.5))') COMMENT_PREFIX &
        // '  theta/PI              =',CoerceNumber(thetaG/PI,10,5)
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The columns below contain the following data:'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  mass         WIMP mass [GeV].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmapSI     The WIMP-proton spin-independent cross-section upper limit [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmanSI     The WIMP-neutron spin-independent cross-section upper limit [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a column header for the table of SI cross-
! section limits tabulated by mass.
! 
SUBROUTINE WriteLimitsSIColumnHeader()
  IMPLICIT NONE
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  ! Mass column
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'mass [GeV]'
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'    mass  '
  END IF
  
  ! WIMP-proton cross-section column
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(1(1X,A12))',ADVANCE='NO') 'sigmapSI[pb]'
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(1(1X,A12))',ADVANCE='NO') '  sigmapSI  '
  END IF
  
  ! WIMP-neutron cross-section column
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(1(1X,A12))',ADVANCE='NO') '  sigmanSI  '
  END IF
  
  IF (VerbosityLevel .GE. 1) WRITE(*,'(A)') ''
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output the mass and SI cross-section limit for the
! current WIMP.  Used for tabulation of limits by WIMP mass.
! 
! Input argument:
!   lnp         Logarithm of the p-value for the limit CL
! 
SUBROUTINE WriteLimitsSIData(lnp)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: lnp
  REAL*8 :: x,sigmapSI,sigmanSI
  
  ! Need scale factor
  x = ScaleToPValue(lnp=lnp)
  
  ! Cross-sections (multiply by x for limit)
  CALL GetWIMP(sigmapSI=sigmapSI,sigmanSI=sigmanSI)
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  ! Mass
  WRITE(*,'(A,1(1X,F10.4))',ADVANCE='NO') DATA_PREFIX,                  &
      CoerceNumber(WIMP%m,10,4)
  
  ! WIMP-proton cross-section
  WRITE(*,'(1(1X,1PG12.5))',ADVANCE='NO')                               &
      CoerceExponent(x*sigmapSI,2,5)
  
  ! WIMP-neutron cross-section
  IF (ABS(VerbosityLevel) .GE. 3) THEN
    WRITE(*,'(1(1X,1PG12.5))',ADVANCE='NO')                             &
        CoerceExponent(x*sigmanSI,2,5)
  END IF
  
  WRITE(*,'(A)') ''
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a header containing information regarding
! the SD cross-section limits table (tabulated by mass) to follow.
! 
! Required input arguments:
!   lnp             The logarithm of the p-value to use for the limit
!                   CL.
!   thetaG          The angle of (Gp,Gn) that is fixed for each limit
!                   determination.
! Optional input arguments:
!   extra_lines     Blank lines (albeit with prefix) to add after
!                   output
! 
SUBROUTINE WriteLimitsSDHeader(lnp,thetaG,extra_lines)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: lnp,thetaG
  INTEGER, INTENT(IN), OPTIONAL :: extra_lines
  
  !WRITE(*,'(A)') COMMENT_LINE
  
  ! Description
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The table below gives the upper limit on spin-dependent (SD) cross-'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'section(s) that are not excluded.  Cross-sections are excluded if their'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'p-value is smaller than the given p-value, where the p-value is'
    IF (DefaultDetector%intervals .AND. (DefaultDetector%Neff .EQ. DefaultDetector%Nevents+1)) THEN
      WRITE(*,'(A)') COMMENT_PREFIX &
        // 'determined using the maximum gap method; see Yellin, Phys. Rev. D 66,'
      WRITE(*,'(A)') COMMENT_PREFIX &
        // '032005 (2002) [physics/0203002].'
    ELSE
      WRITE(*,'(A)') COMMENT_PREFIX &
        // 'determined using the Poisson distribution (signal only: no background'
      WRITE(*,'(A)') COMMENT_PREFIX &
        // ').'
    END IF
    WRITE(*,'(A,1(2X,1PG12.4))') COMMENT_PREFIX &
        // '  p-value               =',CoerceExponent(EXP(lnp),2,4)
    WRITE(*,'(A)') COMMENT_PREFIX
    
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The limit is determined using a fixed ratio of the two WIMP-nucleon'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'couplings Gp & Gn, described by the angle theta s.t. tan(theta) = Gn/Gp.'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The isospin-invariant case Gn=Gp corresponds to theta=PI/4, proton-only'
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'coupling is theta=0, and neutron-only coupling is theta=PI/2.'
    WRITE(*,'(A,1(2X,F9.5))') COMMENT_PREFIX &
        // '  theta/PI              =',CoerceNumber(thetaG/PI,10,5)
    WRITE(*,'(A)') COMMENT_PREFIX
    
    !WRITE(*,'(A)') COMMENT_PREFIX &
    !    // 'NOTE: SPIN-DEPENDENT FORM FACTORS NOT IMPLEMENTED.'
    !WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // 'The columns below contain the following data:'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  mass         WIMP mass [GeV].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmapSD     The WIMP-proton spin-dependent cross-section upper limit [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(A)') COMMENT_PREFIX &
        // '  sigmanSD     The WIMP-neutron spin-dependent cross-section upper limit [pb].'
  END IF
  
  IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A)') COMMENT_PREFIX
  END IF
  
  IF (PRESENT(extra_lines)) CALL WriteEmptyCommentLines(extra_lines)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output a column header for the table of SD cross-
! section limits tabulated by mass.
! 
SUBROUTINE WriteLimitsSDColumnHeader()
  IMPLICIT NONE
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  ! Mass column
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'mass [GeV]'
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(A,1(1X,A10))',ADVANCE='NO') COMMENT_PREFIX,'    mass  '
  END IF
  
  ! WIMP-proton cross-section column
  IF (VerbosityLevel .EQ. 1) THEN
    WRITE(*,'(1(1X,A12))',ADVANCE='NO') 'sigmapSD[pb]'
  ELSE IF (VerbosityLevel .GE. 2) THEN
    WRITE(*,'(1(1X,A12))',ADVANCE='NO') '  sigmapSD  '
  END IF
  
  ! WIMP-neutron cross-section column
  IF (VerbosityLevel .GE. 3) THEN
    WRITE(*,'(1(1X,A12))',ADVANCE='NO') '  sigmanSD  '
  END IF
  
  IF (VerbosityLevel .GE. 1) WRITE(*,'(A)') ''
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Prints to standard output the mass and SD cross-section limit for the
! current WIMP.  Used for tabulation of limits by WIMP mass.
! 
! Input argument:
!   lnp         Logarithm of the p-value for the limit CL
! 
SUBROUTINE WriteLimitsSDData(lnp)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: lnp
  REAL*8 :: x,sigmapSD,sigmanSD
  
  ! Need scale factor
  x = ScaleToPValue(lnp=lnp)
  
  ! Cross-sections (multiply by x for limit)
  CALL GetWIMP(sigmapSD=sigmapSD,sigmanSD=sigmanSD)
  
  ! Columns to print depend on the verbosity level.
  ! For data, only absolute value of verbosity is used.
  
  ! Mass
  WRITE(*,'(A,1(1X,F10.4))',ADVANCE='NO') DATA_PREFIX,                  &
      CoerceNumber(WIMP%m,10,4)
  
  ! WIMP-proton cross-section
  WRITE(*,'(1(1X,1PG12.5))',ADVANCE='NO')                               &
      CoerceExponent(x*sigmapSD,2,5)
  
  ! WIMP-neutron cross-section
  IF (ABS(VerbosityLevel) .GE. 3) THEN
    WRITE(*,'(1(1X,1PG12.5))',ADVANCE='NO')                             &
        CoerceExponent(x*sigmanSD,2,5)
  END IF
  
  WRITE(*,'(A)') ''
  
END SUBROUTINE



!=======================================================================
! WIMP ROUTINES
!=======================================================================

! ----------------------------------------------------------------------
! Get various WIMP quantities.
! 
! Optional output arguments:
!   m           WIMP mass [GeV].
!   GpSI        Spin-independent WIMP-proton coupling [GeV^-2].
!   GnSI        Spin-independent WIMP-neutron coupling [GeV^-2].
!   GpSD        Spin-dependent WIMP-proton coupling [GeV^-2].
!   GnSD        Spin-dependent WIMP-neutron coupling [GeV^-2].
!   GpSI0       Reference spin-independent WIMP-proton coupling [GeV^-2]
!               corresponding to \sigma_{SI,p} = 1 pb.
!   GnSI0       Reference spin-independent WIMP-neutron coupling [GeV^-2]
!               corresponding to \sigma_{SI,n} = 1 pb.
!   GpSD0       Reference spin-dependent WIMP-proton coupling [GeV^-2]
!               corresponding to \sigma_{SD,p} = 1 pb.
!   GnSD0       Reference spin-dependent WIMP-neutron coupling [GeV^-2]
!               corresponding to \sigma_{SD,n} = 1 pb.
!   fp          Spin-independent WIMP-proton coupling [GeV^-2].
!               Related by GpSI = 2 fp.
!   fn          Spin-independent WIMP-neutron coupling [GeV^-2].
!               Related by GnSI = 2 fn.
!   ap          Spin-dependent WIMP-proton coupling [unitless].
!               Related by GpSD = 2\sqrt{2} G_F ap.
!   an          Spin-dependent WIMP-neutron coupling [unitless].
!               Related by GnSD = 2\sqrt{2} G_F an.
!   sigmapSI    Spin-independent WIMP-proton cross-section [pb].
!   sigmanSI    Spin-independent WIMP-neutron cross-section [pb].
!   sigmapSD    Spin-dependent WIMP-proton cross-section [pb].
!   sigmanSD    Spin-dependent WIMP-neutron cross-section [pb].
! 
SUBROUTINE GetWIMP(m,GpSI,GnSI,GpSD,GnSD,GpSI0,GnSI0,GpSD0,GnSD0,       &
                   fp,fn,ap,an,sigmapSI,sigmanSI,sigmapSD,sigmanSD)
  IMPLICIT NONE
  REAL*8, INTENT(OUT), OPTIONAL :: m,GpSI,GnSI,GpSD,GnSD,GpSI0,GnSI0,GpSD0,GnSD0, &
           fp,fn,ap,an,sigmapSI,sigmanSI,sigmapSD,sigmanSD
  IF (PRESENT(m))     m     = WIMP%m
  IF (PRESENT(GpSI))  GpSI  = WIMP%GpSI
  IF (PRESENT(GnSI))  GnSI  = WIMP%GnSI
  IF (PRESENT(GpSD))  GpSD  = WIMP%GpSD
  IF (PRESENT(GnSD))  GnSD  = WIMP%GnSD
  IF (PRESENT(GpSI0)) GpSI0 = WIMP%GpSI0
  IF (PRESENT(GnSI0)) GnSI0 = WIMP%GnSI0
  IF (PRESENT(GpSD0)) GpSD0 = WIMP%GpSD0
  IF (PRESENT(GnSD0)) GnSD0 = WIMP%GnSD0
  IF (PRESENT(fp))    fp    = GToF(WIMP%GpSI)
  IF (PRESENT(fn))    fn    = GToF(WIMP%GnSI)
  IF (PRESENT(ap))    ap    = GToA(WIMP%GpSD)
  IF (PRESENT(an))    an    = GToA(WIMP%GnSD)
  IF (PRESENT(sigmapSI)) sigmapSI = GpToSigmapSI(WIMP%m,WIMP%GpSI)
  IF (PRESENT(sigmanSI)) sigmanSI = GnToSigmanSI(WIMP%m,WIMP%GnSI)
  IF (PRESENT(sigmapSD)) sigmapSD = GpToSigmapSD(WIMP%m,WIMP%GpSD)
  IF (PRESENT(sigmanSD)) sigmanSD = GnToSigmanSD(WIMP%m,WIMP%GnSD)
END SUBROUTINE


! ----------------------------------------------------------------------
! Set various WIMP quantities.
! 
! Optional input arguments:
!   m           WIMP mass [GeV].
!   GpSI        Spin-independent WIMP-proton coupling [GeV^-2].
!   GnSI        Spin-independent WIMP-neutron coupling [GeV^-2].
!   GpSD        Spin-dependent WIMP-proton coupling [GeV^-2].
!   GnSD        Spin-dependent WIMP-neutron coupling [GeV^-2].
!   fp          Spin-independent WIMP-proton coupling [GeV^-2].
!               Related by GpSI = 2 fp.
!   fn          Spin-independent WIMP-neutron coupling [GeV^-2].
!               Related by GnSI = 2 fn.
!   ap          Spin-dependent WIMP-proton coupling [unitless].
!               Related by GpSD = 2\sqrt{2} G_F ap.
!   an          Spin-dependent WIMP-neutron coupling [unitless].
!               Related by GnSD = 2\sqrt{2} G_F an.
! Optional cross-section arguments (give negative value to set
! corresponding coupling negative):
!   sigmapSI    Spin-independent WIMP-proton cross-section [pb].
!   sigmanSI    Spin-independent WIMP-neutron cross-section [pb].
!   sigmapSD    Spin-dependent WIMP-proton cross-section [pb].
!   sigmanSD    Spin-dependent WIMP-neutron cross-section [pb].
!   sigmaSI     Sets both sigmapSI and sigmanSI to the given value [pb].
!   sigmaSD     Sets both sigmapSD and sigmanSD to the given value [pb].
! 
SUBROUTINE SetWIMP(m,GpSI,GnSI,GpSD,GnSD,fp,fn,ap,an,                   &
                   sigmapSI,sigmanSI,sigmapSD,sigmanSD,sigmaSI,sigmaSD)
  IMPLICIT NONE
  REAL*8, INTENT(IN), OPTIONAL :: m,GpSI,GnSI,GpSD,GnSD,fp,fn,ap,an,    &
           sigmapSI,sigmanSI,sigmapSD,sigmanSD,sigmaSI,sigmaSD
  IF (PRESENT(m)) THEN
     WIMP%m = MAX(m,SQRT(TINY(1d0)))
     WIMP%GpSI0 = SigmapSIToGp(WIMP%m,1d0)
     WIMP%GnSI0 = SigmanSIToGn(WIMP%m,1d0)
     WIMP%GpSD0 = SigmapSDToGp(WIMP%m,1d0)
     WIMP%GnSD0 = SigmanSDToGn(WIMP%m,1d0)
  END IF
  IF (PRESENT(GpSI))  WIMP%GpSI = GpSI
  IF (PRESENT(GnSI))  WIMP%GnSI = GnSI
  IF (PRESENT(GpSD))  WIMP%GpSD = GpSD
  IF (PRESENT(GnSD))  WIMP%GnSD = GnSD
  IF (PRESENT(fp))    WIMP%GpSI = FToG(fp)
  IF (PRESENT(fn))    WIMP%GnSI = FToG(fn)
  IF (PRESENT(ap))    WIMP%GpSD = AToG(ap)
  IF (PRESENT(an))    WIMP%GnSD = AToG(an)
  IF (PRESENT(sigmapSI)) WIMP%GpSI = SigmapSIToGp(WIMP%m,sigmapSI)
  IF (PRESENT(sigmanSI)) WIMP%GnSI = SigmanSIToGn(WIMP%m,sigmanSI)
  IF (PRESENT(sigmapSD)) WIMP%GpSD = SigmapSDToGp(WIMP%m,sigmapSD)
  IF (PRESENT(sigmanSD)) WIMP%GnSD = SigmanSDToGn(WIMP%m,sigmanSD)
  IF (PRESENT(sigmaSI)) THEN
    WIMP%GpSI = SigmapSIToGp(WIMP%m,sigmaSI)
    WIMP%GnSI = SigmanSIToGn(WIMP%m,sigmaSI)
  END IF
  IF (PRESENT(sigmaSD)) THEN
    WIMP%GpSD = SigmapSDToGp(WIMP%m,sigmaSD)
    WIMP%GnSD = SigmanSDToGn(WIMP%m,sigmaSD)
  END IF
END SUBROUTINE


!-----------------------------------------------------------------------
! Initializes WIMP.
! Simply sets some default values for the WIMP parameters.
! 
SUBROUTINE InitWIMP()
  IMPLICIT NONE
  
  ! Default mass of 100 GeV
  CALL SetWIMP(m=100d0)
  
  ! Default cross-sections of 1 pb.
  WIMP%GpSI = WIMP%GpSI0
  WIMP%GnSI = WIMP%GnSI0
  WIMP%GpSD = WIMP%GpSD0
  WIMP%GnSD = WIMP%GnSD0
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Initializes WIMP from command-line parameters.
! 
! Possible options:
!   --m=<value>          ! WIMP mass [GeV]
!   --GpSI=<value>       ! Spin-independent WIMP-proton coupling [GeV^-2].
!   --GnSI=<value>       ! Spin-independent WIMP-neutron coupling [GeV^-2].
!   --GpSD=<value>       ! Spin-dependent WIMP-proton coupling [GeV^-2].
!   --GnSD=<value>       ! Spin-dependent WIMP-neutron coupling [GeV^-2].
!   --fp=<value>         ! Spin-independent WIMP-proton coupling [GeV^-2].
!                        ! Related by GpSI = 2 fp.
!   --fn=<value>         ! Spin-independent WIMP-neutron coupling [GeV^-2].
!                        ! Related by GnSI = 2 fn.
!   --ap=<value>         ! Spin-dependent WIMP-proton coupling [unitless].
!                        ! Related by GpSD = 2\sqrt{2} G_F ap.
!   --an=<value>         ! Spin-dependent WIMP-neutron coupling [unitless].
!                        ! Related by GnSD = 2\sqrt{2} G_F an.
! Cross-section options may be given as negative values to indicate the
! corresponding coupling should be negative:
!   --sigmapSI=<value>   ! Spin-independent WIMP-proton cross-section [pb].
!   --sigmanSI=<value>   ! Spin-independent WIMP-neutron cross-section [pb].
!   --sigmapSD=<value>   ! Spin-dependent WIMP-proton cross-section [pb].
!   --sigmanSD=<value>   ! Spin-dependent WIMP-neutron cross-section [pb].
!   --sigmaSI=<value>    ! Sets both sigmapSI and sigmanSI to the given value [pb].
!   --sigmaSD=<value>    ! Sets both sigmapSD and sigmanSD to the given value [pb].
! 
SUBROUTINE InitWIMPCL()
  IMPLICIT NONE
  LOGICAL :: status
  REAL*8 :: x
  
  ! Process mass: default value of 100 GeV
  CALL SetWIMP(m=100d0)
  IF (GetLongArgReal('m',x)) CALL SetWIMP(m=x)
  IF (Arguments%Nparameters .GE. 1) THEN
    x = Arguments%values(1)
    IF (x .GT. 0d0) CALL SetWIMP(m=x)
  END IF
  
  ! Process couplings: defaults of 1 pb.
  WIMP%GpSI = WIMP%GpSI0
  WIMP%GnSI = WIMP%GnSI0
  WIMP%GpSD = WIMP%GpSD0
  WIMP%GnSD = WIMP%GnSD0
  IF (GetLongArgReal('GpSI',x)) CALL SetWIMP(GpSI=x)
  IF (GetLongArgReal('GnSI',x)) CALL SetWIMP(GnSI=x)
  IF (GetLongArgReal('GpSD',x)) CALL SetWIMP(GpSD=x)
  IF (GetLongArgReal('GnSD',x)) CALL SetWIMP(GnSD=x)
  IF (GetLongArgReal('fp',x))   CALL SetWIMP(fp=x)
  IF (GetLongArgReal('fn',x))   CALL SetWIMP(fn=x)
  IF (GetLongArgReal('ap',x))   CALL SetWIMP(ap=x)
  IF (GetLongArgReal('an',x))   CALL SetWIMP(an=x)
  IF (GetLongArgReal('sigmapSI',x)) CALL SetWIMP(sigmapSI=x)
  IF (GetLongArgReal('sigmanSI',x)) CALL SetWIMP(sigmanSI=x)
  IF (GetLongArgReal('sigmapSD',x)) CALL SetWIMP(sigmapSD=x)
  IF (GetLongArgReal('sigmanSD',x)) CALL SetWIMP(sigmanSD=x)
  IF (GetLongArgReal('sigmaSI',x))  CALL SetWIMP(sigmaSI=x)
  IF (GetLongArgReal('sigmaSD',x))  CALL SetWIMP(sigmaSD=x)
  
  ! Process command-line arguments (if more than just mass)
  IF (Arguments%Nparameters .GE. 1) THEN
    status = ParseWIMPParameters(Arguments%Nparameters,Arguments%values)
  END IF
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Routine to read a line containing WIMP parameters from standard input
! (if line not given explicitly), parse it, and set WIMP parameters
! accordingly.  Returns true if a non-empty line was found and was
! parsable (false if no line found (EOF), line was empty, or line
! could not be parsed).
! 
! The following forms are allowed:
!   m
!   m sigmaSI
!   m sigmaSI sigmaSD
!   m sigmaSI sigmanSD sigmanSD
!   m sigmapSI sigmanSI sigmanSD sigmanSD
!   (formerly:  m GpSI GnSI GpSD GnSD)
! In the first case, all couplings are set to 1 pb, while in the second
! case, the SD couplings are set to zero.
! 
! Optional input argument:
!   line        String containing WIMP parameters.  If not given,
!               a line is taken from standard output.
! 
FUNCTION ParseWIMPInput(line) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: line
  CHARACTER(LEN=1024) :: line0
  INTEGER :: ios,Np
  REAL*8 :: params(6)
  
  status = .FALSE.
  
  ! Use given line or read from standard input
  IF (PRESENT(line)) THEN
    line0 = line
  ELSE
    READ(*,'(A)',IOSTAT=ios) line0
    IF (ios .NE. 0) RETURN
  END IF
  
  ! Check if empty string
  IF (TRIM(line0) .EQ. '') RETURN
  
  ! Determine number of parameters and read them in
  Np = MIN(NumberOfFields(line0),6)
  IF (Np .LE. 0) RETURN
  READ(line0,*,IOSTAT=ios) params(1:Np)
  IF (ios .NE. 0) RETURN
  
  ! Now parse parameters
  status = ParseWIMPParameters(Np,params)
  
END FUNCTION


!-----------------------------------------------------------------------
! Routine to take an array containing WIMP parameters and set the
! internal WIMP parameters to that.  The meaning of the given
! parameters depends upon the number of parameters and is the same
! as expected for the commandline.  Returns false if an invalid
! number of parameters or invalid value is found.
! 
! The following array lengths and parameters are allowed:
!   N=1:  m
!   N=2:  m sigmaSI
!   N=3:  m sigmaSI sigmaSD
!   N=4:  m sigmaSI sigmanSD sigmanSD
!   N=5:  m sigmapSI sigmanSI sigmanSD sigmanSD
!         (formerly:  m GpSI GnSI GpSD GnSD)
! In the first case, all couplings are set to 1 pb, while in the second
! case, the SD couplings are set to zero.
! 
! Input arguments:
!   N           Number of parameters
!   p           Array of parameters of size [1:N]
! 
FUNCTION ParseWIMPParameters(N,p) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: p(N)
  
  ! Check for bad cases (e.g. non-positive mass)
  status = .FALSE.
  IF (N .LT. 1) RETURN
  IF (p(1) .LE. 0d0) RETURN
  
  status = .TRUE.
  
  ! Meaning of parameters depends on number of parameters,
  ! but first parameter is always mass
  CALL SetWIMP(m=p(1))
  SELECT CASE (N)
  CASE (1)
    ! Form: m
    ! Set WIMP couplings to 1 pb
    WIMP%GpSI = WIMP%GpSI0
    WIMP%GnSI = WIMP%GnSI0
    WIMP%GpSD = WIMP%GpSD0
    WIMP%GnSD = WIMP%GnSD0
  CASE (2)
    ! Form: m sigmaSI
    ! Set SD couplings to zero
    CALL SetWIMP(sigmaSI=p(2),sigmaSD=0d0)
  CASE (3)
    ! Form: m sigmaSI sigmaSD
    CALL SetWIMP(sigmaSI=p(2),sigmaSD=p(3))
  CASE (4)
    ! Form: m sigmaSI sigmapSD sigmanSD
    CALL SetWIMP(sigmaSI=p(2),sigmapSD=p(3),sigmanSD=p(4))
  CASE (5)
    !! Form: m sigmapSI sigmanSI sigmapSD sigmanSD
    CALL SetWIMP(sigmapSI=p(2),sigmanSI=p(3),sigmapSD=p(4),sigmanSD=p(5))
    ! Form: m GpSI GnSI GpSD GnSD
    !CALL SetWIMP(GpSI=p(2),GnSI=p(3),GpSD=p(4),GnSD=p(5))
  CASE (6:)
    status = .FALSE.
  END SELECT
  
END FUNCTION



!=======================================================================
! COUPLING CONVERSION ROUTINES
!=======================================================================

!-----------------------------------------------------------------------
! Converts spin-independent coupling f to G, related by:
!   Gp = 2*fp    Gn = 2*fn
! 
! Input arguments:
!   f           SI coupling [GeV^-2]
! Returns in [GeV^-2].
! 
ELEMENTAL FUNCTION FToG(f) RESULT(G)
  IMPLICIT NONE
  REAL*8 :: G
  REAL*8, INTENT(IN) :: f
  REAL*8, PARAMETER :: F_SCALE = 2d0
  G = f * F_SCALE
END FUNCTION


!-----------------------------------------------------------------------
! Converts spin-independent coupling G to f, related by:
!   Gp = 2*fp    Gn = 2*fn
! 
! Input arguments:
!   G           SI coupling [GeV^-2]
! Returns in [GeV^-2].
! 
ELEMENTAL FUNCTION GToF(G) RESULT(f)
  IMPLICIT NONE
  REAL*8 :: f
  REAL*8, INTENT(IN) :: G
  REAL*8, PARAMETER :: F_SCALE = 2d0
  f = G / F_SCALE
END FUNCTION


!-----------------------------------------------------------------------
! Converts spin-dependent coupling a to G, related by:
!   Gp = 2*sqrt{2}*G_F*ap    Gn = 2*sqrt{2}*G_F*an
! 
! Input arguments:
!   a           SD coupling [unitless]
! Returns in [GeV^-2].
! 
ELEMENTAL FUNCTION AToG(a) RESULT(G)
  IMPLICIT NONE
  REAL*8 :: G
  REAL*8, INTENT(IN) :: a
  REAL*8, PARAMETER :: A_SCALE = SQRT(2d0)*2d0*FERMI_COUPLING_CONSTANT
  G = a * A_SCALE
END FUNCTION


!-----------------------------------------------------------------------
! Converts spin-dependent coupling G to a, related by:
!   Gp = 2*sqrt{2}*G_F*ap    Gn = 2*sqrt{2}*G_F*an
! 
! Input arguments:
!   G           SD coupling [GeV^-2]
! Returns in [unitless].
! 
ELEMENTAL FUNCTION GToA(G) RESULT(a)
  IMPLICIT NONE
  REAL*8 :: a
  REAL*8, INTENT(IN) :: G
  REAL*8, PARAMETER :: A_SCALE = SQRT(2d0)*2d0*FERMI_COUPLING_CONSTANT
  a = G / A_SCALE
END FUNCTION



!=======================================================================
! CROSS-SECTION/COUPLING CONVERSION ROUTINES
!=======================================================================

!-----------------------------------------------------------------------
! Converts WIMP-proton spin-independent cross-section to G, related by:
!   sigmapSI = \mu^2/\pi Gp^2
! where \mu is WIMP-proton reduced mass.
! 
! Input arguments:
!   m           WIMP mass [GeV]
!   sigma       SI cross-section [pb].  A negative value means
!               the coupling should be set negative.
! Returns in [GeV^-2].
! 
ELEMENTAL FUNCTION SigmapSIToGp(m,sigma) RESULT(G)
  IMPLICIT NONE
  REAL*8 :: G
  REAL*8, INTENT(IN) :: m,sigma
  REAL*8, PARAMETER :: Mp = PROTON_MASS
  ! Factor of 1d-10 takes pb to fm^2
  IF (sigma .GE. 0d0) THEN
    G = SQRTPI * (m+Mp)/(m*Mp) * SQRT(1d-10*sigma) / HBARC
  ELSE
    G = - SQRTPI * (m+Mp)/(m*Mp) * SQRT(-1d-10*sigma) / HBARC
  END IF
END FUNCTION


!-----------------------------------------------------------------------
! Converts WIMP-neutron spin-independent cross-section to G, related by:
!   sigmanSI = \mu^2/\pi Gn^2
! where \mu is WIMP-neutron reduced mass.
! 
! Input arguments:
!   m           WIMP mass [GeV]
!   sigma       SI cross-section [pb].  A negative value means
!               the coupling should be set negative.
! Returns in [GeV^-2].
! 
ELEMENTAL FUNCTION SigmanSIToGn(m,sigma) RESULT(G)
  IMPLICIT NONE
  REAL*8 :: G
  REAL*8, INTENT(IN) :: m,sigma
  REAL*8, PARAMETER :: Mn = NEUTRON_MASS
  ! Factor of 1d-10 takes pb to fm^2
  IF (sigma .GE. 0d0) THEN
    G = SQRTPI * (m+Mn)/(m*Mn) * SQRT(1d-10*sigma) / HBARC
  ELSE
    G = - SQRTPI * (m+Mn)/(m*Mn) * SQRT(-1d-10*sigma) / HBARC
  END IF
END FUNCTION


!-----------------------------------------------------------------------
! Converts G to WIMP-proton spin-independent cross-section, related by:
!   sigmapSI = \mu^2/\pi Gp^2
! where \mu is WIMP-proton reduced mass.
! 
! Input arguments:
!   m           WIMP mass [GeV]
!   G           SI coupling [GeV^-2]
! Returns in [pb].
! 
ELEMENTAL FUNCTION GpToSigmapSI(m,G) RESULT(sigma)
  IMPLICIT NONE
  REAL*8 :: sigma
  REAL*8, INTENT(IN) :: m,G
  REAL*8, PARAMETER :: Mp = PROTON_MASS
  ! Factor of 1d10 takes fm^2 to pb
  sigma = 1d10 * INVPI * ((m*Mp/(m+Mp))*G*HBARC)**2
END FUNCTION


!-----------------------------------------------------------------------
! Converts G to WIMP-neutron spin-independent cross-section, related by:
!   sigmanSI = \mu^2/\pi Gn^2
! where \mu is WIMP-neutron reduced mass.
! 
! Input arguments:
!   m           WIMP mass [GeV]
!   G           SI coupling [GeV^-2]
! Returns in [pb].
! 
ELEMENTAL FUNCTION GnToSigmanSI(m,G) RESULT(sigma)
  IMPLICIT NONE
  REAL*8 :: sigma
  REAL*8, INTENT(IN) :: m,G
  REAL*8, PARAMETER :: Mn = NEUTRON_MASS
  ! Factor of 1d10 takes fm^2 to pb
  sigma = 1d10 * INVPI * ((m*Mn/(m+Mn))*G*HBARC)**2
END FUNCTION


!-----------------------------------------------------------------------
! Converts WIMP-proton spin-dependent cross-section to G, related by:
!   sigmapSD = 3\mu^2/\pi Gp^2
! where \mu is WIMP-proton reduced mass.
! 
! Input arguments:
!   m           WIMP mass [GeV]
!   sigma       SD cross-section [pb].  A negative value means
!               the coupling should be set negative.
! Returns in [GeV^-2].
! 
ELEMENTAL FUNCTION SigmapSDToGp(m,sigma) RESULT(G)
  IMPLICIT NONE
  REAL*8 :: G
  REAL*8, INTENT(IN) :: m,sigma
  REAL*8, PARAMETER :: Mp = PROTON_MASS
  ! Factor of 1d-10 takes pb to fm^2
  IF (sigma .GE. 0d0) THEN
    G = SQRTPI/SQRT3 * (m+Mp)/(m*Mp) * SQRT(1d-10*sigma) / HBARC
  ELSE
    G = - SQRTPI/SQRT3 * (m+Mp)/(m*Mp) * SQRT(-1d-10*sigma) / HBARC
  END IF
END FUNCTION


!-----------------------------------------------------------------------
! Converts WIMP-neutron spin-dependent cross-section to G, related by:
!   sigmanSD = 3\mu^2/\pi Gn^2
! where \mu is WIMP-neutron reduced mass.
! 
! Input arguments:
!   m           WIMP mass [GeV]
!   sigma       SD cross-section [pb].  A negative value means
!               the coupling should be set negative.
! Returns in [GeV^-2].
! 
ELEMENTAL FUNCTION SigmanSDToGn(m,sigma) RESULT(G)
  IMPLICIT NONE
  REAL*8 :: G
  REAL*8, INTENT(IN) :: m,sigma
  REAL*8, PARAMETER :: Mn = NEUTRON_MASS
  ! Factor of 1d-10 takes pb to fm^2
  IF (sigma .GE. 0d0) THEN
    G = SQRTPI/SQRT3 * (m+Mn)/(m*Mn) * SQRT(1d-10*sigma) / HBARC
  ELSE
    G = - SQRTPI/SQRT3 * (m+Mn)/(m*Mn) * SQRT(-1d-10*sigma) / HBARC
  END IF
END FUNCTION


!-----------------------------------------------------------------------
! Converts G to WIMP-proton spin-dependent cross-section, related by:
!   sigmapSD = 3\mu^2/\pi Gp^2
! where \mu is WIMP-proton reduced mass.
! 
! Input arguments:
!   m           WIMP mass [GeV]
!   G           SD coupling [GeV^-2]
! Returns in [pb].
! 
ELEMENTAL FUNCTION GpToSigmapSD(m,G) RESULT(sigma)
  IMPLICIT NONE
  REAL*8 :: sigma
  REAL*8, INTENT(IN) :: m,G
  REAL*8, PARAMETER :: Mp = PROTON_MASS
  ! Factor of 1d10 takes fm^2 to pb
  sigma = 1d10 * 3*INVPI * ((m*Mp/(m+Mp))*G*HBARC)**2
END FUNCTION


!-----------------------------------------------------------------------
! Converts G to WIMP-neutron spin-dependent cross-section, related by:
!   sigmanSD = 3\mu^2/\pi Gn^2
! where \mu is WIMP-neutron reduced mass.
! 
! Input arguments:
!   m           WIMP mass [GeV]
!   G           SD coupling [GeV^-2]
! Returns in [pb].
! 
ELEMENTAL FUNCTION GnToSigmanSD(m,G) RESULT(sigma)
  IMPLICIT NONE
  REAL*8 :: sigma
  REAL*8, INTENT(IN) :: m,G
  REAL*8, PARAMETER :: Mn = NEUTRON_MASS
  ! Factor of 1d10 takes fm^2 to pb
  sigma = 1d10 * 3*INVPI * ((m*Mn/(m+Mn))*G*HBARC)**2
END FUNCTION



!=======================================================================
! HALO ROUTINES
!=======================================================================

! ----------------------------------------------------------------------
! Get various halo quantities.
! 
! Optional output arguments regarding galactic motions:
!   vrot        Local galactic disk rotation speed [km/s].
!   vlsr        Local standard of rest velocity vector (array of size 3)
!               [km/s], defined relative to galactic rest frame.
!   vpec        Sun's peculiar velocity vector (array of size 3) [km/s],
!               defined relative to local standard of rest.
!   vsun        Sun's velocity vector (array of size 3) [km/s], defined
!               relative to galactic rest frame.
! Optional output arguments regarding dark matter density:
!   rho         Local dark matter density [GeV/cm^3].
! Optional output arguments regarding SHM distribution, a truncated
! Maxwell-Boltzmann ("MB"):
!   vbulk       Bulk velocity of dark matter (array of size 3) [km/s],
!               defined relative to galactic rest frame.
!   vobs        Observer/detector's speed (i.e. Sun's speed) [km/s],
!               defined relative to MB rest frame.
!   v0          Most probable speed [km/s] in the MB rest frame.
!   vesc        Galactic/population escape speed [km/s] in the MB rest
!               frame (_galactic_ escape speed only if MB has no bulk
!               motion relative to galactic rest frame).
! Optional output arguments regarding tabulated eta(vmin):
!   tabulated   Indicates if a tabulated eta(vmin) is being used
!   eta_file    The file tabulated eta were taken from
!   Nvmin       Number of tabulation points
!   vmin        Allocatable array of vmin [km/s]
!   eta         Allocatable array of mean inverse speeds at vmin [s/km]
! 
SUBROUTINE GetHalo(vrot,vlsr,vpec,vsun,rho,vbulk,vobs,v0,vesc,          &
                   tabulated,eta_file,Nvmin,vmin,eta)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: eta_file
  LOGICAL, INTENT(OUT), OPTIONAL :: tabulated
  INTEGER, INTENT(OUT), OPTIONAL :: Nvmin
  REAL*8, INTENT(OUT), OPTIONAL :: vrot,vlsr(3),vpec(3),vsun(3),       &
                                   rho,vbulk(3),vobs,v0,vesc
  REAL*8, ALLOCATABLE, INTENT(OUT), OPTIONAL :: vmin(:),eta(:)
  
  IF (PRESENT(vrot))  vrot  = Halo%vrot
  IF (PRESENT(vlsr))  vlsr  = Halo%vlsr
  IF (PRESENT(vpec))  vpec  = Halo%vpec
  IF (PRESENT(vsun))  vsun  = Halo%vsun
  
  IF (PRESENT(rho))   rho   = Halo%rho
  
  IF (PRESENT(vbulk)) vbulk = Halo%vbulk
  IF (PRESENT(vobs))  vobs  = Halo%vobs
  IF (PRESENT(v0))    v0    = Halo%v0
  IF (PRESENT(vesc))  vesc  = Halo%vesc
  
  IF (PRESENT(tabulated)) tabulated = Halo%tabulated
  IF (PRESENT(eta_file)) eta_file = Halo%eta_file
  IF (PRESENT(Nvmin)) Nvmin = Halo%Nvmin
  IF (PRESENT(vmin)) THEN
    ALLOCATE(vmin(Halo%Nvmin))
    vmin = Halo%vmin
  END IF
  IF (PRESENT(eta)) THEN
    ALLOCATE(eta(Halo%Nvmin))
    eta = Halo%eta
  END IF
END SUBROUTINE


! ----------------------------------------------------------------------
! Set various halo quantities.
! 
! Optional input arguments regarding galactic motions:
!   vrot        Local galactic disk rotation speed [km/s].
!   vlsr        Local standard of rest velocity vector (array of size 3)
!               [km/s], defined relative to galactic rest frame.
!   vpec        Sun's peculiar velocity vector (array of size 3) [km/s],
!               defined relative to local standard of rest.
!   vsun        Sun's velocity vector (array of size 3) [km/s], defined
!               relative to galactic rest frame.
! Optional input arguments regarding dark matter density:
!   rho         Local dark matter density [GeV/cm^3].
! Optional input arguments regarding SHM distribution, a truncated
! Maxwell-Boltzmann ("MB"):
!   vbulk       Bulk velocity of dark matter (array of size 3) [km/s],
!               defined relative to galactic rest frame.
!   vobs        Observer/detector's speed (i.e. Sun's speed) [km/s],
!               defined relative to MB rest frame.
!   v0          Most probable speed [km/s] in the MB rest frame.
!   vesc        Galactic/population escape speed [km/s] in the MB rest
!               frame (_galactic_ escape speed only if MB has no bulk
!               motion relative to galactic rest frame).
! Optional tabulated eta(vmin) arguments.  Can be loaded from a given
! file or explicitly provided.  If provided, Nvmin, vmin, and eta must
! all be given to take effect.  When a tabulation is not provided, the
! mean inverse speed will be calculated explicitly (not tabulated!)
! using the SHM as described by the above parameters.
!   tabulated   Indicates if a tabulated eta(vmin) is to be used.  Implied
!               by the use of other tabulation arguments, but can be set
!               false to return to the SHM calculation after a tabulation
!               has been loaded.
!   eta_file    File from which tabulated eta(vmin) should be read;
!               default is to perform explicit calculations for the SHM
!               describe The file tabulated eta were taken from
!   eta_filename Sets the stored file name _without_ loading any data
!               from the file.
!   eta_file_K  The column number in the file to take eta from (default
!               is second)
!   Nvmin       Number of tabulation points
!   vmin        Array of size [1:Nvmin] containing tabulation vmin [km/s]
!   eta         Array of size [1:Nvmin] containing tabulated mean inverse
!               speeds at vmin [s/km]
! 
SUBROUTINE SetHalo(vrot,vlsr,vpec,vsun,vobs,rho,vbulk,v0,vesc,          &
                   tabulated,eta_file,eta_filename,eta_file_K,Nvmin,vmin,eta)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: eta_file,eta_filename
  LOGICAL, INTENT(IN), OPTIONAL :: tabulated
  INTEGER, INTENT(IN), OPTIONAL :: eta_file_K,Nvmin
  REAL*8, INTENT(IN), OPTIONAL :: vrot,vlsr(3),vpec(3),vsun(3),         &
                                  rho,vbulk(3),vobs,v0,vesc
  REAL*8, INTENT(IN), OPTIONAL :: vmin(:),eta(:)
  INTEGER :: K
  
  IF (PRESENT(vrot))  CALL SetDiskRotationSpeed(vrot)
  IF (PRESENT(vlsr))  CALL SetLocalStandardOfRest(vlsr)
  IF (PRESENT(vpec))  CALL SetSunPeculiarVelocity(vpec)
  IF (PRESENT(vsun))  CALL SetSunVelocity(vsun)
  
  IF (PRESENT(rho))   CALL SetLocalDensity(rho)
  
  IF (PRESENT(vbulk)) CALL SetBulkVelocity(vbulk)
  IF (PRESENT(vobs))  CALL SetObserverSpeed(vobs)
  IF (PRESENT(v0))    CALL SetMostProbableSpeed(v0)
  IF (PRESENT(vesc))  CALL SetEscapeSpeed(vesc)
  
  IF (PRESENT(tabulated)) THEN
    IF (Halo%Nvmin .GT. 0) Halo%tabulated = tabulated
  END IF
  IF (PRESENT(Nvmin) .AND. PRESENT(vmin) .AND. PRESENT(eta)) THEN
    IF (Nvmin .GT. 0) THEN
      IF (ALLOCATED(Halo%vmin)) DEALLOCATE(Halo%vmin)
      IF (ALLOCATED(Halo%eta))  DEALLOCATE(Halo%eta)
      Halo%Nvmin = Nvmin
      Halo%vmin  = vmin
      Halo%eta   = eta
      Halo%tabulated = .TRUE.
      Halo%eta_file  = ''
    END IF
  END IF
  IF (PRESENT(eta_file)) THEN
    IF (PRESENT(eta_file_K)) THEN
      K = eta_file_K
    ELSE
      K = 2
    END IF
    CALL LoadArrays(file=eta_file,N=Halo%Nvmin,N1=1,C1=Halo%vmin,       &
                    N2=K,C2=Halo%eta)
    Halo%tabulated = .TRUE.
    Halo%eta_file  = eta_file
  END IF
  
  IF (PRESENT(eta_filename)) Halo%eta_file = eta_filename
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Initializes halo.
! Simply sets halo parameters to default values.
! 
SUBROUTINE InitHalo()
  IMPLICIT NONE
  
  Halo%vrot  = 235d0
  Halo%vlsr  = (/ 0d0, 235d0, 0d0 /)
  Halo%vpec  = (/ 11d0, 12d0, 7d0 /)
  Halo%vsun  = (/ 0d0, 235d0, 0d0 /) + (/ 11d0, 12d0, 7d0 /)
  
  Halo%rho   = 0.4d0
  
  Halo%vbulk = (/ 0d0, 0d0, 0d0 /)
  Halo%vobs  = SQRT(11d0**2 + (235d0+12d0)**2 + 7d0**2)
  Halo%v0    = 235d0
  Halo%vesc  = 550d0
  
  Halo%tabulated = .FALSE.
  Halo%eta_file  = ''
  Halo%Nvmin = 0
  IF (ALLOCATED(Halo%vmin)) DEALLOCATE(Halo%vmin)
  IF (ALLOCATED(Halo%eta))  DEALLOCATE(Halo%eta)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Initializes halo from command-line parameters.
! 
! Possible options regarding galactic motions:
!   --vrot=<value>       ! Local galactic disk rotation speed [km/s].
!   --vlsr=<x>,<y>,<z>   ! Local standard of rest velocity vector (array of size 3)
!                        ! [km/s], defined relative to galactic rest frame.
!   --vpec=<x>,<y>,<z>   ! Sun's peculiar velocity vector (array of size 3) [km/s],
!                        ! defined relative to local standard of rest.
!   --vsun=<x>,<y>,<z>   ! Sun's velocity vector (array of size 3) [km/s], defined
!                        ! relative to galactic rest frame.
! Possible options regarding dark matter density:
!   --rho=<value>        ! Local dark matter density [GeV/cm^3]
! Possible options regarding SHM distribution:
!   --vbulk=<x>,<y>,<z>  ! Bulk velocity of dark matter (array of size 3) [km/s],
!                        ! defined relative to galactic rest frame.
!   --vobs=<value>       ! Observer/detector's speed (i.e. Sun's speed) [km/s],
!                        ! defined relative to MB rest frame.
!   --v0=<value>         ! Most probable speed [km/s] in the galactic rest frame.
!   --vesc=<value>       ! Escape speed of the dark matter population [km/s] in
!                        ! its rest frame.
! Possible options for provided a tabulated eta(vmin) instead of using the above
! SHM distribution.
!   --eta-file=<file>    ! File from which tabulated mean inverse speed eta(vmin)
!                        ! should be read.  First column is vmin [km/s] and second
!                        ! column is eta [s/km].  Default behavior is to do explicit
!                        ! calculations for SHM.
!   --eta-file=<file>,<K>! Same as above, but take the Kth column for eta.
! 
SUBROUTINE InitHaloCL()
  IMPLICIT NONE
  CHARACTER(LEN=1024) :: eta_file
  INTEGER :: I,K,Nval,ios
  REAL*8 :: vrot,vobs,rho,v0,vesc
  REAL*8, ALLOCATABLE :: vlsr(:),vpec(:),vsun(:),vbulk(:)
  ! Older compiler compatibility
  INTEGER, PARAMETER :: NCHAR = 1024
  CHARACTER(LEN=NCHAR), DIMENSION(:), ALLOCATABLE :: aval
  ! ...but this would be better better (needs gfortran 4.6+)
  !CHARACTER(LEN=:), DIMENSION(:), ALLOCATABLE :: aval
  
  CALL InitHalo()
  
  IF (GetLongArgReal('vrot',vrot)) CALL SetDiskRotationSpeed(vrot)
  IF (GetLongArgReals('vlsr',vlsr,I)) THEN
    IF (I .EQ. 3) THEN
      CALL SetLocalStandardOfRest(vlsr)
    ELSE
      WRITE(0,*) 'ERROR: Invalid --vlsr=<vx>,<vy>,<vz> parameter.'
      STOP
    END IF
  END IF
  IF (GetLongArgReals('vpec',vpec,I)) THEN
    IF (I .EQ. 3) THEN
      CALL SetSunPeculiarVelocity(vpec)
    ELSE
      WRITE(0,*) 'ERROR: Invalid --vpec=<vx>,<vy>,<vz> parameter.'
      STOP
    END IF
  END IF
  IF (GetLongArgReals('vsun',vsun,I)) THEN
    IF (I .EQ. 3) THEN
      CALL SetSunVelocity(vsun)
    ELSE
      WRITE(0,*) 'ERROR: Invalid --vsun=<vx>,<vy>,<vz> parameter.'
      STOP
    END IF
  END IF
  
  IF (GetLongArgReal('rho',rho))   CALL SetLocalDensity(rho)
  
  IF (GetLongArgReals('vbulk',vbulk,I)) THEN
    IF (I .EQ. 3) THEN
      CALL SetBulkVelocity(vbulk)
    ELSE
      WRITE(0,*) 'ERROR: Invalid --vbulk=<vx>,<vy>,<vz> parameter.'
      STOP
    END IF
  END IF
  IF (GetLongArgReal('vobs',vobs)) CALL SetObserverSpeed(vobs)
  IF (GetLongArgReal('v0',v0))     CALL SetMostProbableSpeed(v0)
  IF (GetLongArgReal('vesc',vesc)) CALL SetEscapeSpeed(vesc)
  
  !IF (GetLongArgString('eta-file',eta_file)) CALL SetHalo(eta_file=eta_file)
  IF (GetLongArgStrings('eta-file',NCHAR,aval,Nval)) THEN
    IF (Nval .GE. 2) THEN
      READ(aval(2),*,IOSTAT=ios) K
      IF (ios .NE. 0) K = 2
    ELSE
      K = 2
    END IF
    CALL SetHalo(eta_file=aval(1),eta_file_K=K)
  END IF
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Get/set local galactic disk rotation speed [km/s].
! Modifies the Local Standard of Rest (LSR) and Sun's velocity relative
! to halo rest frame as well as the most probable speed of the velocity
! distribution (v0 = vrot).  The observer speed is updated.
! 
PURE FUNCTION GetDiskRotationSpeed() RESULT(vrot)
  IMPLICIT NONE
  REAL*8 :: vrot
  vrot = Halo%vrot
END FUNCTION

SUBROUTINE SetDiskRotationSpeed(vrot)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: vrot
  Halo%vrot = vrot
  Halo%vlsr = (/ 0d0, Halo%vrot, 0d0 /)
  Halo%vsun = Halo%vlsr + Halo%vpec
  Halo%vobs = SQRT(SUM((Halo%vsun - Halo%vbulk)**2))
  Halo%v0   = Halo%vrot
END SUBROUTINE


! ----------------------------------------------------------------------
! Get/set Local Standard Of Rest velocity vector [km/s], defined
! relative to galactic rest frame.  Usually assumed to be (0,vrot,0),
! where vrot is disk rotation speed.  Modifies Sun's velocity relative
! to halo rest frame.  The disk rotation speed and the most probable
! speed of the velocity distribution are set to the y component of this
! velocity vector.  The observer speed is updated.
! 
PURE FUNCTION GetLocalStandardOfRest() RESULT(vlsr)
  IMPLICIT NONE
  REAL*8 :: vlsr(3)
  vlsr = Halo%vlsr
END FUNCTION

SUBROUTINE SetLocalStandardOfRest(vlsr)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: vlsr(3)
  Halo%vlsr = vlsr
  Halo%vsun = Halo%vlsr + Halo%vpec
  Halo%vrot = vlsr(2)
  Halo%vobs = SQRT(SUM((Halo%vsun - Halo%vbulk)**2))
  Halo%v0   = ABS(Halo%vrot)
END SUBROUTINE


! ----------------------------------------------------------------------
! Get/set Sun's peculiar velocity vector [km/s], defined relative to
! local standard of rest.  Modifies Sun's velocity relative to halo
! rest frame.  The observer speed is updated.
! 
PURE FUNCTION GetSunPeculiarVelocity() RESULT(vpec)
  IMPLICIT NONE
  REAL*8 :: vpec(3)
  vpec = Halo%vpec
END FUNCTION

SUBROUTINE SetSunPeculiarVelocity(vpec)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: vpec(3)
  Halo%vpec = vpec
  Halo%vsun = Halo%vlsr + Halo%vpec
  Halo%vobs = SQRT(SUM((Halo%vsun - Halo%vbulk)**2))
END SUBROUTINE


! ----------------------------------------------------------------------
! Get/set Sun's velocity vector [km/s], defined relative to galactic
! rest frame.  Normally taken to be vlsr + vpec, i.e. the sum of the
! local standard of rest and the Sun's peculiar velocity.  The preferred
! way to set speeds is modifying the disk rotation speed or Sun's
! peculiar velocity, not by setting this velocity directly, as the
! contributing velocities become ill-defined.  If the Sun's velocity is
! set here, the routine will attempt to assign a rotation speed vrot
! and local standard of rest vlsr = (0,vrot,0) that matches the given
! velocity vector, using the current value of the peculiar velocity; if
! not possible, the peculiar motion is set to zero first.  The most
! probable speed of the velocity distribution is updated to the
! resulting vrot and the observer speed is updated.
! 
PURE FUNCTION GetSunVelocity() RESULT(vsun)
  IMPLICIT NONE
  REAL*8 :: vsun(3)
  vsun = Halo%vsun
END FUNCTION

SUBROUTINE SetSunVelocity(vsun)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: vsun(3)
  REAL*8 :: vrot2
  Halo%vsun = vsun
  vrot2 = SUM((Halo%vsun - Halo%vpec)**2)
  IF (vrot2 .GE. 0d0) THEN
    Halo%vrot = SQRT(vrot2)
    Halo%vlsr = (/ 0d0, Halo%vrot, 0d0 /)
  ELSE
    Halo%vpec = 0d0
    Halo%vrot = SQRT(SUM(Halo%vsun**2))
    Halo%vlsr = (/ 0d0, Halo%vrot, 0d0 /)
  END IF
  Halo%v0 = ABS(Halo%vrot)
  Halo%vobs = SQRT(SUM((Halo%vsun - Halo%vbulk)**2))
END SUBROUTINE


! ----------------------------------------------------------------------
! Get/set local halo density [GeV/cm^3].
! 
PURE FUNCTION GetLocalDensity() RESULT(rho)
  IMPLICIT NONE
  REAL*8 :: rho
  rho = Halo%rho
END FUNCTION

SUBROUTINE SetLocalDensity(rho)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: rho
  Halo%rho = MAX(rho,0d0)
END SUBROUTINE


! ----------------------------------------------------------------------
! Get/set the dark matter population's bulk velocity vector [km/s],
! defined relative to the galactic rest frame.  Modifies the observer
! speed.
! 
PURE FUNCTION GetBulkVelocity() RESULT(vbulk)
  IMPLICIT NONE
  REAL*8 :: vbulk(3)
  vbulk = Halo%vbulk
END FUNCTION

SUBROUTINE SetBulkVelocity(vbulk)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: vbulk(3)
  Halo%vbulk = vbulk
  Halo%vobs = SQRT(SUM((Halo%vsun - Halo%vbulk)**2))
END SUBROUTINE


! ----------------------------------------------------------------------
! Get/set observer/detector's speed (i.e. Sun's speed) [km/s], defined
! relative to Maxwell-Boltzmann population rest frame.  Normally taken
! to be |vlsr + vpec - vMB|, i.e. the sum of the local standard of rest
! and the Sun's peculiar velocity less the bulk velocity of the dark
! matter population.  The preferred way to set speeds is modifying the
! disk rotation speed, Sun's peculiar velocity, or the bulk dark matter
! motion, not by setting this speed directly, as the various
! velocities become ill-defined.  If the observer's speed is set here,
! the routine will set the bulk motion of the DM to zero (relative to
! the galactic rest frame) and attempt to assign a rotation speed vrot
! and local standard of rest vlsr = (0,vrot,0) that matches the given
! speed, using the current value of the peculiar velocity; if not
! possible, the peculiar motion is set to zero first.  The most
! probable speed of the velocity distribution is updated to the
! resulting vrot.
! 
PURE FUNCTION GetObserverSpeed() RESULT(vobs)
  IMPLICIT NONE
  REAL*8 :: vobs
  vobs = Halo%vobs
END FUNCTION

SUBROUTINE SetObserverSpeed(vobs)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: vobs
  REAL*8 :: vy2
  Halo%vobs  = MAX(vobs,0d0)
  Halo%vbulk = (/ 0d0, 0d0, 0d0 /)
  vy2 = Halo%vobs**2 - Halo%vpec(1)**2 - Halo%vpec(3)**2
  IF (vy2 .GE. 0d0) THEN
    Halo%vrot = SQRT(vy2) - Halo%vpec(2)
    Halo%vlsr = (/ 0d0, Halo%vrot, 0d0 /)
    Halo%vsun = Halo%vlsr + Halo%vpec
  ELSE
    Halo%vpec = 0d0
    Halo%vrot = Halo%vobs
    Halo%vlsr = (/ 0d0, Halo%vrot, 0d0 /)
    Halo%vsun = Halo%vlsr + Halo%vpec
  END IF
  Halo%v0 = ABS(Halo%vrot)
END SUBROUTINE


! ----------------------------------------------------------------------
! Get/set most probable speed v0 [km/s] in the dark matter population's
! rest frame. Related to other speeds characterizing velocity
! distribution by:
!     vrms = sqrt(3/2) v0    [rms velocity]
!     vmp  = v0              [most probably velocity]
!     vave = sqrt(4/pi) v0   [mean velocity]
! 
PURE FUNCTION GetMostProbableSpeed() RESULT(v0)
  IMPLICIT NONE
  REAL*8 :: v0
  v0 = Halo%v0
END FUNCTION

SUBROUTINE SetMostProbableSpeed(v0)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: v0
  Halo%v0 = v0
END SUBROUTINE


! ----------------------------------------------------------------------
! Get/set dark matter population escape speed [km/s].  In the case of
! the SHM with no bulk motion relative to the galactic rest frame, this
! is the galactic escape speed.
! 
PURE FUNCTION GetEscapeSpeed() RESULT(vesc)
  IMPLICIT NONE
  REAL*8 :: vesc
  vesc = Halo%vesc
END FUNCTION

SUBROUTINE SetEscapeSpeed(vesc)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: vesc
  Halo%vesc = MAX(vesc,0d0)
END SUBROUTINE


! ----------------------------------------------------------------------
! INTERFACE NAME: EToVmin
! Calculate the minimum velocity for producing a recoil of energy E,
! given by vmin = sqrt{M E/(2\mu^2)} [km/s].
! 
! This is the scalar version (single mass and energy).
! 
! Input arguments:
!   E           Recoil energy [keV]
!   m           WIMP mass [GeV]
!   Miso        Isotope mass [GeV]
! 
ELEMENTAL FUNCTION EToVmin0(E,m,Miso) RESULT(vmin)
  IMPLICIT NONE
  REAL*8 :: vmin
  REAL*8, INTENT(IN) :: E,m,Miso
  REAL*8 :: mu
  REAL*8, PARAMETER :: c = 1d-3*SPEED_OF_LIGHT  ! Speed of light in km/s
  mu = Miso*m / (Miso + m)
  vmin = c * SQRT(1d-6*Miso*E/(2*mu**2))
END FUNCTION


! ----------------------------------------------------------------------
! INTERFACE NAME: EToVmin
! Calculate the minimum velocity for producing a recoil of energy E,
! given by vmin = sqrt{M E/(2\mu^2)} [km/s].  Returns as array of
! size [1:N].
! 
! This is the 1D array version (single mass and array of energies).
! 
! Input arguments:
!   N           Number of recoil energies
!   E           Array of recoil energies [keV]
!   m           WIMP mass [GeV]
!   Miso        Isotope mass [GeV]
! 
PURE FUNCTION EToVmin1(N,E,m,Miso) RESULT(vmin)
  IMPLICIT NONE
  REAL*8 :: vmin(N)
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: E(N),m,Miso
  REAL*8 :: mu
  REAL*8, PARAMETER :: c = 1d-3*SPEED_OF_LIGHT  ! Speed of light in km/s
  mu = Miso*m / (Miso + m)
  vmin = c * SQRT(1d-6*Miso*E/(2*mu**2))
END FUNCTION


! ----------------------------------------------------------------------
! INTERFACE NAME: EToVmin
! Calculate the minimum velocity for producing a recoil of energy E,
! given by vmin = sqrt{M E/(2\mu^2)} [km/s].  Returns as array of
! size [1:N,1:Niso].
! 
! This is the 2D array version (multiple masses and array of energies).
! 
! Input arguments:
!   N           Number of recoil energies
!   E           Array of recoil energies [keV]
!   m           WIMP mass [GeV]
!   Niso        Number of isotopes
!   Miso        Array of isotope masses [GeV]
! 
PURE FUNCTION EToVmin2(N,E,m,Niso,Miso) RESULT(vmin)
  IMPLICIT NONE
  REAL*8 :: vmin(N,Niso)
  INTEGER, INTENT(IN) :: N,Niso
  REAL*8, INTENT(IN) :: E(N),m,Miso(Niso)
  INTEGER :: I
  REAL*8 :: mu(N)
  REAL*8, PARAMETER :: c = 1d-3*SPEED_OF_LIGHT  ! Speed of light in km/s
  mu = Miso*m / (Miso + m)
  DO I = 1,Niso
    vmin(:,I) = c * SQRT(1d-6*Miso(I)*E/(2*mu(I)**2))
  END DO
END FUNCTION


!-----------------------------------------------------------------------
! INTERFACE NAME: MeanInverseSpeed
! Calculates the mean inverse speed (eta) [s/km] for the given vmin,
! with eta define as:
!     eta(vmin) = \int_{|v|>vmin} d^3v 1/|v| f(v)
! 
! This is the scalar version (single vmin).
! 
! Input arguments:
!   vmin        The minimum speed in the eta integral [km/s]
! 
ELEMENTAL FUNCTION MeanInverseSpeed0(vmin) RESULT(eta)
  IMPLICIT NONE
  REAL*8 :: eta
  REAL*8, INTENT(IN) :: vmin
  REAL*8 :: v0,vobs,vesc,x,y,z,Nesc
  
  ! If have tabulation, use it
  IF (Halo%tabulated) THEN
    eta = MeanInverseSpeedT(vmin)
    RETURN
  END IF
  
  ! Easier to use variable names
  v0   = Halo%v0
  vobs = Halo%vobs
  vesc = Halo%vesc
  
  ! Special case: no dispersion
  ! Distribution is delta function
  IF (v0 .EQ. 0) THEN
    IF (vobs .EQ. 0d0) THEN
      eta = 0d0
    ELSE
      IF (vmin .LE. vobs) THEN
        eta = 1d0 / vobs
      ELSE
        eta = 0d0
      END IF
    END IF
    RETURN
  END IF
  
  x    = vmin / v0
  y    = vobs / v0
  z    = vesc / v0
  Nesc = ERF(z) - 2*INVSQRTPI*z*EXP(-z**2)
  
  ! Special case: no relative motion by observer
  !   eta = 2/(sqrt(pi) Nesc v0) [e^{-x^2} - e^{-z^2}]
  ! Note: EXP2(a,b) = e^b - e^a
  IF (y .EQ. 0d0) THEN
    IF (x .LE. z) THEN
      eta = 2*INVSQRTPI/(Nesc*v0) * EXP2(-z**2,-x**2)
    ELSE
      eta = 0d0
    END IF
    RETURN
  END IF
  
  ! Special case: no finite cutoff (vesc is effectively infinite)
  IF (z .GT. 25d0) THEN
    eta = ERF2(x-y,x+y) / (2*vobs)
    RETURN
  END IF
  
  ! General case.
  ! See e.g. Savage, Freese & Gondolo, PRD 74, 043531 (2006)
  ! [astrop-ph/0607121]; use arxiv version as PRD version has type-
  ! setting issues in the formula.
  ! Note: ERF2(a,b) = ERF(b) - ERF(a)
  IF (x .LE. ABS(y-z)) THEN
    IF (y .LT. z) THEN
      eta = 1d0 / (2*Nesc*vobs) * (ERF2(x-y,x+y) - 4*INVSQRTPI*y*EXP(-z**2))
    ELSE
      eta = 1d0 / vobs
    END IF
  ELSE IF (x .LE. y+z) THEN
    eta = 1d0 / (2*Nesc*vobs) * (ERF2(x-y,z) - 2*INVSQRTPI*(z+y-x)*EXP(-z**2))
  ELSE
    eta = 0d0
  END IF
  
END FUNCTION


!-----------------------------------------------------------------------
! INTERFACE NAME: MeanInverseSpeed
! Calculates the mean inverse speed (eta) [s/km] for the given 1D
! array of vmin, with eta define as:
!     eta(vmin) = \int_{|v|>vmin} d^3v 1/|v| f(v)
! Returns as array of size [1:N].
! 
! This is the 1D array version (1D array of vmin).
! 
! Input arguments:
!   N           Number of vmin
!   vmin        The minimum speed in the eta integral [km/s].
!               Array of size [1:N].
! 
PURE FUNCTION MeanInverseSpeed1(N,vmin) RESULT(eta)
  IMPLICIT NONE
  REAL*8 :: eta(N)
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: vmin(N)
  REAL*8 :: v0,vobs,vesc,x(N),y,z,Nesc
  
  ! If have tabulation, use it
  IF (Halo%tabulated) THEN
    eta = MeanInverseSpeedT(vmin)
    RETURN
  END IF
  
  ! Easier to use variable names
  v0   = Halo%v0
  vobs = Halo%vobs
  vesc = Halo%vesc
  
  ! Special case: no dispersion
  ! Distribution is delta function
  IF (v0 .EQ. 0) THEN
    IF (vobs .EQ. 0d0) THEN
      eta = 0d0
    ELSE
      WHERE (vmin .LE. vobs)
        eta = 1d0 / vobs
      ELSE WHERE
        eta = 0d0
      END WHERE
    END IF
    RETURN
  END IF
  
  x    = vmin / v0
  y    = vobs / v0
  z    = vesc / v0
  Nesc = ERF(z) - 2*INVSQRTPI*z*EXP(-z**2)
  
  ! Special case: no relative motion by observer
  !   eta = 2/(sqrt(pi) Nesc v0) [e^{-x^2} - e^{-z^2}]
  ! Note: EXP2(a,b) = e^b - e^a
  IF (y .EQ. 0d0) THEN
    WHERE (x .LE. z)
      eta = 2*INVSQRTPI/(Nesc*v0) * EXP2(-z**2,-x**2)
    ELSE WHERE
      eta = 0d0
    END WHERE
    RETURN
  END IF
  
  ! Special case: no finite cutoff (vesc is effectively infinite)
  IF (z .GT. 25d0) THEN
    eta = ERF2(x-y,x+y) / (2*vobs)
    RETURN
  END IF
  
  ! General case.
  ! See e.g. Savage, Freese & Gondolo, PRD 74, 043531 (2006)
  ! [astrop-ph/0607121]; use arxiv version as PRD version has type-
  ! setting issues in the formula.
  ! Note: ERF2(a,b) = ERF(b) - ERF(a)
  ! Separate y < z & y > z cases to make easier use of WHERE statements.
  IF (y .LT. z) THEN
    WHERE (x .LT. z-y)
      eta = 1d0 / (2*Nesc*vobs) * (ERF2(x-y,x+y) - 4*INVSQRTPI*y*EXP(-z**2))
    ELSE WHERE (x .LT. z+y)
      eta = 1d0 / (2*Nesc*vobs) * (ERF2(x-y,z) - 2*INVSQRTPI*(z+y-x)*EXP(-z**2))
    ELSE WHERE
      eta = 0d0
    END WHERE
  ELSE
    WHERE (x .LT. y-z)
      eta = 1d0 / vobs
    ELSE WHERE (x .LT. y+z)
      eta = 1d0 / (2*Nesc*vobs) * (ERF2(x-y,z) - 2*INVSQRTPI*(z+y-x)*EXP(-z**2))
    ELSE WHERE
      eta = 0d0
    END WHERE
  END IF
  
END FUNCTION


!-----------------------------------------------------------------------
! INTERFACE NAME: MeanInverseSpeed
! Calculates the mean inverse speed (eta) [s/km] for the given 2D
! array of vmin, with eta define as:
!     eta(vmin) = \int_{|v|>vmin} d^3v 1/|v| f(v)
! Returns as array of size [1:N1,1:N2].
! 
! This is the 2D array version (2D array of vmin).
! 
! Input arguments:
!   N1,N2       Size of vmin and eta arrays, i.e. [1:N1,1:N2]
!   vmin        The minimum speed in the eta integral [km/s].
!               Array of size [1:N].
! 
PURE FUNCTION MeanInverseSpeed2(N1,N2,vmin) RESULT(eta)
  IMPLICIT NONE
  REAL*8 :: eta(N1,N2)
  INTEGER, INTENT(IN) :: N1,N2
  REAL*8, INTENT(IN) :: vmin(N1,N2)
  REAL*8 :: v0,vobs,vesc,x(N1,N2),y,z,Nesc
  
  ! If have tabulation, use it
  IF (Halo%tabulated) THEN
    eta = MeanInverseSpeedT(vmin)
    RETURN
  END IF
  
  ! Easier to use variable names
  v0   = Halo%v0
  vobs = Halo%vobs
  vesc = Halo%vesc
  
  ! Special case: no dispersion
  ! Distribution is delta function
  IF (v0 .EQ. 0) THEN
    IF (vobs .EQ. 0d0) THEN
      eta = 0d0
    ELSE
      WHERE (vmin .LE. vobs)
        eta = 1d0 / vobs
      ELSE WHERE
        eta = 0d0
      END WHERE
    END IF
    RETURN
  END IF
  
  x    = vmin / v0
  y    = vobs / v0
  z    = vesc / v0
  Nesc = ERF(z) - 2*INVSQRTPI*z*EXP(-z**2)
  
  ! Special case: no relative motion by observer
  !   eta = 2/(sqrt(pi) Nesc v0) [e^{-x^2} - e^{-z^2}]
  ! Note: EXP2(a,b) = e^b - e^a
  IF (y .EQ. 0d0) THEN
    WHERE (x .LE. z)
      eta = 2*INVSQRTPI/(Nesc*v0) * EXP2(-z**2,-x**2)
    ELSE WHERE
      eta = 0d0
    END WHERE
    RETURN
  END IF
  
  ! Special case: no finite cutoff (vesc is effectively infinite)
  IF (z .GT. 25d0) THEN
    eta = ERF2(x-y,x+y) / (2*vobs)
    RETURN
  END IF
  
  ! General case.
  ! See e.g. Savage, Freese & Gondolo, PRD 74, 043531 (2006)
  ! [astrop-ph/0607121]; use arxiv version as PRD version has type-
  ! setting issues in the formula.
  ! Note: ERF2(a,b) = ERF(b) - ERF(a)
  ! Separate y < z & y > z cases to make easier use of WHERE statements.
  IF (y .LT. z) THEN
    WHERE (x .LT. z-y)
      eta = 1d0 / (2*Nesc*vobs) * (ERF2(x-y,x+y) - 4*INVSQRTPI*y*EXP(-z**2))
    ELSE WHERE (x .LT. z+y)
      eta = 1d0 / (2*Nesc*vobs) * (ERF2(x-y,z) - 2*INVSQRTPI*(z+y-x)*EXP(-z**2))
    ELSE WHERE
      eta = 0d0
    END WHERE
  ELSE
    WHERE (x .LT. y-z)
      eta = 1d0 / vobs
    ELSE WHERE (x .LT. y+z)
      eta = 1d0 / (2*Nesc*vobs) * (ERF2(x-y,z) - 2*INVSQRTPI*(z+y-x)*EXP(-z**2))
    ELSE WHERE
      eta = 0d0
    END WHERE
  END IF
  
END FUNCTION


!-----------------------------------------------------------------------
! Calculates the mean inverse speed (eta) [s/km] for the given vmin,
! with eta define as:
!     eta(vmin) = \int_{|v|>vmin} d^3v 1/|v| f(v)
! using the stored tabulation rather than the explicit calculation.
! 
! Input arguments:
!   vmin        The minimum speed in the eta integral [km/s]
! 
ELEMENTAL FUNCTION MeanInverseSpeedT(vmin) RESULT(eta)
  IMPLICIT NONE
  REAL*8 :: eta
  REAL*8, INTENT(IN) :: vmin
  INTEGER :: K
  REAL*8 :: f
  
  IF (.NOT. Halo%tabulated .OR. (Halo%Nvmin .LE. 0)) THEN
    eta = 0d0
    RETURN
  END IF
  
  K = BSearch(Halo%Nvmin,Halo%vmin,vmin)
  
  IF (K .LE. 0) THEN
    eta = Halo%eta(1)
  ELSE IF (K .GE. Halo%Nvmin) THEN
    IF (vmin .EQ. Halo%vmin(Halo%Nvmin)) THEN
      eta = Halo%eta(Halo%Nvmin)
    ELSE
      eta = 0d0
    END IF
  ELSE IF (Halo%vmin(K) .EQ. Halo%vmin(K+1)) THEN
    eta = Halo%eta(K)
  ELSE
    f = (vmin-Halo%vmin(K)) / (Halo%vmin(K+1)-Halo%vmin(K))
    eta = (1-f)*Halo%eta(K) + f*Halo%eta(K+1)
  END IF
  
END FUNCTION



!=======================================================================
! ISOTOPES/NUCLEAR ROUTINES (INCLUDING FORM FACTORS)
!=======================================================================

! ----------------------------------------------------------------------
! For the given element, fills in allocatable arrays containing isotopic
! atomic numbers (Z), atomic masses (A), mass fractions (f), and
! masses (M).
! 
! Input argument:
!     Z          Atomic number of element
! Output arguments:
!     Niso       Number of isotopes
!     Ziso       Allocatable array (integer) of isotopes' atomic numbers (Z)
!     Aiso       Allocatable array (integer) of isotopes' atomic masses (A)
!     fiso       Allocatable array of isotopes' mass fractions
!     Miso       Allocatable array of isotopes' nuclear masses [GeV]
! 
SUBROUTINE ElementIsotopeList(Z,Niso,Ziso,Aiso,fiso,Miso)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Z
  INTEGER, INTENT(OUT) :: Niso
  INTEGER, ALLOCATABLE, INTENT(OUT) :: Ziso(:),Aiso(:)
  REAL*8, ALLOCATABLE, INTENT(OUT) :: fiso(:),Miso(:)
  INTEGER :: I1,I2
  ! Isotope data for all elements up to Z=92 (Uranium).
  ! Data is combined into single arrays; the ELEMENT_INDEX indicates
  ! where data for a particular element Z begins in those arrays.
  INTEGER, PARAMETER :: NELEMENTS = 92
  INTEGER, PARAMETER :: NISOTOPES = 286
  ! Number of stable isotopes for given element (indexed by Z)
  INTEGER, PARAMETER :: ELEMENT_NISO(NELEMENTS) =                       &
    (/  2,  2,  2,  1,  2,  2,  2,  3,  1,  3,  1,  3,  1,  3,  1,  4,  &
        2,  3,  3,  6,  1,  5,  2,  4,  1,  4,  1,  5,  2,  5,  2,  5,  &
        1,  6,  2,  6,  2,  4,  1,  5,  1,  7,  0,  7,  1,  6,  2,  8,  &
        2, 10,  2,  8,  1,  9,  1,  7,  2,  4,  1,  7,  0,  7,  2,  7,  &
        1,  7,  1,  6,  1,  7,  2,  6,  1,  5,  2,  7,  2,  6,  1,  7,  &
        2,  4,  1,  0,  0,  0,  0,  0,  0,  1,  0,  3 /)
  ! First data array index for given element (indexed by Z)
  INTEGER, PARAMETER :: ELEMENT_INDEX(NELEMENTS) =                      &
    (/  1,  3,  5,  7,  8, 10, 12, 14, 17, 18, 21, 22, 25, 26, 29, 30,  &
       34, 36, 39, 42, 48, 49, 54, 56, 60, 61, 65, 66, 71, 73, 78, 80,  &
       85, 86, 92, 94,100,102,106,107,112,113,120,120,127,128,134,136,  &
      144,146,156,158,166,167,176,177,184,186,190,191,198,198,205,207,  &
      214,215,222,223,229,230,237,239,245,246,251,253,260,262,268,269,  &
      276,278,282,283,283,283,283,283,283,283,284,284 /)
  ! Atomic number for an isotope
  INTEGER, PARAMETER :: ISOTOPE_Z(NISOTOPES) =                          &
    (/  1,  1,  2,  2,  3,  3,  4,  5,  5,  6,  6,  7,  7,  8,  8,  8,  &
        9, 10, 10, 10, 11, 12, 12, 12, 13, 14, 14, 14, 15, 16, 16, 16,  &
       16, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 20, 20, 20, 21,  &
       22, 22, 22, 22, 22, 23, 23, 24, 24, 24, 24, 25, 26, 26, 26, 26,  &
       27, 28, 28, 28, 28, 28, 29, 29, 30, 30, 30, 30, 30, 31, 31, 32,  &
       32, 32, 32, 32, 33, 34, 34, 34, 34, 34, 34, 35, 35, 36, 36, 36,  &
       36, 36, 36, 37, 37, 38, 38, 38, 38, 39, 40, 40, 40, 40, 40, 41,  &
       42, 42, 42, 42, 42, 42, 42, 44, 44, 44, 44, 44, 44, 44, 45, 46,  &
       46, 46, 46, 46, 46, 47, 47, 48, 48, 48, 48, 48, 48, 48, 48, 49,  &
       49, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 51, 51, 52, 52, 52,  &
       52, 52, 52, 52, 52, 53, 54, 54, 54, 54, 54, 54, 54, 54, 54, 55,  &
       56, 56, 56, 56, 56, 56, 56, 57, 57, 58, 58, 58, 58, 59, 60, 60,  &
       60, 60, 60, 60, 60, 62, 62, 62, 62, 62, 62, 62, 63, 63, 64, 64,  &
       64, 64, 64, 64, 64, 65, 66, 66, 66, 66, 66, 66, 66, 67, 68, 68,  &
       68, 68, 68, 68, 69, 70, 70, 70, 70, 70, 70, 70, 71, 71, 72, 72,  &
       72, 72, 72, 72, 73, 74, 74, 74, 74, 74, 75, 75, 76, 76, 76, 76,  &
       76, 76, 76, 77, 77, 78, 78, 78, 78, 78, 78, 79, 80, 80, 80, 80,  &
       80, 80, 80, 81, 81, 82, 82, 82, 82, 83, 90, 92, 92, 92 /)
  ! Atomic mass number for an isotope
  INTEGER, PARAMETER :: ISOTOPE_A(NISOTOPES) =                          &
    (/  1,  2,  3,  4,  6,  7,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18,  &
       19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,  &
       36, 35, 37, 36, 38, 40, 39, 40, 41, 40, 42, 43, 44, 46, 48, 45,  &
       46, 47, 48, 49, 50, 50, 51, 50, 52, 53, 54, 55, 54, 56, 57, 58,  &
       59, 58, 60, 61, 62, 64, 63, 65, 64, 66, 67, 68, 70, 69, 71, 70,  &
       72, 73, 74, 76, 75, 74, 76, 77, 78, 80, 82, 79, 81, 78, 80, 82,  &
       83, 84, 86, 85, 87, 84, 86, 87, 88, 89, 90, 91, 92, 94, 96, 93,  &
       92, 94, 95, 96, 97, 98,100, 96, 98, 99,100,101,102,104,103,102,  &
      104,105,106,108,110,107,109,106,108,110,111,112,113,114,116,113,  &
      115,112,114,115,116,117,118,119,120,122,124,121,123,120,122,123,  &
      124,125,126,128,130,127,124,126,128,129,130,131,132,134,136,133,  &
      130,132,134,135,136,137,138,138,139,136,138,140,142,141,142,143,  &
      144,145,146,148,150,144,147,148,149,150,152,154,151,153,152,154,  &
      155,156,157,158,160,159,156,158,160,161,162,163,164,165,162,164,  &
      166,167,168,170,169,168,170,171,172,173,174,176,175,176,174,176,  &
      177,178,179,180,181,180,182,183,184,186,185,187,184,186,187,188,  &
      189,190,192,191,193,190,192,194,195,196,198,197,196,198,199,200,  &
      201,202,204,203,205,204,206,207,208,209,232,234,235,238 /)
  ! Atomic spin for an isotope
  REAL*8, PARAMETER :: ISOTOPE_J(NISOTOPES) =                           &
    (/0.5d0,1.0d0,0.5d0,0.0d0,1.0d0,1.5d0,1.5d0,3.0d0,1.5d0,0.0d0,      &
      0.5d0,1.0d0,0.5d0,0.0d0,2.5d0,0.0d0,0.5d0,0.0d0,1.5d0,0.0d0,      &
      1.5d0,0.0d0,2.5d0,0.0d0,2.5d0,0.0d0,0.5d0,0.0d0,0.5d0,0.0d0,      &
      1.5d0,0.0d0,0.0d0,1.5d0,1.5d0,0.0d0,0.0d0,0.0d0,1.5d0,4.0d0,      &
      1.5d0,0.0d0,0.0d0,3.5d0,0.0d0,0.0d0,0.0d0,3.5d0,0.0d0,2.5d0,      &
      0.0d0,3.5d0,0.0d0,6.0d0,3.5d0,0.0d0,0.0d0,1.5d0,0.0d0,2.5d0,      &
      0.0d0,0.0d0,0.5d0,0.0d0,3.5d0,0.0d0,0.0d0,1.5d0,0.0d0,0.0d0,      &
      1.5d0,1.5d0,0.0d0,0.0d0,2.5d0,0.0d0,0.0d0,1.5d0,1.5d0,0.0d0,      &
      0.0d0,4.5d0,0.0d0,0.0d0,1.5d0,0.0d0,0.0d0,0.5d0,0.0d0,0.0d0,      &
      0.0d0,1.5d0,1.5d0,0.0d0,0.0d0,0.0d0,4.5d0,0.0d0,0.0d0,2.5d0,      &
      1.5d0,0.0d0,0.0d0,4.5d0,0.0d0,0.5d0,0.0d0,2.5d0,0.0d0,0.0d0,      &
      0.0d0,4.5d0,0.0d0,0.0d0,2.5d0,0.0d0,2.5d0,0.0d0,0.0d0,0.0d0,      &
      0.0d0,2.5d0,0.0d0,2.5d0,0.0d0,0.0d0,0.5d0,0.0d0,0.0d0,2.5d0,      &
      0.0d0,0.0d0,0.0d0,0.5d0,0.5d0,0.0d0,0.0d0,0.0d0,0.5d0,0.0d0,      &
      0.5d0,0.0d0,0.0d0,4.5d0,4.5d0,0.0d0,0.0d0,0.5d0,0.0d0,0.5d0,      &
      0.0d0,0.5d0,0.0d0,0.0d0,0.0d0,2.5d0,3.5d0,0.0d0,0.0d0,0.5d0,      &
      0.0d0,0.5d0,0.0d0,0.0d0,0.0d0,2.5d0,0.0d0,0.0d0,0.0d0,0.5d0,      &
      0.0d0,1.5d0,0.0d0,0.0d0,0.0d0,3.5d0,0.0d0,0.0d0,0.0d0,1.5d0,      &
      0.0d0,1.5d0,0.0d0,5.0d0,3.5d0,0.0d0,0.0d0,0.0d0,0.0d0,2.5d0,      &
      0.0d0,3.5d0,0.0d0,3.5d0,0.0d0,0.0d0,0.0d0,0.0d0,3.5d0,0.0d0,      &
      3.5d0,0.0d0,0.0d0,0.0d0,2.5d0,2.5d0,0.0d0,0.0d0,1.5d0,0.0d0,      &
      1.5d0,0.0d0,0.0d0,1.5d0,0.0d0,0.0d0,0.0d0,2.5d0,0.0d0,2.5d0,      &
      0.0d0,3.5d0,0.0d0,0.0d0,0.0d0,3.5d0,0.0d0,0.0d0,0.5d0,0.0d0,      &
      0.0d0,0.5d0,0.0d0,2.5d0,0.0d0,0.0d0,3.5d0,7.0d0,0.0d0,0.0d0,      &
      3.5d0,0.0d0,4.5d0,0.0d0,3.5d0,0.0d0,0.0d0,0.5d0,0.0d0,0.0d0,      &
      2.5d0,2.5d0,0.0d0,0.0d0,0.5d0,0.0d0,1.5d0,0.0d0,0.0d0,1.5d0,      &
      1.5d0,0.0d0,0.0d0,0.0d0,0.5d0,0.0d0,0.0d0,1.5d0,0.0d0,0.0d0,      &
      0.5d0,0.0d0,1.5d0,0.0d0,0.0d0,0.5d0,0.5d0,0.0d0,0.0d0,0.5d0,      &
      0.0d0,4.5d0,0.0d0,0.0d0,3.5d0,0.0d0 /)
  ! Elemental mass fraction for an isotope
  REAL*8, PARAMETER :: ISOTOPE_F(NISOTOPES) =                           &
    (/0.9997d0,   0.0002997d0,1.032d-6,   1.000d0,    0.06578d0,  0.9342d0,   1.000d0,    0.1834d0,    &
      0.8166d0,   0.9880d0,   0.01202d0,  0.9961d0,   0.003920d0, 0.9973d0,   0.0004037d0,0.002250d0,  &
      1.000d0,    0.8964d0,   0.002809d0, 0.1008d0,   1.000d0,    0.7795d0,   0.1028d0,   0.1177d0,    &
      1.000d0,    0.9187d0,   0.04832d0,  0.03295d0,  1.000d0,    0.9475d0,   0.007712d0, 0.04460d0,   &
      0.0002243d0,0.7474d0,   0.2526d0,   0.003030d0, 0.0006006d0,0.9964d0,   0.9294d0,   0.0001196d0, &
      0.07051d0,  0.9666d0,   0.006773d0, 0.001447d0, 0.02292d0,  4.586d-5,   0.002237d0, 1.000d0,     &
      0.07920d0,  0.07298d0,  0.7385d0,   0.05532d0,  0.05405d0,  0.002451d0, 0.9975d0,   0.04174d0,   &
      0.8370d0,   0.09674d0,  0.02453d0,  1.000d0,    0.05646d0,  0.9190d0,   0.02160d0,  0.002925d0,  &
      1.000d0,    0.6720d0,   0.2678d0,   0.01183d0,  0.03834d0,  0.01009d0,  0.6850d0,   0.3150d0,    &
      0.4754d0,   0.2813d0,   0.04196d0,  0.1948d0,   0.006629d0, 0.5942d0,   0.4058d0,   0.1961d0,    &
      0.2704d0,   0.07790d0,  0.3738d0,   0.08184d0,  1.000d0,    0.008332d0, 0.09009d0,  0.07433d0,   &
      0.2346d0,   0.5021d0,   0.09057d0,  0.5007d0,   0.4993d0,   0.003254d0, 0.02174d0,  0.1132d0,    &
      0.1137d0,   0.5708d0,   0.1774d0,   0.7170d0,   0.2830d0,   0.005363d0, 0.09668d0,  0.06943d0,   &
      0.8285d0,   1.000d0,    0.5071d0,   0.1118d0,   0.1728d0,   0.1789d0,   0.02944d0,  1.000d0,     &
      0.1422d0,   0.09055d0,  0.1575d0,   0.1668d0,   0.09647d0,  0.2463d0,   0.1003d0,   0.05257d0,   &
      0.01812d0,  0.1249d0,   0.1246d0,   0.1703d0,   0.3181d0,   0.1914d0,   1.000d0,    0.009768d0,  &
      0.1088d0,   0.2201d0,   0.2720d0,   0.2683d0,   0.1210d0,   0.5138d0,   0.4862d0,   0.01178d0,   &
      0.008543d0, 0.1221d0,   0.1263d0,   0.2402d0,   0.1227d0,   0.2911d0,   0.07723d0,  0.04218d0,   &
      0.9578d0,   0.009144d0, 0.006333d0, 0.003291d0, 0.1420d0,   0.07563d0,  0.2406d0,   0.08604d0,   &
      0.3291d0,   0.04755d0,  0.06043d0,  0.5681d0,   0.4319d0,   0.0008457d0,0.02436d0,  0.008572d0,  &
      0.04603d0,  0.06920d0,  0.1859d0,   0.3181d0,   0.3470d0,   1.000d0,    0.0008966d0,0.0008535d0, &
      0.01861d0,  0.2592d0,   0.04028d0,  0.2117d0,   0.2703d0,   0.1064d0,   0.09168d0,  1.000d0,     &
      0.001003d0, 0.0009701d0,0.02357d0,  0.06476d0,  0.07773d0,  0.1120d0,   0.7200d0,   0.0008935d0, &
      0.9991d0,   0.001794d0, 0.002470d0, 0.8832d0,   0.1126d0,   1.000d0,    0.2676d0,   0.1209d0,    &
      0.2375d0,   0.08339d0,  0.1740d0,   0.05845d0,  0.05821d0,  0.02938d0,  0.1465d0,   0.1106d0,    &
      0.1369d0,   0.07358d0,  0.2703d0,   0.2329d0,   0.4748d0,   0.5252d0,   0.001932d0, 0.02134d0,   &
      0.1458d0,   0.2030d0,   0.1562d0,   0.2495d0,   0.2223d0,   1.000d0,    0.0005757d0,0.0009719d0, &
      0.02303d0,  0.1873d0,   0.2542d0,   0.2497d0,   0.2843d0,   1.000d0,    0.001346d0, 0.01569d0,   &
      0.3324d0,   0.2282d0,   0.2709d0,   0.1515d0,   1.000d0,    0.001262d0, 0.02985d0,  0.1411d0,    &
      0.2169d0,   0.1612d0,   0.3200d0,   0.1297d0,   0.9740d0,   0.02604d0,  0.001559d0, 0.05185d0,   &
      0.1844d0,   0.2720d0,   0.1366d0,   0.3537d0,   1.000d0,    0.001175d0, 0.2623d0,   0.1424d0,    &
      0.3066d0,   0.2876d0,   0.3715d0,   0.6285d0,   0.0001934d0,0.01554d0,  0.01572d0,  0.1313d0,    &
      0.1610d0,   0.2632d0,   0.4130d0,   0.3706d0,   0.6294d0,   0.0001363d0,0.007695d0, 0.3278d0,    &
      0.3381d0,   0.2536d0,   0.07269d0,  1.000d0,    0.001466d0, 0.09840d0,  0.1673d0,   0.2303d0,    &
      0.1321d0,   0.3007d0,   0.06976d0,  0.2932d0,   0.7068d0,   0.01378d0,  0.2396d0,   0.2207d0,    &
      0.5259d0,   1.000d0,    1.000d0,    5.310d-5,   0.007114d0, 0.9928d0 /)

  IF (Z .LE. NELEMENTS) THEN
    Niso = ELEMENT_NISO(Z)
    ALLOCATE(Ziso(Niso),Aiso(Niso),fiso(Niso),Miso(Niso))
    I1 = ELEMENT_INDEX(Z)
    I2 = I1 + Niso - 1
    Ziso = ISOTOPE_Z(I1:I2)
    Aiso = ISOTOPE_A(I1:I2)
    !Jiso = ISOTOPE_J(I1:I2)
    fiso = ISOTOPE_F(I1:I2)
    Miso = IsotopeMass(Ziso,Aiso)
  ELSE
    Niso = 0
    ! Zero-length arrays (nothing to fill in)
    ALLOCATE(Ziso(Niso),Aiso(Niso),fiso(Niso),Miso(Niso))
  END IF
  
END SUBROUTINE


! ----------------------------------------------------------------------
! For the given compound, specified by a list of element atomic numbers
! and the stoichiometry, fills in allocatable arrays containing isotopic
! atomic numbers (Z), atomic masses (A), mass fractions (f), and
! masses (M).
! 
! Input arguments:
!     N          Number of compound elements.
!     Z          Atomic number of elements, array of size [1:N].
!     stoich     Stoichiometry of the compound elements, array of size
!                [1:N].  For example, CF3Cl would have Z={6,9,17} and
!                stoich={1,3,1}.
! Output arguments:
!     Niso       Number of isotopes
!     Ziso       Allocatable array (integer) of isotopes' atomic numbers (Z)
!     Aiso       Allocatable array (integer) of isotopes' atomic masses (A)
!     fiso       Allocatable array of isotopes' mass fractions
!     Miso       Allocatable array of isotopes' nuclear masses [GeV]
! 
SUBROUTINE CompoundIsotopeList(N,Z,stoich,Niso,Ziso,Aiso,fiso,Miso)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(IN) :: Z(N),stoich(N)
  INTEGER, INTENT(OUT) :: Niso
  INTEGER, ALLOCATABLE, INTENT(OUT) :: Ziso(:),Aiso(:)
  REAL*8, ALLOCATABLE, INTENT(OUT) :: fiso(:),Miso(:)
  INTEGER :: tempNiso,K,I1,I2
  REAL*8 :: weight(N)
  INTEGER, ALLOCATABLE :: tempZ(:),tempA(:)
  REAL*8, ALLOCATABLE :: tempf(:),tempM(:)
  
  ! Get number of isotopes.
  Niso = 0
  DO K = 1,N
    CALL ElementIsotopeList(Z(K),tempNiso,tempZ,tempA,tempf,tempM)
    Niso = Niso + tempNiso
    weight(K) = stoich(K) * SUM(tempf*tempM)
  END DO
  
  ! Relative weight (by mass) of each element.
  IF (N .GT. 0) THEN
    weight = weight / SUM(weight)
  ELSE
    weight = 1d0
  END IF
  
  ! Allocate and fill in total arrays.
  ! Reweight isotopes' mass fractions by element's mass fraction.
  I1 = 1
  ALLOCATE(Ziso(Niso),Aiso(Niso),fiso(Niso),Miso(Niso))
  DO K = 1,N
    CALL ElementIsotopeList(Z(K),tempNiso,tempZ,tempA,tempf,tempM)
    I2 = I1 + tempNiso - 1
    Ziso(I1:I2) = tempZ
    Aiso(I1:I2) = tempA
    fiso(I1:I2) = weight(K)*tempf
    Miso(I1:I2) = tempM
    I1 = I2 + 1
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Determines the mass of the isotope using the semi-empirical mass
! model.  A few low mass (Z <= 2) isotopes are given explicitly to
! avoid accuracy issues in the mass model in those cases.
! 
! Input arguments:
!     Z          Atomic number
!     A          Mass number
! Returns in [GeV].
! 
ELEMENTAL FUNCTION IsotopeMass(Z,A) RESULT(M)
  IMPLICIT NONE
  REAL*8 :: M
  INTEGER, INTENT(IN) :: Z,A
  REAL*8 :: A0,E
  ! 
  ! The semi-empirical mass formula will be used to calculate nuclear
  ! masses, but a few low-mass cases are given explicitly. [GeV]
  REAL*8, PARAMETER :: H2_MASS  = 1.8756129d0  ! Deuteron
  REAL*8, PARAMETER :: H3_MASS  = 2.8089210d0  ! Tritium ("triton")
  REAL*8, PARAMETER :: HE3_MASS = 2.8083915d0  ! Helium-3 ("helion")
  REAL*8, PARAMETER :: HE4_MASS = 3.7273792d0  ! Alpha
  ! 
  ! Semi-empirical mass formula constants
  !   m = Z*Mp + (A-Z)*Mn                 proton+neutron contributions
  !       - aV*A + aS*A^(2/3)             volume & surface terms
  !       + aC*Z^2/A^(1/3)                Coulomb term
  !       + aA*(A-2Z)^2/A                 Pauli (asymmetry) term
  !       + eps_{Z,A} * aP/A^(1/2)        pairing term
  ! where eps_{Z,A} is 0 if A is odd and (-1)^(Z+1) otherwise.
  ! The values below [GeV] are taken from Rohlf (1994).
  REAL*8, PARAMETER :: SEMF_AV = 0.01575d0
  REAL*8, PARAMETER :: SEMF_AS = 0.0178d0
  REAL*8, PARAMETER :: SEMF_AC = 0.000711d0
  REAL*8, PARAMETER :: SEMF_AA = 0.0237d0
  REAL*8, PARAMETER :: SEMF_AP = 0.01118d0

  ! Bad cases
  IF ((A .LE. 0) .OR. (Z .LT. 0)) THEN
    M = 0d0
    RETURN
  END IF
  
  ! Special cases: all natural and/or commonly used isotopes
  ! with Z <= 2 (i.e. nucleons, hydrogen, helium)
  IF (A .EQ. 1) THEN
    IF (Z .EQ. 0) THEN
      M = NEUTRON_MASS
      RETURN
    ELSE IF (Z .EQ. 1) THEN
      M = PROTON_MASS
      RETURN
    END IF
  ELSE IF ((A .EQ. 2) .AND. (Z .EQ. 1)) THEN
    M = H2_MASS
    RETURN
  ELSE IF (A .EQ. 3) THEN
    IF (Z .EQ. 1) THEN
      M = H3_MASS
      RETURN
    ELSE IF (Z .EQ. 2) THEN
      M = HE3_MASS
      RETURN
    END IF
  ELSE IF ((A .EQ. 4) .AND. (Z .EQ. 2)) THEN
    M = HE4_MASS
    RETURN
  END IF
  
  ! Generic semi-empirical mass formula calculation.
  ! Here, E is binding energy.
  A0 = A      ! type conversion
  E = SEMF_AV*A - SEMF_AS*A0**(2d0/3d0) - SEMF_AC*Z**2/A0**(1d0/3d0)    &
      - SEMF_AA*(A-2*Z)**2/A0
  
  ! Pairing term
  IF (MOD(A,2) .EQ. 0) THEN
    IF (MOD(Z,2) .EQ. 0) THEN
      E = E + SEMF_AP/A0**(1d0/2d0)
    ELSE
      E = E - SEMF_AP/A0**(1d0/2d0)
    END IF
  END IF
  
  M = Z*PROTON_MASS + (A-Z)*NEUTRON_MASS - E
  
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the spin-independent weighted form factor Wsi for the given
! isotope at the given momentum transfer q, defined as:
!   Wsi(+1,:) = (1/pi) Z^2 F^2(:)        ! SI proton
!   Wsi( 0,:) = (1/pi) 2*Z*(A-Z) F^2(:)  ! SI crossterm
!   Wsi(-1,:) = (1/pi) (A-Z)^2 F^2(:)    ! SI neutron
! where F^2 is the standard form factor (':' represents momentum array
! q).  Uses the Helm form factor; see:
!   Lewin & Smith, Astropart. Phys. 6, 87 (1996)  [Eqn 4.7]
! 
! Required input arguments:
!     Z,A        The atomic number and mass number of the isotope.
!     N          Number of momentum q values.
!     q          Array of size [1:N] containing momentum values [GeV].
! Required output argument:
!     W          Array of size [-1,1,1:N] to be filled with weighted
!                form factor values [unitless].
! 
PURE SUBROUTINE CalcWSI(Z,A,N,q,W)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Z,A,N
  REAL*8, INTENT(IN) :: q(1:N)
  REAL*8, INTENT(OUT) :: W(-1:1,1:N)
  INTEGER :: K
  REAL*8 :: c,rn,qrn,qs,F2,weights(-1:1)
  REAL*8, PARAMETER :: Arn = 0.52d0
  REAL*8, PARAMETER :: S   = 0.9d0
  REAL*8, PARAMETER :: C1  = 1.23d0
  REAL*8, PARAMETER :: C2  = -0.60d0
  
  weights(+1) = Z**2      / PI
  weights( 0) = 2*Z*(A-Z) / PI
  weights(-1) = (A-Z)**2  / PI
  
  c  = C1*A**(1d0/3d0) + C2
  rn = SQRT(c**2 + (7d0/3d0)*PI**2*Arn**2 - 5*S**2)
  
  ! Helm FF:  F^2(q) = 9 [j1(q rn)/(q rn)]^2 exp(-q s)
  DO K = 1, N
    qrn = q(K)*rn / HBARC
    qs  = q(K)*S / HBARC
    ! avoid numerical issues for small q by using a Taylor expansion
    IF (qrn .LE. 0.01d0) THEN
      F2 = (1 - qrn**2*((1d0/5d0) - qrn**2*(3d0/175d0))) * EXP(-qs**2)
    ELSE
      F2 = 9 * (SIN(qrn) - qrn*COS(qrn))**2 / qrn**6 * EXP(-qs**2)
    END IF
    W(:,K) = F2 * weights
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Calculates the spin-dependent weighted form factor Wsd for the given
! isotope at the given momentum transfer q, defined as:
!   Wsd(+1,:) = 4/(2J+1) Spp(q)          ! SD proton
!   Wsd( 0,:) = 4/(2J+1) Spn(q)          ! SD crossterm
!   Wsd(-1,:) = 4/(2J+1) Snn(q)          ! SD neutron
! where J is the nuclear spin and Spp/Spn/Snn are the spin structure
! functions.  For a comprehensive review of available spin structure
! functions, see:
!   Bednyakov & Simkovic, Phys. Part. Nucl. 37, S106 (2006)
!     [hep-ph/0608097]
! Note that, in the above review, the quantities "ap" & "an" are
! actually the quantities Gp & Gn as used here, not the quantities
! ap & an as used here or in much of the other direct detection
! literature.
! 
! Xenon form factors implemented by Andre Scaffidi.
! 
! NOTE: ONLY A LIMITED SELECTION OF SD FORM FACTORS ARE CURRENTLY
!       IMPLEMENTED.  THEY ARE SIMPLY SET TO ZERO.
! 
! Required input arguments:
!     Z,A        The atomic number and mass number of the isotope.
!     N          Number of momentum q values.
!     q          Array of size [1:N] containing momentum values [GeV].
! Required output argument:
!     W          Array of size [-1,1,1:N] to be filled with weighted
!                form factor values [unitless].
! 
PURE SUBROUTINE CalcWSD(Z,A,N,q,W)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Z,A,N
  REAL*8, INTENT(IN) :: q(1:N)
  REAL*8, INTENT(OUT) :: W(-1:1,1:N)
  REAL*8 :: J,b,umax,Spp(1:N),Spn(1:N),Snn(1:N)
  
  ! Initialization
  ! Most isotopes have zero spin (=> zero form factor)
  J   = 0d0
  Spp = 0d0
  Snn = 0d0
  Spn = 0d0
  W   = 0d0
  
  ! Fluorine ----------------------------------
  IF (Z .EQ. 9) THEN
    ! Fluorine 19 --------------
    ! From Klos et al. [1304.7684]; see Table 6.
    ! The S11 and S01 (1b+2b) band means are used.
    ! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    ! where u = q^2 b^2 / 2.
    IF (A .EQ. 19) THEN
      ! Harmonic oscillator length [GeV^-1]
      ! NOTE: The length is listed incorrectly in 1304.7684 Table 6.
      ! Correct number found in caption of Figure 12.
      b = 1.7608d0 / HBARC
      ! Nuclear spin
      J = 0.5d0
      ! Fits over 0 < u < 5 (probably).  Limit to this range.
      umax = 5d0
      ! Helper routine for exponential-polynomial Sij forms.
      CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,14,   &
               (/  0.108058d0, -0.143789d0,  0.0680848d0, 4.07415d-4,   &
                  -0.0314817d0, 0.0385933d0,-0.0293716d0, 0.0152264d0,  &
                  -5.52655d-3,  1.41965d-3, -2.56989d-4,  3.20688d-5,   &
                  -2.62562d-6,  1.26950d-7, -2.74719d-9  /),            &
               (/  0.0660281d0,-0.137668d0,  0.161957d0, -0.166004d0,   &
                   0.152569d0, -0.111464d0,  0.0609363d0,-0.0246265d0,  &
                   7.35189d-3, -1.61496d-3,  2.57660d-4, -2.90407d-5,   &
                   2.19205d-6, -9.94286d-8,  2.04873d-9  /),            &
               (/  0.167759d0, -0.286581d0,  0.244497d0, -0.176999d0,   &
                   0.134461d0, -0.0921638d0, 0.0494464d0,-0.0198425d0,  &
                   5.87500d-3, -1.26970d-3,  1.96948d-4, -2.12790d-5,   &
                   1.51602d-6, -6.38405d-8,  1.20003d-9  /) )
      
    ! --------------------------
    ! Other isotopes have zero spin
    ELSE
      W = 0d0
      RETURN
    END IF
  
  ! Sodium ------------------------------------
  ELSE IF (Z .EQ. 11) THEN
    ! Sodium 23 ----------------
    ! From Klos et al. [1304.7684]; see Table 7.
    ! The S11 and S01 (1b+2b) band means are used.
    ! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    ! where u = q^2 b^2 / 2.
    IF (A .EQ. 23) THEN
      ! Harmonic oscillator length [GeV^-1]
      b = 1.8032d0 / HBARC
      ! Nuclear spin
      J = 1.5d0
      ! Fits over 0 < u < 5 (probably).  Limit to this range.
      umax = 5d0
      ! Helper routine for exponential-polynomial Sij forms.
      CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
               (/  0.0325305d0,-0.0433531d0, 0.0319487d0,-5.68858d-3,   &
                   2.67783d-4,  2.44643d-5, -4.79620d-6,  5.39846d-7,   &
                  -3.24691d-8,  8.09358d-10 /),                         &
               (/  0.0127243d0,-0.0248722d0, 0.0275805d0,-0.0135587d0,  &
                   3.93910d-3, -7.28827d-4,  8.82003d-5, -6.80402d-6,   &
                   3.02566d-7, -5.81500d-9  /),                         &
               (/  0.0404609d0,-0.0677623d0, 0.0660458d0,-0.0269059d0,  &
                   7.22875d-3, -1.45367d-3,  2.08818d-4, -1.89848d-5,   &
                   9.43875d-7, -1.93865d-8  /) )
      
    ! --------------------------
    ! Other isotopes have zero spin
    ELSE
      W = 0d0
      RETURN
    END IF
  
  ! Aluminum ----------------------------------
  ELSE IF (Z .EQ. 13) THEN
    ! Aluminum 27 --------------
    ! From Klos et al. [1304.7684]; see Table 7.
    ! The S11 and S01 (1b+2b) band means are used.
    ! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    ! where u = q^2 b^2 / 2.
    IF (A .EQ. 27) THEN
      ! Harmonic oscillator length [GeV^-1]
      b = 1.8405d0 / HBARC
      ! Nuclear spin
      J = 2.5d0
      ! Fits over 0 < u < 5 (probably).  Limit to this range.
      umax = 5d0
      ! Helper routine for exponential-polynomial Sij forms.
      CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
               (/  0.0888149d0,-0.117822d0,  0.0631336d0,-9.19554d-3,   &
                   5.84421d-4,  5.54484d-4, -1.15453d-4,  1.40388d-5,   &
                  -9.21830d-7,  2.52336d-8  /),                         &
               (/  0.0334384d0,-0.0710220d0, 0.0805917d0,-0.0514533d0,  &
                   0.0221406d0,-6.14292d-3,  1.08899d-3, -1.17175d-4,   &
                   6.93171d-6, -1.71376d-7  /),                         &
               (/  0.108031d0, -0.182354d0,  0.149131d0, -0.0623283d0,  &
                   0.0196187d0,-4.07283d-3,  6.05456d-4, -5.96107d-5,   &
                   3.28620d-6, -7.69545d-8  /) )
      
    ! --------------------------
    ! Other isotopes have zero spin
    ELSE
      W = 0d0
      RETURN
    END IF
  
  ! Silicon -----------------------------------
  ELSE IF (Z .EQ. 14) THEN
    ! Silicon 29 ---------------
    ! From Klos et al. [1304.7684]; see Table 8.
    ! The S11 and S01 (1b+2b) band means are used.
    ! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    ! where u = q^2 b^2 / 2.
    IF (A .EQ. 29) THEN
      ! Harmonic oscillator length [GeV^-1]
      b = 1.8575d0 / HBARC
      ! Nuclear spin
      J = 0.5d0
      ! Fits over 0 < u < 5 (probably).  Limit to this range.
      umax = 5d0
      ! Helper routine for exponential-polynomial Sij forms.
      CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
               (/  0.0140647d0,-0.0188522d0, 0.0149891d0,-0.00542122d0, &
                   1.17173d-3, -1.15932d-4,  2.47182d-5, -3.04480d-6,   &
                   2.00549d-7, -5.46011d-9  /),                         &
               (/  5.63416d-3, -0.0121901d0, 0.0156006d0,-0.0110712d0,  &
                   4.85653d-3, -1.28086d-3,  2.08618d-4, -2.04711d-5,   &
                   1.10510d-6, -2.47894d-8  /),                         &
               (/ -0.0176295d0, 0.0301066d0,-0.0308628d0, 0.0175759d0,  &
                  -6.78176d-3,  1.69859d-3, -3.00595d-4,  3.39713d-5,   &
                  -2.08093d-6,  5.28653d-8  /) )
      
    ! --------------------------
    ! Other isotopes have zero spin
    ELSE
      W = 0d0
      RETURN
    END IF
  
  ! Germanium ---------------------------------
  ELSE IF (Z .EQ. 32) THEN
    ! Germanium 73 -------------
    ! From Klos et al. [1304.7684]; see Table 5.
    ! The S11 and S01 (1b+2b) band means are used.
    ! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    ! where u = q^2 b^2 / 2.
    IF (A .EQ. 73) THEN
      ! Harmonic oscillator length [GeV^-1]
      b = 2.1058d0 / HBARC
      ! Nuclear spin
      J = 4.5d0
      ! Fits over 0 < u < 5 (probably).  Limit to this range.
      umax = 5d0
      ! Helper routine for exponential-polynomial Sij forms.
      CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
               (/  0.215608d0, -0.578786d0,  0.698020d0, -0.372000d0,   &
                   0.107576d0, -0.0182408d0, 2.17108d-3, -2.07981d-4,   &
                   1.65907d-5, -5.95664d-7  /),                         &
               (/  0.0972089d0,-0.308986d0,  0.450727d0, -0.337355d0,   &
                   0.154809d0, -0.0469625d0, 9.71560d-3, -1.33058d-3,   &
                   1.09084d-4, -4.02514d-6  /),                         &
               (/ -0.287562d0,  0.844765d0, -1.133659d0,  0.745494d0,   &
                  -0.296646d0,  0.0788570d0,-0.0147852d0, 1.87401d-3,   &
                  -1.42195d-4,  4.72898d-6  /) )
      
    ! --------------------------
    ! Other isotopes have zero spin
    ELSE
      W = 0d0
      RETURN
    END IF
  
  ! Iodine ------------------------------------
  ELSE IF (Z .EQ. 53) THEN
    ! Iodine 127 ---------------
    ! From Klos et al. [1304.7684]; see Table 5.
    ! The S11 and S01 (1b+2b) band means are used.
    ! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    ! where u = q^2 b^2 / 2.
    IF (A .EQ. 127) THEN
      ! Harmonic oscillator length [GeV^-1]
      b = 2.2801d0 / HBARC
      ! Nuclear spin
      J = 2.5d0
      ! Fits over 0 < u < 5 (probably).  Limit to this range.
      umax = 5d0
      ! Helper routine for exponential-polynomial Sij forms.
      CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
               (/  0.0928480d0,-0.252496d0,  0.351982d0, -0.260427d0,   &
                   0.118280d0, -0.0319614d0, 4.92618d-3, -4.06546d-4,   &
                   1.55818d-5, -1.64934d-7  /),                         &
               (/  0.0389166d0,-0.119307d0,  0.189835d0, -0.168819d0,   &
                   0.0952229d0,-0.0343338d0, 7.86014d-3, -1.11341d-3,   &
                   8.98377d-5, -3.17792d-6  /),                         &
               (/  0.119382d0, -0.345408d0,  0.515816d0, -0.421111d0,   &
                   0.215622d0, -0.0691557d0, 0.0137850d0,-1.68267d-3,   &
                   1.18375d-4, -3.78243d-6  /) )
      
    ! --------------------------
    ! Other isotopes have zero spin
    ELSE
      W = 0d0
      RETURN
    END IF
  
  ! Xenon -------------------------------------
  ELSE IF (Z .EQ. 54) THEN
    ! Xenon 129 ----------------
    ! From Klos et al. [1304.7684]; see Table 1.
    ! The S11 and S01 (1b+2b) band means are used.
    ! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    ! where u = q^2 b^2 / 2.
    IF (A .EQ. 129) THEN
      ! Harmonic oscillator length [GeV^-1]
      b = 2.2853d0 / HBARC
      ! Nuclear spin
      J = 0.5d0
      ! Fits over 0 < u < 10 (probably).  Limit to this range.
      umax = 10d0
      ! Helper routine for exponential-polynomial Sij forms.
      CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
               (/  0.0547144d0,-0.146407d0,  0.180603d0, -0.125526d0,   &
                   0.0521484d0,-0.0126363d0, 1.76284d-3, -1.32501d-4,   &
                   4.23423d-6, -1.68052d-9  /),                         &
               (/  0.0289650d0,-0.0867525d0, 0.115723d0, -0.0858610d0,  &
                   0.0384596d0,-0.0105918d0, 1.80025d-3, -1.83841d-4,   &
                   1.03293d-5, -2.44338d-7  /),                         &
               (/ -0.0791167d0, 0.225715d0, -0.293581d0,  0.215439d0,   &
                  -0.0959137d0, 0.0260514d0,-4.33883d-3,  4.32823d-4,   &
                  -2.37266d-5,  2.994545d0  /) )
      
    ! Xenon 131 ----------------
    ! From Klos et al. [1304.7684]; see Table 1.
    ! The S11 and S01 (1b+2b) band means are used.
    ! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    ! where u = q^2 b^2 / 2.
    ELSE IF (A .EQ. 131) THEN
      ! Harmonic oscillator length [GeV^-1]
      b = 2.2905d0 / HBARC
      ! Nuclear spin
      J = 1.5d0
      ! Fits over 0 < u < 10 (probably).  Limit to this range.
      umax = 10d0
      ! Helper routine for exponential-polynomial Sij forms.
      CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
               (/  0.0417857d0,-0.111132d0,  0.171306d0, -0.132481d0,   &
                   0.0630161d0,-0.0177684d0, 2.82192d-3, -2.32247d-4,   &
                   7.81471d-6,  1.25984d-9  /),                         &
               (/  0.0219206d0,-0.0642919d0, 0.0957262d0,-0.0727452d0,  &
                   0.0338802d0,-9.86454d-3,  1.75888d-3, -1.83629d-4,   &
                   1.01679d-5, -2.25247d-7  /),                         &
               (/ -0.0602463d0, 0.171349d0, -0.265846d0,  0.211589d0,   &
                  -0.103611d0,  0.0311471d0,-5.60929d-3,  5.81416d-4,   &
                  -3.16217d-5,  6.82201d-7  /) )
      
    !! Xenon 129 ----------------
    !! From Menendez et al. [1208.1094]; see Table 1.
    !! The 1b+2b results are taken where possible.
    !! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    !! where u = q^2 b^2 / 2.
    !IF (A .EQ. 129) THEN
    !  ! Harmonic oscillator length [GeV^-1]
    !  b = 2.2853d0 / HBARC
    !  ! Nuclear spin
    !  J = 0.5d0
    !  ! Fits over 0 < u < 3 (probably).  Limit to this range.
    !  umax = 3d0
    !  ! Helper routine for exponential-polynomial Sij forms.
    !  CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
    !           (/  0.054731d0, -0.146897d0,  0.182479d0,  -0.128112d0,  &
    !               0.0539978d0,-0.0133335d0, 0.00190579d0,-1.48373d-4,  &
    !               5.11732d-6, -2.06597d-8 /),                          &
    !           (/  0.02933d0,  -0.0905396d0, 0.122783d0,  -0.0912046d0, &
    !               0.0401076d0,-0.010598d0,  0.00168737d0,-1.56768d-4,  &
    !               7.69202d-6, -1.48874d-7 /),                          &
    !           (/ -0.0796645d0, 0.231997d0, -0.304198d0,   0.222024d0,  &
    !              -0.096693d0,  0.0251835d0,-0.00392356d0, 3.53343d-4,  &
    !              -1.65058d-5,  2.88576d-7 /) )
      
    !! Xenon 131 ----------------
    !! From Menendez et al. [1208.1094]; see Table 1.
    !! The 1b+2b results are taken where possible.
    !! Functional form is Sij(u) = e^{-u} \sum_{k=0} C_{ij,k} u^k
    !! where u = q^2 b^2 / 2.
    !ELSE IF (A .EQ. 131) THEN
    !  ! Harmonic oscillator length [GeV^-1]
    !  b = 2.2905d0 / HBARC
    !  ! Nuclear spin
    !  J=1.5d0
    !  ! Fits over 0 < u < 3 (probably).  Limit to this range.
    !  umax = 3d0
    !  ! Helper routine for exponential-polynomial Sij forms.
    !  CALL ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,1d0,1d0,1d0,9,    &
    !           (/  0.0417889d0,-0.111171d0,  0.171966d0,  -0.133219d0,  &
    !               0.0633805d0,-0.0178388d0, 0.00282476d0,-2.31681d-4,  &
    !               7.78223d-6, -4.49287d-10 /),                         &
    !           (/  0.022446d0, -0.0733931d0, 0.110509d0,  -0.0868752d0, &
    !               0.0405399d0,-0.0113544d0, 0.00187572d0,-1.75285d-4,  &
    !               8.40043d-6, -1.53632d-7  /),                         &
    !           (/ -0.0608808d0, 0.181473d0, -0.272533d0,   0.211776d0,  &
    !              -0.0985956d0, 0.027438d0, -0.0044424d0,  3.97619d-4,  &
    !              -1.74758d-5,  2.55979d-7  /) )
      
    ! --------------------------
    ! Other isotopes have zero spin
    ELSE
      W = 0d0
      RETURN
    END IF
  
  ! Zero spin or unimplemented isotope ---------
  ELSE
    W = 0d0
    RETURN
  ENDIF
  
  ! Weighted form factors
  W(+1,:) = (4d0 / (2d0*J + 1d0)) * Spp    ! SD proton
  W( 0,:) = (4d0 / (2d0*J + 1d0)) * Spn    ! SD crossterm
  W(-1,:) = (4d0 / (2d0*J + 1d0)) * Snn    ! SD neutron
  
  
  CONTAINS
  
  ! --------------------------------------------
  ! Calculates form factors of the form:
  !   Sij(u) = e^{-u} A_ij \sum_{k=0} C_{ij,k} u^k
  ! where u = q^2 b^2 / 2 and A = (A00,A11,A01).  Sets to zero
  ! when u > umax.  Fills the arrays S00, S11, & S01 in the
  ! parent routine.
  PURE SUBROUTINE ExponentialPolynomial(N,Spp,Spn,Snn,b,umax,           &
                                        A00,A11,A01,NC,C00,C11,C01)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N,NC
    REAL*8, INTENT(IN) :: b,umax,A00,A11,A01,C00(0:NC),C11(0:NC),C01(0:NC)
    REAL*8, INTENT(OUT) :: Spp(1:N),Spn(1:N),Snn(1:N)
    INTEGER :: I,K
    REAL*8 :: u,uk,expu,S00,S11,S01
    DO I=1,N
      u  = 0.5d0 * (q(I)*b)**2
      uk = 1
      S00 = 0d0
      S11 = 0d0
      S01 = 0d0
      ! Conservatively set to zero outside range of validity.
      IF (u .LE. umax) THEN
        DO K=0,NC
          S00 = S00 + uk * C00(K)
          S11 = S11 + uk * C11(K)
          S01 = S01 + uk * C01(K)
          uk = u*uk
        END DO
        expu = EXP(-u)
        S00 = expu * A00 * S00
        S11 = expu * A11 * S11
        S01 = expu * A01 * S01
      END IF
      ! Basis transformation: a0=ap+an, a1=ap-an
      Spp(I) = S00 + S11 + S01
      Snn(I) = S00 + S11 - S01
      Spn(I) = 2 * (S00 - S11)
    END DO
    
  END SUBROUTINE ExponentialPolynomial
  
  
END SUBROUTINE


! ----------------------------------------------------------------------
! INTERFACE NAME: EToQ
! Calculates the momentum transfer q [GeV] corresponding to a nuclear
! recoil energy E [keV] for the given isotope mass m [GeV].
! 
! This is the scalar version (single mass and energy).
! 
! Input arguments:
!   E           Recoil energy [keV]
!   Miso        Isotope mass [GeV]
! 
ELEMENTAL FUNCTION EToQ(E,Miso) RESULT(q)
  IMPLICIT NONE
  REAL*8 :: Q
  REAL*8, INTENT(IN) :: E,Miso
  ! Factor of 1d-6 to take keV -> GeV
  q = SQRT(2*Miso*(1d-6*E))
END FUNCTION



!=======================================================================
! DETECTOR EFFICIENCY/RESPONSE ROUTINES
!=======================================================================

! ----------------------------------------------------------------------
! Sets the efficiencies to those found in the given file, presumably
! an efficiency or efficiencies file generated by TPCMC.
! 
! Required input argument:
!     file       The file to load efficiencies from.
! Required output arguments:
!     NE         Number of recoil energy tabulation points loaded
!     E          Allocatable array to contain energies [keV].
!                Allocated to size [1:NE].
!     Neff       Number of interval/bin efficiencies loaded.
!     eff        Allocatable array to contain efficiency curves for
!                the total range and each interval/bin.  Allocated to
!                size [1:NE,0:Neff].
! 
SUBROUTINE LoadEfficiencyFile(file,NE,E,Neff,eff)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: file
  INTEGER, INTENT(OUT) :: NE,Neff
  REAL*8, ALLOCATABLE, INTENT(OUT) :: E(:),eff(:,:)
  LOGICAL :: status
  INTEGER :: Kcol,Keff,Nrow,Ncol,Nvalid
  REAL*8, ALLOCATABLE :: data(:,:)
  
  ! Load table from file
  CALL LoadTable(file=file,Nrow=Nrow,Ncol=Ncol,data=data,status=status)
  IF ((.NOT. status) .OR. (Ncol .LT. 2)) THEN
    WRITE(0,*) 'ERROR: Failed to load data from file ' // TRIM(file) // '.'
    STOP
  END IF
  
  ! Energies
  NE = Nrow
  ALLOCATE(E(NE))
  E = data(:,1)
  
  ! Find number of valid efficiency columns
  Nvalid = 0
  DO Kcol = 2,Ncol
    IF (ALL(data(:,Kcol) .GE. 0d0) .AND. ALL(data(:,Kcol) .LE. 1.00001d0)) Nvalid = Nvalid + 1
  END DO
  IF (Nvalid .LE. 0) THEN
    WRITE(0,*) 'ERROR: Failed to find valid data in file ' // TRIM(file) // '.'
    STOP
  END IF
  
  ! Now get efficiencies
  Neff = Nvalid - 1
  ALLOCATE(eff(NE,0:Neff))
  Keff = 0
  DO Kcol = 2,Ncol
    IF (ALL(data(:,Kcol) .GE. 0d0) .AND. ALL(data(:,Kcol) .LE. 1.00001d0) &
        .AND. (Keff .LE. Neff)) THEN
      eff(:,Keff) = data(:,Kcol)
      Keff = Keff + 1
    END IF
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! FOR FUTURE USE
! Retrieves the efficiencies found in the given file, presumably an
! efficiency or efficiencies file generated by TPCMC.
! 
! Required input argument:
!     file       The file to load efficiencies from.
! Required output argument:
!     effS       DetectorEfficiencyStruct to contain tabulated
!                efficiency data.
! 
SUBROUTINE NEWLoadEfficiencyFile(file,effS)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: file
  TYPE(DetectorEfficiencyStruct), INTENT(OUT) :: effS
  LOGICAL :: status
  INTEGER :: Kcol,Keff,Nrow,Ncol,Nvalid
  REAL*8, ALLOCATABLE :: data(:,:)
  
  effS%file = file
  
  ! Load table from file
  CALL LoadTable(file=file,Nrow=Nrow,Ncol=Ncol,data=data,status=status)
  IF ((.NOT. status) .OR. (Ncol .LT. 2)) THEN
    WRITE(0,*) 'ERROR: Failed to load data from file ' // TRIM(file) // '.'
    STOP
  END IF
  
  ! Energies
  effS%NE = Nrow
  ALLOCATE(effS%E(effS%NE))
  effS%E = data(:,1)
  
  ! Find number of valid efficiency columns
  Nvalid = 0
  DO Kcol = 2,Ncol
    IF (ALL(data(:,Kcol) .GE. 0d0) .AND. ALL(data(:,Kcol) .LE. 1.00001d0)) Nvalid = Nvalid + 1
  END DO
  IF (Nvalid .LE. 0) THEN
    WRITE(0,*) 'ERROR: Failed to find valid data in file ' // TRIM(file) // '.'
    STOP
  END IF
  
  ! Now get efficiencies
  effS%Neff = Nvalid
  ALLOCATE(effS%eff(effS%NE,0:effS%Neff))
  Keff = 0
  DO Kcol = 2,Ncol
    IF (ALL(data(:,Kcol) .GE. 0d0) .AND. ALL(data(:,Kcol) .LE. 1.00001d0) &
        .AND. (Keff .LE. effS%Neff)) THEN
      effS%eff(:,Keff) = data(:,Kcol)
      Keff = Keff + 1
    END IF
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Retabulates the given efficiency data at the given array of energies.
! The efficiency is taken to be zero beyond the tabulation range of
! the provided efficiency data.
! 
! Input arguments:
!     NEeff      Number of recoil energy tabulation points for raw
!                efficiency data.
!     Eeff       Array of size [1:NEeff] containing tabulation energies
!                for raw efficiency data [keV].
!     Neff       Number of interval/bin efficiencies in raw efficiency
!                data (0 if total only).
!     eff        Array of size [1:NEeff,0:Neff] containing raw
!                efficiency data.  The second index is over the
!                interval/bin number (0 for total range).
!     NE         Number of recoil energies in desired tabulation.
!     E          Array of size [1:NE] containing desired tabulation
!                energies [keV].
!     intervals  LOGICAL indicating if calculations for intervals/bins
!                will be performed (necessary for max gap).  If .FALSE.,
!                efficiency data for intervals/bins will be dropped.
! Output argument:
!     eff        Allocatable array to contain efficiency curves for
!                the total range and each interval/bin at the energies
!                in E.  Allocated to size [1:NE,0:Neff] or [1:NE,0:0],
!                depending on the intervals flag.
! 
SUBROUTINE RetabulateEfficiency(NEeff,Eeff,Neff,eff,NE,E,intervals,eff0)
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: intervals
  INTEGER, INTENT(IN) :: NEeff,Neff,NE
  REAL*8, INTENT(IN) :: Eeff(1:NEeff),eff(1:NEeff,0:Neff),E(1:NE)
  REAL*8, ALLOCATABLE, INTENT(OUT) :: eff0(:,:)
  INTEGER :: Neff0,Keff,K0,K0start
  REAL*8 :: f
  
  ! Allocate output array
  IF (intervals) THEN
    Neff0 = Neff
  ELSE
    Neff0 = 0
  END IF
  ALLOCATE(eff0(1:NE,0:Neff0))
  
  ! Special cases
  IF (E(NE) .LT. Eeff(1)) THEN
    eff0 = 0d0
    RETURN
  ELSE IF (E(1) .GT. Eeff(NEeff)) THEN
    eff0 = 0d0
    RETURN
  END IF
  
  ! Energies below tabulation range
  K0 = 1
  DO WHILE (E(K0) .LT. Eeff(1))
    eff0(K0,:) = 0d0
    K0 = K0 + 1
  END DO
  
  ! Cycle over remaining energies
  Keff    = -1
  K0start = K0
  DO K0 = K0start,NE
    Keff = BSearch(NEeff,Eeff,E(K0),Keff)
    ! Below tabulation range: set to zero (already accounted for)
    !IF (Keff .LE. 0) THEN
    !  eff0(K0,:) = 0d0
    ! Within tabulation range: linear interpolation
    IF (Keff .LT. NEeff) THEN
      ! Determine fractional distance between tabulation points
      IF (Eeff(Keff) .EQ. Eeff(Keff+1)) THEN
        f = 0.5d0
      ELSE
        f = (E(K0)-Eeff(Keff)) / (Eeff(Keff+1)-Eeff(Keff))
      END IF
      eff0(K0,:) = (1d0-f)*eff(Keff,0:Neff0) + f*eff(Keff+1,0:Neff0)
    ! Above tabulation range: set to zero
    ELSE
      IF (E(K0) .EQ. Eeff(NEeff)) THEN
        eff0(K0,:) = eff(NEeff,0:Neff0)
      ELSE
        eff0(K0,:) = 0d0
      END IF
      ! Can skip remaining points
      IF (K0 .LT. NE) THEN
        eff0(K0+1:NE,:) = 0d0
      END IF
      EXIT
    END IF
  END DO
  
END SUBROUTINE



!=======================================================================
! DETECTOR SETUP ROUTINES
!=======================================================================

! ----------------------------------------------------------------------
! Get various detector quantities.
! Note this only allows access to fixed quantities, not quantities that
! are recalculated for each WIMP (see GetRates() for that).
! 
! Optional input argument:
!   D           The DetectorStruct structure to extract detector
!               quantities from.  If not given, a default structure
!               (internally stored) will be used.
! Optional output arguments:
!   mass        Detector fiducial mass [kg]
!   time        Detector exposure time [day]
!   exposure    Detector exposure [kg day]
!   Nevents     Number of observed events
!   background  Average expected background events
!   Niso        Number of isotopes
!   Ziso        Allocatable integer array to be filled with isotope
!               atomic numbers. Allocated to size [1:Niso].
!   Aiso        Allocatable integer array to be filled with isotope
!               atomic mass numbers. Allocated to size [1:Niso].
!   fiso        Allocatable real array to be filled with isotope
!               mass fractions. Allocated to size [1:Niso].
!   Miso        Allocatable real array to be filled with isotope
!               nuclear masses [GeV]. Allocated to size [1:Niso].
!   NE          Number of tabulated recoil energies E
!   E           Allocatable real array to be filled with recoil energies
!               [keV]. Allocated to size [1:NE].
!   eff_file    File from which efficiencies were read
!   NEeff       Number of recoil energies E for efficiencies tabulation
!   Eeff        Allocatable real array to be filled with recoil energies
!               used for efficiencies tabulation [keV]. Allocated to
!               size [1:NEeff].
!   Neff        Number of subintervals for which efficiencies are
!               available (0 for only total interval)
!   eff         Allocatable dimension=2 array containing efficiencies
!               as a function of recoil energy. Allocated to size
!               [1:NEeff,0:Neff], where the first index is over recoil
!               energies and the second index is over the sub-interval
!               number (0 for the total interval).
!   Wsi,Wsd     Allocatable dimension=3 array containing weighted form
!               factors for spin-independent (SI) and spin-dependent
!               (SD) couplings.  Allocated to size [-1:1,1:NE,1:Niso].
!   intervals   LOGICAL indicating if rates for intervals/bins are
!               to be calculated (used for max gap).
! 
SUBROUTINE GetDetector(D,mass,time,exposure,Nevents,background,         &
                       Niso,Ziso,Aiso,fiso,Miso,NE,E,                   &
                       eff_file,NEeff,Eeff,Neff,eff,                    &
                       Wsi,Wsd,intervals)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN), TARGET, OPTIONAL :: D
  LOGICAL, INTENT(OUT), OPTIONAL :: intervals
  CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: eff_file
  INTEGER, INTENT(OUT), OPTIONAL :: Nevents,Niso,NE,NEeff,Neff
  INTEGER, ALLOCATABLE, INTENT(OUT), OPTIONAL :: Ziso(:),Aiso(:)
  REAL*8, INTENT(OUT), OPTIONAL :: mass,time,exposure,background
  REAL*8, ALLOCATABLE, INTENT(OUT), OPTIONAL :: fiso(:),Miso(:),E(:),   &
          Eeff(:),eff(:,:),Wsi(:,:,:),Wsd(:,:,:)
  TYPE(DetectorStruct), POINTER :: DP
  
  ! Get pointer to detector structure to be used: either an explicitly
  ! provided structure or the default, internally-stored structure.
  IF (PRESENT(D)) THEN
    DP => D
  ELSE
    DP => DefaultDetector
  END IF
  
  ! Exposures
  IF (PRESENT(mass))     mass     = DP%mass
  IF (PRESENT(time))     time     = DP%time
  IF (PRESENT(exposure)) exposure = DP%exposure
  
  ! Observed events and expected background events
  IF (PRESENT(Nevents))    Nevents    = DP%Nevents
  IF (PRESENT(background)) background = DP%MuBackground
  
  ! Isotope data
  IF (PRESENT(Niso)) Niso = DP%Niso
  IF (PRESENT(Ziso)) THEN
    ALLOCATE(Ziso(DP%Niso))
    Ziso = DP%Ziso
  END IF
  IF (PRESENT(Aiso)) THEN
    ALLOCATE(Aiso(DP%Niso))
    Aiso = DP%Aiso
  END IF
  IF (PRESENT(fiso)) THEN
    ALLOCATE(fiso(DP%Niso))
    fiso = DP%fiso
  END IF
  IF (PRESENT(Miso)) THEN
    ALLOCATE(Miso(DP%Niso))
    Miso = DP%Miso
  END IF
  
  ! Recoil energies
  IF (PRESENT(NE)) NE = DP%NE
  IF (PRESENT(E)) THEN
    ALLOCATE(E(DP%NE))
    E = DP%E
  END IF
  
  ! Efficiencies
  IF (PRESENT(eff_file)) eff_file = DP%eff_file
  IF (PRESENT(NEeff))    NEeff    = DP%NEeff
  IF (PRESENT(Eeff)) THEN
    ALLOCATE(E(DP%NEeff))
    Eeff = DP%Eeff
  END IF
  IF (PRESENT(Neff))     Neff     = DP%Neff
  IF (PRESENT(eff)) THEN
    ALLOCATE(eff(DP%NEeff,0:DP%Neff))
    eff = DP%eff
  END IF
  
  ! Weighted form factors
  IF (PRESENT(Wsi)) THEN
    ALLOCATE(Wsi(-1:1,DP%NE,0:DP%Niso))
    Wsi = DP%Wsi
  END IF
  IF (PRESENT(Wsd)) THEN
    ALLOCATE(Wsd(-1:1,DP%NE,0:DP%Niso))
    Wsd = DP%Wsd
  END IF
  
  ! Calculate rates for intervals/bins?
  IF (PRESENT(intervals)) intervals = DP%intervals
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Set various detector quantities.  Making changes will reset all
! per-WIMP calculations like rates.
! Note this only allows access to fixed quantities, or quantities that
! need be calculated only once, not quantities that are recalculated
! for each WIMP (see SetRates() for that).
! 
! Optional input/output argument:
!   D           The DetectorStruct structure containing detector
!               data to be modified.  If not given, a default structure
!               (internally stored) will be used.
! Optional input arguments:
!   mass        Detector fiducial mass [kg]
!   time        Detector exposure time [day]
!   exposure    Detector exposure [kg day]
!   Nevents     Number of observed events
!   background  Average expected background events
! Optional isotope-related input arguments.  Niso must be given for
! any of the arrays to be used, regardless if the number has changed.
! If Niso changes, then all of the A/Z/f arrays must be specified for
! changes to take effect (otherwise these inputs are ignored).
!   Niso        Number of isotopes
!   Ziso        Integer array of size [1:Niso] containing atomic
!               numbers.
!   Aiso        Integer array of size [1:Niso] containing atomic
!               mass numbers.
!   fiso        Array of size [1:Niso] containing isotope mass
!               fractions.
! Pre-populated isotopic abundances for compounds with given
! stoichiometry.  The first two of these must be given to take effect
! (otherwise these inputs are ignored).
!   Nelem       Number of elements in compound.
!   Zelem       Integer array of size [1:Nelem] containing atomic
!               numbers of compound elements.
!   stoich      Integer array of size [1:Nelem] containing compound
!               stoichiometry.  For example, CF3Cl would be {1,3,1}.
!               Default is 1 for each element.
! Optional recoil energy tabulation arguments.  Both are required for
! changes to take effect.
!   NE          Number of tabulated recoil energies E
!   E           Array of size [1:NE] containing recoil energies [keV].
!               This defines the energies used for calculating dR/dE
!               and its integral.
! Optional efficiency curve(s) arguments.  Efficiencies can be provided
! by file or by providing tabulation data.  If tabulation data is given,
! all of NEeff, Eeff, Neff, and eff must be provided.
!   eff_file    File from which efficiencies shoud be read.
!   eff_filename Sets the stored file name _without_ loading any data
!               from the file.
!   NEeff       Number of recoil energy tabulation points for
!               efficiency data.
!   Eeff        Array of size [1:NEeff] containing tabulation energies
!               for efficiency data [keV].
!   Neff        Number of sub-interval/bin efficiencies in efficiency
!               data (0 if total only).
!   eff         Array of size [1:NEeff,0:Neff] containing efficiencies
!               as a function of recoil energy, where the first index
!               is over recoil energies and the second index is over
!               the sub-interval number (0 for the total interval).
! Optional analysis-related arguments.
!   intervals   Specify if rates for sub-intervals should be calculated;
!               otherwise, only the full interval is used for
!               calculations (default: true).  This basically specifies
!               if maximum gap or a Poisson is used in exclusion
!               constraints.
!   Emin        If given, sets all efficiencies below the given energy
!               [keV] to zero, removing all contributions from recoils
!               at lower energies.
! Optional form factor arguments.  The NE and Niso arguments must also
! be given for changes to take effect.
!   Wsi,Wsd     Arrays of size [-1:1,1:NE,1:Niso] containing weighted
!               form factors for spin-independent (SI) and spin-
!               dependent (SD) couplings.  Note: these are replaced with
!               internal calculations if the energy tabulation changes
!               in subsequent calls.  Should not be used in combination
!               with Emin.
! 
SUBROUTINE SetDetector(D,mass,time,exposure,Nevents,background,         &
                       Niso,Ziso,Aiso,fiso,Nelem,Zelem,stoich,          &
                       NE,E,                                            &
                       eff_file,eff_filename,NEeff,Eeff,Neff,eff,       &
                       intervals,Emin,Wsi,Wsd)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(INOUT), TARGET, OPTIONAL :: D
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: eff_file,eff_filename
  LOGICAL, INTENT(IN), OPTIONAL :: intervals
  INTEGER, INTENT(IN), OPTIONAL :: Nevents,Niso,Nelem,NE,NEeff,Neff
  !INTEGER, INTENT(IN), OPTIONAL :: Ziso(Niso),Aiso(Niso),Zelem(Nelem), &
  !                                 stoich(Nstoich)
  INTEGER, INTENT(IN), OPTIONAL :: Ziso(:),Aiso(:),Zelem(:),stoich(:)
  REAL*8, INTENT(IN), OPTIONAL :: mass,time,exposure,background,Emin
  !REAL*8, INTENT(IN), OPTIONAL :: fiso(Niso),E(NE),eff(NE,0:Neff),      &
  !        Wsi(-1:1,NE,Niso),Wsd(-1:1,NE,Niso)
  REAL*8, INTENT(IN), OPTIONAL :: fiso(:),E(:),Eeff(:),eff(:,0:),       &
          Wsi(-1:,:,:),Wsd(-1:,:,:)
  TYPE(DetectorStruct), POINTER :: DP
  LOGICAL :: iso_change,E_change,eff_change
  INTEGER :: KE,Kiso,Neff0
  INTEGER, ALLOCATABLE :: stoich0(:)
  REAL*8 :: f,dx
  ! For setting default energy tabulation (logarithmic spacing)
  INTEGER, PARAMETER :: E_PER_DECADE = 100
  REAL*8, PARAMETER :: DEFAULT_EMIN = 0.1d0
  REAL*8, PARAMETER :: DEFAULT_EMAX = 1000d0
  
  ! Get pointer to detector structure to be used: either an explicitly
  ! provided structure or the default, internally-stored structure.
  IF (PRESENT(D)) THEN
    DP => D
  ELSE
    DP => DefaultDetector
  END IF
  
  ! Indicate if quantities changed, indicating need for
  ! array resizing and initialization
  iso_change = .FALSE.
  E_change   = .FALSE.
  eff_change = .FALSE.
  
  ! Exposures
  IF (PRESENT(mass)) THEN
    DP%mass     = mass
    DP%exposure = DP%mass * DP%time
  END IF
  IF (PRESENT(time)) THEN
    DP%time     = time
    DP%exposure = DP%mass * DP%time
  END IF
  IF (PRESENT(exposure)) THEN
    DP%mass     = -1d0
    DP%time     = -1d0
    DP%exposure = exposure
  END IF
  
  ! Observed events and expected background events
  IF (PRESENT(Nevents))    DP%Nevents      = Nevents
  IF (PRESENT(background)) DP%MuBackground = background
  
  ! Isotope data
  IF (PRESENT(Niso)) THEN
    IF (Niso .EQ. DP%Niso) THEN
      IF (PRESENT(Ziso)) THEN
        DP%Ziso = Ziso(1:Niso)
        iso_change = .TRUE.
      END IF
      IF (PRESENT(Aiso)) THEN
        DP%Aiso = Aiso(1:Niso)
        iso_change = .TRUE.
      END IF
      IF (PRESENT(fiso)) THEN
        DP%fiso = fiso(1:Niso)
        iso_change = .TRUE.
      END IF
      IF (iso_change) THEN
        DP%Miso = IsotopeMass(DP%Ziso,DP%Aiso)
      END IF
    ELSE IF (PRESENT(Ziso) .AND. PRESENT(Aiso) .AND. PRESENT(fiso)) THEN
      IF (ALLOCATED(DP%Ziso)) DEALLOCATE(DP%Ziso)
      IF (ALLOCATED(DP%Aiso)) DEALLOCATE(DP%Aiso)
      IF (ALLOCATED(DP%fiso)) DEALLOCATE(DP%fiso)
      IF (ALLOCATED(DP%Miso)) DEALLOCATE(DP%Miso)
      ALLOCATE(DP%Ziso(Niso),DP%Aiso(Niso),DP%fiso(Niso),DP%Miso(Niso))
      DP%Niso = Niso
      DP%Ziso = Ziso(1:Niso)
      DP%Aiso = Aiso(1:Niso)
      DP%fiso = fiso(1:Niso)
      DP%Miso = IsotopeMass(DP%Ziso,DP%Aiso)
      iso_change = .TRUE.
    END IF
  END IF
  IF (PRESENT(Nelem) .AND. PRESENT(Zelem)) THEN
    ALLOCATE(stoich0(Nelem))
    IF (PRESENT(stoich)) THEN
      stoich0 = stoich
    ELSE
      stoich0 = 1
    END IF
    CALL CompoundIsotopeList(Nelem,Zelem,stoich0,                       &
                             DP%Niso,DP%Ziso,DP%Aiso,DP%fiso,DP%Miso)
    iso_change = .TRUE.
  END IF
  
  ! Energies
  IF (PRESENT(NE) .AND. PRESENT(E)) THEN
    IF (ALLOCATED(DP%E) .AND. (NE .NE. DP%NE)) DEALLOCATE(DP%E)
    IF (.NOT. ALLOCATED(DP%E)) ALLOCATE(DP%E(NE))
    DP%NE = NE
    DP%E = E(1:NE)
    IF (ALLOCATED(DP%E_cache)) DEALLOCATE(DP%E_cache)
    ALLOCATE(DP%E_cache(DP%NE))
    DP%E_cache = DP%E
    E_change = .TRUE.
  END IF
  
  ! Set default energies if not initialized or NE <= 0.
  ! Will use logarithmic spacing, but including E=0.
  IF (PRESENT(NE)) THEN
    IF (NE .LE. 0) DP%NE = NE
  END IF
  IF (DP%NE .LE. 0) THEN
    IF (ALLOCATED(DP%E)) DEALLOCATE(DP%E)
    DP%NE = NINT(E_PER_DECADE*(LOG10(DEFAULT_EMAX/DEFAULT_EMIN))) + 2
    ALLOCATE(DP%E(DP%NE))
    dx = 1d0/E_PER_DECADE
    DP%E(1) = 0d0
    DP%E(2) = DEFAULT_EMIN
    DO KE = 2, DP%NE-1
      DP%E(KE) = DEFAULT_EMIN * 10d0**((KE-2)*dx)
    END DO
    DP%E(DP%NE) = DEFAULT_EMAX
    IF (ALLOCATED(DP%E_cache)) DEALLOCATE(DP%E_cache)
    ALLOCATE(DP%E_cache(DP%NE))
    DP%E_cache = DP%E
    E_change = .TRUE.
  END IF
  
  ! Set efficiencies
  ! ...from file
  IF (PRESENT(eff_file)) THEN
    IF (eff_file .NE. '') THEN
      DP%eff_file = eff_file
      CALL LoadEfficiencyFile(eff_file,DP%NEeff,DP%Eeff,DP%Neff,DP%eff)
      eff_change = .TRUE.
    END IF
  ! ...by arguments
  ELSE IF (PRESENT(NEeff) .AND. PRESENT(Eeff) .AND. PRESENT(Neff)       &
           .AND. PRESENT(eff)) THEN
    IF (ALLOCATED(DP%Eeff)) DEALLOCATE(DP%Eeff)
    IF (ALLOCATED(DP%eff))  DEALLOCATE(DP%eff)
    ALLOCATE(DP%Eeff(NEeff),DP%eff(NEeff,0:Neff))
    DP%NEeff = NEeff
    DP%Eeff  = Eeff(1:NEeff)
    DP%Neff  = Neff
    DP%eff   = eff(1:NEeff,0:Neff)
    eff_change = .TRUE.
  END IF
  
  ! Set efficiency to 100% if not initialized or NEeff <= 0.
  ! No intervals.
  IF (PRESENT(NEeff)) THEN
    IF (NEeff .LE. 0) DP%NEeff = NEeff
  END IF
  IF (DP%NEeff .LE. 0) THEN
    IF (ALLOCATED(DP%Eeff)) DEALLOCATE(DP%Eeff)
    IF (ALLOCATED(DP%eff))  DEALLOCATE(DP%eff)
    ALLOCATE(DP%Eeff(2),DP%eff(2,0:0))
    DP%NEeff = 2
    DP%Eeff  = (/ 0d0, HUGE(DP%Eeff) /)
    DP%Neff  = 0
    DP%eff   = 1d0
    eff_change = .TRUE.
  END IF
  
  ! Saved efficiency file name
  IF (PRESENT(eff_filename)) DP%eff_file = eff_filename
  
  ! Include sub-intervals?
  IF (PRESENT(intervals)) THEN
    IF (intervals .NEQV. DP%intervals) eff_change = .TRUE.
    DP%intervals = intervals
  END IF
  
  ! Apply threshold cut.
  ! We move all E < Emin tabulation points to Emin.
  IF (PRESENT(Emin)) THEN
    ! First reset to original tabulation
    DP%E = MAX(DP%E_cache,Emin)
    E_change = .TRUE.
  END IF
  
  ! Update efficiencies
  IF (eff_change .OR. E_change) THEN
    CALL RetabulateEfficiency(DP%NEeff,DP%Eeff,DP%Neff,DP%eff,          &
                              DP%NE,DP%E,DP%intervals,DP%eff0)
  END IF
  
  ! Weighted form factors (SI)
  IF (PRESENT(NE) .AND. PRESENT(Niso) .AND. PRESENT(Wsi)) THEN
    IF (ALLOCATED(DP%Wsi)) DEALLOCATE(DP%Wsi)
    ALLOCATE(DP%Wsi(-1:1,NE,Niso))
    DP%Wsi = Wsi(-1:1,1:NE,1:Niso)
  ELSE IF ((iso_change .OR. E_change) .AND. (DP%Niso .GE. 0)            &
           .AND. (DP%NE .GE. 0)) THEN
    IF (ALLOCATED(DP%Wsi)) DEALLOCATE(DP%Wsi)
    ALLOCATE(DP%Wsi(-1:1,DP%NE,DP%Niso))
    DO Kiso = 1,DP%Niso
      CALL CalcWSI(DP%Ziso(Kiso),DP%Aiso(Kiso),DP%NE,                   &
                   EToQ(DP%E,DP%Miso(Kiso)),DP%Wsi(:,:,Kiso))
    END DO
  END IF
  
  ! Weighted form factors (SD)
  IF (PRESENT(NE) .AND. PRESENT(Niso) .AND. PRESENT(Wsd)) THEN
    IF (ALLOCATED(DP%Wsd)) DEALLOCATE(DP%Wsd)
    ALLOCATE(DP%Wsd(-1:1,NE,Niso))
    DP%Wsd = Wsd(-1:1,1:NE,1:Niso)
  ELSE IF ((iso_change .OR. E_change) .AND. (DP%Niso .GE. 0)            &
           .AND. (DP%NE .GE. 0)) THEN
    IF (ALLOCATED(DP%Wsd)) DEALLOCATE(DP%Wsd)
    ALLOCATE(DP%Wsd(-1:1,DP%NE,DP%Niso))
    DO Kiso = 1,DP%Niso
      CALL CalcWSD(DP%Ziso(Kiso),DP%Aiso(Kiso),DP%NE,                   &
                   EToQ(DP%E,DP%Miso(Kiso)),DP%Wsd(:,:,Kiso))
    END DO
  END IF
  
  ! Resize halo velocity arrays if necessary
  IF ((iso_change .OR. E_change) .AND. (DP%Niso .GE. 0)                 &
      .AND. (DP%NE .GE. 0)) THEN
    IF (ALLOCATED(DP%vmin)) DEALLOCATE(DP%vmin)
    ALLOCATE(DP%vmin(DP%NE,DP%Niso))
    IF (ALLOCATED(DP%eta))  DEALLOCATE(DP%eta)
    ALLOCATE(DP%eta(DP%NE,DP%Niso))
  END IF
  
  ! Number of intervals/bins to do calculations for.
  ! Used for array sizing below.
  IF (DP%intervals) THEN
    Neff0 = DP%Neff
  ELSE
    Neff0 = 0
  END IF
  
  ! Resize rate arrays if necessary
  IF ((E_change .OR. eff_change) .AND. (DP%NE .GE. 0)                   &
      .AND. (DP%Neff .GE. 0)) THEN
    IF (ALLOCATED(DP%dRdEsi0)) DEALLOCATE(DP%dRdEsi0)
    ALLOCATE(DP%dRdEsi0(-1:1,DP%NE))
    IF (ALLOCATED(DP%dRdEsd0)) DEALLOCATE(DP%dRdEsd0)
    ALLOCATE(DP%dRdEsd0(-1:1,DP%NE))
    IF (ALLOCATED(DP%Rsi0)) DEALLOCATE(DP%Rsi0)
    ALLOCATE(DP%Rsi0(-1:1,0:Neff0))
    IF (ALLOCATED(DP%Rsd0)) DEALLOCATE(DP%Rsd0)
    ALLOCATE(DP%Rsd0(-1:1,0:Neff0))
    IF (ALLOCATED(DP%dRdEsi)) DEALLOCATE(DP%dRdEsi)
    ALLOCATE(DP%dRdEsi(DP%NE))
    IF (ALLOCATED(DP%dRdEsd)) DEALLOCATE(DP%dRdEsd)
    ALLOCATE(DP%dRdEsd(DP%NE))
    IF (ALLOCATED(DP%dRdE)) DEALLOCATE(DP%dRdE)
    ALLOCATE(DP%dRdE(DP%NE))
    IF (ALLOCATED(DP%Rsi)) DEALLOCATE(DP%Rsi)
    ALLOCATE(DP%Rsi(0:Neff0))
    IF (ALLOCATED(DP%Rsd)) DEALLOCATE(DP%Rsd)
    ALLOCATE(DP%Rsd(0:Neff0))
    IF (ALLOCATED(DP%R)) DEALLOCATE(DP%R)
    ALLOCATE(DP%R(0:Neff0))
  END IF
  
  ! Resize event arrays if necessary
  IF (eff_change .AND. (DP%Neff .GE. 0)) THEN
    IF (ALLOCATED(DP%MuSignalSI0)) DEALLOCATE(DP%MuSignalSI0)
    ALLOCATE(DP%MuSignalSI0(-1:1,0:Neff0))
    IF (ALLOCATED(DP%MuSignalSD0)) DEALLOCATE(DP%MuSignalSD0)
    ALLOCATE(DP%MuSignalSD0(-1:1,0:Neff0))
    IF (ALLOCATED(DP%MuSignalSI)) DEALLOCATE(DP%MuSignalSI)
    ALLOCATE(DP%MuSignalSI(0:Neff0))
    IF (ALLOCATED(DP%MuSignalSD)) DEALLOCATE(DP%MuSignalSD)
    ALLOCATE(DP%MuSignalSD(0:Neff0))
    IF (ALLOCATED(DP%MuSignal)) DEALLOCATE(DP%MuSignal)
    ALLOCATE(DP%MuSignal(0:Neff0))
  END IF
  
  ! Set all calculable quantities to zero
  IF ((iso_change .OR. E_change .OR. eff_change)                        &
       .AND. (DP%Niso .GE. 0) .AND. (DP%NE .GE. 0)                      &
       .AND. (DP%Neff .GE. 0)) THEN
    DP%vmin        = 0d0
    DP%eta         = 0d0
    DP%dRdEsi0     = 0d0
    DP%dRdEsd0     = 0d0
    DP%Rsi0        = 0d0
    DP%Rsd0        = 0d0
    DP%dRdEsi      = 0d0
    DP%dRdEsd      = 0d0
    DP%dRdE        = 0d0
    DP%Rsi         = 0d0
    DP%Rsd         = 0d0
    DP%R           = 0d0
    DP%MuSignalSI0 = 0d0
    DP%MuSignalSD0 = 0d0
    DP%MuSignalSI  = 0d0
    DP%MuSignalSD  = 0d0
    DP%MuSignal    = 0d0
  END IF
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Initializes the detector.
! Simply sets the detector to match that of the LUX 2013 analysis.
! If something other than this is desired, use the SetDetector()
! function after this.
! 
! Optional output argument:
!   D           The DetectorStruct structure to be initialized.
!               If not given, the default structure (internally stored)
!               will be initialized.
! Optional input argument:
!   intervals   Specify if sub-intervals should be loaded; otherwise,
!               only the full interval is used for calculations
!               (default: true).
! 
SUBROUTINE InitDetector(D,intervals)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(OUT), TARGET, OPTIONAL :: D
  LOGICAL, INTENT(IN), OPTIONAL :: intervals
  TYPE(DetectorStruct), POINTER :: DP
  LOGICAL :: intervals0
  
  ! Get pointer to detector structure to be used: either an explicitly
  ! provided structure or the default, internally-stored structure.
  IF (PRESENT(D)) THEN
    DP => D
  ELSE
    DP => DefaultDetector
  END IF
  
  ! Include sub-intervals?
  intervals0 = .TRUE.
  IF (PRESENT(intervals)) intervals0 = intervals
  
  ! Set to LUX 2013 analysis
  CALL LUX_2013_InitTo(DP,intervals0)
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Initializes the detector setup from command-line parameters.
! 
! Possible options:
!   --mass=<value>       ! Detector fiducial mass [kg]
!   --time=<value>       ! Detector exposure time [day]
!   --exposure=<value>   ! Detector exposure [kg day]
!   --events=<N>         ! Number of observed events
!   --background=<value> ! Average expected background events
! Isotope-related options specifying the atomic number, atomic mass
! number, and mass fraction of the isotopes.  Must specify all three
! and they must have the same number of isotopes for changes to take
! effect:
!   --isotope-Z=<Z1>,<Z2>,<Z3>,...
!   --isotope-A=<A1>,<A2>,<A3>,...
!   --isotope-f=<f1>,<f2>,<f3>,...
! Use naturally occurring isotopes for the elements in the given
! compound specified by atomic numbers and the stoichiometry.  If
! stoichiometry is not given, it is assumed to be 1 for each element
! (equal abundances):
!   --element-Z=<Z1>,<Z2>,<Z3>,...
!   --stoichiometry=<stoich1>,<stoich2>,<stoich3>,...
! Shortcuts for specific materials by name (using natural abundances):
!   --argon
!   --germanium
!   --sodium-iodide
!   --silicon
!   --xenon
! Efficiency related options:
!   --file=<file>        ! File from which tabulated efficiencies should be read.
!                        ! First column is E [keV], any following consecutive
!                        ! columns containing values above 1.0 are ignored, and
!                        ! the remaining columns are taken to be efficiency curves.
!                        ! First efficiency is for full range, remaining efficiencies
!                        ! are for sub-intervals.  If not given, a default is used.
!   --no-intervals       ! Disable use of sub-intervals, even if available in file.
!   --Emin=<value>       ! If given, sets all efficiencies below the given energy
!                        ! [keV] to zero, removing all contributions from recoils
!                        ! at lower energies.  Not reversible.
! 
! Optional output argument:
!   D           The DetectorStruct structure to be initialized.
!               If not given, the default structure (internally stored)
!               will be initialized.
! Optional input arguments:
!   eff_file    File from which efficiencies shoud be read.  If not
!               given or set to '', an internal default will be used
!               (the LUX 2013 result).  Non-empty value takes precedence
!               over --file option.
!   intervals   Specify if sub-intervals should be loaded from the
!               efficiency file (if available); otherwise, only the
!               full interval is used for calculations (default: true).
!               Superceded by --no-intervals setting.
! 
SUBROUTINE InitDetectorCL(D,eff_file,intervals)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(OUT), TARGET, OPTIONAL :: D
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: eff_file
  LOGICAL, INTENT(IN), OPTIONAL :: intervals
  TYPE(DetectorStruct), POINTER :: DP
  CHARACTER(LEN=1024) :: eff_file0
  LOGICAL :: intervals0
  INTEGER :: Nevents,Niso,Nelem,N1,N2,N3
  INTEGER, ALLOCATABLE :: Ziso(:),Aiso(:),Zelem(:),stoich(:)
  REAL*8 :: mass,time,exposure,background,Emin
  REAL*8, ALLOCATABLE :: fiso(:)
  
  ! Get pointer to detector structure to be used: either an explicitly
  ! provided structure or the default, internally-stored structure.
  IF (PRESENT(D)) THEN
    DP => D
  ELSE
    DP => DefaultDetector
  END IF
  
  Niso = -1
  
  ! File name ('' for internal default)
  IF (.NOT. GetLongArgString('file',eff_file0)) eff_file0 = ''
  IF (PRESENT(eff_file)) THEN
    IF (TRIM(eff_file) .NE. '') eff_file0 = eff_file
  END IF
  
  ! Include sub-intervals?
  intervals0 = .TRUE.
  IF (PRESENT(intervals)) intervals0 = intervals
  IF (GetLongArg('no-intervals')) intervals0 = .FALSE.

  ! Experiment-specific settings
  ! XENON100 2012 result
  IF (GetLongArg('XENON100-2012')) THEN
    CALL XENON100_2012_InitTo(DP,intervals0)
  ! LUX 2013 result
  ELSE IF (GetLongArg('LUX-2013')) THEN
    CALL LUX_2013_InitTo(DP,intervals0)
  ! SuperCDMS 2014 result (low-energy analysis)
  ELSE IF (GetLongArg('SuperCDMS-2014')) THEN
    CALL SuperCDMS_2014_InitTo(DP,intervals0)
  ! SIMPLE 2014 result
  ELSE IF (GetLongArg('SIMPLE-2014')) THEN
    CALL SIMPLE_2014_InitTo(DP,intervals0)
  ! DARWIN proposal, xenon-based (as of 2015)
  ELSE IF (GetLongArg('DARWIN-Xe-2015')) THEN
    CALL DARWIN_Xe_2015_InitTo(DP,intervals0)
  ! DARWIN proposal, argon-based (as of 2015)
  ELSE IF (GetLongArg('DARWIN-Ar-2015')) THEN
    CALL DARWIN_Ar_2015_InitTo(DP,intervals0)
  ! Defaults (LUX 2013)
  ELSE
    CALL LUX_2013_InitTo(DP,intervals0)
    CALL SetDetector(eff_filename='[LUX 2013 (default)]')
  END IF
  
  ! If given, load efficiency file
  IF (TRIM(eff_file0) .NE. '') THEN
    CALL SetDetector(DP,eff_file=eff_file0,intervals=intervals0)
  END IF
  
  ! Update isotopes
  IF (GetLongArgIntegers('isotope-Z',Ziso,N1)                           &
      .AND. GetLongArgIntegers('isotope-A',Aiso,N2)                     &
      .AND. GetLongArgReals('isotope-f',fiso,N3)) THEN
    IF ((N1 .EQ. N2) .AND. (N1 .EQ. N3)) THEN
      Niso = N1
      CALL SetDetector(DP,Niso=Niso,Ziso=Ziso,Aiso=Aiso,fiso=fiso)
    END IF
  ELSE IF (GetLongArgIntegers('element-Z',Zelem,Nelem)) THEN
    IF (GetLongArgIntegers('stoichiometry',stoich,N2)) THEN
      IF (N2 .EQ. Nelem) THEN
        CALL SetDetector(DP,Nelem=Nelem,Zelem=Zelem,stoich=stoich)
      END IF
    ELSE
      IF (ALLOCATED(stoich)) DEALLOCATE(stoich)
      ALLOCATE(stoich(Nelem))
      stoich = 1
      CALL SetDetector(DP,Nelem=Nelem,Zelem=Zelem,stoich=stoich)
    END IF
  ELSE IF (GetLongArg('argon')) THEN
    CALL SetDetector(DP,Nelem=1,Zelem=(/18/))
  ELSE IF (GetLongArg('germanium')) THEN
    CALL SetDetector(DP,Nelem=1,Zelem=(/32/))
  ELSE IF (GetLongArg('silicon')) THEN
    CALL SetDetector(DP,Nelem=1,Zelem=(/14/))
  ELSE IF (GetLongArg('xenon')) THEN
    CALL SetDetector(DP,Nelem=1,Zelem=(/54/))
  ELSE IF (GetLongArg('sodium-iodide')) THEN
    CALL SetDetector(DP,Nelem=2,Zelem=(/11,53/),stoich=(/1,1/))
  END IF
  
  ! Update exposures
  IF (GetLongArgReal('mass',mass))         CALL SetDetector(DP,mass=mass)
  IF (GetLongArgReal('time',time))         CALL SetDetector(DP,time=time)
  IF (GetLongArgReal('exposure',exposure)) CALL SetDetector(DP,exposure=exposure)
  
  ! Update observed events and expected background events
  IF (GetLongArgInteger('events',Nevents))     CALL SetDetector(DP,Nevents=Nevents)
  IF (GetLongArgReal('background',background)) CALL SetDetector(DP,background=background)
  
  ! Set an energy threshold
  IF (GetLongArgReal('Emin',Emin)) THEN
    CALL SetDetector(DP,Emin=Emin)
  END IF
  
END SUBROUTINE



!=======================================================================
! RATE ROUTINES
!=======================================================================

! ----------------------------------------------------------------------
! Get various rate quantities.  These require CalcRates() to have been
! run first.
! 
! Optional input argument:
!   D           The DetectorStruct structure to extract rate
!               quantities from.  If not given, a default structure
!               (internally stored) will be used.
! Optional output arguments:
!   Nevents     Number of observed events
!   background  Average expected background events
!   signal      Average expected signal events
!   signal_si   Average expected spin-independent signal events
!   signal_sd   Average expected spin-dependent signal events
!   rate        Signal rate [cpd/kg]
!   rate_si     Spin-independent signal rate [cpd/kg]
!   rate_sd     Spin-dependent signal rate [cpd/kg]
! Optional sub-interval/bin arguments.  Arrays will be allocated to
! size [1:Nbins]:
!   Nbins       Number of bins/intervals
!   binsignal   Allocatable array to be filled with average expected
!               signal in each bin.  Allocated to size [1:Nbins].
!               The sum of all bins is equal to 'signal'.
!   binsignal_si,binsignal_sd
!               Same as 'binsignal' but for only spin-independent or
!               spin-dependent events, respectively.
!   binrate      Allocatable array to be filled with recoil rate in
!               each bin.  Allocated to size [1:Nbins].
!               The sum of all bins is equal to 'rate'.
!   binrate_si,binrate_sd
!               Same as 'binrate' but for only spin-independent or
!               spin-dependent events, respectively.
!   
SUBROUTINE GetRates(D,Nevents,background,signal,signal_si,signal_sd,    &
                    rate,rate_si,rate_sd,Nbins,binsignal,binsignal_si,  &
                    binsignal_sd,binrate,binrate_si,binrate_sd)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(IN), TARGET, OPTIONAL :: D
  INTEGER, INTENT(OUT), OPTIONAL :: Nevents,Nbins
  REAL*8, INTENT(OUT), OPTIONAL :: background,signal,signal_si,signal_sd,&
          rate,rate_si,rate_sd
  REAL*8, ALLOCATABLE, INTENT(OUT), OPTIONAL :: binsignal(:),           &
          binsignal_si(:),binsignal_sd(:),binrate(:),binrate_si(:),     &
          binrate_sd(:)
  TYPE(DetectorStruct), POINTER :: DP
  INTEGER :: Neff
  
  ! Get pointer to detector structure to be used: either an explicitly
  ! provided structure or the default, internally-stored structure.
  IF (PRESENT(D)) THEN
    DP => D
  ELSE
    DP => DefaultDetector
  END IF
  
  Neff = DP%Neff
  
  ! Observed events and expected background events
  IF (PRESENT(Nevents))    Nevents    = DP%Nevents
  IF (PRESENT(background)) background = DP%MuBackground
  
  ! Signal events
  IF (PRESENT(signal))    signal    = DP%MuSignal(0)
  IF (PRESENT(signal_si)) signal_si = DP%MuSignalSI(0)
  IF (PRESENT(signal_sd)) signal_sd = DP%MuSignalSD(0)
  
  ! Signal rates
  IF (PRESENT(rate))    rate    = DP%R(0)
  IF (PRESENT(rate_si)) rate_si = DP%Rsi(0)
  IF (PRESENT(rate_sd)) rate_sd = DP%Rsd(0)
  
  ! Bins
  IF (PRESENT(Nbins)) Nbins = Neff
  
  ! Signal events by bin
  IF (PRESENT(binsignal)) THEN
    ALLOCATE(binsignal(Neff))
    binsignal = DP%MuSignal(1:Neff)
  END IF
  IF (PRESENT(binsignal_si)) THEN
    ALLOCATE(binsignal_si(Neff))
    binsignal_si = DP%MuSignalSI(1:Neff)
  END IF
  IF (PRESENT(binsignal_sd)) THEN
    ALLOCATE(binsignal_sd(Neff))
    binsignal_sd = DP%MuSignalSD(1:Neff)
  END IF
  
  ! Signal rates by bin
  IF (PRESENT(binrate)) THEN
    ALLOCATE(binrate(Neff))
    binrate = DP%MuSignal(1:Neff)
  END IF
  IF (PRESENT(binrate_si)) THEN
    ALLOCATE(binrate_si(Neff))
    binrate_si = DP%MuSignalSI(1:Neff)
  END IF
  IF (PRESENT(binrate_sd)) THEN
    ALLOCATE(binrate_sd(Neff))
    binrate_sd = DP%MuSignalSD(1:Neff)
  END IF
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Calculates various rate quantities using the current WIMP.
! This is the main routine intended for putting together the WIMP
! rate calculations and must be called any time the WIMP mass
! and/or couplings are modified.
! 
! Optional input/output argument:
!   D           The DetectorStruct structure containing detector
!               data, also where calculated rate data will be placed.
!               If not given, a default structure (internally stored)
!               will be used.
! 
SUBROUTINE CalcRates(D)
  IMPLICIT NONE
  TYPE(DetectorStruct), INTENT(INOUT), TARGET, OPTIONAL :: D
  TYPE(DetectorStruct), POINTER :: DP
  INTEGER :: Kiso,KE,Keff,Neff0
  REAL*8 :: alphasi(-1:1),alphasd(-1:1)
  ! Constant used to convert units:
  !   s / (cm^3 km GeV^4)  -->  cpd/kg/keV
  ! Includes factor of hbar^2 c^4.
  REAL*8, PARAMETER :: TO_CPD_KG_KEV = 1.695e14
  
  !---------------------------------------------
  ! Differential rate is given by:
  !   dR/dE(E) = 1/(2m\mu^2) \sigma(q) \rho \eta(vmin)
  ! where q(E) = \sqrt{2ME}, vmin(E) = \sqrt{ME/2\mu^2}, and:
  !   \sigma(E) = \mu^2 (hbar c)^2 [W(1,E)*Gp^2 + W(0,E)*Gp*Gn + W(-1,E)*Gn^2]
  ! This gives:
  !   dR/dE(E) = 1/2m \rho \eta(E) [W(1)*Gp^2 + W(0)*Gp*Gn + W(-1)*Gn^2]
  ! The total differential rate is a mass-fraction weighted sum of
  ! dR/dE over different isotopes and is summed over spin-independent
  ! (SI) and spin-dependent (SD) contributions.
  ! 
  ! The weighted form factors W are defined as:
  !   Wsi(+1,E) = (1/pi) Z^2 F^2(E)        ! SI proton
  !   Wsi( 0,E) = (1/pi) 2*Z*(A-Z) F^2(E)  ! SI crossterm
  !   Wsi(-1,E) = (1/pi) (A-Z)^2 F^2(E)    ! SI neutron
  !   Wsd(+1,E) = 4/(2J+1) Spp(E)          ! SD proton
  !   Wsd( 0,E) = 4/(2J+1) Spn(E)          ! SD crossterm
  !   Wsd(-1,E) = 4/(2J+1) Snn(E)          ! SD neutron
  ! The above definitions give for the conventional SI and SD
  ! cross-sections:
  !   \sigma(E) = \mu^2 (hbar c)^2 [W(+1,E)*Gp^2 + W(0,E)*Gp*Gn + W(-1,E)*Gn^2]
  ! where Gp and Gn are the effective proton and neutron couplings
  ! in units of [GeV^-2] and \mu is the reduced mass.  In terms of
  ! more commonly used notation:
  !   SI (scalar):        G = 2f
  !   SD (axial-vector):  G = 2\sqrt{2} G_F a
  ! where G, f, and a have 'p' and 'n' subscripts.
  ! 
  ! Thes spin-independent and spin-dependent weighted form factors are
  ! tabulated over E for each of the isotopes, so the Wsi and Wsd arrays
  ! are of size of size [-1:1,1:NE,1:Niso] where NE is the number of
  ! tabulation points (energies) and Niso is the number of isotopes.
  ! As these factors are independent of the WIMP, they are calculated
  ! only once, during initialization (nothing needs to be done here).
  ! SD FORM FACTORS ARE ONLY IMPLEMENTED FOR XENON.  ALL OTHER SD FORM
  ! FACTORS ARE SET TO ZERO.
  ! 
  ! The mean inverse speed
  !   eta(E) = \int_{v>vmin} d^3v (1/v) f(v),
  ! with
  !   vmin(E) = \sqrt{M E/(2\mu^2)}
  ! the minimum WIMP speed capable of causing a nucleus to recoil with
  ! energy E, depends upon the WIMP mass, so it will need to be
  ! recalculated for every WIMP.  Both vmin and eta are tabulated over
  ! E for each of the isotopes, so the vmin and eta arrays are of size
  ! [1:NE,1:Niso].
  ! 
  ! When factoring in detector efficiencies and energy resolution, the
  ! total rate of observed events is given by:
  !   R = \int_0^\infty dE \phi(E) dR/dE(E)
  ! where \phi(E) is the efficiency/response factor defined as the
  ! fraction of events at a given energy E that will be both observed
  ! and fall within some particular S1 analysis range.  The actual
  ! observed spectrum is dR/dS1(S1), NOT the quantity \phi(E) dR/dE(E);
  ! the former cannot be directly calculated using \phi(E).  However,
  ! the definition of \phi(E) allows the integral of each spectrum to
  ! be related as we can also write the rate as:
  !   R = \int_{S1min}^{S1max} dS1 dR/dS1(S1)
  ! The \phi(E) are calculated for _specific_ ranges of S1.  As these
  ! efficiencies are independent of the WIMP, they do not need to be
  ! recalculated, which is good, since this code cannot calculate
  ! them.  Instead, the the efficiency curve(s) are loaded from a file
  ! (or uses the internal default) during the initialization stage.
  ! 
  ! The efficiency file should contain a tabulated efficiency for the
  ! full desired S1 range, but can optionally include additional
  ! tabulated efficiencies for sub-intervals of that full range.
  ! For example, LUX 2013 analyzed the S1 range [2,30] and observed one
  ! event at S1=3.1 below their nuclear recoil mean.  For the maximum
  ! gap method of determining exclusion constraints in the presence of
  ! unknown backgrounds, the expected rate in intervals between events
  ! is used, so the rate in the S1 ranges [2,3.1] and [3.1,20] must be
  ! calculated.  The default efficiencies are the \phi(E) efficiencies
  ! for these three ranges: the full range and the two sub-intervals.
  ! The tabulated efficiency array (called 'eff') is an array of size
  ! [1:NE,0:Neff] where Neff is the number of sub-intervals (0 if
  ! separate interval efficiencies are not to be used or were not
  ! available).  The second index is that of the sub-interval, with 0
  ! used for the total range.
  ! 
  ! To avoid having to deal with a variety of interpolation issues,
  ! everything is tabulated at the same energies E, given in the array
  ! E of size [1:NE].  The tabulation is taken from that of the
  ! provided efficiency file (or internal default) since efficiencies
  ! cannot be recalculated with this program (everything else can).
  ! 
  ! In the calculations below, we calculate the differential rates
  ! dR/dE(E) for each of the isotopes and each of the SI & SD couplings.
  ! Then we perform an efficiency-weighted integral over E and sum
  ! contributions together to get the total rates and expected events.
  ! The separate coupling spectra are saved to allow for possible
  ! further inspection.
  ! 
  ! To be clear, the following are required to determine the full rate:
  !   * Sum over SI & SD contributions (with a Gp^2, Gn^2, and Gp*Gn
  !     term for each)
  !   * Mass fraction weighted sum over isotope contributions
  !   * Integration over energy
  !   * All of the above possibly for multiple analysis intervals
  ! Though not necesssary, various intermediate quantities are kept
  ! to allow for possible further inspection (aside from individual
  ! isotope contributions, which are not kept).
  ! Note there is no required order for performing the sums and
  ! integrations; we take advantage of this to perform the SI/SD
  ! component sum _last_ so that the rate is very easily recalculated
  ! for different couplings.  Notably, the two cross-sections are
  ! treated as:
  !    \sigma(E) = (Gp/Gp0)^2 \sigmapp0(E) + (Gn/Gn0)^2 \sigmann0(E)
  !                + (Gp*Gn/Gp0*Gn0) \sigmapn0(E)
  ! where
  !    \sigmapp0(E) = \mu^2 (hbar c)^2 W(+1,E) * Gp0^2
  !    \sigmann0(E) = \mu^2 (hbar c)^2 W(-1,E) * Gn0^2
  !    \sigmapn0(E) = \mu^2 (hbar c)^2 W( 0,E) * Gp0*Gn0
  ! and Gp0 & Gn0 are reference couplings, which we take to be those
  ! that yield WIMP-nucleon cross-sections of 1 pb.  In that case,
  ! we can also write:
  !    \sigma(E) = (sigmap/[pb])*\sigmapp0(E) + (sigman/[pb])*\sigmann0(E)
  !                +/- \sqrt{(sigmap/[pb])(sigman/[pb])} \sigmapn0(E)
  ! with the sign of the cross-term equal to the sign of Gp*Gn.  The
  ! quantities below with names ending with 'si0' and 'sd0' have an
  ! index [-1:1] identical to that described for the W terms above,
  ! corresponds to the that quantity calculated for sigma(E) =
  ! sigmapp0(E) [+1], sigmann0(E) [-1], or sigmapn0(E) [0].
  !---------------------------------------------
  
  ! Get pointer to detector structure to be used: either an explicitly
  ! provided structure or the default, internally-stored structure.
  IF (PRESENT(D)) THEN
    DP => D
  ELSE
    DP => DefaultDetector
  END IF
  
  ! Update mean inverse speed.
  DP%vmin = EToVmin(DP%NE,DP%E,WIMP%m,DP%Niso,DP%Miso)
  DP%eta  = MeanInverseSpeed(DP%NE,DP%Niso,DP%vmin)
  
  ! Below, we will calculate reference rates for reference couplings
  ! GpSI0, GnSI0, GpSD0, & GnSD0 defined such that WIMP-nucleon
  ! cross-sections are 1 pb.  To convert between reference rates and
  ! actual rates, we multiply by alpha = G^2/G0^2.  Indices of alpha
  ! match that of the W arrays.
  alphasi(+1) = (WIMP%GpSI/WIMP%GpSI0)**2                      ! SI proton
  alphasi( 0) = (WIMP%GpSI/WIMP%GpSI0)*(WIMP%GnSI/WIMP%GnSI0)  ! SI crossterm
  alphasi(-1) = (WIMP%GnSI/WIMP%GnSI0)**2                      ! SI neutron
  alphasd(+1) = (WIMP%GpSD/WIMP%GpSD0)**2                      ! SD proton
  alphasd( 0) = (WIMP%GpSD/WIMP%GpSD0)*(WIMP%GnSD/WIMP%GnSD0)  ! SD crossterm
  alphasd(-1) = (WIMP%GnSD/WIMP%GnSD0)**2                      ! SD neutron
  
  ! Reference differential rates.
  ! Keep separate contributions from Gp^2, Gn^2, and Gp*Gn terms.
  DP%dRdEsi0 = 0d0
  DP%dRdEsd0 = 0d0
  ! Mass-fraction weighted sum over isotopes.
  DO Kiso = 1,DP%Niso
    ! SI proton term
    DP%dRdEsi0(+1,:) = DP%dRdEsi0(+1,:)                                 &
        +  DP%fiso(Kiso) / (2 * WIMP%m)                                 &
           * Halo%rho * DP%eta(:,Kiso)                                  &
           * DP%Wsi(+1,:,Kiso) * WIMP%GpSI0**2                          &
           * TO_CPD_KG_KEV
    ! SI cross term
    DP%dRdEsi0( 0,:) = DP%dRdEsi0( 0,:)                                 &
        +  DP%fiso(Kiso) / (2 * WIMP%m)                                 &
           * Halo%rho * DP%eta(:,Kiso)                                  &
           * DP%Wsi( 0,:,Kiso) * WIMP%GpSI0*WIMP%GnSI0                  &
           * TO_CPD_KG_KEV
    ! SI neutron term
    DP%dRdEsi0(-1,:) = DP%dRdEsi0(-1,:)                                 &
        +  DP%fiso(Kiso) / (2 * WIMP%m)                                 &
           * Halo%rho * DP%eta(:,Kiso)                                  &
           * DP%Wsi(-1,:,Kiso) * WIMP%GnSI0**2                          &
           * TO_CPD_KG_KEV
    ! SD proton term
    DP%dRdEsd0(+1,:) = DP%dRdEsd0(+1,:)                                 &
        +  DP%fiso(Kiso) / (2 * WIMP%m)                                 &
           * Halo%rho * DP%eta(:,Kiso)                                  &
           * DP%Wsd(+1,:,Kiso) * WIMP%GpSD0**2                          &
           * TO_CPD_KG_KEV
    ! SD cross term
    DP%dRdEsd0( 0,:) = DP%dRdEsd0( 0,:)                                 &
        +  DP%fiso(Kiso) / (2 * WIMP%m)                                 &
           * Halo%rho * DP%eta(:,Kiso)                                  &
           * DP%Wsd( 0,:,Kiso) * WIMP%GpSD0*WIMP%GnSD0                  &
           * TO_CPD_KG_KEV
    ! SD neutron term
    DP%dRdEsd0(-1,:) = DP%dRdEsd0(-1,:)                                 &
        +  DP%fiso(Kiso) / (2 * WIMP%m)                                 &
           * Halo%rho * DP%eta(:,Kiso)                                  &
           * DP%Wsd(-1,:,Kiso) * WIMP%GnSD0**2                          &
           * TO_CPD_KG_KEV
  END DO
  
  ! Number of intervals/bins to do calculations for.
  ! If intervals=.FALSE. then efficiency index is over [0:0]
  ! (total only).
  IF (DP%intervals) THEN
    Neff0 = DP%Neff
  ELSE
    Neff0 = 0
  END IF
  
  ! Integrate (efficiency-weighted) to find total rates.
  ! Uses a simple trapezoidal integration.
  DP%Rsi0 = 0d0
  DP%Rsd0 = 0d0
  ! Cycle over E bins and efficiency curves.
  DO KE = 1,DP%NE-1
    DO Keff = 0,Neff0
      DP%Rsi0(:,Keff) = DP%Rsi0(:,Keff)                                 &
          + 0.5d0 * (DP%E(KE+1) - DP%E(KE))                             &
            * (DP%eff0(KE,Keff)*DP%dRdEsi0(:,KE)                        &
               + DP%eff0(KE+1,Keff)*DP%dRdEsi0(:,KE+1))
      DP%Rsd0(:,Keff) = DP%Rsd0(:,Keff)                                 &
          + 0.5d0 * (DP%E(KE+1) - DP%E(KE))                             &
            * (DP%eff0(KE,Keff)*DP%dRdEsd0(:,KE)                        &
               + DP%eff0(KE+1,Keff)*DP%dRdEsd0(:,KE+1))
    END DO
  END DO
  
  ! Rates for actual couplings.  Here do separate SI and SD rates
  ! as well as total.
  DP%dRdEsi(:) =  alphasi(+1) * DP%dRdEsi0(+1,:)                        &
                + alphasi( 0) * DP%dRdEsi0( 0,:)                        &
                + alphasi(-1) * DP%dRdEsi0(-1,:)
  DP%dRdEsd(:) =  alphasd(+1) * DP%dRdEsd0(+1,:)                        &
                + alphasd( 0) * DP%dRdEsd0( 0,:)                        &
                + alphasd(-1) * DP%dRdEsd0(-1,:)
  DP%dRdE(:)   = DP%dRdEsi(:) + DP%dRdEsd(:)
  DP%Rsi(:) =  alphasi(+1) * DP%Rsi0(+1,:)                              &
             + alphasi( 0) * DP%Rsi0( 0,:)                              &
             + alphasi(-1) * DP%Rsi0(-1,:)
  DP%Rsd(:) =  alphasd(+1) * DP%Rsd0(+1,:)                              &
             + alphasd( 0) * DP%Rsd0( 0,:)                              &
             + alphasd(-1) * DP%Rsd0(-1,:)
  DP%R(:)   = DP%Rsi(:) + DP%Rsd(:)
  
  ! Average expected event (components at reference couplings)
  DP%MuSignalSI0(:,:) = DP%exposure * DP%Rsi0(:,:)
  DP%MuSignalSD0(:,:) = DP%exposure * DP%Rsd0(:,:)
  
  ! Average expected events.
  ! Could also achieve this by e.g. exposure * Rsi(:).
  DP%MuSignalSI(:) =  alphasi(+1) * DP%MuSignalSI0(+1,:)                &
                    + alphasi( 0) * DP%MuSignalSI0( 0,:)                &
                    + alphasi(-1) * DP%MuSignalSI0(-1,:)
  DP%MuSignalSD(:) =  alphasd(+1) * DP%MuSignalSD0(+1,:)                &
                    + alphasd( 0) * DP%MuSignalSD0( 0,:)                &
                    + alphasd(-1) * DP%MuSignalSD0(-1,:)
  DP%MuSignal(:)   = DP%MuSignalSI(:) + DP%MuSignalSD(:)
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Returns the observed number of events.
! 
! Required input argument:
!   D           A DetectorStruct containing detector info.
! 
FUNCTION GetEvents(D) RESULT(N)
  IMPLICIT NONE
  INTEGER :: N
  TYPE(DetectorStruct), INTENT(IN) :: D
  CALL GetRates(D,Nevents=N)
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of background events.
! 
! Required input argument:
!   D           A DetectorStruct containing detector info.
! 
FUNCTION GetBackground(D) RESULT(b)
  IMPLICIT NONE
  REAL*8 :: b
  TYPE(DetectorStruct), INTENT(IN) :: D
  CALL GetRates(D,background=b)
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of signal events for the
! current WIMP.
! 
! Required input argument:
!   D           A DetectorStruct containing rates (CalcRates(D) must
!               have already been called).
! 
FUNCTION GetSignal(D) RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  TYPE(DetectorStruct), INTENT(IN) :: D
  CALL GetRates(D,signal=s)
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of spin-independent signal events
! for the current WIMP.
! 
! Required input argument:
!   D           A DetectorStruct containing rates (CalcRates(D) must
!               have already been called).
! 
FUNCTION GetSignalSI(D) RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  TYPE(DetectorStruct), INTENT(IN) :: D
  CALL GetRates(D,signal_si=s)
END FUNCTION


! ----------------------------------------------------------------------
! Returns the average expected number of spin-dependent signal events
! for the current WIMP.
! 
! Required input argument:
!   D           A DetectorStruct containing rates (CalcRates(D) must
!               have already been called).
! 
FUNCTION GetSignalSD(D) RESULT(s)
  IMPLICIT NONE
  REAL*8 :: s
  TYPE(DetectorStruct), INTENT(IN) :: D
  CALL GetRates(D,signal_sd=s)
END FUNCTION



!=======================================================================
! LIKELIHOOD/P-VALUE ROUTINES
!=======================================================================

! ----------------------------------------------------------------------
! Calculates the log-likelihood for the current WIMP mass and couplings.
! Uses a Poisson distribution in the number of observed events N:
!    P(N|s+b)
! where s is the average expected signal and b is the average expected
! background.
! 
! Optional input argument:
!   D           The DetectorStruct structure containing detector
!               and rate date to calculate the likelihood for.  If
!               not given, a default structure (internally stored)
!               will be used.
! 
FUNCTION LogLikelihood(D) RESULT(lnlike)
  IMPLICIT NONE
  REAL*8 :: lnlike
  TYPE(DetectorStruct), INTENT(IN), TARGET, OPTIONAL :: D
  TYPE(DetectorStruct), POINTER :: DP
  INTEGER :: N
  REAL*8 :: b,s,mu
  
  ! Get pointer to detector structure to be used: either an explicitly
  ! provided structure or the default, internally-stored structure.
  IF (PRESENT(D)) THEN
    DP => D
  ELSE
    DP => DefaultDetector
  END IF
  
  ! Get observed events and expected events
  N = DP%Nevents
  b = DP%MuBackground
  s = DP%MuSignal(0)
  mu = s + b
  
  ! Poisson likelihood (routine handles special cases)
  lnlike = PoissonLogPDF(N,mu)
  
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the log of the p-value for the current WIMP mass and
! couplings (NO BACKGROUND SUBTRACTION).  Uses the maximum gap method
! if efficiencies were provided in the efficiencies file for the
! intervals, otherwise uses a Poisson distribution in the number of
! observed events N:
!    P(N|s)
! where s is the average expected signal (background contributions are
! ignored).
! 
! Optional input argument:
!   D           The DetectorStruct structure containing detector
!               and rate date to calculate the p-value for.  If not
!               given, a default structure (internally stored) will
!               be used.
! 
FUNCTION LogPValue(D) RESULT(lnp)
  IMPLICIT NONE
  REAL*8 :: lnp
  TYPE(DetectorStruct), INTENT(IN), TARGET, OPTIONAL :: D
  TYPE(DetectorStruct), POINTER :: DP
  INTEGER :: N,I
  REAL*8 :: s,mu
  
  ! Get pointer to detector structure to be used: either an explicitly
  ! provided structure or the default, internally-stored structure.
  IF (PRESENT(D)) THEN
    DP => D
  ELSE
    DP => DefaultDetector
  END IF
  
  ! Get observed events and expected events
  N  = DP%Nevents
  mu = DP%MuSignal(0)
  
  ! Check if rates are available for each interval
  IF (DP%intervals .AND. (DP%Neff .EQ. DP%Nevents+1)) THEN
    lnp = LogMaximumGapP(mu,MAXVAL(DP%MuSignal(1:N+1)))
  ELSE
    lnp = LogPoissonP(N,mu)
  END IF
  
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the factor x by which the cross-sections must be scaled
! (sigma -> x*sigma) to achieve the desired p-value (given as log(p)).
! See LogPValue() above for a description of the statistics.
! 
! Optional input arguments:
!   D           The DetectorStruct structure containing detector
!               and rate date to calculate the scaling for.  If not
!               given, a default structure (internally stored) will
!               be used.
!   lnp         The logarithm of the desired p-value.  Default is
!               p=0.1 (90% CL).
! 
FUNCTION ScaleToPValue(D,lnp) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  TYPE(DetectorStruct), INTENT(IN), TARGET, OPTIONAL :: D
  REAL*8, INTENT(IN), OPTIONAL :: lnp
  TYPE(DetectorStruct), POINTER :: DP
  INTEGER :: N
  REAL*8 :: lnp0,mu,f
  
  ! Get pointer to detector structure to be used: either an explicitly
  ! provided structure or the default, internally-stored structure.
  IF (PRESENT(D)) THEN
    DP => D
  ELSE
    DP => DefaultDetector
  END IF
  
  ! logarithm of p-value to use
  IF (PRESENT(lnp)) THEN
    lnp0 = lnp
  ELSE
    lnp0 = LOG(0.1d0)
  END IF
  
  ! Get observed events and expected events
  N  = DP%Nevents
  mu = DP%MuSignal(0)
  
  IF (mu .LE. 0d0) THEN
    x = HUGE(1d0)
    RETURN
  END IF
  
  IF (lnp0 .GE. 0d0) THEN
    x = 0d0
    RETURN
  END IF
  
  ! Use maximum gap if intervals available, otherwise Poisson
  IF (DP%intervals .AND. (DP%Neff .EQ. DP%Nevents+1)) THEN
    f = MAXVAL(DP%MuSignal(1:N+1)) / mu
    x = MaximumGapScaleToPValue(lnp0,mu,f)
  ELSE
    x = PoissonScaleToPValue(lnp0,N,mu)
  END IF
  
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the factor x by which the signal must be scaled (s -> x*s)
! to achieve the desired p-value (given as log(p)) using the maximum
! gap method.
! 
! NOTE: Uses quickly implemented numerical methods.  Could be much
!       better optimized if the need arises.
! 
! Input arguments:
!   lnp         The logarithm of the desired p-value.
!   mu          The total expected number of events.
!   f           The largest fraction of events expected in any interval
!               (i.e. the "maximum" gap).
! 
PURE FUNCTION MaximumGapScaleToPValue(lnp,mu,f) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN) :: lnp,mu,f
  REAL*8 :: x1,lnp1,x2,lnp2,xm,lnpm
  
  IF (mu .LE. 0d0) THEN
    x = HUGE(1d0)
    RETURN
  END IF
  
  ! Should not happen...
  IF (lnp .GE. 0d0) THEN
    x = 0d0
    RETURN
  END IF
  
  ! Should not happen...
  IF (f .LE. 0d0) THEN
    x = HUGE(1d0)
    RETURN
  END IF
  
  ! Special case
  IF (f .GE. 1d0) THEN
    x = -lnp / mu
    RETURN
  END IF
  
  ! Starting point
  x1   = 1d0 / mu
  lnp1 = LogMaximumGapP(x1*mu,x1*f*mu)
  ! Bracket
  IF (lnp1 .GT. lnp) THEN
    x2   = 2d0*x1
    lnp2 = LogMaximumGapP(x2*mu,x2*f*mu)
    DO WHILE (lnp2 .GE. lnp)
      x1   = x2
      lnp1 = lnp2
      x2   = 2d0*x2
      lnp2 = LogMaximumGapP(x2*mu,x2*f*mu)
    END DO
  ELSE
    DO WHILE (lnp1 .LE. lnp)
      x2   = x1
      lnp2 = lnp1
      x1   = 0.5d0*x1
      lnp1 = LogMaximumGapP(x1*mu,x1*f*mu)
    END DO
  END IF
  
  ! Bisection (geometric)
  DO WHILE ((ABS(lnp2-lnp1) .GT. 1d-5) .AND. (ABS(LOG(x2/x1)) .GT. 1d-5))
    xm   = SQRT(x1*x2)
    lnpm = LogMaximumGapP(xm*mu,xm*f*mu)
    IF (lnpm .GE. lnp) THEN
      x1   = xm
      lnp1 = lnpm
    ELSE
      x2   = xm
      lnp2 = lnpm
    END IF
  END DO
  x = 0.5d0*(x1+x2)
  
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the factor x by which the signal must be scaled (s -> x*s)
! to achieve the desired p-value (given as log(p)) using a Poisson
! distribution.
! 
! NOTE: Uses quickly implemented numerical methods.  Could be much
!       better optimized if the need arises.
! 
! Input arguments:
!   lnp         The logarithm of the desired p-value.
!   N           Number of observed events.
!   mu          The total expected number of events.
! 
PURE FUNCTION PoissonScaleToPValue(lnp,N,mu) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: lnp,mu
  REAL*8 :: x1,lnp1,x2,lnp2,xm,lnpm
  
  IF (mu .LE. 0d0) THEN
    x = HUGE(1d0)
    RETURN
  END IF
  
  ! Should not happen...
  IF (lnp .GE. 0d0) THEN
    x = 0d0
    RETURN
  END IF
  
  ! Analytic formula for N=0
  IF (N .EQ. 0) THEN
    x = -lnp / mu
    RETURN
  END IF
  
  ! Starting point
  x1   = N / mu
  lnp1 = LogPoissonP(N,x1*mu)
  
  ! Bracket
  IF (lnp1 .GT. lnp) THEN
    x2   = 2d0*x1
    lnp2 = LogPoissonP(N,x2*mu)
    DO WHILE (lnp2 .GE. lnp)
      x1   = x2
      lnp1 = lnp2
      x2   = 2d0*x2
      lnp2 = LogPoissonP(N,x2*mu)
    END DO
  ELSE
    DO WHILE (lnp1 .LE. lnp)
      x2   = x1
      lnp2 = lnp1
      x1   = 0.5d0*x1
      lnp1 = LogPoissonP(N,x1*mu)
    END DO
  END IF
  
  ! Bisection (geometric)
  DO WHILE ((ABS(lnp2-lnp1) .GT. 1d-5) .AND. (ABS(x2-x1) .GT. 1d-5))
    xm   = SQRT(x1*x2)
    lnpm = LogPoissonP(N,xm*mu)
    IF (lnpm .GE. lnp) THEN
      x1   = xm
      lnp1 = lnpm
    ELSE
      x2   = xm
      lnp2 = lnpm
    END IF
  END DO
  x = 0.5d0*(x1+x2)
  
END FUNCTION


! ----------------------------------------------------------------------
! Calculates the confidence interval [s1,s2] for the signal
! contribution s in a Poisson distribution with background expectation
! b and observed events N.  Calculated at the CL corresponding to the
! given p-value (p = 1-CL).  Uses the Feldman-Cousins method as
! described in:
!    G. Feldman & R. Cousins, Phys. Rev. D 57, 3873 (1998) [physics/9711021]
! 
! Input arguments:
!   lnp         The logarithm of the p-value (CL = 1-p)
!   N           Number of observed events
!   b           Mean number of expected background events
! Output arguments:
!   s1,s2       The signal contribution confidence interval [s1,s2]
! 
PURE SUBROUTINE FeldmanCousinsPoissonCI(lnp,N,b,s1,s2)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: lnp,b
  REAL*8, INTENT(OUT) :: s1,s2
  REAL*8 :: step,sa,sb,x,lnsum1,lnsum2
  ! Specifies the desired precision for boundary searches.
  ! This should not be smaller than EPSILON(1d0)!
  REAL*8, PARAMETER :: RELATIVE_PRECISION = 100*EPSILON(1d0)
  
  ! It is possible to get a measure zero confidence interval,
  ! though that is only possible if the CL is low (61% or lower,
  ! safely below the 1-sigma CL = 68.3%).  The measure zero
  ! interval will occur if the following is true:
  !   \Sum_{k=N+1}^{floor(b)} P(k|b) > 1-p = CL
  ! in which case the interval is [0,0].  The l.h.s. is bounded from
  ! above by Q(floor(b),b) - P(0|b), the sum for N=0 (Q is the
  ! regularized incomplete gamma function) which has a global
  ! maximum of 0.611 at b=4.0.  Other maxima occur at integer values
  ! of b, with maxima only slowly changing (goes to 0.5 as b becomes
  ! large).  Catch this case here:
  IF (lnp .GT. LOG(0.38d0)) THEN
    CALL LogPoissonSums(N,b,lnsum1,x)
    CALL LogPoissonSums(INT(b+EPSILON(1d0)),b,lnsum2,x)
    IF (EXP(lnsum2)-EXP(lnsum1) .GE. (1-EXP(lnp))) THEN
      s1 = 0d0
      s2 = 0d0
      RETURN
    END IF
  END IF
  
  ! Get lower boundary
  IF (accept(N,b,0d0,lnp)) THEN
    ! This should include b > N
    s1 = 0d0
  ELSE
    ! This _should_ be an accepted point
    s1 = N - b
    ! Take s1 -> s1/2 until s1 no longer accepted
    DO WHILE (accept(N,b,s1,lnp))
      s1 = 0.5d0*s1
    END DO
    step = s1
    ! Boundary is now between [s1,s1+step]
    DO WHILE (step .GE. RELATIVE_PRECISION*s1)
      step = 0.5d0*step
      IF (.NOT. accept(N,b,s1+step,lnp)) s1 = s1 + step
    END DO
  END IF
  
  ! Get upper boundary
  ! Need at least one good starting point
  IF (N .GT. b) THEN
    ! This _should_ be an accepted point
    s2 = N - b
  ELSE
    s2 = 1d0
    DO WHILE (.NOT. accept(N,b,s2,lnp))
      s2 = 0.5d0*s2
    END DO
  END IF
  ! Take s2 -> 2*s2 until s2 no longer accepted
  DO WHILE (accept(N,b,s2,lnp))
    s2 = 2*s2
  END DO
  step = s2
  ! Boundary is now between [s2-step,s2]
  DO WHILE (step .GE. RELATIVE_PRECISION*s2)
    step = 0.5d0*step
    IF (.NOT. accept(N,b,s2-step,lnp)) s2 = s2 - step
  END DO
  
  
  CONTAINS
  
  ! ------------------------------------------
  ! Determines if the given N is within the acceptance
  ! region for a Poisson of mean b+s with given CL=1-p.
  PURE FUNCTION accept(N,b,s,lnp)
    IMPLICIT NONE
    LOGICAL :: accept
    INTEGER, INTENT(IN) :: N
    REAL*8, INTENT(IN) :: b,s,lnp
    INTEGER :: Ktail,K1,K2
    REAL*8 :: x,lnP1,lnP2,lnR1,lnR2,lnPsum
    
    accept = .TRUE.
    
    ! Below, we are going to determine the probability outside
    ! some interval [K1,K2] (Psum) by starting with an interval
    ! of the form [0,K2] and working our way inward until the
    ! probability exceeds the p-value or N falls outside the
    ! interval.  Note this means we are adding terms in order of
    !  _increasing_ R (the Feldman-Cousins ordering parameter).
    
    ! To avoid ambiguity, terms of equal R are either accepted
    ! or rejected from the interval as one.  To avoid under-
    ! coverage, we accept such a group if necessary to ensure the
    ! interval covers at least 1-p, leading to (conservative)
    ! overcoverage.  There is a case (s=0) where the R has a
    ! plateau over a range of K that must be handled carefully.
    ! The only other possibility of equal R is by coincidence for
    ! two terms on opposite sides of the Poisson peak; this case
    ! is handled appropriately.
    
    ! Special case: s = 0 with N <= b
    ! For the s = 0 case, R=1 for all K <= b with R < 1 for
    ! K > b.  By our construction, that means that [0,floor(b)]
    ! is _always_ part of the accepted interval.  We catch this
    ! case here because the later algorithms will not handle
    ! it correctly (the N > b case is handled correctly below).
    IF ((s .EQ. 0d0) .AND. (N .LE. b)) THEN
      accept = .TRUE.
      RETURN
    END IF
    
    ! Special case: N = b+s
    ! R is maximized at K=b+s, so this is the first accepted term.
    ! Though the code below handles this correctly, we avoid that
    ! possibly long calculation for what is a common case.
    IF (ABS(b+s-N) .LT. 0.4d0) THEN
      accept = .TRUE.
      RETURN
    END IF
    
    ! For K >= Ktail, ordering parameter R is smaller than
    ! for any K < Ktail and is asymptotically decreasing.
    Ktail = FClnRsearch0(b,s) + 1
    
    ! Special case: N is in the tail region.
    ! Just have to check if probability for K >= N exceeds p-value.
    IF (N .GE. Ktail) THEN
      CALL LogPoissonSums(N,b+s,x,lnPsum)
      accept = (lnPsum .GT. lnp)
      RETURN
    END IF
    
    ! Get area of tail region
    CALL LogPoissonSums(Ktail,b+s,x,lnPsum)
    
    ! Special case: acceptance region contains at least [0,Ntail-1].
    ! By our above check, N is in [0,Ntail-1]
    IF (lnPsum .GT. lnp) THEN
      accept = .TRUE.
      RETURN
    END IF
    
    ! Now we narrow the interval from both sides using the
    ! ordering parameter.  We narrow until N falls outside
    ! the interval or the interval no longer covers enough
    ! probability (i.e. area outside the interval exceeds
    ! p-value).  We treats terms of equal R as one in terms
    ! of acceptance: this is handled appropriately below
    ! except for the case where R(K) is flat, which only
    ! occurs for the s=0 case already handled above.
    !
    K1   = 0
    lnR1 = FClnR(K1,b,s)
    K2   = Ktail-1
    lnR2 = FClnR(K2,b,s)
    ! Psum is probability outside of [K1,K2]
    DO WHILE ((lnPsum .LT. lnp))
      ! See if N is no longer in the valid region
      IF ((N .LT. K1) .OR. (N .GT. K2)) THEN
        accept = .FALSE.
        RETURN
      END IF
      ! Find term with _smaller_ R value: add its probability
      ! to Psum and drop from interval.  Note special case where
      ! both terms have same R value (drop both).
      IF (lnR1 .LT. lnR2) THEN
        lnPsum = LOG_SUM(lnPsum,LogPoisson(K1,b+s))
        K1     = K1+1
        lnR1   = FClnR(K1,b,s)
      ELSE IF (lnR1 .GT. lnR2) THEN
        lnPsum = LOG_SUM(lnPsum,LogPoisson(K2,b+s))
        K2     = K2-1
        lnR2   = FClnR(K2,b,s)
      ELSE
        lnPsum = LOG_SUM(lnPsum,LogPoisson(K1,b+s))
        K1     = K1+1
        lnR1   = FClnR(K1,b,s)
        lnPsum = LOG_SUM(lnPsum,LogPoisson(K2,b+s))
        K2     = K2-1
        lnR2   = FClnR(K2,b,s)
      END IF
    END DO
    
    ! If we made it here, we narrowed the interval too far before
    ! losing the K=N term.  N must be within the proper interval.
    accept = .TRUE.
    
  END FUNCTION
  
  ! ------------------------------------------
  ! Calculates the logarithm of the Poisson probability
  !   P(k|mu) = e^{-mu} mu^k / k!
  PURE FUNCTION LogPoisson(K,mu) RESULT(lnP)
    IMPLICIT NONE
    REAL*8 :: lnP
    INTEGER, INTENT(IN) :: K
    REAL*8, INTENT(IN) :: mu
    IF (K .EQ. 0) THEN
      lnP = -mu
    ELSE IF (mu .EQ. 0d0) THEN
      lnP = -HUGE(1d0)
    ELSE
      lnP = -mu + K*LOG(mu) - LOG_GAMMA(K+1d0)
    END IF
  END FUNCTION
  
  ! ------------------------------------------
  ! Calculates ln(r), where the ratio
  !   R = P(K|b+s)/P(K|b+s0(K))
  ! is the ordering term used by Feldman-Cousins, with
  !   s0(K) = max(0,K-b)
  ! the best-fit signal.
  PURE FUNCTION FClnR(K,b,s) RESULT(lnR)
    IMPLICIT NONE
    REAL*8 :: lnR
    INTEGER, INTENT(IN) :: K
    REAL*8, INTENT(IN) :: b,s
    IF (K .EQ. 0) THEN
      lnR = -s
    ELSE IF (b+s .LE. 0d0) THEN
      lnR = -HUGE(1d0)
    ELSE IF ( K .LE. b) THEN
      lnR = K*LOGp1(s/b) - s
    ELSE
      lnR = K*LOG((b+s)/K) + K - (b+s)
    END IF
  END FUNCTION
  
  ! ------------------------------------------
  ! Finds the index K such that the Feldman-Cousins
  ! ratio R is maximized for the given b & s.
  PURE FUNCTION FCpeakRloc(b,s) RESULT(K)
    IMPLICIT NONE
    INTEGER :: K
    REAL*8, INTENT(IN) :: b,s
    ! When analytically continuing R(k) for continuous k,
    ! there is a single maximum located at k=b+s.  For the
    ! discrete case, it must be at either the floor or
    ! ceiling of b+s.
    K = MAX(INT(b+s),0)
    IF (FClnR(K+1,b,s) .GT. FClnR(K,b,s)) K = K+1
  END FUNCTION
  
  ! ------------------------------------------
  ! Finds the largest index K such that the Feldman-Cousins
  ! ratio R is as large as R(K=0) for the given b & s.
  PURE FUNCTION FClnRsearch0(b,s) RESULT(K)
    IMPLICIT NONE
    INTEGER :: K
    REAL*8, INTENT(IN) :: b,s
    INTEGER :: K1,K2,Km,Kstep
    REAL*8 :: lnR0,lnR1,lnR2,lnRm
    Kstep = 1
    K1 = FCpeakRloc(b,s)
    K2 = K1+Kstep
    lnR0 = FClnR(0,b,s)
    lnR1 = FClnR(K1,b,s)
    lnR2 = FClnR(K2,b,s)
    ! Bracket
    DO WHILE (lnR2 .GE. lnR0)
      Kstep = 2*Kstep
      K1    = K2
      lnR1  = lnR2
      K2    = K2+Kstep
      lnR2 = FClnR(K2,b,s)
    END DO
    ! Bisection
    DO WHILE (K2-K1 .GT. 1)
      Km   = K1 + (K2-K1)/2
      lnRm = FClnR(Km,b,s)
      IF (lnRm .GE. lnR0) THEN
        K1   = Km
        lnR1 = lnRm
      ELSE
        K2   = Km
        lnR2 = lnRm
      END IF
    END DO
    K = K1
  END FUNCTION
  
END SUBROUTINE


! ----------------------------------------------------------------------
! For the Poisson distribution P(k|mu), calculates the log of the
! p-value, where p is defined as:
!   p = \Sum_{k\le N} P(k|mu)
! 
PURE FUNCTION LogPoissonP(N,mu) RESULT(lnp)
  IMPLICIT NONE
  REAL*8 :: lnp
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: mu
  REAL*8 :: lnlesum,lngesum
  
  ! Utility routine calculates P(k|mu) for k\le N and k\ge N.
  CALL LogPoissonSums(N,mu,lnlesum,lngesum)
  lnp = lnlesum
  
END FUNCTION


! ----------------------------------------------------------------------
! For the Poisson distribution P(k|mu), calculates the log of the sums:
!   lesum = \Sum_{k\le N} P(k|mu)
!   gesum = \Sum_{k\ge N} P(k|mu)
! This routine calculates the smaller of the two and uses the relation
!   lesum + gesum = 1 + P(N|mu)
! to find the other (note the k=N term appears in both sums).
! 
! Input arguments:
!   N           Maximum/minimum k in sum
!   mu          Mean of the Poisson
! Output arguments:
!   lnlesum     The logarithm of the sum of Poisson terms k = 0,1,...,N.
!   lngesum     The logarithm of the sum of Poisson terms k = N,N+1,....
! 
PURE SUBROUTINE LogPoissonSums(N,mu,lnlesum,lngesum)
  IMPLICIT NONE
  REAL*8 :: lnp
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: mu
  REAL*8, INTENT(OUT) :: lnlesum,lngesum
  INTEGER :: K,Kmax
  REAL*8 :: lnpmf0,lnpmf
  REAL*8, PARAMETER :: LN_PRECISION = -35d0  ! Precision of ~ 1d-15
  
  ! Special cases
  IF (mu .LE. 0d0) THEN
    lnlesum = 0d0
    IF (N .EQ. 0) THEN
      lngesum = 0d0
    ELSE
      lngesum = -HUGE(1d0)
    END IF
    RETURN
  ELSE IF (N .EQ. 0d0) THEN
    lnlesum = -mu
    lngesum = 0d0
    RETURN
  END IF
  
  ! General case
  lnpmf = N*LOG(mu) - mu - LOG_GAMMA(N+1d0)
  lnpmf0 = lnpmf
  
  ! NOTE: This algorithm calculates the log of the sum of both
  !       all K <= N (lnlesum) and all K >= N (lngesum).
  ! The distribution peaks around n = mu; to avoid a loss of
  ! precision, we will do the lower sum for n < mu and the
  ! upper sum for n > mu.
  IF (N .LE. mu) THEN
    K = N
    lnlesum = lnpmf
    DO WHILE (k .GT. 0)
      K = K-1
      lnpmf = lnpmf + LOG((K+1)/mu)
      lnlesum = LOG_SUM(lnlesum,lnpmf)
      IF (lnpmf - lnlesum .LE. LN_PRECISION) EXIT
    END DO
    lngesum = LOG(1d0 - (EXP(lnlesum) - EXP(lnpmf0)))
  ELSE
    K = N
    lngesum = lnpmf
    ! Determine a conservative upper limit on sum terms.
    ! Will hopefully hit precision condition first, but include
    ! this upper limit in case we didn't think of some pathological
    ! cases.
    IF ((N-mu)**2 .GT. mu) THEN
      Kmax = N + NINT(-LN_PRECISION * (N/(N-mu))) + 10
    ELSE
      Kmax = N + NINT(-LN_PRECISION * (1d0+SQRT(mu))) + 10
    END IF
    DO WHILE (K .LT. Kmax)
      K = K+1
      lnpmf = lnpmf + LOG(mu/K)
      lngesum = LOG_SUM(lngesum,lnpmf)
      IF (lnp - lngesum .LE. LN_PRECISION) EXIT
    END DO
    lnlesum = LOG(1d0 - (EXP(lngesum) - EXP(lnpmf0)))
  END IF
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Calculates the log of the maximum gap statistic.  See:
!   S. Yellin, Phys. Rev. D 66, 032005 (2002) [arXiv:physics/0203002]
! This statistic is an upper bound on the p-value in the event the
! background spectrum is unknown.  The p-value is equal to 1-C_0 and
! is given by:
!   p = \sum_{k=1}^{floor(mu/x)} (k*x-mu)^(k-1) e^(-k*x) (mu-k(x-1)) / k!
! 
! NOTE: IMPLEMENTATION NEEDS TESTING.
! 
! Input arguments:
!   mu          Average expected number of signal events
!   x           Largest expected number of signal events in any
!               interval/gap.
! 
PURE FUNCTION LogMaximumGapP(mu,x) RESULT(lnp)
  IMPLICIT NONE
  REAL*8 :: lnp
  REAL*8, INTENT(IN) :: mu,x
  INTEGER :: K,Kmax
  REAL*8 :: p,psum,z
  
  ! Bad cases
  IF ((mu .LT. 0d0) .OR. (x .LT. 0d0)) THEN
    lnp = HUGE(1d0)
    RETURN
  END IF
  
  ! Special cases
  IF (mu .EQ. 0d0) THEN
    lnp = 0d0
    RETURN
  ELSE IF (x .GE. mu) THEN
    lnp = -mu
    RETURN
  ELSE IF (x .EQ. 0d0) THEN
    lnp = 0d0
    RETURN
  END IF
  
  ! Not optimized for log case or any extreme case
  ! WARNING: CONDITIONS AND ALGORITHMS NOT FULLY TESTED.
  
  IF (mu*EXP(-x) .GT. 12.5d0) THEN
    lnp = 0d0
    RETURN
  END IF
  
  ! Calculation involves sum from 1 to Kmax = floor(mu/x).  By
  ! definition, Kmax should never exceed the number of intervals.
  ! We handle the lowest Kmax cases explicitly and do the sum for
  ! larger Kmax.
  
  ! Small epsilon added to avoid bad truncation results
  ! when mu/x is an integer.
  Kmax = INT(mu/x+SQRT(EPSILON(1d0)))
  SELECT CASE (Kmax)
  CASE (:0)
    lnp = HUGE(1d0)
  CASE (1)
    ! p = (mu-x+1) e^-x
    ! Note LOGp1(z) = ln(1+z)
    lnp = -x + LOGp1(mu-x)
  CASE (2)
    ! p = (mu-x+1) e^-x [1 - 1/2 (mu-2x) e^(-x) (1 - (x-1)/(mu-(x-1)))]
    ! maximum z in range of validity is 0.1925 (safe here)
    z = 0.5d0 * (mu-2*x) * EXP(-x) * (mu-2*(x-1)) / (mu-x+1)
    lnp = -x + LOGp1(mu-x) + LOGp1(-z)
  CASE (3:)
    ! Do explicit sum.
    ! Note finite alternating series makes accuracy/precision
    ! and convergence determinations difficult.
    psum = 0d0
    DO K = 1,Kmax
      p = (K*x-mu)**(K-1) * EXP(-K*x) * (mu - K*(x-1)) / GAMMAI(K+1)
      psum = psum + p
    END DO
    lnp = LOG(psum)
  END SELECT
  
END FUNCTION



!=======================================================================
! USAGE
!=======================================================================

! ----------------------------------------------------------------------
! Checks if the command line contains any requests to show the program
! usage.
! 
FUNCTION ShowUsageRequested() RESULT(flag)
  IMPLICIT NONE
  LOGICAL :: flag
  INTEGER :: I,Narg
  CHARACTER*32 :: arg
  
  flag = .FALSE.
  Narg = IARGC()
  
  DO I=1,Narg
    CALL GETARG(I,arg)
    IF (arg(1:1) .EQ. '?') THEN
      flag = .TRUE.
      RETURN
    ELSE IF ( (arg(1:2) .EQ. '-?') .OR. (arg(1:2) .EQ. '-h') ) THEN
      flag = .TRUE.
      RETURN
    ELSE IF (arg(1:6) .EQ. '--help') THEN
      flag = .TRUE.
      RETURN
    END IF
  END DO
  
END FUNCTION ShowUsageRequested


! ----------------------------------------------------------------------
! Shows usage.
! Prints out contents of given file.
! 
! Optional input argument:
!   file        File containing usage contents.  Default is the
!               calling program's name with '.use' appended.
!   
SUBROUTINE ShowUsage(file)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN), OPTIONAL :: file
  CHARACTER*1024 :: file0
  INTEGER, PARAMETER :: USAGE_FID = 19
  CHARACTER*1024 :: line
  INTEGER :: ios,K
  
  IF (PRESENT(file)) THEN
    file0 = file
  ELSE
    CALL GETARG(0,file0)
    K = INDEX(file0,'/',.TRUE.)
    IF (K .NE. 0) file0 = file0(K+1:)
    file0 = TRIM(file0) // '.use'
  END IF
  
  OPEN(FILE=file0,UNIT=USAGE_FID,STATUS='OLD',                          &
       FORM='FORMATTED',ACTION='READ',IOSTAT=ios)
  IF (ios .NE. 0) THEN
    WRITE(0,*) 'Unable to display usage:'
    WRITE(0,*) 'Could not open usage file ' // TRIM(file0)
    WRITE(0,*)
    RETURN
  END IF
  
  DO WHILE (.TRUE.)
    READ(UNIT=USAGE_FID,FMT='(A)',IOSTAT=ios) line
    IF (ios .LT. 0) EXIT
    WRITE(*,'(A)') TRIM(line)
  END DO
  
  CLOSE(UNIT=USAGE_FID,IOSTAT=ios)
  
END SUBROUTINE ShowUsage



!=======================================================================
! COMMAND LINE PROCESSING
! Several types:
!   *) GetXXXArg(index,status)
!      Parses the command line argument at the given index and returns
!      as type XXX (logical,integer,real,complex).  A negative index
!      will count from the end of the argument list.  The optional
!      'status' parameter indicates if the argument exists and was
!      properly parsed.
!   *) GetShortArg(akey)
!      Indicates if the given akey appears in any argument of the form
!      -<akey> or -<akey1><akey2><akey3>....  In this case, akey
!      is meant to be a single character.
!   *) GetLongArg(akey)
!      Indicates if the given akey appears in any argument of the form
!      --<akey> or -<akey>=<aval>.
!   *) GetLongArgXXX(akey,aval)
!      Indicates if the given akey appears in any argument of the form
!      --<akey> or -<akey>=<aval>.  If present, <aval> will be parsed
!      as type XXX (string,logical,integer,real,complex) and placed in
!      aval.  The first instance of a particular key will always be
!      parsed; if additional arguments use the same key, they will be
!      ignored.  This routine returns .FALSE. if there is a failure to
!      parse <aval>.
!=======================================================================

! ----------------------------------------------------------------------
! Obtains the command line argument at the given index as a LOGICAL.
! Negative index counts from the end of the argument list.
! The optional status argument indicates if the given argument exists
! and was properly parsed.
! 
FUNCTION GetLogicalArg(index,status) RESULT(aval)
  IMPLICIT NONE
  LOGICAL :: aval
  INTEGER, INTENT(IN) :: index
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  CHARACTER*256 :: arg
  INTEGER :: index0,Narg,ios
  
  aval = .FALSE.
  IF (PRESENT(status)) status = .FALSE.
  
  ! Check index validity, count from end if necessary
  Narg   = IARGC()
  IF (index .LT. 0) THEN
    index0 = Narg + 1 + index
  ELSE
    index0 = index
  END IF
  IF ((index0 .LT. 1) .OR. (index0 .GT. Narg)) RETURN
  
  ! Parse argument
  CALL GETARG(index0,arg)
  READ(UNIT=arg,FMT=*,IOSTAT=ios) aval
  IF (ios .EQ. 0) THEN
    IF (PRESENT(status)) status = .TRUE.
  END IF
  
END FUNCTION GetLogicalArg


! ----------------------------------------------------------------------
! Obtains the command line argument at the given index as an INTEGER.
! Negative index counts from the end of the argument list.
! The optional status argument indicates if the given argument exists
! and was properly parsed.
! Accepts real valued strings, which will be rounded to the nearest
! integer.
! 
FUNCTION GetIntegerArg(index,status) RESULT(aval)
  IMPLICIT NONE
  INTEGER :: aval
  INTEGER, INTENT(IN) :: index
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  CHARACTER*256 :: arg
  INTEGER :: index0,Narg,ios
  REAL*8 :: x
  
  aval = -HUGE(aval)
  IF (PRESENT(status)) status = .FALSE.
  
  ! Check index validity, count from end if necessary
  Narg   = IARGC()
  IF (index .LT. 0) THEN
    index0 = Narg + 1 + index
  ELSE
    index0 = index
  END IF
  IF ((index0 .LT. 1) .OR. (index0 .GT. Narg)) RETURN
  
  ! Parse argument
  ! Try reading as integer first, then as real and convert
  ! to integer
  CALL GETARG(index0,arg)
  READ(UNIT=arg,FMT=*,IOSTAT=ios) aval
  IF (ios .EQ. 0) THEN
    IF (PRESENT(status)) status = .TRUE.
  ELSE
    READ(UNIT=arg,FMT=*,IOSTAT=ios) x
    IF (ios .EQ. 0) THEN
      IF (PRESENT(status)) status = .TRUE.
      aval = NINT(x)
    END IF
  END IF
  
END FUNCTION GetIntegerArg


! ----------------------------------------------------------------------
! Obtains the command line argument at the given index as a REAL*8.
! Negative index counts from the end of the argument list.
! The optional status argument indicates if the given argument exists
! and was properly parsed.
! 
FUNCTION GetRealArg(index,status) RESULT(aval)
  IMPLICIT NONE
  REAL*8 :: aval
  INTEGER, INTENT(IN) :: index
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  CHARACTER*256 :: arg
  INTEGER :: index0,Narg,ios
  
  aval = -HUGE(aval)
  IF (PRESENT(status)) status = .FALSE.
  
  ! Check index validity, count from end if necessary
  Narg   = IARGC()
  IF (index .LT. 0) THEN
    index0 = Narg + 1 + index
  ELSE
    index0 = index
  END IF
  IF ((index0 .LT. 1) .OR. (index0 .GT. Narg)) RETURN
  
  ! Parse argument
  CALL GETARG(index0,arg)
  READ(UNIT=arg,FMT=*,IOSTAT=ios) aval
  IF (ios .EQ. 0) THEN
    IF (PRESENT(status)) status = .TRUE.
  END IF
  
END FUNCTION GetRealArg


! ----------------------------------------------------------------------
! Obtains the command line argument at the given index as a COMPLEX*16.
! Negative index counts from the end of the argument list.
! The optional status argument indicates if the given argument exists
! and was properly parsed.
! 
FUNCTION GetComplexArg(index,status) RESULT(aval)
  IMPLICIT NONE
  COMPLEX*16 :: aval
  INTEGER, INTENT(IN) :: index
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  CHARACTER*256 :: arg
  INTEGER :: index0,Narg,ios
  
  aval = -HUGE(1d0)*(1,1)
  IF (PRESENT(status)) status = .FALSE.
  
  ! Check index validity, count from end if necessary
  Narg   = IARGC()
  IF (index .LT. 0) THEN
    index0 = Narg + 1 + index
  ELSE
    index0 = index
  END IF
  IF ((index0 .LT. 1) .OR. (index0 .GT. Narg)) RETURN
  
  ! Parse argument
  CALL GETARG(index0,arg)
  READ(UNIT=arg,FMT=*,IOSTAT=ios) aval
  IF (ios .EQ. 0) THEN
    IF (PRESENT(status)) status = .TRUE.
  END IF
  
END FUNCTION GetComplexArg


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form -akey or -<akey1><akey2><akey3>...
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
FUNCTION GetShortArg(akey) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  CHARACTER*256 :: arg
  INTEGER :: I,Narg,pos
  
  Narg   = IARGC()
  status = .FALSE.
  
  DO I=1,Narg
    CALL GETARG(I,arg)
    IF (LEN_TRIM(arg) .LT. 2)  CYCLE
    IF (arg(1:1) .NE. '-')  CYCLE
    IF (arg(2:2) .EQ. '-')  CYCLE
    pos = INDEX(arg,akey)
    IF (pos .GE. 2) THEN
      status = .TRUE.
      RETURN
    END IF
  END DO
  
END FUNCTION GetShortArg


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey or --akey=<val>.
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
FUNCTION GetLongArg(akey) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  CHARACTER*256 :: arg
  INTEGER :: I,Narg,len
  
  Narg   = IARGC()
  len    = LEN_TRIM(akey)
  status = .FALSE.
  
  DO I=1,Narg
    CALL GETARG(I,arg)
    IF (arg .EQ. '--' // akey) THEN
      status = .TRUE.
      RETURN
    ELSE IF (arg(:3+len) .EQ. '--' // akey // '=') THEN
      status = .TRUE.
      RETURN
    END IF
  END DO
  
END FUNCTION GetLongArg


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey[<suffix>] or --akey[<suffix>]=<val>.
! This differs from above in that it indicates if there are any command
! line arguments where the key name begins with the given string.
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
FUNCTION GetLongArgPrefix(prefix) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: prefix
  CHARACTER*256 :: arg
  INTEGER :: I,Narg,len
  
  Narg   = IARGC()
  len    = LEN_TRIM(prefix)
  status = .FALSE.
  
  DO I=1,Narg
    CALL GETARG(I,arg)
    IF (arg(:2+len) .EQ. '--' // prefix) THEN
      status = .TRUE.
      RETURN
    END IF
  END DO
  
END FUNCTION GetLongArgPrefix


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval> and inserts the appropriate string into
! aval.
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
FUNCTION GetLongArgString(akey,aval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  CHARACTER*(*), INTENT(OUT) :: aval
  CHARACTER*1024 :: arg
  INTEGER :: I,Narg,len,pos
  
  Narg   = IARGC()
  len    = LEN_TRIM(akey)
  status = .FALSE.
  
  DO I=1,Narg
    CALL GETARG(I,arg)
    IF (arg(:2+len) .EQ. '--' // akey) THEN
      pos = INDEX(arg,'=')
      IF (pos .EQ. 3+len) THEN
        status = .TRUE.
        aval = TRIM(arg(pos+1:))
        RETURN
      END IF
    END IF
  END DO
  
END FUNCTION GetLongArgString


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval1>,<aval2>,... and inserts the appropriate
! CHARACTER*N values into aval (ALLOCATABLE array of length Nval).
! 
! Currently requires allocatable array of fixed length strings,
! declared as:
!   CHARACTER(LEN=N), DIMENSION(:), ALLOCATABLE :: aval
! 
! A better implementation (but requires gfortran 4.6+) uses deferred-
! length strings, where aval is declared as:
!   CHARACTER(LEN=:), DIMENSION(:), ALLOCATABLE :: aval
! 
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
! NOTE: This function must be defined within a module or the use of
!       interfaces will be required.
! 
FUNCTION GetLongArgStrings(akey,N,aval,Nval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  INTEGER, INTENT(IN) :: N
  ! Fixed length strings
  CHARACTER(LEN=N), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: aval
  ! Better: use deferred length strings (but needs gfortran 4.6+)
  !CHARACTER(LEN=:), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: aval
  INTEGER, INTENT(OUT) :: Nval
  CHARACTER*1024 :: sval
  INTEGER :: I1,I2,Ia,ios
  
  IF (.NOT. GetLongArgString(akey,sval)) THEN
    status = .FALSE.
    RETURN
  END IF
  
  status = .TRUE.
  
  ! Find number of values
  I1   = 1
  Nval = 1
  DO WHILE (INDEX(sval(I1:),',') .GT. 0)
    I1 = I1 + INDEX(sval(I1:),',')
    Nval = Nval + 1
  END DO
  IF (ALLOCATED(aval)) DEALLOCATE(aval)
  ALLOCATE(CHARACTER(LEN=N) :: aval(1:Nval))
  aval = ''
  
  ! Parse values
  READ(UNIT=sval(I1:),FMT='(A)',IOSTAT=ios) aval(Nval)
  IF (ios .NE. 0) THEN
    status = .FALSE.
    Nval   = 0
    IF (ALLOCATED(aval)) DEALLOCATE(aval)
    ALLOCATE(CHARACTER(LEN=N) :: aval(0))
    RETURN
  END IF
  I2 = 1
  DO Ia = 1, Nval-1
    I1 = I2
    I2 = I1 + INDEX(sval(I1:),',')
    READ(UNIT=sval(I1:I2-2),FMT='(A)',IOSTAT=ios) aval(Ia)
    IF (ios .NE. 0) THEN
      status = .FALSE.
      Nval   = 0
      IF (ALLOCATED(aval)) DEALLOCATE(aval)
      ALLOCATE(CHARACTER(LEN=N) :: aval(0))
      RETURN
    END IF
  END DO
  
END FUNCTION GetLongArgStrings


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval> and inserts the appropriate LOGICAL value
! into aval.
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
FUNCTION GetLongArgLogical(akey,aval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  LOGICAL, INTENT(OUT) :: aval
  CHARACTER*256 :: sval
  INTEGER :: ios
  IF (.NOT. GetLongArgString(akey,sval)) THEN
    status = .FALSE.
    RETURN
  END IF
  READ(UNIT=sval,FMT=*,IOSTAT=ios) aval
  status = (ios .EQ. 0)
END FUNCTION GetLongArgLogical


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval1>,<aval2>,... and inserts the appropriate
! LOGICAL values into aval (ALLOCATABLE array of length Nval).
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
! NOTE: This function must be defined within a module or the use of
!       interfaces will be required.
! 
FUNCTION GetLongArgLogicals(akey,aval,Nval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  LOGICAL, ALLOCATABLE, INTENT(OUT) :: aval(:)
  INTEGER, INTENT(OUT) :: Nval
  ! Fixed length strings
  INTEGER, PARAMETER :: N = 64
  CHARACTER(LEN=N), ALLOCATABLE :: sval(:)
  ! Better: use deferred length strings (but needs gfortran 4.6+)
  !CHARACTER(LEN=:), ALLOCATABLE :: sval(:)
  INTEGER :: I,ios
  
  IF (.NOT. GetLongArgStrings(akey,N,sval,Nval)) THEN
    status = .FALSE.
    Nval   = 0
    IF (ALLOCATED(aval)) DEALLOCATE(aval)
    ALLOCATE(aval(0))
    RETURN
  END IF
  
  status = .TRUE.
  IF (ALLOCATED(aval)) DEALLOCATE(aval)
  ALLOCATE(aval(1:Nval))
  
  DO I=1,Nval
    READ(UNIT=sval(I),FMT=*,IOSTAT=ios) aval(I)
    IF (ios .NE. 0) THEN
      status = .FALSE.
      Nval   = 0
      IF (ALLOCATED(aval)) DEALLOCATE(aval)
      ALLOCATE(aval(0))
      RETURN
    END IF
  END DO
  
END FUNCTION GetLongArgLogicals


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval> and inserts the appropriate INTEGER value
! into aval.
! Accepts real valued strings, which will be rounded to the nearest
! integer.
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
FUNCTION GetLongArgInteger(akey,aval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  INTEGER, INTENT(OUT) :: aval
  CHARACTER*256 :: sval
  INTEGER :: ios
  REAL*8 :: x
  IF (.NOT. GetLongArgString(akey,sval)) THEN
    status = .FALSE.
    RETURN
  END IF
  ! Try reading as integer first, then as real and convert to integer
  READ(UNIT=sval,FMT=*,IOSTAT=ios) aval
  IF (ios .NE. 0) THEN
    READ(UNIT=sval,FMT=*,IOSTAT=ios) x
    IF (ios .EQ. 0) aval = NINT(x)
  END IF
  status = (ios .EQ. 0)
END FUNCTION GetLongArgInteger


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval1>,<aval2>,... and inserts the appropriate
! INTEGER values into aval (ALLOCATABLE array of length Nval).
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
! NOTE: This function must be defined within a module or the use of
!       interfaces will be required.
! 
FUNCTION GetLongArgIntegers(akey,aval,Nval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  INTEGER, ALLOCATABLE, INTENT(OUT) :: aval(:)
  INTEGER, INTENT(OUT) :: Nval
  ! Fixed length strings
  INTEGER, PARAMETER :: N = 64
  CHARACTER(LEN=N), ALLOCATABLE :: sval(:)
  ! Better: use deferred length strings (but needs gfortran 4.6+)
  !CHARACTER(LEN=:), ALLOCATABLE :: sval(:)
  INTEGER :: I,ios
  
  IF (.NOT. GetLongArgStrings(akey,N,sval,Nval)) THEN
    status = .FALSE.
    Nval   = 0
    IF (ALLOCATED(aval)) DEALLOCATE(aval)
    ALLOCATE(aval(0))
    RETURN
  END IF
  
  status = .TRUE.
  IF (ALLOCATED(aval)) DEALLOCATE(aval)
  ALLOCATE(aval(1:Nval))
  
  DO I=1,Nval
    READ(UNIT=sval(I),FMT=*,IOSTAT=ios) aval(I)
    IF (ios .NE. 0) THEN
      status = .FALSE.
      Nval   = 0
      IF (ALLOCATED(aval)) DEALLOCATE(aval)
      ALLOCATE(aval(0))
      RETURN
    END IF
  END DO
  
END FUNCTION GetLongArgIntegers


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval> and inserts the appropriate REAL*8 value
! into aval.
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
FUNCTION GetLongArgReal(akey,aval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  REAL*8, INTENT(OUT) :: aval
  CHARACTER*256 :: sval
  INTEGER :: ios
  IF (.NOT. GetLongArgString(akey,sval)) THEN
    status = .FALSE.
    RETURN
  END IF
  READ(UNIT=sval,FMT=*,IOSTAT=ios) aval
  status = (ios .EQ. 0)
END FUNCTION GetLongArgReal


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval1>,<aval2>,... and inserts the appropriate
! REAL*8 values into aval (ALLOCATABLE array of length Nval).
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
! NOTE: This function must be defined within a module or the use of
!       interfaces will be required.
! 
FUNCTION GetLongArgReals(akey,aval,Nval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  REAL*8, ALLOCATABLE, INTENT(OUT) :: aval(:)
  INTEGER, INTENT(OUT) :: Nval
  ! Fixed length strings
  INTEGER, PARAMETER :: N = 64
  CHARACTER(LEN=N), ALLOCATABLE :: sval(:)
  ! Better: use deferred length strings (but needs gfortran 4.6+)
  !CHARACTER(LEN=:), ALLOCATABLE :: sval(:)
  INTEGER :: I,ios
  
  IF (.NOT. GetLongArgStrings(akey,N,sval,Nval)) THEN
    status = .FALSE.
    Nval   = 0
    IF (ALLOCATED(aval)) DEALLOCATE(aval)
    ALLOCATE(aval(0))
    RETURN
  END IF
  
  status = .TRUE.
  IF (ALLOCATED(aval)) DEALLOCATE(aval)
  ALLOCATE(aval(1:Nval))
  
  DO I=1,Nval
    READ(UNIT=sval(I),FMT=*,IOSTAT=ios) aval(I)
    IF (ios .NE. 0) THEN
      status = .FALSE.
      Nval   = 0
      IF (ALLOCATED(aval)) DEALLOCATE(aval)
      ALLOCATE(aval(0))
      RETURN
    END IF
  END DO
  
END FUNCTION GetLongArgReals


! ----------------------------------------------------------------------
! Checks if the given argument key is given as a command line argument
! of the form --akey=<aval> and inserts the appropriate COMPLEX*16 value
! into aval.
! Returns .TRUE. if valid key was found, .FALSE. otherwise.
! 
FUNCTION GetLongArgComplex(akey,aval) RESULT(status)
  IMPLICIT NONE
  LOGICAL :: status
  CHARACTER*(*), INTENT(IN) :: akey
  COMPLEX*16, INTENT(OUT) :: aval
  CHARACTER*256 :: sval
  INTEGER :: ios
  IF (.NOT. GetLongArgString(akey,sval)) THEN
    status = .FALSE.
    RETURN
  END IF
  READ(UNIT=sval,FMT=*,IOSTAT=ios) aval
  status = (ios .EQ. 0)
END FUNCTION GetLongArgComplex


! ----------------------------------------------------------------------
! Returns the full command line used to run the program.
! This routine is nearly identical to the Fortran 2003 routine
! GET_COMMAND, but is provided here because some older compilers
! do not provide that routine.  White space between arguments will
! not be conserved: arguments will be separated by a single space,
! regardless of actual white space on the command line.
! 
SUBROUTINE GetFullCommand(cmd)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(OUT) :: cmd
  INTEGER :: pos,I,Narg,alen,clen
  CHARACTER*256 :: arg
  
  ! Fortran 2003
  !CALL GET_COMMAND(cmd)
  
  cmd = ''
  pos = 1
  clen = LEN(cmd)
  
  ! Command
  !CALL GETARG(0,arg,alen)
  CALL GETARG(0,arg)
  alen = LEN_TRIM(arg)
  cmd = arg
  pos = pos + alen
  
  ! Arguments
  Narg = IARGC()
  DO I=1, Narg
    !CALL GETARG(I,arg,alen)
    CALL GETARG(I,arg)
    alen = LEN_TRIM(arg)
    cmd(pos:pos) = ' '
    pos = pos + 1
    IF (pos+alen-1 .GT. clen) THEN
      cmd(pos:) = arg
    ELSE
      cmd(pos:pos+alen-1) = arg(1:alen)
    END IF
    pos = pos + alen
    IF (pos .GE. clen) EXIT
  END DO
  
  ! Following necessary to terminate/clear string
  cmd(pos:) = ''
  
END SUBROUTINE GetFullCommand



!=======================================================================
! EXPORT/IMPORT DATA
!=======================================================================

! ----------------------------------------------------------------------
! Returns an available (unused) file I/O unit number or -1 if no
! available unit number could be found.
! NOTE: A free unit number is a unit number not currently connected
!       at the time this function is called.  Numbers are not reserved,
!       so do not assume this number will remain free if it is not used
!       to open a file immediately.
! 
FUNCTION FreeIOUnit() RESULT(Nio)
  IMPLICIT NONE
  INTEGER :: Nio
  INTEGER, PARAMETER :: START_UNIT = 50
  INTEGER, PARAMETER :: MAX_UNIT   = 99
  LOGICAL opened
  
  ! Cycle through unit numbers until we find one that is not connected
  Nio = START_UNIT
  INQUIRE(UNIT=Nio,OPENED=opened)
  DO WHILE (opened .AND. (Nio .LT. MAX_UNIT))
    Nio = Nio + 1
    INQUIRE(UNIT=Nio,OPENED=opened)
  END DO
  
  ! Did not find unopen unit number in allowed range
  IF (opened) THEN
    WRITE(UNIT=5,FMT=*) 'ERROR: an error occurred while trying '        &
        // 'to find an available unit number.'
    Nio = -1
    RETURN
  END IF
  
END FUNCTION


! ----------------------------------------------------------------------
! Determines if the given line is a comment (non-data) line.
! A comment line is any blank line or line beginning with the comment
! character (default: '#').
!   line            the line to check
!   commentchar     (Optional) comment character.  Default is '#'.
! 
LOGICAL FUNCTION IsCommentLine(line,commentchar)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN) :: line
  CHARACTER, INTENT(IN), OPTIONAL :: commentchar
  
  IsCommentLine = .FALSE.
  IF (LEN_TRIM(line) .EQ. 0) THEN
    IsCommentLine = .TRUE.
  ELSE IF (PRESENT(commentchar)) THEN
    IF (line(1:1) .EQ. commentchar) IsCommentLine = .TRUE.
  ELSE
    IF (line(1:1) .EQ. '#') IsCommentLine = .TRUE.
  END IF
  
END FUNCTION


! ----------------------------------------------------------------------
! Determines the number of lines in the given file.
!   fid             The identifier for the file to be read
!   status          (Optional) .TRUE. if file was successfully read
! 
FUNCTION FileLines(fid,status) RESULT(lines)
  IMPLICIT NONE
  INTEGER :: lines
  INTEGER, INTENT(IN) :: fid
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  INTEGER :: ios
  
  lines = 0
  IF (PRESENT(status)) status = .FALSE.
  
  REWIND(UNIT=fid,IOSTAT=ios)
  IF (ios .LT. 0) RETURN
  
  DO WHILE (ios .GE. 0)
    lines = lines + 1
    READ(UNIT=fid,FMT=*,IOSTAT=ios)
  END DO
  
  ! Count included failed final read, so subtract 1
  lines = lines - 1
  
  REWIND(UNIT=fid,IOSTAT=ios)
  IF (PRESENT(status)) status = .TRUE.
  
END FUNCTION


! ----------------------------------------------------------------------
! Determines the number of data (non-comment) lines in the given file.
! Blank lines and lines beginning with the comment character (default:
! '#') are considered comment lines.
!   fid             The identifier for the file to be read
!   commentchar     (Optional) comment character.  Default is '#'.
!   status          (Optional) .TRUE. if file was successfully read
! 
FUNCTION FileDataLines(fid,commentchar,status) RESULT(lines)
  IMPLICIT NONE
  INTEGER :: lines
  INTEGER, INTENT(IN) :: fid
  CHARACTER, INTENT(IN), OPTIONAL :: commentchar
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  INTEGER :: ios
  CHARACTER :: commentchar0
  CHARACTER*1024 :: line
  
  IF (PRESENT(commentchar)) THEN
    commentchar0 = commentchar
  ELSE
    commentchar0 = '#'
  END IF
  
  lines = 0
  IF (PRESENT(status)) status = .FALSE.
  
  REWIND(UNIT=fid,IOSTAT=ios)
  IF (ios .LT. 0) RETURN
  
  DO WHILE (ios .GE. 0)
    READ(UNIT=fid,FMT='(A)',IOSTAT=ios) line
    IF (ios .GE. 0) THEN
      IF (.NOT. IsCommentLine(line)) lines = lines + 1
    END IF
  END DO
  
  REWIND(UNIT=fid,IOSTAT=ios)
  IF (PRESENT(status)) status = .TRUE.
  
END FUNCTION


! ----------------------------------------------------------------------
! Determines number of fields in the given string.
! 
PURE FUNCTION NumberOfFields(in) RESULT(N)
  IMPLICIT NONE
  INTEGER :: N
  CHARACTER*(*), INTENT(IN) :: in
  INTEGER :: I,len,pos
  LOGICAL :: is_blank
  
  len = LEN_TRIM(in)
  
  pos      = 0
  N        = 0
  is_blank = .TRUE.
  
  DO I = 1, len
    IF (is_blank .AND. (in(I:I) .NE. ' ')) THEN
      N = N + 1
      is_blank = .FALSE.
    ELSE
      is_blank = (in(I:I) .EQ. ' ')
    END IF
  END DO
  
END FUNCTION


! ----------------------------------------------------------------------
! Determines the number of data columns in the given file.
! This is based upon the number of fields separated by spaces on a
! data line.
! NOTE: Currently limited to first 1024 characters in a line.
!       This can be increased below.
!   fid             The identifier for the file to be read
!   commentchar     (Optional) comment character.  Default is '#'.
!   status          (Optional) .TRUE. if file was successfully read
! 
FUNCTION FileDataColumns(fid,commentchar,status) RESULT(columns)
  IMPLICIT NONE
  INTEGER :: columns
  INTEGER, INTENT(IN) :: fid
  CHARACTER, INTENT(IN), OPTIONAL :: commentchar
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  INTEGER :: ios
  CHARACTER :: commentchar0
  CHARACTER*1024 :: dataline    ! Increase for wide files
  
  IF (PRESENT(commentchar)) THEN
    commentchar0 = commentchar
  ELSE
    commentchar0 = '#'
  END IF
  
  columns = 0
  IF (PRESENT(status)) status = .FALSE.
  
  REWIND(UNIT=fid,IOSTAT=ios)
  IF (ios .LT. 0) RETURN
  
  ! Find a data line
  ios = 0
  dataline = ''
  DO WHILE ((ios .GE. 0) .AND. IsCommentLine(dataline,commentchar0))
    READ(UNIT=fid,FMT='(A)',IOSTAT=ios) dataline
  END DO
  ! No data lines?
  IF (ios .LT. 0) RETURN
  
  columns = NumberOfFields(dataline)
  
  REWIND(UNIT=fid,IOSTAT=ios)
  IF (PRESENT(status)) status = .TRUE.
  
END FUNCTION


! ----------------------------------------------------------------------
! Loads data from a formatted (text) file.
! Take an allocatable 2D array as an argument, which will be allocated
! to the correct size.  Any line beginning with a letter or symbol is
! treated as a comment line and ignored.
! Optional input arguments (at least one required):
!     file       The name of the file to be read
!     fid        The identifier for the file to be read.
!                File must already be open if given.
! Other optional input arguments:
!     commentchar  Comment character.  Default is '#'.
! Output arguments:
!     Nrow,Ncol  Number of data rows and columns; this indicates the
!                size of the allocated array
!     data       Allocatable REAL*8 array that will be allocated to
!                size data(Nrow,Ncol) and filled with data from the
!                file
! Optional output arguments:
!     success    Set to .TRUE. if the data was successfully loaded;
!                otherwise .FALSE.
! 
SUBROUTINE LoadTable(file,fid,commentchar,Nrow,Ncol,data,status)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN), OPTIONAL :: file
  INTEGER, INTENT(IN), OPTIONAL :: fid
  CHARACTER, INTENT(IN), OPTIONAL :: commentchar
  INTEGER, INTENT(OUT), OPTIONAL :: Nrow,Ncol
  REAL*8, ALLOCATABLE, INTENT(OUT), OPTIONAL :: data(:,:)
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  INTEGER :: fid0,Nrow0,Ncol0,ios,I
  LOGICAL :: status0
  CHARACTER :: commentchar0
  CHARACTER*256 :: buf
  
  ! Argument checking
  IF (.NOT. (PRESENT(file) .OR. PRESENT(fid))) THEN
    WRITE(0,*) 'ERROR: LoadTable requires a file or fid argument.'
    STOP
  END IF
  
  IF (PRESENT(commentchar)) THEN
    commentchar0 = commentchar
  ELSE
    commentchar0 = '#'
  END IF
  
  IF (PRESENT(status)) status = .FALSE.
  
  ! Open file, if necessary
  IF (PRESENT(fid)) THEN
    fid0 = fid
  ELSE
    fid0 = FreeIOUnit()
    OPEN(UNIT=fid0,FILE=file,STATUS='OLD',FORM='FORMATTED',             &
         ACTION='READ',IOSTAT=ios)
    IF (ios .NE. 0) RETURN
  END IF
  
  ! Rewind file, in case it has already been read from
  REWIND(UNIT=fid0,IOSTAT=ios)
  IF (ios .LT. 0) RETURN
  
  ! Get number of data lines
  Nrow0 = FileDataLines(fid0,commentchar0,status0)
  IF (.NOT. status0) RETURN
  IF (PRESENT(Nrow)) Nrow = Nrow0
  
  ! Get number of data columns
  Ncol0 = FileDataColumns(fid0,commentchar0,status0)
  IF (.NOT. status0) RETURN
  IF (PRESENT(Ncol)) Ncol = Ncol0
  
  ! No data array argument: nothing left to do
  IF (.NOT. PRESENT(data)) THEN
    IF (.NOT. PRESENT(fid)) THEN
      CLOSE(UNIT=fid0,IOSTAT=ios)
    END IF
    IF (PRESENT(status)) status = .TRUE.
    RETURN
  END IF
  
  ! Load data into array
  IF (ALLOCATED(data)) DEALLOCATE(data)
  ALLOCATE(data(Nrow0,Ncol0))
  I = 1
  DO WHILE (I .LE. Nrow0)
    READ(UNIT=fid0,FMT='(A)',IOSTAT=ios) buf
    ! Check if end of file or read error (ios < 0 for eof and
    ! ios > 0 for an I/O error)
    IF (ios .NE. 0) RETURN
    ! Check if comment line; if so, go back and read another line
    IF (IsCommentLine(buf)) CYCLE
    ! Parse line
    BACKSPACE(UNIT=fid0)
    READ(UNIT=fid0,FMT=*) data(I,:)
    I = I+1
  END DO
  
  ! Cleanup
  IF (.NOT. PRESENT(fid)) THEN
    CLOSE(UNIT=fid0,IOSTAT=ios)
  END IF
  
  IF (PRESENT(status)) status = .TRUE.
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Loads data columns from a formatted (text) file.
! Currently allows up to 5 selectable columns to be read.
! Takes allocatable arrays as arguments, which will be allocated to
! the correct size.  Any line beginning with a letter or symbol is
! treated as a comment line and ignored.
! Optional input arguments (at least one required):
!     file       The name of the file to be read
!     fid        The identifier for the file to be read.
!                File must already be open if given.
! Other optional input arguments:
!     commentchar  Comment character.  Default is '#'.
! Output arguments:
!     N          Length of data columns; this indicates the length
!                of the allocated arrays
! Optional input arguments:
!     N1,N2,N3,N4,N5
!                File columns to read into corresponding arrays.
!                If not specified, columns 1-5 will be read.
!                A negative number counts from the end (-1 is last
!                column).
! Optional output arguments:
!     C1,C2,C3,C4,C5
!                Allocatable REAL*8 arrays that will be allocated and
!                filled with data from the columns in the file
!                corresponding to N1 - N5.
!     success    Set to .TRUE. if the data was successfully loaded;
!                otherwise .FALSE.
! 
SUBROUTINE LoadArrays(file,fid,commentchar,N,N1,C1,N2,C2,N3,C3,N4,C4,   &
                      N5,C5,status)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN), OPTIONAL :: file
  INTEGER, INTENT(IN), OPTIONAL :: fid
  CHARACTER, INTENT(IN), OPTIONAL :: commentchar
  INTEGER, INTENT(OUT), OPTIONAL :: N
  INTEGER, INTENT(IN), OPTIONAL :: N1,N2,N3,N4,N5
  REAL*8, ALLOCATABLE, INTENT(OUT), OPTIONAL :: C1(:),C2(:),C3(:),      &
                                                C4(:),C5(:)
  LOGICAL, INTENT(OUT), OPTIONAL :: status
  INTEGER :: fid0,N0,Ncol0,ios,I,J,J1,J2,J3,J4,J5
  LOGICAL :: status0
  REAL*8, ALLOCATABLE :: data(:,:)
  CHARACTER :: commentchar0
  CHARACTER*256 :: buf
  
  ! Argument checking
  IF (.NOT. (PRESENT(file) .OR. PRESENT(fid))) THEN
    WRITE(0,*) 'ERROR: LoadArrays requires a file or fid argument.'
    STOP
  END IF

  IF (PRESENT(commentchar)) THEN
    commentchar0 = commentchar
  ELSE
    commentchar0 = '#'
  END IF
  
  IF (PRESENT(status)) status = .FALSE.
  
  ! Open file, if necessary
  IF (PRESENT(fid)) THEN
    fid0 = fid
  ELSE
    fid0 = FreeIOUnit()
    OPEN(UNIT=fid0,FILE=file,STATUS='OLD',FORM='FORMATTED',             &
         ACTION='READ',IOSTAT=ios)
    IF (ios .NE. 0) RETURN
  END IF
  
  ! Rewind file, in case it has already been read from
  REWIND(UNIT=fid0,IOSTAT=ios)
  IF (ios .LT. 0) RETURN
  
  ! Get number of data lines
  N0 = FileDataLines(fid0,commentchar0,status0)
  IF (.NOT. status0) RETURN
  IF (PRESENT(N)) N = N0
  
  ! Get number of data columns
  Ncol0 = FileDataColumns(fid0,commentchar0,status0)
  IF (.NOT. status0) RETURN
  
  ! Find how many columns should be read
  ! Following sets indices to 0, default, or appropriate
  ! optional argument
  IF (UpdateColNum(PRESENT(C1),PRESENT(N1),1,J1)) J1 = N1
  IF (PRESENT(C1)) CALL UpdateColNum2(Ncol0,J1,status0)
  IF (UpdateColNum(PRESENT(C2),PRESENT(N2),2,J2)) J2 = N2
  IF (PRESENT(C2)) CALL UpdateColNum2(Ncol0,J2,status0)
  IF (UpdateColNum(PRESENT(C3),PRESENT(N3),3,J3)) J3 = N3
  IF (PRESENT(C3)) CALL UpdateColNum2(Ncol0,J3,status0)
  IF (UpdateColNum(PRESENT(C4),PRESENT(N4),4,J4)) J4 = N4
  IF (PRESENT(C4)) CALL UpdateColNum2(Ncol0,J4,status0)
  IF (UpdateColNum(PRESENT(C5),PRESENT(N5),5,J5)) J5 = N5
  IF (PRESENT(C5)) CALL UpdateColNum2(Ncol0,J5,status0)
  IF (.NOT. status0) RETURN
  
  ! Only read in to last necessary column
  Ncol0 = MAX(J1,J2,J3,J4,J5)
  
  ! Load data into temporary array
  ALLOCATE(data(N0,Ncol0))
  I = 1
  DO WHILE (I .LE. N0)
    READ(UNIT=fid0,FMT='(A)',IOSTAT=ios) buf
    ! Check if end of file or read error (ios < 0 for eof and
    ! ios > 0 for an I/O error)
    IF (ios .NE. 0) RETURN
    ! Check if comment line; if so, go back and read another line
    IF (IsCommentLine(buf)) CYCLE
    ! Parse line
    BACKSPACE(UNIT=fid0)
    READ(UNIT=fid0,FMT=*) data(I,:)
    I = I+1
  END DO
  
  ! Allocate and fill arrays
  IF (PRESENT(C1)) CALL UpdateArray(J1,C1)
  IF (PRESENT(C2)) CALL UpdateArray(J2,C2)
  IF (PRESENT(C3)) CALL UpdateArray(J3,C3)
  IF (PRESENT(C4)) CALL UpdateArray(J4,C4)
  IF (PRESENT(C5)) CALL UpdateArray(J5,C5)
  
  ! Cleanup
  DEALLOCATE(data)
  IF (.NOT. PRESENT(fid)) THEN
    CLOSE(UNIT=fid0,IOSTAT=ios)
  END IF
  
  IF (PRESENT(status)) status = .TRUE.
  
  
  CONTAINS
    
  ! ------------------------------------------
  ! Utility to determine column number
  ! Will set column number to 0 if no column will be read or
  ! the default defCN if a column will be read.  If givenCN
  ! is .TRUE. (indicating if the optional NJ argument was given)
  ! and the column is to be read, this function will return
  ! .TRUE. indicating CN should be read from NJ instead.
  FUNCTION UpdateColNum(usecol,givenCN,defCN,CN)
    IMPLICIT NONE
    LOGICAL :: UpdateColNum
    LOGICAL, INTENT(IN) :: usecol,givenCN
    INTEGER, INTENT(IN) :: defCN
    INTEGER, INTENT(OUT) :: CN
    UpdateColNum = .FALSE.
    IF (usecol) THEN
      CN = defCN
      IF (givenCN) UpdateColNum = .TRUE.
    ELSE
      CN = 0
    END IF
  END FUNCTION UpdateColNum
  
  ! ------------------------------------------
  ! Utility to update column number
  ! Negative indices (indicating column number from end) are
  ! converted to true indices.  Validity of the column number
  ! is checked.
  SUBROUTINE UpdateColNum2(Ncol,CN,status0)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Ncol
    INTEGER, INTENT(INOUT) :: CN
    LOGICAL, INTENT(INOUT) :: status0
    IF (Ncol .LT. 0) THEN
      CN = Ncol + 1 - CN
    END IF
    status0 = status0 .AND. (CN .GT. 0) .AND. (CN .LE. Ncol)
  END SUBROUTINE UpdateColNum2
  
  ! ------------------------------------------
  ! Allocates and fills an array arr with data from column J
  SUBROUTINE UpdateArray(J,arr)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: J
    REAL*8, ALLOCATABLE, INTENT(OUT) :: arr(:)
    IF (ALLOCATED(arr)) DEALLOCATE(arr)
    ALLOCATE(arr(1:N0))
    arr = data(:,J)
  END SUBROUTINE UpdateArray
  
END SUBROUTINE



!=======================================================================
! TABULATION ROUTINES
!=======================================================================

! ----------------------------------------------------------------------
! Extracts parameters from arguments of the form:
!   --<akey>=<min>,<max>,<N>,<use_log>
! Does not change the existing value if the parameter is not given
! or cannot be parsed.  The integer argument N can be given as a
! floating-point number and the logical argument use_log can be
! as 'T'/'F' or numerically (zero is false, all other values true).
! 
SUBROUTINE GetTabulationArgs(akey,xmin,xmax,N,use_log)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN) :: akey
  REAL*8, INTENT(INOUT) :: xmin,xmax
  INTEGER, INTENT(INOUT) :: N
  LOGICAL, INTENT(INOUT) :: use_log
  LOGICAL :: status,Ltmp
  INTEGER :: Nval,ios,Itmp
  REAL*8 :: Rtmp
  ! Older compiler compatibility
  INTEGER, PARAMETER :: NCHAR = 32
  CHARACTER(LEN=NCHAR), DIMENSION(:), ALLOCATABLE :: aval
  ! ...but this would be better better (needs gfortran 4.6+)
  !CHARACTER(LEN=:), DIMENSION(:), ALLOCATABLE :: aval
  
  IF (.NOT. GetLongArgStrings(akey,NCHAR,aval,Nval)) RETURN
  
  IF (Nval .GE. 1) THEN
    READ(UNIT=aval(1),FMT=*,IOSTAT=ios) Rtmp
    IF (ios .EQ. 0) xmin = Rtmp
  END IF
  
  IF (Nval .GE. 2) THEN
    READ(UNIT=aval(2),FMT=*,IOSTAT=ios) Rtmp
    IF (ios .EQ. 0) xmax = Rtmp
  END IF
  
  IF (Nval .GE. 3) THEN
    READ(UNIT=aval(3),FMT=*,IOSTAT=ios) Rtmp  ! read as real
    IF (ios .EQ. 0) N = Rtmp
  END IF
  
  IF (Nval .GE. 4) THEN
    READ(UNIT=aval(4),FMT=*,IOSTAT=ios) Ltmp
    IF (ios .EQ. 0) THEN
      use_log = Ltmp
    ELSE
      READ(UNIT=aval(4),FMT=*,IOSTAT=ios) Rtmp
      IF (ios .EQ. 0) use_log = .NOT. (Rtmp .EQ. 0)
    END IF
  END IF
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Initializes the given TabulationStruct to use a tabulation specified
! by the given parameters.  Some sanity checks are performed.
! 
! Output argument:
!   TS              TabulationStruct to be filled with tabulation
!                   parameterization data.
! Input arguments:
!   xmin,xmax       Tabulation range.  Might be adjusted by up to 1/2
!                   bin size to ensure even sized bins in logarithmic
!                   case with N < 0.
!   N               Number of intervals between tabulation points.
!                   A negative number indicates intervals per decade.
!   use_log         If logarithmic rather than linear spacing is to
!                   be used.
! 
SUBROUTINE InitTabulation(TS,xmin,xmax,N,use_log)
  IMPLICIT NONE
  TYPE(TabulationStruct), INTENT(OUT) :: TS
  REAL*8, INTENT(IN) :: xmin,xmax
  INTEGER, INTENT(IN) :: N
  LOGICAL, INTENT(IN) :: use_log
  
  TS%logarithmic = use_log
  
  IF (use_log) THEN
    IF (xmin .LE. 0d0) THEN
      IF (xmax .LE. 0d0) THEN
        TS%xmin =   1d0
        TS%xmax = 100d0
      ELSE
        TS%xmin = xmax/100
        TS%xmax = xmax
      END IF
    ELSE
      IF (xmax .LE. 0d0) THEN
        TS%xmin = xmin
        TS%xmax = 100*xmin
      ELSE
        TS%xmin = xmin
        TS%xmax = xmax
      END IF
    END IF
    IF (N .EQ. 0) THEN
      TS%N = -10
    ELSE
      TS%N = N
    END IF
    TS%lnxmin = LOG(TS%xmin)
    TS%lnxmax = LOG(TS%xmax)
  ELSE
    TS%xmin = MIN(xmin,xmax)
    TS%xmax = MAX(xmin,xmax)
    IF (N .EQ. 0) THEN
      TS%N = 100
    ELSE
      TS%N = ABS(N)
    END IF
    TS%lnxmin = 0d0
    TS%lnxmax = 0d0
  END IF
  
  ! Special case (single point)
  IF (TS%xmax .EQ. TS%xmin) THEN
    TS%N     = 0
    TS%delta = 0d0
    RETURN
  END IF
  
  ! Logarithmic case
  IF (use_log) THEN
    IF (TS%N .LT. 0) THEN
      TS%delta = LOG(10d0) / ABS(TS%N)
      TS%N     = NINT(ABS(TS%N) * LOG10(TS%xmax/TS%xmin))
      TS%xmax  = TS%xmin*EXP(TS%N*TS%delta)
      TS%lnxmax = LOG(TS%xmax)
    ELSE
      TS%delta = (TS%lnxmax - TS%lnxmin) / TS%N
    END IF
  ! Linear case
  ELSE
    TS%delta = (TS%xmax - TS%xmin) / TS%N
  END IF
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Returns the Kth tabulation value, with K=0 (TS%N) corresponding to
! the minimum (maximum) of the tabulation range.
! 
PURE FUNCTION TabulationValue(TS,K) RESULT (x)
  IMPLICIT NONE
  REAL*8 :: x
  TYPE(TabulationStruct), INTENT(IN) :: TS
  INTEGER, INTENT(IN) :: K
  
  IF (TS%logarithmic) THEN
    x = EXP(TS%lnxmin + K*TS%delta)
  ELSE
    x = TS%xmin + K*TS%delta
  END IF
  
END FUNCTION


! ----------------------------------------------------------------------
! Returns the tabulation index K such that x_K <= x < x_{K+1}, or -1 if
! x does not fall within the tabulation range.  Note tabulation points
! are x_K with K=0,...,TS%N.
! 
PURE FUNCTION TabulationInterval(TS,x) RESULT (K)
  IMPLICIT NONE
  INTEGER :: K
  TYPE(TabulationStruct), INTENT(IN) :: TS
  REAL*8, INTENT(IN) :: x
  
  IF (TS%N .EQ. 0) THEN
    IF (x .EQ. TS%xmin) THEN
      K = 0
    ELSE
      K = -1
    END IF
    RETURN
  END IF
  
  IF (x .LT. TS%xmin) THEN
    K = -1
  ELSE IF (TS%logarithmic) THEN
    K = INT((LOG(x)-TS%lnxmin)/TS%delta)
  ELSE
    K = INT((x-TS%xmin)/TS%delta)
  END IF
  IF ((K .LT. 0) .OR. (K .GT. TS%N-1)) K = -1
  
END FUNCTION


! ----------------------------------------------------------------------
! Returns the tabulation index K such that x falls within a bin centered
! on x_K and with a width equal to the tabulation spacing.  A value of
! -1 is returned if x falls outside the tabulation bins.  For
! logarithmic tabulation, the bin center and size are taken in
! logarithmic space.  Note tabulation points are x_K with K=0,...,TS%N.
! 
PURE FUNCTION TabulationBin(TS,x) RESULT (K)
  IMPLICIT NONE
  INTEGER :: K
  TYPE(TabulationStruct), INTENT(IN) :: TS
  REAL*8, INTENT(IN) :: x
  
  IF (TS%N .EQ. 0) THEN
    IF (x .EQ. TS%xmin) THEN
      K = 0
    ELSE
      K = -1
    END IF
    RETURN
  END IF
  
  IF (TS%logarithmic) THEN
    K = NINT((LOG(x)-TS%lnxmin)/TS%delta)
  ELSE
    K = NINT((x-TS%xmin)/TS%delta)
  END IF
  IF ((K .LT. 0) .OR. (K .GT. TS%N)) K = -1
  
END FUNCTION



!=======================================================================
! INTERPOLATION
! Interpolation routines are given in two forms: scalar and array
! interpolation location(s).  An interface is defined at the top of
! the file to allow for a common name.
!=======================================================================

! ----------------------------------------------------------------------
! NOTE: Will use interface 'LinearInterpolate'.
! Determines y(x=x0) using linear interpolation between the points
! specified by the arrays x(N) and y(N).
! 
! Input arguments:
!   N               Length of x & y arrays
!   x,y             Arrays of x and y points defining y(x).
!                   Array in x must be increasing.
!   x0              Value to interpolate at
! Optional input arguments:
!   extrapolate     Indicate if y(x0) should be extrapolated from the
!                   given points if x0 falls outside of the domain of x.
!                   If set to .FALSE., this routine will return 0 if
!                   x0 is outside the domain.  Default is .TRUE.
! 
PURE FUNCTION LinearInterpolate_S(N,x,y,x0,extrapolate) RESULT(y0)
  IMPLICIT NONE
  REAL*8 :: y0
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: x(N),y(N),x0
  LOGICAL, INTENT(IN), OPTIONAL :: extrapolate
  INTEGER :: I
  REAL*8 :: a,b
  
  y0 = 0d0
  IF (N .EQ. 0) RETURN
  
  IF (PRESENT(extrapolate)) THEN
    IF ((.NOT. extrapolate)                                             &
        .AND. ((x0 .LT. x(1)) .OR. (x0 .GT. x(N)))) THEN
      RETURN
    END IF
  END IF
  
  IF (N .EQ. 1) THEN
    y0 = y(1)
    RETURN
  END IF
  
  ! Find interval to use for linear interpolation.
  ! Index here is for lower point in interval.
  I = BSearch(N,x,x0)
  I = MIN(I,N-1)
  I = MAX(I,1)
  
  IF (x(I) .EQ. x(I+1)) THEN
    y0 = 0.5d0 * (y(I) + y(I+1))
  ELSE
    a  = (y(I+1) - y(I)) / (x(I+1) - x(I))
    b  = y(I) - a*x(I)
    y0 = a*x0 + b
  END IF
  
END FUNCTION


! ----------------------------------------------------------------------
! NOTE: Will use interface 'LinearInterpolate'.
! Determines y(x=x0) using linear interpolation between the points
! specified by the arrays x(N) and y(N).
! 
! Input arguments:
!   N               Length of x & y arrays
!   x,y             Arrays of x and y points defining y(x).
!                   Array in x must be increasing.
!   N0              Length of x0 array
!   x0              Array of values to interpolate at
! Optional input arguments:
!   extrapolate     Indicate if y(x0) should be extrapolated from the
!                   given points if x0 falls outside of the domain of x.
!                   If set to .FALSE., this routine will return 0 if
!                   x0 is outside the domain.  Default is .TRUE.
! 
PURE FUNCTION LinearInterpolate_A(N,x,y,N0,x0,extrapolate) RESULT(y0)
  IMPLICIT NONE
  REAL*8 :: y0(N0)
  INTEGER, INTENT(IN) :: N,N0
  REAL*8, INTENT(IN) :: x(N),y(N),x0(N0)
  LOGICAL, INTENT(IN), OPTIONAL :: extrapolate
  LOGICAL :: extrapolate0
  INTEGER :: I,I0,Is
  REAL*8 :: a,b
  
  extrapolate0 = .TRUE.
  IF (PRESENT(extrapolate)) extrapolate0 = extrapolate
  
  ! Nothing to interpolate
  IF (N .EQ. 0) THEN
    y0 = 0d0
    RETURN
  END IF
  
  ! Single point: no searching necessary
  IF (N .EQ. 1) THEN
    IF (extrapolate0) THEN
      y0 = y(1)
    ELSE
      WHERE(x0 .EQ. x(1))
        y0 = y(1)
      ELSE WHERE
        y0 = 0d0
      END WHERE
    END IF
    RETURN
  END IF
  
  ! Cycle through x0 points
  Is = 0
  IF (extrapolate0) THEN
    DO I0 = 1,N0
      Is = BSearch(N,x,x0(I0),Is)
      Is = MAX(MIN(I,N-1),1)
      IF (x(I0) .EQ. x(I0+1)) THEN
        y0 = 0.5d0 * (y(I0) + y(I0+1))
      ELSE
        a      = (y(I0+1) - y(I0)) / (x(I0+1) - x(I0))
        b      = y(I0) - a*x(I0)
        y0(I0) = a*x0(I0) + b
      END IF
    END DO
  ELSE
    DO I0 = 1,N0
      IF ((x0(I0) .LT. x(1)) .OR. (x0(I0) .GE. x(N))) THEN
        y0(I0) = 0d0
      ELSE
        Is = BSearch(N,x,x0(I0),Is)
        IF (x(I0) .EQ. x(I0+1)) THEN
          y0 = 0.5d0 * (y(I0) + y(I0+1))
        ELSE
          a      = (y(I0+1) - y(I0)) / (x(I0+1) - x(I0))
          b      = y(I0) - a*x(I0)
          y0(I0) = a*x0(I0) + b
        END IF
      END IF
    END DO
  END IF
  
END FUNCTION



!=======================================================================
! RANDOM NUMBER GENERATOR
! The generator here is based on the algorithm of Marsaglia & Zaman
! (and later Tsang).  Some of the routines below are based upon an
! implementation by Christian Walck.  Due to the use of an OpenMP
! THREADPRIVATE directive when declaring the internal state structure,
! the various random number routines below should not have any
! performance issues when using multi-threading through OpenMP (as
! opposed to the RANDOM_NUMBER intrinsic, which may lead to _slower_
! code in multi-threaded cases due to its serial-only implementation
! leading to possible thread pile-up).  Upon entering the first OpenMP
! parallel region, each thread is given its own random seed unless
! explicit seeds are set by calling rand_init() within that parallel
! region; per-thread RNG states are then maintained for later parallel
! regions.  Note the RNG state from the serial region may be replaced
! with that of the master thread.
! 
! See:
!   G. Marsaglia and A. Zaman, "Toward a Universal Random Number
!     Generator," Florida State University Report: FSU-SCRI-87-50
!     (1987).
!   G. Marsaglia, A. Zaman and W.W. Tsang, "Toward a universal
!     random number generator," Statistics & Probability Letters,
!     9, 35, (1990).
!   F. James, "A Review of Pseudorandom Number Generators,"
!     Comput. Phys. Commun. 60, 329 (1990).
! 
!=======================================================================

! ----------------------------------------------------------------------
! Initializes state data for the random number generator of Marsaglia,
! Zaman & Tsang, optionally using the given seed.
! 
! Optional input argument:
!     seed       INTEGER seed used to generate state data.
!                If not given, a random seed will be used.
! Optional output argument:
!     state      A RandState structure containing RNG state
!                data that will be initialized.  If not given,
!                the internal RNG state will be initialized.
! 
SUBROUTINE rand_init(seed,state)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: seed
  TYPE(RandState), INTENT(OUT), OPTIONAL :: state
  INTEGER :: seed0,I,J,IF1,IF2,IF3,IC1,M
  TYPE(RandState) :: state0
  REAL*8 :: x,S,T
  INTEGER, PARAMETER :: MAX_SEED = 921350144
  
  IF (PRESENT(seed)) THEN
    seed0 = MODULO(seed,MAX_SEED)
  ELSE
    CALL RANDOM_NUMBER(x)
    seed0 = INT(x*MAX_SEED)
  END IF
  
  state0%initialized = .TRUE.
  state0%I = 97
  state0%J = 33
  
  ! Three Fibonacci generator seeds (2-177) and one congruential
  ! generator seed (0-168)
  IF1 = MOD(seed0/169/176/176,176) + 2
  IF2 = MOD(seed0/169/176,176) + 2
  IF3 = MOD(seed0/169,176) + 2
  IC1 = MOD(seed0,169)
  
  DO I = 1,97
     S = 0.0d0
     T = 0.5d0
     DO J = 1,24
        M = MOD(MOD(IF1*IF2,179)*IF3,179)
        IF1 = IF2
        IF2 = IF3
        IF3 = M
        IC1 = MOD(53*IC1+1,169)
        IF ( MOD(IC1*M,64) .GE. 32 ) S = S + T
        T = 0.5d0 * T
     END DO
     state0%U(I) = S
  END DO
  
  state0%C  =   362436.0d0 / 16777216.0d0
  state0%CD =  7654321.0d0 / 16777216.0d0
  state0%CM = 16777213.0d0 / 16777216.0d0
  
  IF (PRESENT(state)) THEN
    state = state0
  ELSE
    DEFAULT_RandState = state0
  END IF
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Returns a unit random number using the algorithm of Marsaglia, Zaman
! & Tsang, optionally using/updating the given state data.
! 
! Optional input/output argument:
!     state      A RandState structure containing RNG state
!                data that will be used for generating the random
!                number; the state will be updated.  If not given,
!                an internal RNG state will be used.
! 
FUNCTION rand(state) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  TYPE(RandState), INTENT(INOUT), OPTIONAL :: state
  IF (PRESENT(state)) THEN
    IF (.NOT. state%initialized) CALL rand_init(state=state)
    CALL rand_number(state,x)
  ELSE
    IF (.NOT. DEFAULT_RandState%initialized) CALL rand_init(state=DEFAULT_RandState)
    CALL rand_number(DEFAULT_RandState,x)
  END IF
END FUNCTION


! ----------------------------------------------------------------------
! Generates a unit random number using the algorithm of Marsaglia, Zaman
! & Tsang.  Requires a state structure to be given, which _must_ be
! initialized.  The rand() function above is the intended to be the
! main routine for external use.
! 
! Required input/output argument:
!     state      A RandState structure containing RNG state
!                data that will be used for generating the random
!                number; the state will be updated.
! Required output argument:
!     x          The uniform random number.
! 
ELEMENTAL SUBROUTINE rand_number(state,x)
  IMPLICIT NONE
  TYPE(RandState), INTENT(INOUT) :: state
  REAL*8, INTENT(OUT) :: x
  LOGICAL :: valid
  
  valid = .FALSE.
  DO WHILE (.NOT. valid)
    x = state%U(state%I) - state%U(state%J)
    IF (x .LT. 0d0) x = x + 1d0
    state%U(state%I) = x
    
    state%I = state%I - 1
    IF (state%I .EQ. 0) state%I = 97
    state%J = state%J - 1
    IF (state%J .EQ. 0) state%J = 97
    
    state%C = state%C - state%CD
    IF (state%C .LT. 0d0) state%C = state%C + state%CM
    x = x - state%C
    IF (x .LT. 0d0) x = x + 1d0
    
    ! Avoid returning zero
    valid = (x .GT. 0d0) .AND. (x .LT. 1d0)
  END DO
  
END SUBROUTINE



!=======================================================================
! NUMBER FORMATTING
!=======================================================================

! ----------------------------------------------------------------------
! Converts the given number to a string of width w, with nl characters
! to the left of the decimal point, np digits of precision, and ne
! digits in the exponential (set ne=0 for floating point only).  The
! number will be printed in floating point format if possible,
! otherwise exponential notation is used.  Numbers are coerced into
! a representable range.
! 
! Requirements:
!   nl + np + 2 + ne <= w  [ne > 0]
! 
PURE FUNCTION NumberToString(x,w,nl,np,ne) RESULT(s)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: x
  INTEGER, INTENT(IN) :: w
  INTEGER, INTENT(IN), OPTIONAL :: nl,np,ne
  CHARACTER*(w) :: s
  LOGICAL :: no_precision,use_fixed
  INTEGER :: nl0,nr0,np0,ne0,xsgn,exponent,expmax,decloc
  REAL*8 :: x0,eps,mantissa
  CHARACTER*64 :: XFMT,MFMT,EFMT,t
  
  ! Determine number of characters
  nl0 = 2
  np0 = 0
  ne0 = 0
  IF (PRESENT(nl)) nl0 = nl
  IF (PRESENT(np)) np0 = np
  IF (PRESENT(ne)) ne0 = ne
  nl0 = MAX(nl0,2)
  nr0 = w - nl0 - 1
  ne0 = MIN(MAX(ne0,0),9)
  no_precision = .FALSE.
  IF (Np0 .LE. 0) THEN
    no_precision = .TRUE.
    IF (Ne0 .GT. 0) THEN
      Np0 = Nr0 - 1 - Ne
    ELSE
      Np0 = Nr0 + 1
    END IF
  END IF
  
  ! Check for bad cases
  IF ((nr0 .LT. 0) .OR. (Np0 .LE. 0)) THEN
    s = REPEAT('*',w)
    RETURN
  END IF
  
  ! Use modifiable variable (x is fixed)
  x0 = x
  
  IF (x0 .GE. 0) THEN
    xsgn = +1
  ELSE
    xsgn = -1
  END IF
  eps = 10d0**(-Np0)
  
  ! Exponential notation components
  IF (x0 .NE. 0d0) THEN
    exponent = FLOOR(LOG10(ABS(x0*(1+0.51d0*eps))))
  ELSE
    exponent = 0d0
  END IF
  mantissa = x * 10d0**(-exponent)
  
  ! Coerce number into representable range
  IF (Ne0 .GT. 0) THEN
    expmax = 10**Ne0
    IF (ABS(exponent) .GE. expmax) THEN
      IF (exponent .GT. 0) THEN
        mantissa = xsgn * 10 * (1-eps)
        exponent = expmax - 1
      ELSE IF (Np0 - 1 + exponent .GT. -expmax) THEN
        mantissa = mantissa * 10d0**(exponent+expmax-1)
        exponent = -expmax + 1
      ELSE
        mantissa = 0d0
        exponent = 0
        x0       = 0d0
      END IF
    END IF
  ELSE
    IF ((x0 .GT. 0) .AND. (exponent + 1 .GT. Nl0)) THEN
      IF (no_precision) eps = MAX(10d0**(-w+2),EPSILON(eps))
      mantissa = 10 * (1-eps)
      exponent = Nl0 - 1
      x0       = mantissa * 10d0**exponent
    ELSE IF ((x0 .LT. 0) .AND. (exponent + 2 .GT. Nl0)) THEN
      IF (no_precision) eps = MAX(10d0**(-w+3),EPSILON(eps))
      mantissa = -10 * (1-eps)
      exponent = Nl0 - 2
      x0       = mantissa * 10d0**exponent
    END IF
  END IF
  
  ! Determine if fixed notation should be used
  use_fixed = .TRUE.
  IF (Ne0 .GT. 0) THEN
    IF (Np0 - 1 - exponent .GT. Nr0) THEN
      use_fixed = .FALSE.
    ELSE IF ((x0 .GT. 0) .AND. (exponent + 1 .GT. Nl0)) THEN
      use_fixed = .FALSE.
    ELSE IF ((x0 .LT. 0) .AND. (exponent + 2 .GT. Nl0)) THEN
      use_fixed = .FALSE.
    END IF
  END IF
  
  ! Construct string
  ! Fixed format
  IF (use_fixed) THEN
    IF (no_precision) THEN
      WRITE(XFMT,'(A,I2,A,I2,A)') '(F',w,'.',Nr0,')'
    ELSE
      WRITE(XFMT,'(A,I2,A,I2,A)') '(F',w,'.',MAX(Np0-exponent-1,0),')'
    END IF
    WRITE(t,XFMT) x0
    t = ADJUSTL(t)
    decloc = INDEX(t,'.')
    IF (decloc .EQ. 0) THEN
      decloc = LEN_TRIM(t) + 1
      !t = TRIM(t) // '.'
    END IF
    IF (decloc .LT. nl0 + 1) THEN
      s = REPEAT(' ',nl0 + 1 - decloc) // t
    ELSE IF (decloc .GT. nl0 + 1) THEN
      ! This should not happen....
      s = t(decloc-nl0:)
    ELSE
      s = t
    END IF
  ! Exponential format
  ELSE
    ! Mantissa part
    WRITE(MFMT,'(A,I2,A,I2,A)') '(F',Np0+2,'.',Np0-1,')'
    WRITE(t,MFMT) mantissa
    t = ADJUSTL(t)
    decloc = INDEX(t,'.')
    IF (decloc .EQ. 0) THEN
      decloc = LEN_TRIM(t) + 1
      !t = TRIM(t) // '.'
    END IF
    IF (decloc .LT. nl0 + 1) THEN
      s = REPEAT(' ',nl0 + 1 - decloc) // t
    ELSE IF (decloc .GT. nl0 + 1) THEN
      ! This should not happen....
      s = t(decloc-nl0:)
    ELSE
      s = t
    END IF
    ! Exponential part
    WRITE(EFMT,'(A,I1,A)') '(I',Ne0,')'
    WRITE(t,EFMT) ABS(exponent)
    t = ADJUSTL(t)
    IF (LEN_TRIM(t) .LT. Ne0) t = REPEAT('0',Ne0-LEN_TRIM(t)) // TRIM(t)
    IF (exponent .GE. 0) THEN
      t = 'E+' // TRIM(t)
    ELSE
      t = 'E-' // TRIM(t)
    END IF
    ! Combine parts
    s = TRIM(s) // t
  END IF
  
END FUNCTION


! ----------------------------------------------------------------------
! Returns the given number coerced into a range where it can be printed
! in at most w characters, optionally with a fixed number of decimal
! digits d (not using scientific notation).  Large magnitude numbers
! are reduced and very small magnitude numbers are set to zero.
! 
ELEMENTAL FUNCTION CoerceNumber(x,w,d) RESULT(y)
  IMPLICIT NONE
  REAL*8 :: y
  REAL*8, INTENT(IN) :: x
  INTEGER, INTENT(IN) :: w
  INTEGER, INTENT(IN), OPTIONAL :: d
  REAL*8 :: xmin,xmax
  y = x
  IF (PRESENT(d)) THEN
    xmin = (0.5d0+EPSILON(1d0)) * 10d0**(-d)
    IF (ABS(x) .LT. xmin) THEN
      y = 0d0
    ELSE IF (x .GE. 0) THEN
      xmax = 10d0**(w-d-1) - 2*xmin
      IF (x .GT. xmax) y = xmax
    ELSE
      xmax = 10d0**(w-d-2) - 2*xmin
      IF (ABS(x) .GT. xmax) y = -xmax
    END IF
  ELSE IF (x .GE. 0) THEN
    xmin = (0.5d0+EPSILON(1d0)) * 10d0**(2-w)
    xmax = 10d0**w - 1d0
    IF (x .LT. xmin) THEN
      y = 0d0
    ELSE IF (x .GT. xmax) THEN
      y = xmax
    END IF
  ELSE
    xmin = (0.5d0+EPSILON(1d0)) * 10d0**(3-w)
    xmax = 10d0**(w-1) - 1d0
    IF (ABS(x) .LT. xmin) THEN
      y = 0d0
    ELSE IF (ABS(x) .GT. xmax) THEN
      y = -xmax
    END IF
  END IF
  RETURN
END FUNCTION


! ----------------------------------------------------------------------
! Returns the given number coerced into a range with an exponent of at
! most N digits.  Large magnitude numbers are reduced and very small
! magnitude numbers are set to zero.  The number of significant digits
! Ns can optionally be given (controls how many 9's appear in upper
! cutoff).
! 
ELEMENTAL FUNCTION CoerceExponent(x,N,Ns) RESULT(y)
  IMPLICIT NONE
  REAL*8 :: y
  REAL*8, INTENT(IN) :: x
  INTEGER, INTENT(IN) :: N
  INTEGER, INTENT(IN), OPTIONAL :: Ns
  REAL*8 :: xmin,xmax
  xmin = 10d0**(-10**N+1)
  IF (PRESENT(Ns)) THEN
    xmax = (1-0.1d0**Ns)*10d0**(10**N)
  ELSE
    xmax = 0.9*10d0**(10**N)
  END IF
  IF (ABS(x) .LT. xmin) THEN
    y = 0d0
  ELSE IF (ABS(x) .GE. xmax) THEN
    y = SIGN(1d0,x)*xmax
  ELSE
    y = x
  END IF
  RETURN
END FUNCTION



!=======================================================================
! MATH FUNCTIONS
!=======================================================================

!-----------------------------------------------------------------------
! Error function between x1 and x2, i.e. erf(x2)-erf(x1).
! This routine accounts for the case when x1 and x2 are similar
! and loss of precision will result from canceling when an explicit
! subtraction of two calls to the implicit ERF function is performed.
! 
! Implementation here appears to be valid to ~30*epsilon, where
! epsilon is the smallest unit of precision (e.g. 1e-15 for double
! precision).  Precision can get somewhat worse than this when
! x1,x2 > 15 or x1,x2 < -15, but only in cases where erf(x1,x2)
! < 1e-100.  The precision has been determined only by testing various
! cases using Mathematica and has not been formally proven.
! 
! Possibly useful identities:
!   erf(x)  = GammaP(1/2,x^2)
!   erfc(x) = GammaQ(1/2,x^2)
! 
ELEMENTAL FUNCTION ERF2(x1,x2) RESULT(z)
  IMPLICIT NONE
  REAL*8 :: z
  REAL*8, INTENT(IN) :: x1,x2
  REAL*8 :: xc,delx
  REAL*8, PARAMETER :: SQRTPI = 1.7724538509055160d0
  ! 1-eps is approximate level of cancelation at which to handle
  ! the difference more carefully
  REAL*8, PARAMETER :: EPS = 0.03d0
  
  ! Opposite sign: no canceling to worry about here
  IF (x1*x2 .LE. 0d0) THEN
    z = erf(x2) - erf(x1)
    RETURN
  END IF
  
  xc   = 0.5d0 * (x1+x2)
  delx = 0.5d0 * (x2-x1)
  
  ! Smaller arguments:
  !   |xc| < 1    --> |x1|,|x2| < 2
  ! Canceling is significant if |delx| < eps*|xc|
  IF ((ABS(xc) .LE. 1d0) .AND. (ABS(delx) .GT. EPS*ABS(xc))) THEN
    ! erf(x2) - erf(x1)
    z = erf(x2) - erf(x1)
    RETURN
    
  ! At least one argument is "large": 
  !   |xc| > 1    --> |x1| > 1 or |x2| > 1
  ! Canceling is significant if |4*xc*delx| < eps
  ELSE IF ((ABS(xc) .GT. 1d0) .AND. (ABS(4*xc*delx) .GT. EPS)) THEN
    IF (xc .GT. 0d0) THEN
      ! Difference of complementary error function gives better
      ! precision here:
      ! erf(x2) - erf(x1) = erfc(x1) - erfc(x2)
      z = -(erfc(x2) - erfc(x1))
    ELSE
      ! Difference of complementary error function gives better
      ! precision here (use symmetry of erf(x)):
      ! erf(x2) - erf(x1) = erf(-x1) - erf(-x2) = erfc(-x2) - erfc(-x1)
      z = erfc(-x2) - erfc(-x1)
    END IF
    RETURN
  END IF
  
  ! If we reached this point, there is significant canceling.
  ! For these cases, the integrand in the error function does not
  ! change much over x1 < x < x2.  Taylor expand the integrand
  ! about xc = (x1+x2)/2 and integrate.  The following keeps terms
  ! up through tenth order.
  z = 4 * delx * EXP(-xc**2) / SQRTPI                                   &
        * (1 + (2*xc**2 - 1)*delx**2 / 3                                &
             + (4*xc**4 - 12*xc**2 + 3)*delx**4 / 30                    &
             + (8*xc**6 - 60*xc**4 + 90*xc**2 - 15)*delx**6 / 630       &
             + (16*xc**8 - 224*xc**6 + 840*xc**4 - 840*xc**2            &
                         + 105)*delx**8 / 22680                         &
          )
  
END FUNCTION


!-----------------------------------------------------------------------
! The logarithm of the error function.  Requires x > 0.
! Accounts for large x case to avoid loss of precision.
! 
ELEMENTAL FUNCTION LOG_ERF(x) RESULT(z)
  IMPLICIT NONE
  REAL*8 :: z
  REAL*8, INTENT(IN) :: x
  
  ! Invalid input
  IF (x .LE. 0d0) THEN
    z = -HUGE(z)
    RETURN
  END IF
  
  ! Negligible loss of precision for smaller x
  IF (x .LE. 1d0) THEN
    z = LOG(ERF(x))
  ! Work with smaller complementary error function
  ELSE
    z = LOGp1(-ERFC(x))
  END IF
  
END FUNCTION


!-----------------------------------------------------------------------
! The logarithm of the complementary error function.
! Accounts for large x case to avoid loss of precision.
! 
ELEMENTAL FUNCTION LOG_ERFC(x) RESULT(z)
  IMPLICIT NONE
  REAL*8 :: z
  REAL*8, INTENT(IN) :: x
  REAL*8 :: y,w
  REAL*8, PARAMETER :: SQRTPI = 1.7724538509055160d0  ! Sqrt(Pi)
  
  ! Work with smaller error function for smaller |x|
  IF (ABS(x) .LE. 0.1d0) THEN
    z = LOGp1(-ERF(x))
  ! Negligible loss of precision for x not too large
  ELSE IF (x .LE. 25d0) THEN
    z = LOG(ERFC(x))
  ! Use asymptotic expansion:
  !   w = sqrt(pi) x e^{x^2} erfc(x) - 1
  !       ->  \sum_{k=1} (-1)^k (2k-1)!!/(2x^2)^k
  ! For x > 25, double precision in eight terms (k <= 8)
  ELSE
    y = 1 / (2*x**2)
    w = -y*(1-3*y*(1-5*y*(1-7*y*(1-9*y*(1-11*y*(1-13*y*(1-15*y)))))))
    z = -x**2 - LOGp1(w) - LOG(SQRTPI*x)
  END IF
  
END FUNCTION


!-----------------------------------------------------------------------
! The inverse error function, finding x given y such that y = ERF(x).
! 
! This implementation uses Halley's method.  It seems to be reasonable
! in terms of speed and precision and has had a cursory check for
! robustness, but there may be more optimal algorithms than this.
! Precision seems to be within ~10*EPSILON (usually even closer).
! 
ELEMENTAL FUNCTION ERFINV(y) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN) :: y
  INTEGER :: K
  REAL*8 :: y0,y2,w,fk,delta
  REAL*8, PARAMETER :: QUARTERPI  = 0.78539816339744831d0  ! Pi/4
  REAL*8, PARAMETER :: INVSQRTPI  = 0.56418958354775629d0  ! 1/Sqrt(Pi)
  
  y0 = ABS(y)
  y2 = y*y
  
  ! Initial guess.  Taken from B.A. Popov, ACM SIGSAM 34, 25 (2000).
  IF (y0 .EQ. 0d0) THEN
    x = 0d0
    RETURN
  ELSE IF (y0 .LE. 0.97314979d0) THEN
    x = y*(-0.95493118d0+y2*(0.53160534d0+y2*0.23343441d0)) / (-1.0977154d0+y2)
  ELSE IF (y0 .LE. 0.99767065d0) THEN
    x = y*(1.6200516d0+y2*(-4.9295187d0+y2*3.2890636d0)) / (-1.0083317d0+y2)
  ELSE IF (y0 .LE. 0.99978842d0) THEN
    x = y*(29.849915d0+y2*(-61.896833d0+y2*32.044810d0)) / (-1.0007339d0+y2)
  ELSE IF (y0 .LT. 1d0) THEN
    w = -LOG(1-y2)
    x = SIGN(SQRT(w - 0.5d0*LOG(QUARTERPI*w)),y)
  ELSE
    x = SIGN(HUGE(x),y)
    RETURN
  END IF
  
  ! Halley's method.  Checked by hand that double precision is achieved
  ! by three iterations for all tested cases (epsilon, 1-epsilon, +/-,
  ! etc).  No safety checks as the algorithm appears to be convergent
  ! in all cases (though not proved).
  DO K=1,3
    fk    = ERF(x) - y
    ! Newton's
    !delta = - fk / (2*INVSQRTPI*EXP(-x**2))
    ! Halley's
    delta = fk / (x*fk - 2*INVSQRTPI*EXP(-x**2))
    x     = x + delta
  END DO
  
END FUNCTION


!-----------------------------------------------------------------------
! The inverse complementary error function, finding x given y such that
! y = ERFC(x).
! 
! This is a quick implementation that has not been optimized for speed
! and has not been tested for precision (though it should be close to
! double precision in most cases).
! 
ELEMENTAL FUNCTION ERFCINV(y) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN) :: y
  INTEGER :: K
  REAL*8 :: fk,delta
  REAL*8, PARAMETER :: INVPI      = 0.31830988618379067d0  ! 1/Pi
  REAL*8, PARAMETER :: INVSQRTPI  = 0.56418958354775629d0  ! 1/Sqrt(Pi)

  ! Use ERFINV via the relation erfc^-1(1-z) = erf^-1(z).
  ! Below algorith (Halley's method, 5 iterations) works fine for
  ! y < 0.5 and possibly to higher y values, but ERFINV is faster
  ! (by x2-5), so we use that where loss of precision is minimal.
  IF (y .GT. 0.01d0) THEN
    x = ERFINV(1d0-y)
    RETURN
  ELSE IF (y .LE. 0d0) THEN
    x = HUGE(x)
    RETURN
  END IF

  ! Only small, positive y at this point.

  ! Initial guess: invert first term of asymptotic expansion
  !   y = erfc(x) = e^{-x^2} / (x\sqrt{\pi}) [1 + O(1/x^2)]
  x = SQRT(0.5d0*LAMBERTW(2*INVPI/y**2))

  ! Halley's method.  Checked by hand that double precision is achieved
  ! by five iterations for tested cases (epsilon, ~ 0.5).  No safety
  ! checks as the algorithm appears to be convergent in all cases
  ! (though not proved).
  DO K=1,5
    fk    = ERFC(x) - y
    ! Newton's
    !delta = fk / (2*INVSQRTPI*EXP(-x**2))
    ! Halley's
    delta = fk / (x*fk + 2*INVSQRTPI*EXP(-x**2))
    x     = x + delta
  END DO

END FUNCTION


!-----------------------------------------------------------------------
! The value of exp(x2)-exp(x1).
! Accounts for the case when x2-x1 is small and loss of precision
! due to canceling of the two terms can occur.
ELEMENTAL FUNCTION EXP2(x1,x2) RESULT(z)
  IMPLICIT NONE
  REAL*8 :: z
  REAL*8, INTENT(IN) :: x1,x2
  REAL*8 :: x,yk
  INTEGER :: K
  x = x2-x1
  ! Good to ~ 100*EPSILON precision (matched to sum terms below)
  IF (ABS(x) .GT. 0.007d0) THEN
    z = EXP(x2) - EXP(x1)
    RETURN
  END IF
  IF (x1 .EQ. x2) THEN
    z = 0d0
    RETURN
  END IF
  ! Can write: e^x2 - e^x1 = e^x1 [e^(x2-x1) - 1]
  ! Here, we find e^(x2-x1) - 1
  K = 1
  yk = x
  z = yk
  ! Nearly full precision, but does check limit optimization?
  !DO WHILE (ABS(yk) .GT. EPSILON(1d0)*ABS(z))
  !  K = K+1
  !  yk = yk*x/K
  !  z = z + yk
  !END DO
  ! Optimization: might be faster to calculate a fixed number
  !               of terms rather than check after each term if
  !               further terms are required.
  ! Possibilities:
  !     Good to ~   20*EPSILON in 7 terms for |x| < 0.04.
  !     Good to ~   50*EPSILON in 6 terms for |x| < 0.02.
  !     Good to ~  100*EPSILON in 5 terms for |x| < 0.007.
  !     Good to ~  500*EPSILON in 4 terms for |x| < 0.002.
  !     Good to ~ 4000*EPSILON in 3 terms for |x| < 0.00025.
  DO K = 2,5
    yk = yk*x/K
    z = z + yk
  END DO
  z = EXP(x1)*z
END FUNCTION


!-----------------------------------------------------------------------
! The value of exp(x)-1.
! Accounts for the case when x is small and loss of precision due
! to canceling of the two terms can occur.
ELEMENTAL FUNCTION EXPm1(x) RESULT(z)
  IMPLICIT NONE
  REAL*8 :: z
  REAL*8, INTENT(IN) :: x
  z = EXP2(0d0,x)
END FUNCTION


! ----------------------------------------------------------------------
! Function to calculate the quantity:
!   ln(1+x)
! Accounts for the case where x is small.
! 
! Precision is ~ 100*EPSILON.
! 
ELEMENTAL FUNCTION LOGp1(x) RESULT(z)
  IMPLICIT NONE
  REAL*8 :: z
  REAL*8, INTENT(IN) :: x
  REAL*8 :: xabs,xk
  INTEGER :: k,klast
  
  xabs = ABS(x)
  
  ! If x is not too small, we can just evaluate explicitly
  IF (xabs .GT. 0.01d0) THEN
    z = LOG(1+x)
    RETURN
  END IF
  
  ! We will use an expansion about zero, keeping terms in the
  ! expansion up to x^klast/klast.
  ! klast is chosen below to give almost double precision.
  ! Precision can be as poor as ~ 100*EPSILON, but is better
  ! for smaller x.  The expansion is:
  !     \sum_{k=1}^{\infty} (-1)^{k+1} x^k / k
  ! The precision is approximately:
  !     x^klast / (klast+1)
  ! where klast is the last term in the sum.
  
  IF (xabs .LT. 1d-4) THEN
    klast = 4
  ELSE
    klast = 8
  END IF
  ! Go to klast=12 for |x| < 0.05
  ! Go to klast=16 for |x| < 0.1
  
  ! Use expansion about zero
  xk = x
  z  = xk
  DO k = 2,klast
    xk      = -x*xk
    z = z + xk/k
  END DO
  
END FUNCTION LOGp1


! ----------------------------------------------------------------------
! Function to calculate the quantity ln(a+b) given ln(a) and ln(b).
! Precision is ~ EPSILON.
! 
ELEMENTAL FUNCTION LOG_SUM(lna,lnb) RESULT(z)
  IMPLICIT NONE
  REAL*8 :: z
  REAL*8, INTENT(IN) :: lna,lnb
  REAL*8 :: lnx,r
  
  ! Write a+b as x(1+r) with x = max(a,b) and r = min(b/a,a/b).
  IF (lna .GE. lnb) THEN
    lnx = lna
    r   = EXP(lnb-lna)
  ELSE
    lnx = lnb
    r   = EXP(lna-lnb)
  END IF
  
  ! Below, we use an expansion of ln(1+r) if r is small.
  ! The sum is terminated at approximately double precision.
  IF (r .EQ. 0d0) THEN
    z = lnx
  ELSE IF (r .GT. 0.01d0) THEN
    z = lnx + LOG(1d0 + r)
  ELSE IF (r .GT. 0.001d0) THEN
    z = lnx + r*(1d0                                                    &
                 -r*((1d0/2d0)                                          &
                     -r*((1d0/3d0)                                      &
                         -r*((1d0/4d0)                                  &
                             -r*((1d0/5d0)                              &
                                 -r*((1d0/6d0)                          &
                                     -r*((1d0/7d0)                      &
                                         -r*((1d0/8d0)                  &
                                             -r*(1d0/9d0)               &
                ))))))))
  ELSE IF (r .GT. 0.00001d0) THEN
    z = lnx + r*(1d0                                                    &
                 -r*((1d0/2d0)                                          &
                     -r*((1d0/3d0)                                      &
                         -r*((1d0/4d0)                                  &
                             -r*((1d0/5d0)                              &
                                 -r*(1d0/6d0)                           &
                )))))
  ELSE
    z = lnx + r*(1d0                                                    &
                 -r*((1d0/2d0)                                          &
                     -r*((1d0/3d0)                                      &
                         -r*(1d0/4d0)                                   &
                )))
  END IF
  
END FUNCTION


! ----------------------------------------------------------------------
! Function to calculate the quantity ln(a-b) given ln(a) and ln(b).
! Requires a > b.  Precision is ~ EPSILON.
! 
ELEMENTAL FUNCTION LOG_DIFF(lna,lnb) RESULT(z)
  IMPLICIT NONE
  REAL*8 :: z
  REAL*8, INTENT(IN) :: lna,lnb
  REAL*8 :: lnx,r
  
  ! Bad case
  IF (lnb .GE. lna) THEN
    z = -HUGE(z)
    RETURN
  END IF
  
  ! Write a-b as x(1-r) with x = a and r = b/a.
  lnx = lna
  r   = EXP(lnb-lna)
  
  ! Below, we use an expansion of ln(1-r) if r is small.
  ! The sum is terminated at approximately double precision.
  IF (r .EQ. 0d0) THEN
    z = lnx
  ELSE IF (r .GT. 0.01d0) THEN
    z = lnx + LOG(1d0 - r)
  ELSE IF (r .GT. 0.001d0) THEN
    z = lnx - r*(1d0                                                    &
                 +r*((1d0/2d0)                                          &
                     +r*((1d0/3d0)                                      &
                         +r*((1d0/4d0)                                  &
                             +r*((1d0/5d0)                              &
                                 +r*((1d0/6d0)                          &
                                     +r*((1d0/7d0)                      &
                                         +r*((1d0/8d0)                  &
                                             +r*(1d0/9d0)               &
                ))))))))
  ELSE IF (r .GT. 0.00001d0) THEN
    z = lnx - r*(1d0                                                    &
                 +r*((1d0/2d0)                                          &
                     +r*((1d0/3d0)                                      &
                         +r*((1d0/4d0)                                  &
                             +r*((1d0/5d0)                              &
                                 +r*(1d0/6d0)                           &
                )))))
  ELSE
    z = lnx - r*(1d0                                                    &
                 +r*((1d0/2d0)                                          &
                     +r*((1d0/3d0)                                      &
                         +r*(1d0/4d0)                                   &
                )))
  END IF
  
END FUNCTION


!-----------------------------------------------------------------------
! Factorial
! 
! NOTE: Valid for INTEGER*4 or INTEGER*8 as the default INTEGER type.
! In the former case, some of the constants below overflow and must
! be coerced to the representable range, which is done automatically
! by ifort, but gfortran requires the additional compilation flag
! '-fno-range-check' (this flag only affects compile-time constants).
! 
ELEMENTAL FUNCTION FACTORIAL(n) RESULT(z)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  INTEGER :: z
  ! Tabulated values
  INTEGER, PARAMETER :: NMAX = 20
  INTEGER, PARAMETER, DIMENSION(0:NMAX) :: FVALS = &
    (/ 1,                                           & ! n = 0
       1,                  2,                    & ! n = 1-2
       6,                  24,                   & ! n = 3-4
       120,                720,                  & ! n = 5-6
       5040,               40320,                & ! n = 7-8
       362880,             3628800,              & ! n = 9-10
       39916800,           479001600,            & ! n = 11-12
       6227020800,         87178291200,          & ! n = 13-14
       1307674368000,      20922789888000,       & ! n = 15-16
       355687428096000,    6402373705728000,     & ! n = 17-18
       121645100408832000, 2432902008176640000  /) ! n = 19-20
  IF (n .LT. 0) THEN
    z = -HUGE(n)
    !STOP 'ERROR: factorial cannot be called with negative argument'
  ELSE IF (n .LE. NMAX) THEN
    z = FVALS(n)
  ELSE
    z = HUGE(n)
  END IF
END FUNCTION


!-----------------------------------------------------------------------
! Gamma function of integer argument [double precision]
! NOTE: Intrinsic routine of real argument available in Fortran 2008;
!       the speed of the intrinsic routine may be comparable to this
!       tabulated case.
! 
ELEMENTAL FUNCTION GAMMAI(n) RESULT(z)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL*8 :: z
  REAL*8 :: x,xinv
  INTEGER, PARAMETER :: NMAX  = 171    ! n > NMAX exceeds HUGE(1d0)
  ! Tabulated values
  ! NOTE: number of continuation lines required to place all tabulated
  !       values into one table exceeds F95 spec (>39).
  !       Both Intel and GNU Fortran compilers allow more continuation
  !       lines, but we use multiple tables for portability.
  INTEGER, PARAMETER :: NMAX1 = 100
  INTEGER, PARAMETER :: NMAX2 = 171
  ! Table for n=1-100
  REAL*8, PARAMETER, DIMENSION(1:NMAX1) :: GVALS1 = &
    (/  1.0d0,                  1.0d0,                  2.0d0,                  6.0d0,                  & ! 1-4
        24.0d0,                 120.0d0,                720.0d0,                5040.0d0,               & ! 5-8
        40320.0d0,              362880.0d0,             3.6288000000000000d6,   3.9916800000000000d7,   & ! 9-12
        4.7900160000000000d8,   6.2270208000000000d9,   8.7178291200000000d10,  1.3076743680000000d12,  & ! 13-16
        2.0922789888000000d13,  3.5568742809600000d14,  6.4023737057280000d15,  1.2164510040883200d17,  & ! 17-20
        2.4329020081766400d18,  5.1090942171709440d19,  1.1240007277776077d21,  2.5852016738884977d22,  & ! 21-24
        6.2044840173323944d23,  1.5511210043330986d25,  4.0329146112660564d26,  1.0888869450418352d28,  & ! 25-28
        3.0488834461171386d29,  8.8417619937397020d30,  2.6525285981219106d32,  8.2228386541779228d33,  & ! 29-32
        2.6313083693369353d35,  8.6833176188118865d36,  2.9523279903960414d38,  1.0333147966386145d40,  & ! 33-36
        3.7199332678990122d41,  1.3763753091226345d43,  5.2302261746660111d44,  2.0397882081197443d46,  & ! 37-40
        8.1591528324789773d47,  3.3452526613163807d49,  1.4050061177528799d51,  6.0415263063373836d52,  & ! 41-44
        2.6582715747884488d54,  1.1962222086548019d56,  5.5026221598120889d57,  2.5862324151116818d59,  & ! 45-48
        1.2413915592536073d61,  6.0828186403426756d62,  3.0414093201713378d64,  1.5511187532873823d66,  & ! 49-52
        8.0658175170943879d67,  4.2748832840600256d69,  2.3084369733924138d71,  1.2696403353658276d73,  & ! 53-56
        7.1099858780486345d74,  4.0526919504877217d76,  2.3505613312828786d78,  1.3868311854568984d80,  & ! 57-60
        8.3209871127413901d81,  5.0758021387722480d83,  3.1469973260387938d85,  1.9826083154044401d87,  & ! 61-64
        1.2688693218588416d89,  8.2476505920824707d90,  5.4434493907744306d92,  3.6471110918188685d94,  & ! 65-68
        2.4800355424368306d96,  1.7112245242814131d98,  1.1978571669969892d100, 8.5047858856786232d101, & ! 69-72
        6.1234458376886087d103, 4.4701154615126843d105, 3.3078854415193864d107, 2.4809140811395398d109, & ! 73-76
        1.8854947016660503d111, 1.4518309202828587d113, 1.1324281178206298d115, 8.9461821307829753d116, & ! 77-80
        7.1569457046263802d118, 5.7971260207473680d120, 4.7536433370128417d122, 3.9455239697206587d124, & ! 81-84
        3.3142401345653533d126, 2.8171041143805503d128, 2.4227095383672732d130, 2.1077572983795277d132, & ! 85-88
        1.8548264225739844d134, 1.6507955160908461d136, 1.4857159644817615d138, 1.3520015276784030d140, & ! 89-92
        1.2438414054641307d142, 1.1567725070816416d144, 1.0873661566567431d146, 1.0329978488239059d148, & ! 93-96
        9.9167793487094969d149, 9.6192759682482120d151, 9.4268904488832477d153, 9.3326215443944153d155 /) ! 97-100
  ! Table for n=101-171
  REAL*8, PARAMETER, DIMENSION(NMAX1+1:NMAX2) :: GVALS2 = &
    (/  9.3326215443944153d157, 9.4259477598383594d159, 9.6144667150351266d161, 9.9029007164861804d163, & ! 101-104
        1.0299016745145628d166, 1.0813967582402909d168, 1.1462805637347084d170, 1.2265202031961379d172, & ! 105-108
        1.3246418194518290d174, 1.4438595832024936d176, 1.5882455415227429d178, 1.7629525510902447d180, & ! 109-112
        1.9745068572210740d182, 2.2311927486598136d184, 2.5435597334721876d186, 2.9250936934930157d188, & ! 113-116
        3.3931086844518982d190, 3.9699371608087209d192, 4.6845258497542907d194, 5.5745857612076059d196, & ! 117-120
        6.6895029134491271d198, 8.0942985252734437d200, 9.8750442008336014d202, 1.2146304367025330d205, & ! 121-124
        1.5061417415111409d207, 1.8826771768889261d209, 2.3721732428800469d211, 3.0126600184576595d213, & ! 125-128
        3.8562048236258042d215, 4.9745042224772874d217, 6.4668554892204737d219, 8.4715806908788205d221, & ! 129-132
        1.1182486511960043d224, 1.4872707060906857d226, 1.9929427461615189d228, 2.6904727073180505d230, & ! 133-136
        3.6590428819525487d232, 5.0128887482749917d234, 6.9177864726194885d236, 9.6157231969410890d238, & ! 137-140
        1.3462012475717525d241, 1.8981437590761710d243, 2.6953641378881628d245, 3.8543707171800728d247, & ! 141-144
        5.5502938327393048d249, 8.0479260574719919d251, 1.1749972043909108d254, 1.7272458904546389d256, & ! 145-148
        2.5563239178728656d258, 3.8089226376305697d260, 5.7133839564458546d262, 8.6272097742332404d264, & ! 149-152
        1.3113358856834525d267, 2.0063439050956824d269, 3.0897696138473509d271, 4.7891429014633939d273, & ! 153-156
        7.4710629262828944d275, 1.1729568794264144d278, 1.8532718694937348d280, 2.9467022724950383d282, & ! 157-160
        4.7147236359920613d284, 7.5907050539472187d286, 1.2296942187394494d289, 2.0044015765453026d291, & ! 161-164
        3.2872185855342962d293, 5.4239106661315888d295, 9.0036917057784374d297, 1.5036165148649990d300, & ! 165-168
        2.5260757449731984d302, 4.2690680090047053d304, 7.2574156153079990d306  /)                        ! 169-171
  ! Coefficients for asymptotic expansion
  ! Not all terms are necessary or used here
  INTEGER, PARAMETER :: NC = 20
  REAL*8, PARAMETER, DIMENSION(0:NC) :: C = &
    (/  1.0d0,                                                  &
        0.083333333333333333d0,     0.0034722222222222222d0,    &
       -0.0026813271604938272d0,   -0.00022947209362139918d0,   &
        0.00078403922172006663d0,   0.000069728137583658578d0,  &
       -0.00059216643735369388d0,  -0.000051717909082605922d0,  &
        0.00083949872067208728d0,   0.000072048954160200106d0,  &
       -0.0019144384985654775d0,   -0.00016251626278391582d0,   &
        0.0064033628338080698d0,    0.00054016476789260452d0,   &
       -0.029527880945699121d0,    -0.0024817436002649977d0,    &
        0.17954011706123486d0,      0.015056113040026424d0,     &
       -1.3918010932653375d0,      -0.1165462765994632d0       /)
  
  IF (n .LE. 0) THEN
    z = -HUGE(z)
    !STOP 'ERROR: gamma cannot be called with non-positive argument'
  ELSE IF (n .GT. NMAX) THEN
    z = HUGE(z)
  ELSE IF (n .LE. NMAX1) THEN
    z = GVALS1(n)
  ELSE IF (n .LE. NMAX2) THEN
    z = GVALS2(n)
  ELSE
    ! This case should not occur
    z = HUGE(z)
  END IF
  ! Algorithm here unnecessary, but shown for reference
  !! Asymptotic expansion (Stirling's)
  !! Error will be ~ C[5]/x^5
  !! Good to ~ 1e-15 for n > 171
  !! Result would be HUGE(1d0) for n > 171
  !x = n
  !xinv = 1d0/x
  !z = x**(x-0.5d0) * EXP(-x) * SQRT2PI                               &
  !    * (C(0) +                                                      &
  !       xinv*(C(1) +                                                &
  !             xinv*(C(2) +                                          &
  !                   xinv*(C(3) +                                    &
  !                         xinv*(C(4) +                              &
  !                               xinv*C(5) ) ) ) ) )
  !! Alternatively, could use LogGamma function for large n
  !!x = n
  !!z = EXP(LOG_GAMMA(x))
END FUNCTION


!-----------------------------------------------------------------------
! Gamma function of real argument [double precision]
! NOTE: Intrinsic routine of real argument available in Fortran 2008.
! 
! We make use of the identity:
!   Gamma(1-z) = pi*z / Gamma(1+z) / sin(pi*z)
! 
!ELEMENTAL FUNCTION GAMMA(x) RESULT(z)
!  IMPLICIT NONE
!  REAL*8, INTENT(IN) :: x
!  REAL*8 :: z
!  REAL*8 :: x0
!  REAL*8, PARAMETER :: PI = 3.1415926535897932d0
!  ! We make use of the identity:
!  !   Gamma(1-z) = pi*z / Gamma(1+z) / sin(pi*z)
!  IF (x .GE. 1d0) THEN
!    z = EXP(LOG_GAMMA(x))
!  ! Gamma is undefined when x=0,-1,-2,...
!  ELSE IF (MOD(x,1d0) .EQ. 0d0) THEN
!    z = -HUGE(z)
!  ELSE
!    x0 = 1 - x
!    z = PI*x0 / (SIN(PI*x0) * EXP(LOG_GAMMA(1+x0)))
!  END IF
!END FUNCTION


!-----------------------------------------------------------------------
! Logarithm of the gamma function [double precision]
! NOTE: Intrinsic routine available in Fortran 2008.
! 
! Calculated using Lanczos' approximation.  Valid only for x > 0.
! 
! Lanczos's approximation:
!   Gamma(z+1) = (z+r+1/2)^{z+1/2} * EXP[-(z+r+1/2)] * sqrt(2*pi)
!                * (B_0 + \Sum_{k=0}^{N} B_k/(z+k) + eps)
! where r is some constant s.t. z+r+1/2 > 0.
! Coefficients in the series expansion are dependent on r and the
! number of terms at which the series is truncated.  The determination
! of these coefficients is too complicated to show here.
! 
! Note the formula above is valid for complex cases as well, as long as
! Re(x) > 1.  The following identity also applies for the complex case:
!   Gamma(1-z) = pi*z / Gamma(1+z) / sin(pi*z)
! 
! The coefficients are calculated as described in Section 6.7 of
! Glendon Pugh's thesis (2004), which provides an extensive discussion
! of the Lanczos approximation.  The coefficients are given for the
! following modified form of the formula:
!   Gamma(z+1) = 2 * sqrt(e/pi) * ((z+r+1/2)/e)^(z+1/2)
!                * (D_0 + \Sum_{k=0}^{N} D_k/(z+k) + eps)
! The quantities N and r are chosen from Table C.1 to achieve the
! desired precision:
!   Double precision:  N=10, r=10.900511  (see Table 8.5)
! 
!ELEMENTAL FUNCTION LOG_GAMMA(x) RESULT(z)
!  IMPLICIT NONE
!  REAL*8, INTENT(IN) :: x
!  REAL*8 :: z
!  REAL*8 :: series_sum,x0
!  INTEGER :: K
!  REAL*8, PARAMETER :: PI = 3.1415926535897932d0
!  REAL*8, PARAMETER :: R = 10.900511d0
!  INTEGER, PARAMETER :: ND = 10
!  REAL*8, PARAMETER, DIMENSION(0:ND) :: Dk = &
!    (/ +2.4857408913875357d-5,                          &
!       +1.0514237858172197d0,  -3.4568709722201624d0,   &
!       +4.5122770946689482d0,  -2.9828522532357666d0,   &
!       +1.0563971157712671d0,  -1.9542877319164587d-1,  &
!       +1.7097054340444122d-2, -5.7192611740430578d-4,  &
!       +4.6339947335990564d-6, -2.7199490848860770d-9  /)
!  
!  IF (x .LE. 0d0) THEN
!    z = -HUGE(z)
!    RETURN
!    !STOP 'ERROR: loggamma cannot be called with non-positive argument'
!  END IF
!  ! Error for x>0 should theoretically be smaller than 1e-15,
!  ! but loss of precision in sum below might make it slightly worse
!  ! (maybe ~ 1e-14?)
!  IF (x .GE. 1) THEN
!    ! Write as Gamma(x) = Gamma(1+x0)
!    x0 = x-1
!  ELSE
!    ! Special case: for 0<x<1, find Gamma(1+(1-x)) = Gamma(1+x0)
!    ! and use identity Gamma(1-y) = pi*y / Gamma(1+y) / sin(pi*y)
!    x0 = 1 - x
!  END IF
!  series_sum = Dk(0)
!  DO K = 1,ND
!    series_sum = series_sum + Dk(K)/(x0+K)
!  END DO
!  z = (x0+0.5d0)*(LOG(x0+R+0.5d0) - 1)                                  &
!      + 0.5d0*(1 + LOG(4d0/PI))                                         &
!      + LOG(series_sum)
!  IF (x .LT. 1) THEN
!    ! Special case: for 0<x<1, we found LogGamma(1+(1-x)).
!    ! Now use identity Gamma(1-y) = pi*y / Gamma(1+y) / sin(pi*y)
!    x0 = 1 - x
!    z = LOG(PI*x0 / SIN(PI*x0)) - z
!  END IF
!END FUNCTION


!-----------------------------------------------------------------------
! Lower gamma function of real arguments [double precision]
!   \gamma(s,x) = \int_0^x dt t^{s-1} e^{-t}
! 
ELEMENTAL FUNCTION LOWER_GAMMA(s,x) RESULT(z)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s,x
  REAL*8 :: z
  REAL*8 :: P,Q
  CALL GAMMA_PQ(s,x,P,Q)
  z = GAMMA(s) * P
END FUNCTION


!-----------------------------------------------------------------------
! Upper gamma function of real arguments [double precision]
!   \Gamma(s,x) = \int_x^\infty dt t^{s-1} e^{-t}
! 
ELEMENTAL FUNCTION UPPER_GAMMA(s,x) RESULT(z)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s,x
  REAL*8 :: z
  REAL*8 :: P,Q
  CALL GAMMA_PQ(s,x,P,Q)
  z = GAMMA(s) * Q
END FUNCTION


!-----------------------------------------------------------------------
! Regularized gamma function P of real arguments [double precision]
!   P(s,x) = \gamma(s,x) / \Gamma(s)
! 
ELEMENTAL FUNCTION GAMMA_P(s,x) RESULT(z)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s,x
  REAL*8 :: z
  REAL*8 :: P,Q
  CALL GAMMA_PQ(s,x,P,Q)
  z = P
END FUNCTION


!-----------------------------------------------------------------------
! Regularized gamma function Q of real arguments [double precision]
!   Q(s,x) = \Gamma(s,x) / \Gamma(s)
! 
ELEMENTAL FUNCTION GAMMA_Q(s,x) RESULT(z)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s,x
  REAL*8 :: z
  REAL*8 :: P,Q
  CALL GAMMA_PQ(s,x,P,Q)
  z = Q
END FUNCTION


!-----------------------------------------------------------------------
! Regularized incomplete gamma functions P, Q of real arguments
! [double precision].  These are defined as:
!   P(s,x) = \gamma(s,x) / \Gamma(s)
!   Q(s,x) = \Gamma(s,x) / \Gamma(s)
! where \gamma(s,x) & \Gamma(s,x) are the lower & upper incomplete
! gamma functions:
!   \gamma(s,x) = \int_0^x dt t^{s-1} e^{-t}
!   \Gamma(s,x) = \int_x^\infty dt t^{s-1} e^{-t}
! The following relations hold:
!   P(s,x) + Q(s,x) = 1
!   \gamma(s,x) + \Gamma(s,x) = \Gamma(s)
! The below routine calculates the approximately smaller of P or Q
! and uses P+Q=1 to determine the other.
! 
! NOTE: Valid only for s > 0 and x >= 0.
! 
ELEMENTAL SUBROUTINE GAMMA_PQ(s,x,P,Q)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s,x
  REAL*8, INTENT(OUT) :: P,Q
  
  ! Special cases
  IF (x .EQ. 0d0) THEN
    P = 0d0
    Q = 1d0
    RETURN
  END IF
  
  ! Bad cases
  IF ((x .LT. 0d0) .OR. (s .LT. 0d0)) THEN
    P = -HUGE(1d0)
    Q = HUGE(1d0)
    RETURN
  END IF
  
  ! Calculate P,Q using uniform asymptotic expansion
  IF ((s .GE. 10d0) .AND. (x-10d0 .GT. 0.302d0*(s-10d0)) &
      .AND. (x-10d0 .LT. 2.357d0*(s-10d0))) THEN
    CALL GAMMA_PQ_UA(s,x,P,Q)
  ! Calculate P using Taylor series
  ELSE IF (x .LE. s) THEN
    P = GAMMA_P_TS(s,x)
    Q = 1-P
  ! Calculate Q using continued fraction
  ELSE
    Q = GAMMA_Q_CF(s,x)
    P = 1-Q
  END IF
  
  CONTAINS
  
  !---------------------------------------------
  ! Calculates gamma function quantity P(s,x) using a Taylor
  ! series expansion.
  !   P(s,x) = \gamma(s,x) / \Gamma(s)
  !   \gamma(s,x) = \int_0^x dt t^{s-1} e^{-t}
  ! NOTE: Requires 0 < s, 0 <= x.  Should only be used for x
  !       not (much) larger than s or convergence may be slow.
  !       Also best for s not large (timing scales as maybe
  !       sqrt(s)).
  PURE FUNCTION GAMMA_P_TS(s,x) RESULT(P)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: s,x
    REAL*8 :: P
    INTEGER :: K
    REAL*8 :: zk,sum
    ! P(s,x) = x^s e^{-x} \Sum_{k=0}^\infty x^k / \Gamma(s+k+1)
    K   = 0
    zk  = 1d0
    sum = zk
    DO WHILE (zk .GT. EPSILON(zk))
      K   = K+1
      zk  = zk * x / (s+K)
      sum = sum + zk
      IF (K .GE. 10000) EXIT
    END DO
    P = sum * EXP(s*LOG(x) - x - LOG_GAMMA(s+1))
  END FUNCTION
  
  !---------------------------------------------
  ! Calculates gamma function quantity Q(s,x) using a
  ! continued fraction.
  !   Q(s,x) = \Gamma(s,x) / \Gamma(s)
  !   \Gamma(s,x) = \int_x^{\infty} dt t^{s-1} e^{-t}
  ! NOTE: Requires 0 < x,s.  Should only be used for x not
  !       (much) smaller than s or convergence may be slow
  !       (or simply fail).  Also best for s not large
  !       (timing scales as maybe sqrt(s)).
  PURE FUNCTION GAMMA_Q_CF(s,x) RESULT(Q)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: s,x
    REAL*8 :: Q
    INTEGER :: K
    REAL*8 :: xs1,fk,Ck,Dk,Deltak
    ! Continued fraction (see Numerical Recipes):
    !   Q(s,x) = x^s e^{-x} / \Gamma(s)
    !            / (x-s+1 + K_k(-k(k-s),-s+2k+x+1)_{1}^{\infty})
    ! where
    !   K_k(a_k,b_k)_1 = a_1/b1+ a2/b2+ a3/b3+ ...
    ! is the Gauss notation for a continued fraction.
    xs1 = x-s+1
    K  = 0
    Ck = xs1
    Dk = 0
    fk = xs1
    Deltak = HUGE(x)
    DO WHILE (ABS(Deltak-1) .GE. EPSILON(x))
      K  = K+1
      Ck = (xs1+2*K) - K*(K-s)/Ck
      IF (Ck .EQ. 0d0) Ck = 1d-30
      Dk = (xs1+2*K) - K*(K-s)*Dk
      IF (Dk .EQ. 0d0) Dk = 1d-30
      Dk = 1/Dk
      Deltak = Ck*Dk
      fk = Deltak*fk
      IF (K .GE. 10000) EXIT
    END DO
    Q = EXP(s*LOG(x) - x - LOG_GAMMA(s)) / fk
  END FUNCTION
  
  !---------------------------------------------
  ! Calculates gamma function quantities P(s,x) and Q(s,x)
  ! using a uniform asymptotic expansion (see 1306.1754).
  ! Nearly constant evaluation time for any s & x, unlike
  ! the Taylor series or continued fraction algorithms.
  ! NOTE: Intended for s not small (larger than ~10) and
  !       0.30 < x/s < 2.35; accuracy starts to drop
  !       outside this region.
  PURE SUBROUTINE GAMMA_PQ_UA(s,x,P,Q)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: s,x
    REAL*8, INTENT(OUT) :: P,Q
    INTEGER :: K
    REAL*8 :: eta,lambda,u,Sa,Ra,sum
    REAL*8 :: beta(0:32)
    INTEGER, PARAMETER :: KMAX = 25
    REAL*8, PARAMETER :: D(0:31) = &
      (/  1.0000000000000000d+00,-3.3333333333333333d-01, 8.3333333333333333d-02, &
         -1.4814814814814815d-02, 1.1574074074074074d-03, 3.5273368606701940d-04, &
         -1.7875514403292181d-04, 3.9192631785224378d-05,-2.1854485106799922d-06, &
         -1.8540622107151600d-06, 8.2967113409530860d-07,-1.7665952736826079d-07, &
          6.7078535434014986d-09, 1.0261809784240308d-08,-4.3820360184533532d-09, &
          9.1476995822367902d-10,-2.5514193994946250d-11,-5.8307721325504251d-11, &
          2.4361948020667416d-11,-5.0276692801141756d-12, 1.1004392031956135d-13, &
          3.3717632624009854d-13,-1.3923887224181621d-13, 2.8534893807047443d-14, &
         -5.1391118342425726d-16,-1.9752288294349443d-15, 8.0995211567045613d-16, &
         -1.6522531216398162d-16, 2.5305430097478884d-18, 1.1686939738559577d-17, &
         -4.7700370498204848d-18, 9.6991260590562371d-19 /)
    REAL*8, PARAMETER :: TWOPI = 6.2831853071795865d0   ! 2*Pi
    
    ! Define ratio of x & s
    lambda = x/s
    
    ! eta is defined as 1/2 eta^2 = lambda - 1 - ln(lambda) with
    ! the same sign as lambda-1.  Use an expansion for lambda ~ 1
    ! to avoid cancellation issues.  Rapid convergent of below
    ! series in eta requires |eta| < 1, which corresponds to
    ! 0.3017 < lambda < 2.3577.
    u = lambda - 1
    IF (ABS(u) .GT. 0.01d0) THEN
      eta = SQRT(2*(lambda-1-LOG(lambda)))
    ELSE
      eta = u*(1+u*(-1d0/3d0+u*(7d0/36d0+u*(-73d0/540d0+u*(1331d0/12960d0)))))
    END IF
    IF (x .LT. s) eta = -eta
    
    ! Use erfc in eta*sqrt{s/2} as first approximation and calculate
    ! correction term:
    !   R_a(\eta) = e^{-1/2 s \eta^2} S_a(\eta) / sqrt{2\pi s}
    ! where Sa(eta) can be given by an expansion in eta:
    !   S_a(\eta) \approx a/(a+\beta_1) \Sum_{k=0}^N \beta_k \eta^k
    ! The beta terms are dependent upon s and can be calculated using
    ! a recursion relation.
    beta(KMAX+2) = 0d0
    beta(KMAX+1) = 0d0
    sum = 0d0
    DO K = KMAX,0,-1
      beta(K) = (K+2)/s * beta(K+2) + D(K+1)
      sum = beta(K) + eta*sum
    END DO
    Sa = s/(s+beta(1)) * sum
    Ra = EXP(-0.5d0*s*eta**2) * Sa / SQRT(TWOPI*s)
    
    ! Formulas for P & Q always valid, but just calculate the
    ! smaller one and use P+Q=1 for the other.
    IF (x .LE. s) THEN
      P = 0.5d0*ERFC(-eta*SQRT(0.5d0*s)) - Ra
      Q = 1-P
    ELSE
      Q = 0.5d0*ERFC(+eta*SQRT(0.5d0*s)) + Ra
      P = 1-Q
    END IF
    
  END SUBROUTINE
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Logarithm of regularized gamma function P of real arguments
! [double precision]
!   P(s,x) = \gamma(s,x) / \Gamma(s)
! 
ELEMENTAL FUNCTION LOG_GAMMA_P(s,x) RESULT(z)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s,x
  REAL*8 :: z
  REAL*8 :: lnP,lnQ
  CALL LOG_GAMMA_PQ(s,x,lnP,lnQ)
  z = lnP
END FUNCTION


!-----------------------------------------------------------------------
! Logarithm of regularized gamma function Q of real arguments
! [double precision]
!   Q(s,x) = \Gamma(s,x) / \Gamma(s)
! 
ELEMENTAL FUNCTION LOG_GAMMA_Q(s,x) RESULT(z)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s,x
  REAL*8 :: z
  REAL*8 :: lnP,lnQ
  CALL LOG_GAMMA_PQ(s,x,lnP,lnQ)
  z = lnQ
END FUNCTION


!-----------------------------------------------------------------------
! Logarithm of regularized incomplete gamma functions P, Q of real
! arguments [double precision].  These are defined as:
!   P(s,x) = \gamma(s,x) / \Gamma(s)
!   Q(s,x) = \Gamma(s,x) / \Gamma(s)
! where \gamma(s,x) & \Gamma(s,x) are the lower & upper incomplete
! gamma functions:
!   \gamma(s,x) = \int_0^x dt t^{s-1} e^{-t}
!   \Gamma(s,x) = \int_x^\infty dt t^{s-1} e^{-t}
! The following relations hold:
!   P(s,x) + Q(s,x) = 1
!   \gamma(s,x) + \Gamma(s,x) = \Gamma(s)
! The below routine calculates the approximately smaller of P or Q
! and uses P+Q=1 to determine the other.
! 
! NOTE: Valid only for s > 0 and x >= 0.
! 
ELEMENTAL SUBROUTINE LOG_GAMMA_PQ(s,x,lnP,lnQ)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: s,x
  REAL*8, INTENT(OUT) :: lnP,lnQ
  
  ! Special cases
  IF (x .EQ. 0d0) THEN
    lnP = -HUGE(1d0)
    lnQ = 0d0
    RETURN
  END IF
  
  ! Bad cases
  IF ((x .LT. 0d0) .OR. (s .LE. 0d0)) THEN
    lnP = -HUGE(1d0)
    lnQ = HUGE(1d0)
    RETURN
  END IF
  
  ! Calculate P,Q using uniform asymptotic expansion
  IF ((s .GE. 10d0) .AND. (x-10d0 .GT. 0.302d0*(s-10d0)) &
      .AND. (x-10d0 .LT. 2.357d0*(s-10d0))) THEN
    CALL LOG_GAMMA_PQ_UA(s,x,lnP,lnQ)
  ! Calculate P using Taylor series
  ELSE IF (x .LE. s) THEN
    lnP = LOG_GAMMA_P_TS(s,x)
    lnQ = LOGp1(-EXP(lnP))
  ! Calculate Q using continued fraction
  ELSE
    lnQ = LOG_GAMMA_Q_CF(s,x)
    lnP = LOGp1(-EXP(lnQ))
  END IF
  
  CONTAINS
  
  !---------------------------------------------
  ! Calculates gamma function quantity P(s,x) using a Taylor
  ! series expansion.
  !   P(s,x) = \gamma(s,x) / \Gamma(s)
  !   \gamma(s,x) = \int_0^x dt t^{s-1} e^{-t}
  ! NOTE: Requires 0 < s, 0 <= x.  Should only be used for x
  !       not (much) larger than s or convergence may be slow.
  !       Also best for s not large (timing scales as maybe
  !       sqrt(s)).
  PURE FUNCTION LOG_GAMMA_P_TS(s,x) RESULT(lnP)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: s,x
    REAL*8 :: lnP
    INTEGER :: K
    REAL*8 :: zk,sum
    ! P(s,x) = x^s e^{-x} \Sum_{k=0}^\infty x^k / \Gamma(s+k+1)
    ! sum here excludes first term (which is 1)
    K   = 1
    zk  = x / (s+1)
    sum = zk
    DO WHILE (zk .GT. EPSILON(zk))
      K   = K+1
      zk  = zk * x / (s+K)
      sum = sum + zk
      IF (K .GE. 10000) EXIT
    END DO
    lnP = s*LOG(x) - x - LOG_GAMMA(s+1) + LOGp1(sum)
  END FUNCTION
  
  !---------------------------------------------
  ! Calculates gamma function quantity Q(s,x) using a
  ! continued fraction.
  !   Q(s,x) = \Gamma(s,x) / \Gamma(s)
  !   \Gamma(s,x) = \int_x^{\infty} dt t^{s-1} e^{-t}
  ! NOTE: Requires 0 < x,s.  Should only be used for x not
  !       (much) smaller than s or convergence may be slow
  !       (or simply fail).  Also best for s not large
  !       (timing scales as maybe sqrt(s)).
  PURE FUNCTION LOG_GAMMA_Q_CF(s,x) RESULT(lnQ)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: s,x
    REAL*8 :: lnQ
    INTEGER :: K
    REAL*8 :: xs1,fk,Ck,Dk,Deltak
    ! Continued fraction (see Numerical Recipes):
    !   Q(s,x) = x^s e^{-x} / \Gamma(s)
    !            / (x-s+1 + K_k(-k(k-s),-s+2k+x+1)_{1}^{\infty})
    ! where
    !   K_k(a_k,b_k)_1 = a_1/b1+ a2/b2+ a3/b3+ ...
    ! is the Gauss notation for a continued fraction.
    xs1 = x-s+1
    K  = 0
    Ck = xs1
    Dk = 0
    fk = xs1
    Deltak = HUGE(x)
    DO WHILE (ABS(Deltak-1) .GE. EPSILON(x))
      K  = K+1
      Ck = (xs1+2*K) - K*(K-s)/Ck
      IF (Ck .EQ. 0d0) Ck = 1d-30
      Dk = (xs1+2*K) - K*(K-s)*Dk
      IF (Dk .EQ. 0d0) Dk = 1d-30
      Dk = 1/Dk
      Deltak = Ck*Dk
      fk = Deltak*fk
      IF (K .GE. 10000) EXIT
    END DO
    lnQ = s*LOG(x) - x - LOG_GAMMA(s) - LOG(fk)
  END FUNCTION
  
  !---------------------------------------------
  ! Calculates gamma function quantities P(s,x) and Q(s,x)
  ! using a uniform asymptotic expansion (see 1306.1754).
  ! Nearly constant evaluation time for any s & x, unlike
  ! the Taylor series or continued fraction algorithms.
  ! NOTE: Intended for s not small (larger than ~10) and
  !       0.30 < x/s < 2.35; accuracy starts to drop
  !       outside this region.
  ! This routine is about ~ 2 slower than the non-
  ! logarithmic version.
  PURE SUBROUTINE LOG_GAMMA_PQ_UA(s,x,lnP,lnQ)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: s,x
    REAL*8, INTENT(OUT) :: lnP,lnQ
    INTEGER :: K,sgnRa
    REAL*8 :: eta,lambda,u,lnSa,lnRa,sum
    REAL*8 :: beta(0:32)
    INTEGER, PARAMETER :: KMAX = 25
    REAL*8, PARAMETER :: D(0:31) = &
      (/  1.0000000000000000d+00,-3.3333333333333333d-01, 8.3333333333333333d-02, &
         -1.4814814814814815d-02, 1.1574074074074074d-03, 3.5273368606701940d-04, &
         -1.7875514403292181d-04, 3.9192631785224378d-05,-2.1854485106799922d-06, &
         -1.8540622107151600d-06, 8.2967113409530860d-07,-1.7665952736826079d-07, &
          6.7078535434014986d-09, 1.0261809784240308d-08,-4.3820360184533532d-09, &
          9.1476995822367902d-10,-2.5514193994946250d-11,-5.8307721325504251d-11, &
          2.4361948020667416d-11,-5.0276692801141756d-12, 1.1004392031956135d-13, &
          3.3717632624009854d-13,-1.3923887224181621d-13, 2.8534893807047443d-14, &
         -5.1391118342425726d-16,-1.9752288294349443d-15, 8.0995211567045613d-16, &
         -1.6522531216398162d-16, 2.5305430097478884d-18, 1.1686939738559577d-17, &
         -4.7700370498204848d-18, 9.6991260590562371d-19 /)
    REAL*8, PARAMETER :: TWOPI = 6.2831853071795865d0   ! 2*Pi
    REAL*8, PARAMETER :: LN2   = 0.69314718055994531d0  ! ln(2)
    
    ! Define ratio of x & s
    lambda = x/s
    
    ! eta is defined as 1/2 eta^2 = lambda - 1 - ln(lambda) with
    ! the same sign as lambda-1.  Use an expansion for lambda ~ 1
    ! to avoid cancellation issues.  Rapid convergent of below
    ! series in eta requires |eta| < 1, which corresponds to
    ! 0.3017 < lambda < 2.3577.
    u = lambda - 1
    IF (ABS(u) .GT. 0.01d0) THEN
      eta = SQRT(2*(lambda-1-LOG(lambda)))
    ELSE
      eta = u*(1+u*(-1d0/3d0+u*(7d0/36d0+u*(-73d0/540d0+u*(1331d0/12960d0)))))
    END IF
    IF (x .LT. s) eta = -eta
    
    ! Use erfc in eta*sqrt{s/2} as first approximation and calculate
    ! correction term:
    !   R_a(\eta) = e^{-1/2 s \eta^2} S_a(\eta) / sqrt{2\pi s}
    ! where Sa(eta) can be given by an expansion in eta:
    !   S_a(\eta) \approx a/(a+\beta_1) \Sum_{k=0}^N \beta_k \eta^k
    ! The beta terms are dependent upon s and can be calculated using
    ! a recursion relation.
    beta(KMAX+2) = 0d0
    beta(KMAX+1) = 0d0
    sum = 0d0
    DO K = KMAX,0,-1
      beta(K) = (K+2)/s * beta(K+2) + D(K+1)
      sum = beta(K) + eta*sum
    END DO
    IF (sum .GE. 0d0) THEN
      sgnRa = +1
    ELSE
      sgnRa = -1
    END IF
    ! Assuming b_1 > -s, which seems to be the case.
    ! Have not checked if LOG(sum) loses accuracy here....
    lnSa = LOG(ABS(sum)) - LOGp1(beta(1)/s)
    lnRa = -0.5d0*s*eta**2 + lnSa - 0.5d0*LOG(TWOPI*s)
    
    ! Formulas for P & Q always valid, but just calculate the
    ! smaller one and use P+Q=1 for the other.
    IF (x .LE. s) THEN
      IF (sgnRa .GT. 0) THEN
        lnP = LOG_DIFF(LOG_ERFC(-eta*SQRT(0.5d0*s))-LN2,lnRa)
      ELSE
        lnP = LOG_SUM(LOG_ERFC(-eta*SQRT(0.5d0*s))-LN2,lnRa)
      END IF
      lnQ = LOGp1(-EXP(lnP))
    ELSE
      IF (sgnRa .GT. 0) THEN
        lnQ = LOG_SUM(LOG_ERFC(+eta*SQRT(0.5d0*s))-LN2,lnRa)
      ELSE
        lnQ = LOG_DIFF(LOG_ERFC(+eta*SQRT(0.5d0*s))-LN2,lnRa)
      END IF
      lnP = LOGp1(-EXP(lnQ))
    END IF
    
  END SUBROUTINE
  
END SUBROUTINE


!-----------------------------------------------------------------------
! Binomial coefficient (n k) = n!/k!(n-k)! [integer]
! Various techniques are used to provide a fast calculation in
! different cases, but no check for integer overflow is performed
! for the final result.
! 
ELEMENTAL FUNCTION BINOMIAL_COEFF(N,K) RESULT(Cnk)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N,K
  INTEGER :: Cnk
  INTEGER :: K0,J
  ! Array of binomial coefficients C(n,k).
  ! For k <= n/2, coefficient is at index:
  !    J = floor[(n+1)/2]*floor[(n+2)/2] + k + 1
  ! For k > n/2, we will use C(n,k) = C(n,n-k).
  INTEGER, PARAMETER :: NMAX_TAB = 30    ! Coefficients for N in [0,NMAX_TAB]
  INTEGER, PARAMETER :: NC       = 256
  INTEGER, PARAMETER :: BC_ARRAY(1:NC) =                                &
    (/  1,1,1,2,1,3,1,4,6,1,5,10,1,6,15,20,1,7,21,35,1,8,28,56,70,1,9,  &  ! n=1-9
        36,84,126,1,10,45,120,210,252,1,11,55,165,330,462,1,12,66,220,  &  ! n=9-12
        495,792,924,1,13,78,286,715,1287,1716,1,14,91,364,1001,2002,    &  ! n=12-14
        3003,3432,1,15,105,455,1365,3003,5005,6435,1,16,120,560,1820,   &  ! n=14-16
        4368,8008,11440,12870,1,17,136,680,2380,6188,12376,19448,24310, &  ! n=16-17
        1,18,153,816,3060,8568,18564,31824,43758,48620,1,19,171,969,    &  ! n=18-19
        3876,11628,27132,50388,75582,92378,1,20,190,1140,4845,15504,    &  ! n=19-20
        38760,77520,125970,167960,184756,1,21,210,1330,5985,20349,      &  ! n=20-21
        54264,116280,203490,293930,352716,1,22,231,1540,7315,26334,     &  ! n=21-22
        74613,170544,319770,497420,646646,705432,1,23,253,1771,8855,    &  ! n=22-23
        33649,100947,245157,490314,817190,1144066,1352078,1,24,276,     &  ! n=23-24
        2024,10626,42504,134596,346104,735471,1307504,1961256,2496144,  &  ! n=24
        2704156,1,25,300,2300,12650,53130,177100,480700,1081575,        &  ! n=24-25
        2042975,3268760,4457400,5200300,1,26,325,2600,14950,65780,      &  ! n=25-26
        230230,657800,1562275,3124550,5311735,7726160,9657700,10400600, &  ! n=26
        1,27,351,2925,17550,80730,296010,888030,2220075,4686825,        &  ! n=27
        8436285,13037895,17383860,20058300,1,28,378,3276,20475,98280,   &  ! n=27-28
        376740,1184040,3108105,6906900,13123110,21474180,30421755,      &  ! n=28
        37442160,40116600,1,29,406,3654,23751,118755,475020,1560780,    &  ! n=28-29
        4292145,10015005,20030010,34597290,51895935,67863915,77558760,  &  ! n=29
        1,30,435,4060,27405,142506,593775,2035800,5852925,14307150,     &  ! n=30
        30045015,54627300,86493225,119759850,145422675,155117520 /)        ! n=30
  ! Largest N at each K0 = min[K,N-K] for which recursive method
  ! below will not lead to overflow for 4-byte integers.
  INTEGER, PARAMETER :: K0MAX = 14
  INTEGER, PARAMETER :: NMAX_AT_K0(0:K0MAX) =                           &
    (/  HUGE(1),HUGE(1),46341,1626,338,140,82,58,46,39,35,33,31,30,30 /)
  ! Bad cases
  IF ((N .LT. K) .OR. (K .LT. 0)) THEN
    Cnk = -HUGE(N)
    !STOP 'ERROR: invalid binomial coefficient arguments (0 <= k <= n)'
    RETURN
  END IF
  ! Use (n k) = (n n-k)
  K0 = MIN(K,N-K)
  ! Use tabulated values for small n
  IF (N .LE. NMAX_TAB) THEN
    ! Parentheses are important: act as floor function
    J = ((N+1)/2) * ((N+2)/2) + K0 + 1
    Cnk = BC_ARRAY(J)
    RETURN
  ! Use recursion: (n k) = n/k (n-1 k-1)
  ! Note overflow for large enough k and/or n.
  ! Below check should avoid all potential overflow cases for
  ! 4-byte integers.
  ELSE IF (K0 .LE. K0MAX) THEN
    IF (N .LE. NMAX_AT_K0(K0)) THEN
      Cnk = 1
      DO J=1,K0
        Cnk = (Cnk*(N-K0+J))/J
      END DO
      RETURN
    END IF
  END IF
  ! Use LogGamma for large arguments.
  ! Seems to give full INTEGER*4 precision when not in an
  ! overflow condition.
  Cnk = NINT(EXP(LOG_GAMMA(N+1d0) - LOG_GAMMA(K+1d0) - LOG_GAMMA(N-K+1d0)))
END FUNCTION


!-----------------------------------------------------------------------
! Logarithm of binomial coefficient (n k) = n!/k!(n-k)! [double precision]
! 
ELEMENTAL FUNCTION LOG_BINOMIAL_COEFF(N,K) RESULT(lnCnk)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N,K
  REAL*8 :: lnCnk
  INTEGER :: Cnk,K0,J
  INTEGER, PARAMETER :: NMAX_TAB = 30
  ! Largest N at each K0 = min[K,N-K] for which recursive method
  ! below will not lead to overflow for 4-byte integers.
  INTEGER, PARAMETER :: K0MAX = 14
  INTEGER, PARAMETER :: NMAX_AT_K0(0:K0MAX) =                           &
    (/  HUGE(1),HUGE(1),46341,1626,338,140,82,58,46,39,35,33,31,30,30 /)
  ! Bad cases
  IF ((N .LT. K) .OR. (K .LT. 0)) THEN
    lnCnk = -HUGE(N)
    !STOP 'ERROR: invalid binomial coefficient arguments (0 <= k <= n)'
    RETURN
  END IF
  ! Use (n k) = (n n-k)
  K0 = MIN(K,N-K)
  ! Get tabulated value from other routine
  IF (N .LE. NMAX_TAB) THEN
    lnCnk = LOG(1d0*BINOMIAL_COEFF(N,K))
    RETURN
  ! Use recursion: (n k) = n/k (n-1 k-1)
  ! Note overflow for large enough k and/or n.
  ! Below check should avoid all potential overflow cases for
  ! 4-byte integers.
  ELSE IF (K0 .LE. K0MAX) THEN
    IF (N .LE. NMAX_AT_K0(K0)) THEN
      Cnk = 1
      DO J=1,K0
        Cnk = (Cnk*(N-K0+J))/J
      END DO
      lnCnk = LOG(1d0*Cnk)
      RETURN
    END IF
  END IF
  ! Use LogGamma for large arguments.
  lnCnk = LOG_GAMMA(N+1d0) - LOG_GAMMA(K+1d0) - LOG_GAMMA(N-K+1d0)
END FUNCTION


!-----------------------------------------------------------------------
! Calculates W(x) to double precision, where W is the Lambert W
! function defined as the solution w = W(x) to the equation:
!     x = w e^w
! Returns the principal branch value (the larger of two solutions
! over -1/e < x < 0; there is only one solution for x > 0).  The
! W function is undefined for x < -1/e.
! 
! Valid input:  -1/e < x < \infty
! 
ELEMENTAL FUNCTION LAMBERTW(x) RESULT(w)
  IMPLICIT NONE
  REAL*8 :: w
  REAL*8, INTENT(IN) :: x
  REAL*8 :: epsk,zk,qk,wk,wk1,p,num,den,a,b,ia
  REAL*8, PARAMETER :: e    = 2.7182818284590452d0
  REAL*8, PARAMETER :: einv = 0.36787944117144232d0
  
  ! No solution for x < 1/e
  IF (x .LT. -einv) THEN
    !STOP 'Error in lambertw: argument must be larger than -EXP(-1)'
    w = -HUGE(w)
    RETURN
  END IF
  
  ! Could use Newton's or Halley's iteration method to find the
  ! solution, but the Lambert function has a faster method,
  ! Fritsch's iteration:
  !    W_{k+1} = W_k * (1 + eps_k)
  ! with
  !    eps_k   = z_k/(1 + W_k) * (q_k - z_k)/(q_k - 2*z_k)
  !    z_k     = ln(x/W_k) - W_k
  !    q_k     = 2 (1 + W_k) (1 + W_k + (2/3)*z_k)
  ! If eps_k is the error in W_k, then the error in W_{k+1} is of
  ! order (eps_k)^4, a faster convergent that Halley's method.
  ! For a first guess accurate to order O(10^-4), double precision
  ! can be achieved in only a single iteration.  We use estimates
  ! for the initial guess as determined by Veberic.
  ! 
  ! For further information, see:
  !   D. Veberic, arXiv:1003.1628.
  !   F.N. Fritsch, R.E. Shafer and W.P. Crowley,
  !       Commun. ACM 16, 123 (1973).
  
  ! Initial estimate by Veberic
  IF (x .LT. -0.32358170806015724d0) THEN
    ! branch point expansion
    p = SQRT(2*(1+e*x))
    wk = -1 + p*(1 + p*(-1d0/3 + p*(11d0/72 + p*(-43d0/540              &
              + p*(769d0/17280 + p*(-221d0/8505                         &
              + p*(680863d0/43545600 + p*(-1963d0/204120                &
              + p*(226287557d0/37623398400d0) ))))))))
  ELSE IF (x .LT. 0.14546954290661823d0) THEN
    ! rational fit
    num = x * (1 + x *                                                  &
              (5.931375839364438d0 + x *                                &
              (11.392205505329132d0+ x *                                &
              (7.338883399111118d0 + x * 0.6534490169919599d0) )))
    den = 1 + x *                                                       &
             (6.931373689597704d0 + x *                                 &
             (16.82349461388016d0 + x *                                 &
             (16.43072324143226d0 + x * 5.115235195211697d0) ))
    wk = num / den
  ELSE IF (x .LT. 8.706658967856612d0) THEN
    ! rational fit
    num = x * (1 + x *                                                  &
              (2.4450530707265568d0 + x *                               &
              (1.3436642259582265d0 + x *                               &
              (0.14844005539759195d0+ x * 8.047501729129999d-4) )))
    den = 1 + x *                                                       &
             (3.4447089864860025d0 + x *                                &
             (3.2924898573719523d0 + x *                                &
             (0.9164600188031222d0 + x * 0.05306864044833221d0) ))
    wk = num / den
  ELSE
    ! asymptotic expansion
    a = LOG(x)
    b = LOG(a)
    ia = 1/a
    wk = a - b + b * ia *                                               &
                (1 + ia *                                               &
                (0.5d0*(-2 + b) + ia *                                  &
                (1/6d0*(6 + b*(-9 + b*2)) + ia *                        &
                (1/12d0*(-12 + b*(36 + b*(-22 + b*3))) + ia *           &
                 1/60d0*(60 + b*(-300 + b*(350 + b*(-125 + b*12))))     &
                ))))
  END IF
  
  ! Special cases:
  ! For x equal to 0 or -1/e, the Fritsch iteration does
  ! not work as some of the terms go to infinity.  However,
  ! for x sufficiently near 0 or -1/e, the above first
  ! approximation is already nearly double precision.
  IF ((ABS(x) .LT. 1d-7) .OR. (x .LT. -einv+1d-6) ) THEN
    w = wk
    RETURN
  END IF
  
  ! Now apply Fritsch iteration
  wk1  = wk + 1
  zk   = LOG(x/wk) - wk
  qk   = 2 * wk1 * (wk1 + (2d0/3)*zk)
  epsk = zk * (qk - zk) / (wk1 * (qk - 2*zk))
  wk   = wk * (1 + epsk)
  ! In most cases, no further iterations will be necessary
  DO WHILE (ABS(epsk) .GT. 1d-5)
    wk1  = wk + 1
    zk   = LOG(x/wk) - wk
    qk   = 2 * wk1 * (wk1 + (2d0/3)*zk)
    epsk = zk * (qk - zk) / (wk1 * (qk - 2*zk))
    wk   = wk * (1 + epsk)
  END DO
  
  w = wk
  
END FUNCTION


!-----------------------------------------------------------------------
! Calculates W(x) to double precision, where W is the Lambert W
! function defined as the solution w = W(x) to the equation:
!     x = w e^w
! Returns the secondary branch value (the smaller of two solutions
! over -1/e < x < 0).  This branch is defined only for
! -1/e < x < 0.
! 
! Valid input:  -1/e < x < 0
! 
ELEMENTAL FUNCTION LAMBERTW2(x) RESULT(w2)
  IMPLICIT NONE
  REAL*8 :: w2
  REAL*8, INTENT(IN) :: x
  INTEGER :: k
  REAL*8 :: epsk,zk,qk,wk,wk1,p,num,den,a
  REAL*8, PARAMETER :: e    = 2.7182818284590452d0
  REAL*8, PARAMETER :: einv = 0.36787944117144232d0
  
  ! No solution for x < 1/e
  IF (x .LT. -einv) THEN
    !STOP 'Error in lambertw2: argument must be larger than -EXP(-1)'
    w2 = HUGE(w2)
    RETURN
  END IF
  
  ! No second branch for x > 0
  IF (x .GT. 0) THEN
    !STOP 'Error in lambertw2: argument must be smaller than 0'
    w2 = HUGE(w2)
    RETURN
  END IF
  
  ! Limit as x ->0 is -\infty
  IF (x .EQ. 0) THEN
    w2 = -HUGE(w2)
    RETURN
  END IF
  
  ! Could use Newton's or Halley's iteration method to find the
  ! solution, but the Lambert function has a faster method,
  ! Fritsch's iteration:
  !    W_{k+1} = W_k * (1 + eps_k)
  ! with
  !    eps_k   = z_k/(1 + W_k) * (q_k - z_k)/(q_k - 2*z_k)
  !    z_k     = ln(x/W_k) - W_k
  !    q_k     = 2 (1 + W_k) (1 + W_k + (2/3)*z_k)
  ! If eps_k is the error in W_k, then the error in W_{k+1} is of
  ! order (eps_k)^4, a faster convergent that Halley's method.
  ! For a first guess accurate to order O(10^-4), double precision
  ! can be achieved in only a single iteration.  We use estimates
  ! for the initial guess as determined by Veberic.
  ! 
  ! For further information, see:
  !   D. Veberic, arXiv:1003.1628.
  !   F.N. Fritsch, R.E. Shafer and W.P. Crowley,
  !       Commun. ACM 16, 123 (1973).
  
  ! Initial estimate by Veberic
  IF (x .LT. -0.30298541769d0) THEN
    ! branch point expansion
    p = -SQRT(2*(1+e*x))
    wk = -1 + p*(1 + p*(-1d0/3 + p*(11d0/72 + p*(-43d0/540              &
              + p*(769d0/17280 + p*(-221d0/8505                         &
              + p*(680863d0/43545600 + p*(-1963d0/204120                &
              + p*(226287557d0/37623398400d0) ))))))))
  ELSE IF (x .LT. -0.051012917658221676d0) THEN
    ! rational fit
    num = -7.814176723907436d0 + x *                                    &
            (253.88810188892484d0 + x * 657.9493176902304d0)
    den = 1 + x *                                                       &
             (-60.43958713690808d0+ x *                                 &
             (99.98567083107612d0 + x *                                 &
             (682.6073999909428d0 + x *                                 &
             (962.1784396969866d0 + x * 1477.9341280760887d0) )))
    wk = num / den
  ELSE
    ! continued logarithm
    a  = LOG(-x)
    wk = a
    DO k=1,9
      wk = a - LOG(-wk)
    END DO
  END IF
  
  ! Special cases:
  ! For x equal to -1/e, the Fritsch iteration does not
  ! not work as some of the terms go to infinity.  However,
  ! for x sufficiently near -1/e, the above first
  ! approximation is already nearly double precision.
  IF (x .LT. -einv+1d-6) THEN
    w2 = wk
    RETURN
  END IF
  
  ! Now apply Fritsch iteration
  wk1  = wk + 1
  zk   = LOG(x/wk) - wk
  qk   = 2 * wk1 * (wk1 + (2d0/3)*zk)
  epsk = zk * (qk - zk) / (wk1 * (qk - 2*zk))
  wk   = wk * (1 + epsk)
  ! In most cases, no further iterations will be necessary
  DO WHILE (ABS(epsk) .GT. 1d-5)
    wk1  = wk + 1
    zk   = LOG(x/wk) - wk
    qk   = 2 * wk1 * (wk1 + (2d0/3)*zk)
    epsk = zk * (qk - zk) / (wk1 * (qk - 2*zk))
    wk   = wk * (1 + epsk)
  END DO
  
  w2 = wk
  
END FUNCTION


!-----------------------------------------------------------------------
! Calculates W(x=e^ln(x)) to double precision, where W is the Lambert W
! function defined as the solution w = W(x) to the equation:
!     x = w e^w
! Returns the principal branch value (the larger of two solutions
! over -1/e < x < 0; there is only one solution for x > 0).  The
! W function is undefined for x < -1/e.
! 
! This version takes ln(x) as the input to allow for cases where x
! is large.
! 
! Valid input:  -\infty < lnx < \infty
! 
ELEMENTAL FUNCTION LAMBERTWLN(lnx) RESULT(w)
  IMPLICIT NONE
  REAL*8 :: w
  REAL*8, INTENT(IN) :: lnx
  REAL*8 :: epsk,zk,qk,wk,wk1,a,b,ia
  REAL*8, PARAMETER :: e    = 2.7182818284590452d0
  REAL*8, PARAMETER :: einv = 0.36787944117144232d0
  
  ! Here, we only calculate W(x) for very large x.  If x is a
  ! not very large, we use the lambertw routine.
  IF (lnx .LT. 300d0) THEN
    w = LAMBERTW(EXP(lnx))
    RETURN
  END IF
  
  ! Could use Newton's or Halley's iteration method to find the
  ! solution, but the Lambert function has a faster method,
  ! Fritsch's iteration:
  !    W_{k+1} = W_k * (1 + eps_k)
  ! with
  !    eps_k   = z_k/(1 + W_k) * (q_k - z_k)/(q_k - 2*z_k)
  !    z_k     = ln(x/W_k) - W_k
  !    q_k     = 2 (1 + W_k) (1 + W_k + (2/3)*z_k)
  ! If eps_k is the error in W_k, then the error in W_{k+1} is of
  ! order (eps_k)^4, a faster convergent that Halley's method.
  ! For a first guess accurate to order O(10^-4), double precision
  ! can be achieved in only a single iteration.  We use estimates
  ! for the initial guess as determined by Veberic.
  ! 
  ! For further information, see:
  !   D. Veberic, arXiv:1003.1628.
  !   F.N. Fritsch, R.E. Shafer and W.P. Crowley,
  !       Commun. ACM 16, 123 (1973).
  
  ! Initial estimate by Veberic
  ! asymptotic expansion
  a = lnx
  b = LOG(a)
  ia = 1/a
  wk = a - b + b * ia *                                                 &
              (1 + ia *                                                 &
              (0.5d0*(-2 + b) + ia *                                    &
              (1/6d0*(6 + b*(-9 + b*2)) + ia *                          &
              (1/12d0*(-12 + b*(36 + b*(-22 + b*3))) + ia *             &
               1/60d0*(60 + b*(-300 + b*(350 + b*(-125 + b*12))))       &
              ))))
  
  ! Now apply Fritsch iteration
  wk1  = wk + 1
  zk   = lnx - LOG(wk) - wk
  qk   = 2 * wk1 * (wk1 + (2d0/3)*zk)
  epsk = zk * (qk - zk) / (wk1 * (qk - 2*zk))
  wk   = wk * (1 + epsk)
  ! In most cases, no further iterations will be necessary
  DO WHILE (ABS(epsk) .GT. 1d-5)
    wk1  = wk + 1
    zk   = lnx - LOG(wk) - wk
    qk   = 2 * wk1 * (wk1 + (2d0/3)*zk)
    epsk = zk * (qk - zk) / (wk1 * (qk - 2*zk))
    wk   = wk * (1 + epsk)
  END DO
  
  w = wk
  
END FUNCTION



!=======================================================================
! PROBABILITY DISTRIBUTIONS
!=======================================================================

!----------UNIFORM DISTRIBUTION-----------------------------------------
! Uniform distribution has probability:
!   P(x|xmin,xmax) = 1/(xmax-xmin)    xmin < x < xmax
!                  = 0                otherwise
! The range of the distribution [xmin,xmax] has default [0,1].
!-----------------------------------------------------------------------

!----------------------------------------
! Probability for uniform distribution
!   x           Value distributed normally
!   xmin        Minimum of uniform values (optional, default is 0)
!   xmax        Maximum of uniform values (optional, default is 1)
ELEMENTAL FUNCTION UniformPDF(x,xmin,xmax) RESULT(pdf)
  IMPLICIT NONE
  REAL*8 :: pdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax
  REAL*8 :: xmin0,xmax0
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = 1d0
  END IF
  IF ((x .GE. xmin0) .AND. (x .LE. xmax0)) THEN
    pdf = 1d0 / (xmax0-xmin0)
  ELSE
    pdf = 0d0
  END IF
END FUNCTION UniformPDF

!----------------------------------------
! Log of probability for uniform distribution
!   x           Value distributed normally
!   xmin        Minimum of uniform values (optional, default is 0)
!   xmax        Maximum of uniform values (optional, default is 1)
ELEMENTAL FUNCTION UniformLogPDF(x,xmin,xmax) RESULT(lnpdf)
  IMPLICIT NONE
  REAL*8 :: lnpdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax
  REAL*8 :: xmin0,xmax0
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = 1d0
  END IF
  IF ((x .GE. xmin0) .AND. (x .LE. xmax0)) THEN
    lnpdf = LOG(1d0 / (xmax0-xmin0))
  ELSE
    lnpdf = -HUGE(1d0)
  END IF
END FUNCTION UniformLogPDF

!----------------------------------------
! CDF for uniform distribution
!   x           Value distributed normally
!   xmin        Minimum of uniform values (optional, default is 0)
!   xmax        Maximum of uniform values (optional, default is 1)
ELEMENTAL FUNCTION UniformCDF(x,xmin,xmax) RESULT(cdf)
  IMPLICIT NONE
  REAL*8 :: cdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax
  REAL*8 :: xmin0,xmax0
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = 1d0
  END IF
  IF (x .LE. xmin0) THEN
    cdf = 0d0
  ELSE IF (x .GE. xmax0) THEN
    cdf = 1d0
  ELSE
    cdf = (x-xmin0) / (xmax0-xmin0)
  END IF
END FUNCTION UniformCDF

!----------------------------------------
! 1-CDF for uniform distribution
!   x           Value distributed normally
!   xmin        Minimum of uniform values (optional, default is 0)
!   xmax        Maximum of uniform values (optional, default is 1)
ELEMENTAL FUNCTION Uniform1mCDF(x,xmin,xmax) RESULT(p)
  IMPLICIT NONE
  REAL*8 :: p
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax
  REAL*8 :: xmin0,xmax0
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = 1d0
  END IF
  IF (x .LE. xmin0) THEN
    p = 1d0
  ELSE IF (x .GE. xmax0) THEN
    p = 0d0
  ELSE
    p = (xmax-x) / (xmax0-xmin0)
  END IF
END FUNCTION Uniform1mCDF

!----------------------------------------
! Random number for uniform distribution
!   xmin        Minimum of uniform values (optional, default is 0)
!   xmax        Maximum of uniform values (optional, default is 1)
FUNCTION RandUniform(xmin,xmax) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax
  REAL*8 :: urand
  urand = rand()
  IF (PRESENT(xmin)) THEN
    IF (PRESENT(xmax)) THEN
      x = xmin + urand*(xmax-xmin)
    ELSE
      x = xmin + urand*(1d0-xmin)
    END IF
  ELSE
    IF (PRESENT(xmax)) THEN
      x = urand*xmax
    ELSE
      x = urand
    END IF
  END IF
END FUNCTION RandUniform


!----------NORMAL DISTRIBUTION------------------------------------------
! Normal distribution has probability:
!   P(x|mu,sigma) = EXP(-(x-mu)^2/(2 sigma^2)) / SQRT(2 PI sigma^2)
! The average mu has default 0 and the standard deviation sigma
! has default 1.
!-----------------------------------------------------------------------

!----------------------------------------
! Probability for normal distribution
!   x           Value distributed normally
!   mu          Average (optional, default is 0)
!   sigma       Standard deviation (optional, default is 1)
ELEMENTAL FUNCTION NormalPDF(x,mu,sigma) RESULT(pdf)
  IMPLICIT NONE
  REAL*8 :: pdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: z
  REAL*8, PARAMETER :: SQRT2PI = 2.5066282746310005d0
  IF (PRESENT(mu)) THEN
    z = x - mu
  ELSE
    z = x
  END IF
  IF (PRESENT(sigma)) THEN
    pdf = EXP(-z**2/(2*sigma**2)) / (SQRT2PI*sigma)
  ELSE
    pdf = EXP(-z**2/2) / SQRT2PI
  END IF
END FUNCTION NormalPDF

!----------------------------------------
! Log of probability for normal distribution
!   x           Value distributed normally
!   mu          Average (optional, default is 0)
!   sigma       Standard deviation (optional, default is 1)
ELEMENTAL FUNCTION NormalLogPDF(x,mu,sigma) RESULT(lnpdf)
  IMPLICIT NONE
  REAL*8 :: lnpdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: z
  REAL*8, PARAMETER :: LOGSQRT2PI = 0.91893853320467274d0
  IF (PRESENT(mu)) THEN
    z = x - mu
  ELSE
    z = x
  END IF
  IF (PRESENT(sigma)) THEN
    lnpdf = -z**2/(2*sigma**2) - LOGSQRT2PI - LOG(sigma)
  ELSE
    lnpdf = -z**2/2 - LOGSQRT2PI
  END IF
END FUNCTION NormalLogPDF

!----------------------------------------
! CDF for normal distribution
!   x           Value distributed normally
!   mu          Average (optional, default is 0)
!   sigma       Standard deviation (optional, default is 1)
ELEMENTAL FUNCTION NormalCDF(x,mu,sigma) RESULT(cdf)
  IMPLICIT NONE
  REAL*8 :: cdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: z
  REAL*8, PARAMETER :: SQRT2 = 1.4142135623730950d0
  IF (PRESENT(mu)) THEN
    z = x - mu
  ELSE
    z = x
  END IF
  IF (PRESENT(sigma)) THEN
    z = z / (SQRT2*sigma)
  ELSE
    z = z / SQRT2
  END IF
  cdf = 0.5d0 * ERFC(-z)
END FUNCTION NormalCDF

!----------------------------------------
! 1-CDF for normal distribution
!   x           Value distributed normally
!   mu          Average (optional, default is 0)
!   sigma       Standard deviation (optional, default is 1)
ELEMENTAL FUNCTION Normal1mCDF(x,mu,sigma) RESULT(p)
  IMPLICIT NONE
  REAL*8 :: p
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: z
  REAL*8, PARAMETER :: SQRT2 = 1.4142135623730950d0
  IF (PRESENT(mu)) THEN
    z = x - mu
  ELSE
    z = x
  END IF
  IF (PRESENT(sigma)) THEN
    z = z / (SQRT2*sigma)
  ELSE
    z = z / SQRT2
  END IF
  p = 0.5d0 * ERFC(z)
END FUNCTION Normal1mCDF

!----------------------------------------
! Random number for normal distribution
!   mu          Average (optional, default is 0)
!   sigma       Standard deviation (optional, default is 1)
FUNCTION RandNormal(mu,sigma) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: u,v,r2,w,urand(2),rand1
  REAL*8, SAVE :: rand2
  LOGICAL, SAVE :: needs_calc = .TRUE.
  ! Algorithm calculates two random normal numbers at a time.
  ! Store second for next call to this function.
  IF (needs_calc) THEN
    r2 = 1d0
    DO WHILE (r2 .GE. 1d0)
      urand(1) = rand()
      urand(2) = rand()
      u = 2d0*urand(1) - 1d0
      v = 2d0*urand(2) - 1d0
      r2 = u**2 + v**2
    END DO
    w = SQRT(-2d0 * LOG(r2) / r2)
    rand1 = u * w
    rand2 = v * w
    x = rand1
    needs_calc = .FALSE.
  ELSE
    x = rand2
    needs_calc = .TRUE.
  END IF
  IF (PRESENT(sigma)) x = sigma*x
  IF (PRESENT(mu)) x = mu + x
END FUNCTION RandNormal


!----------LOG NORMAL DISTRIBUTION--------------------------------------
! Normal distribution has probability:
!   P(x|mu,sigma) = EXP(-(ln(x)-mu)^2/(2 sigma^2)) / SQRT(2 PI sigma^2) / x
! The location mu has default 0 and the scale sigma has default 1.
! Note the parameters mu and sigma are not the mean and standard
! deviation.
!   mean:      muL      = exp(mu + sigma^2/2)
!   variance:  sigmaL^2 = (exp(sigma^2)-1) exp(2*mu+sigma^2)
!       =>     mu      = -1/2 log[(1+sigmaL^2/muL^2)/muL^2]
!              sigma^2 = log[1+sigmaL^2/muL^2]
! A normal distribution with parameters sigmaN << muN is asymptotically
! similar to a log normal distribution with parameters:
!   mu    = ln(muN)
!   sigma = sigmaN/muN
! Specifically, setting muL = muN and sigmaL = sigmaN:
!   mu      = log(muN) - 1/2 log(1+sigmaN^2/muN^2)
!           = log(muN) + O[sigmaN^2/muN^2]
!   sigma^2 = log(1+sigmaN^2/muN^2)
!           = sigmaN^2/muN^2 (1 + O[sigmaN^2/muN^2])
!-----------------------------------------------------------------------

!----------------------------------------
! Probability for log normal distribution
!   x           Value distributed log-normally
!   mu          Location (optional, default is 0)
!   sigma       Scale (optional, default is 1)
ELEMENTAL FUNCTION LogNormalPDF(x,mu,sigma) RESULT(pdf)
  IMPLICIT NONE
  REAL*8 :: pdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: z
  REAL*8, PARAMETER :: SQRT2PI = 2.5066282746310005d0
  ! Bad case
  IF (x .LE. 0d0) THEN
    pdf = 0d0
    RETURN
  END IF
  IF (PRESENT(mu)) THEN
    z = LOG(x) - mu
  ELSE
    z = LOG(x)
  END IF
  IF (PRESENT(sigma)) THEN
    pdf = EXP(-z**2/(2*sigma**2)) / (SQRT2PI*sigma*x)
  ELSE
    pdf = EXP(-z**2/2) / (SQRT2PI*x)
  END IF
END FUNCTION LogNormalPDF

!----------------------------------------
! Log of probability for log normal distribution
!   x           Value distributed normally
!   mu          Average (optional, default is 0)
!   sigma       Standard deviation (optional, default is 1)
ELEMENTAL FUNCTION LogNormalLogPDF(x,mu,sigma) RESULT(lnpdf)
  IMPLICIT NONE
  REAL*8 :: lnpdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: z
  REAL*8, PARAMETER :: LOGSQRT2PI = 0.91893853320467274d0
  ! Bad case
  IF (x .LE. 0d0) THEN
    lnpdf = -HUGE(1d0)
    RETURN
  END IF
  IF (PRESENT(mu)) THEN
    z = LOG(x) - mu
  ELSE
    z = LOG(x)
  END IF
  IF (PRESENT(sigma)) THEN
    lnpdf = -z**2/(2*sigma**2) - LOGSQRT2PI - LOG(x*sigma)
  ELSE
    lnpdf = -z**2/2 - LOGSQRT2PI - LOG(x)
  END IF
END FUNCTION LogNormalLogPDF

!----------------------------------------
! CDF for log normal distribution
!   x           Value distributed log-normally
!   mu          Location (optional, default is 0)
!   sigma       Scale (optional, default is 1)
ELEMENTAL FUNCTION LogNormalCDF(x,mu,sigma) RESULT(cdf)
  IMPLICIT NONE
  REAL*8 :: cdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: z
  REAL*8, PARAMETER :: SQRT2 = 1.4142135623730950d0
  ! Bad case
  IF (x .LE. 0d0) THEN
    cdf = 0d0
    RETURN
  END IF
  IF (PRESENT(mu)) THEN
    z = LOG(x) - mu
  ELSE
    z = LOG(x)
  END IF
  IF (PRESENT(sigma)) THEN
    z = LOG(z) / (SQRT2*sigma)
  ELSE
    z = LOG(z) / SQRT2
  END IF
  cdf = 0.5d0 * ERFC(-z)
END FUNCTION LogNormalCDF

!----------------------------------------
! 1-CDF for log normal distribution
!   x           Value distributed log-normally
!   mu          Location (optional, default is 0)
!   sigma       Scale (optional, default is 1)
ELEMENTAL FUNCTION LogNormal1mCDF(x,mu,sigma) RESULT(p)
  IMPLICIT NONE
  REAL*8 :: p
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: z
  REAL*8, PARAMETER :: SQRT2 = 1.4142135623730950d0
  ! Bad case
  IF (x .LE. 0d0) THEN
    p = 1d0
    RETURN
  END IF
  IF (PRESENT(mu)) THEN
    z = LOG(x) - mu
  ELSE
    z = LOG(x)
  END IF
  IF (PRESENT(sigma)) THEN
    z = LOG(z) / (SQRT2*sigma)
  ELSE
    z = LOG(z) / SQRT2
  END IF
  p = 0.5d0 * ERFC(z)
END FUNCTION LogNormal1mCDF

!----------------------------------------
! Random number for log normal distribution
!   mu          Location (optional, default is 0)
!   sigma       Scale (optional, default is 1)
FUNCTION RandLogNormal(mu,sigma) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN), OPTIONAL :: mu,sigma
  REAL*8 :: xn
  xn = RandNormal()
  IF (PRESENT(sigma)) xn = sigma*xn
  IF (PRESENT(mu)) xn = mu + xn
  x = EXP(xn)
END FUNCTION RandLogNormal


!----------EXPONENTIAL DISTRIBUTION-------------------------------------
! Exponential distribution has probability:
!   P(x|xmin,xmax,xs) = A EXP(-x/xs)    xmin < x < xmax
!                     = 0               otherwise
! where the constant 'A' is a normalization factor.
! The range of the distribution [xmin,xmax] has default [0,infinity].
! The scale of the exponential xs has default 1.
!-----------------------------------------------------------------------

!----------------------------------------
! Probability for exponential distribution
!   x           Value distributed exponentially
!   xmin        Minimum value (optional, default is 0)
!   xmax        Maximum value (optional, default is infinity)
!   xs          The exponential scale (optional, default is 1)
ELEMENTAL FUNCTION ExponentialPDF(x,xmin,xmax,xs) RESULT(pdf)
  IMPLICIT NONE
  REAL*8 :: pdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax,xs
  REAL*8 :: xmin0,xmax0
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = HUGE(1d0)
  END IF
  IF ((x .LT. xmin0) .OR. (x .GT. xmax0)) THEN
    pdf = 0d0
    RETURN
  END IF
  IF (PRESENT(xs)) THEN
    IF (xs .EQ. 0d0) THEN
      pdf = 1d0 / (xmax0-xmin0)
    ELSE
      pdf = EXP((xmin0-x)/xs) / (1d0 - EXP((xmin0-xmax0)/xs)) / xs
    END IF
  ELSE
    pdf = EXP(xmin0-x) / (1d0 - EXP(xmin0-xmax0))
  END IF
END FUNCTION ExponentialPDF

!----------------------------------------
! Log of probability for exponential distribution
!   x           Value distributed exponentially
!   xmin        Minimum value (optional, default is 0)
!   xmax        Maximum value (optional, default is infinity)
!   xs          The exponential scale (optional, default is 1)
ELEMENTAL FUNCTION ExponentialLogPDF(x,xmin,xmax,xs) RESULT(lnpdf)
  IMPLICIT NONE
  REAL*8 :: lnpdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax,xs
  REAL*8 :: xmin0,xmax0
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = HUGE(1d0)
  END IF
  IF ((x .LT. xmin0) .OR. (x .GT. xmax0)) THEN
    lnpdf = -HUGE(1d0)
    RETURN
  END IF
  IF (PRESENT(xs)) THEN
    IF (xs .EQ. 0d0) THEN
      lnpdf = LOG(1d0 / (xmax0-xmin0))
    ELSE
      lnpdf = (xmin0-x)/xs - LOGp1(-EXP((xmin0-xmax0)/xs)) - LOG(xs)
    END IF
  ELSE
    lnpdf = (xmin0-x) - LOGp1(-EXP(xmin0-xmax0))
  END IF
END FUNCTION ExponentialLogPDF

!----------------------------------------
! CDF for exponential distribution
!   x           Value distributed exponentially
!   xmin        Minimum value (optional, default is 0)
!   xmax        Maximum value (optional, default is infinity)
!   xs          The exponential scale (optional, default is 1)
ELEMENTAL FUNCTION ExponentialCDF(x,xmin,xmax,xs) RESULT(cdf)
  IMPLICIT NONE
  REAL*8 :: cdf
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax,xs
  REAL*8 :: xmin0,xmax0
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = HUGE(1d0)
  END IF
  IF (x .LE. xmin0) THEN
    cdf = 0d0
  ELSE IF (x .GE. xmax0) THEN
    cdf = 1d0
  ELSE
    IF (PRESENT(xs)) THEN
      IF (xs .EQ. 0d0) THEN
        cdf = (x-xmin0) / (xmax0-xmin0)
      ELSE
        !cdf = (1d0 - EXP((xmin0-x)/xs)) / (1d0 - EXP((xmin0-xmax0)/xs))
        cdf = EXPm1((xmin0-x)/xs) / EXPm1((xmin0-xmax0)/xs)
      END IF
    ELSE
      !cdf = (1d0 - EXP(xmin0-x)) / (1d0 - EXP(xmin0-xmax0))
      cdf = EXPm1(xmin0-x) / EXPm1(xmin0-xmax0)
    END IF
  END IF
END FUNCTION ExponentialCDF

!----------------------------------------
! 1-CDF for exponential distribution
!   x           Value distributed exponentially
!   xmin        Minimum value (optional, default is 0)
!   xmax        Maximum value (optional, default is infinity)
!   xs          The exponential scale (optional, default is 1)
ELEMENTAL FUNCTION Exponential1mCDF(x,xmin,xmax,xs) RESULT(p)
  IMPLICIT NONE
  REAL*8 :: p
  REAL*8, INTENT(IN) :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax,xs
  REAL*8 :: xmin0,xmax0
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = HUGE(1d0)
  END IF
  IF (x .LE. xmin0) THEN
    p = 1d0
  ELSE IF (x .GE. xmax0) THEN
    p = 0d0
  ELSE
    IF (PRESENT(xs)) THEN
      IF (xs .EQ. 0d0) THEN
        p = (xmax0-x) / (xmax0-xmin0)
      ELSE
        !p = (EXP((xmax0-x)/xs) - 1d0) / (EXP((xmax0-xmin0)/xs) - 1d0)
        p = EXPm1((xmax0-x)/xs) / EXPm1((xmax0-xmin0)/xs)
      END IF
    ELSE
      !p = (EXP(xmax0-x) - 1d0) / (EXP(xmax0-xmin0) - 1d0)
      p = EXPm1(xmax0-x) / EXPm1(xmax0-xmin0)
    END IF
  END IF
END FUNCTION Exponential1mCDF

!----------------------------------------
! Random number for exponential distribution
!   xmin        Minimum value (optional, default is 0)
!   xmax        Maximum value (optional, default is infinity)
!   xs          The exponential scale (optional, default is 1)
FUNCTION RandExponential(xmin,xmax,xs) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8, INTENT(IN), OPTIONAL :: xmin,xmax,xs
  REAL*8 :: xmin0,xmax0
  REAL*8 :: urand
  urand = rand()
  IF (PRESENT(xmin)) THEN
    xmin0 = xmin
  ELSE
    xmin0 = 0d0
  END IF
  IF (PRESENT(xmax)) THEN
    xmax0 = xmax
  ELSE
    xmax0 = HUGE(1d0)
  END IF
  IF (PRESENT(xs)) THEN
    IF (xs .EQ. 0d0) THEN
      x = xmin0 + urand*(xmax0-xmin0)
    ELSE
      x = xmin0 - xs*LOG(1d0 - urand*(1d0-EXP((xmin0-xmax0)/xs)))
    END IF
  ELSE
    x = xmin0 - LOG(1d0 - urand*(1d0-EXP(xmin0-xmax0)))
  END IF
END FUNCTION RandExponential


!----------CHI-SQUARE DISTRIBUTION--------------------------------------
! Chi-square distribution has probability:
!   P(x|k) = x^(k/2 - 1) EXP(-x/2) / 2^(k/2) / GAMMA(k/2)
!-----------------------------------------------------------------------

!----------------------------------------
! Probability for chi-square distribution
!   x           Value chi-square distributed
!   k           Degrees of freedom
ELEMENTAL FUNCTION ChiSquarePDF(x,k) RESULT(pdf)
  IMPLICIT NONE
  REAL*8 :: pdf
  REAL*8, INTENT(IN) :: x
  INTEGER, INTENT(IN) :: k
  IF (x .GE. 0d0) THEN
    pdf = x**(0.5d0*k-1) * EXP(-0.5d0*x) / (2**(0.5d0*k)*GAMMA(0.5d0*k))
  ELSE
    pdf = 0d0
  END IF
END FUNCTION ChiSquarePDF

!----------------------------------------
! Log of probability for chi-square distribution
!   x           Value chi-square distributed
!   k           Degrees of freedom
ELEMENTAL FUNCTION ChiSquareLogPDF(x,k) RESULT(lnpdf)
  IMPLICIT NONE
  REAL*8 :: lnpdf
  REAL*8, INTENT(IN) :: x
  INTEGER, INTENT(IN) :: k
  REAL*8, PARAMETER :: LOG2 = 0.69314718055994531d0
  IF (x .GE. 0d0) THEN
    lnpdf = (0.5d0*k-1)*LOG(x) - 0.5d0*x                                &
            - (0.5d0*k)*LOG2 - LOG_GAMMA(0.5d0*k)
  ELSE
    lnpdf = 0d0
  END IF
END FUNCTION ChiSquareLogPDF

!----------------------------------------
! CDF for chi-square distribution
!   x           Value chi-square distributed
!   k           Degrees of freedom
ELEMENTAL FUNCTION ChiSquareCDF(x,k) RESULT(cdf)
  IMPLICIT NONE
  REAL*8 :: cdf
  REAL*8, INTENT(IN) :: x
  INTEGER, INTENT(IN) :: k
  IF (x .GE. 0d0) THEN
    ! Regularized gamma function P(k/2,x/2)
    cdf = GAMMA_P(0.5d0*k,0.5d0*x)
  ELSE
    cdf = 0d0
  END IF
END FUNCTION ChiSquareCDF

!----------------------------------------
! 1-CDF for chi-square distribution
!   x           Value chi-square distributed
!   k           Degrees of freedom
ELEMENTAL FUNCTION ChiSquare1mCDF(x,k) RESULT(p)
  IMPLICIT NONE
  REAL*8 :: p
  REAL*8, INTENT(IN) :: x
  INTEGER, INTENT(IN) :: k
  IF (x .GE. 0d0) THEN
    ! Regularized gamma function Q(k/2,x/2)
    p = GAMMA_Q(0.5d0*k,0.5d0*x)
  ELSE
    p = 1d0
  END IF
END FUNCTION ChiSquare1mCDF

!----------------------------------------
! Random number for normal distribution
!   k           Degrees of freedom
FUNCTION RandChiSquare(k) RESULT(x)
  IMPLICIT NONE
  REAL*8 :: x
  INTEGER, INTENT(IN) :: k
  REAL*8 :: u
  u = rand()
  !x = REDUCEDGAMMAINV(0.5d0*k,u)
  x = -1d0
  WRITE(*,*) 'ERROR: RandChiSquare not implemented (REDUCEDGAMMAINV)'
  STOP
END FUNCTION RandChiSquare


!----------POISSON DISTRIBUTION-----------------------------------------
! Poisson distribution has probability:
!   P(N|mu) = EXP(-mu) mu^N / N!
!-----------------------------------------------------------------------

!----------------------------------------
! Probability for Poisson distribution
!   N           Number Poisson distributed
!   mu          Average
ELEMENTAL FUNCTION PoissonPDF(N,mu) RESULT(pdf)
  IMPLICIT NONE
  REAL*8 :: pdf
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: mu
  ! Bad cases
  IF (N .LT. 0) THEN
    pdf = 0d0
    RETURN
  ELSE IF (mu .LT. 0d0) THEN
    pdf = -HUGE(1d0)
    RETURN
  END IF
  ! Catch special case mu = 0
  IF (mu .EQ. 0d0) THEN
    IF (N .EQ. 0) THEN
      pdf = 1d0
    ELSE
      pdf = 0d0
    END IF
    RETURN
  END IF
  ! General case
  pdf = EXP(-mu + N*LOG(mu) - LOG_GAMMA(N+1d0))
END FUNCTION PoissonPDF

!----------------------------------------
! Log of probability for Poisson distribution
!   N           Number Poisson distributed
!   mu          Average
ELEMENTAL FUNCTION PoissonLogPDF(N,mu) RESULT(lnpdf)
  IMPLICIT NONE
  REAL*8 :: lnpdf
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: mu
  ! Bad cases
  IF ((N .LT. 0) .OR. (mu .LT. 0d0)) THEN
    lnpdf = -HUGE(1d0)
    RETURN
  END IF
  ! Catch special case mu = 0
  IF (mu .EQ. 0d0) THEN
    IF (N .EQ. 0) THEN
      lnpdf = 0d0
    ELSE
      lnpdf = -HUGE(1d0)
    END IF
    RETURN
  END IF
  ! General case
  lnpdf = -mu + N*LOG(mu) - LOG_GAMMA(N+1d0)
END FUNCTION PoissonLogPDF

!----------------------------------------
! CDF for Poisson distribution
!   N           Number Poisson distributed
!   mu          Average
ELEMENTAL FUNCTION PoissonCDF(N,mu) RESULT(cdf)
  IMPLICIT NONE
  REAL*8 :: cdf
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: mu
  ! Bad cases
  IF (N .LT. 0) THEN
    cdf = 0d0
    RETURN
  ELSE IF (mu .LT. 0d0) THEN
    cdf = -HUGE(1d0)
    RETURN
  END IF
  ! Catch special case mu = 0
  IF (mu .EQ. 0d0) THEN
    cdf = 1d0
    RETURN
  END IF
  ! General case
  ! Regularized gamma function Q(N+1,mu)
  cdf = GAMMA_Q(N+1d0,mu)
END FUNCTION PoissonCDF

!----------------------------------------
! 1-CDF for Poisson distribution
!   N           Number Poisson distributed
!   mu          Average
ELEMENTAL FUNCTION Poisson1mCDF(N,mu) RESULT(p)
  IMPLICIT NONE
  REAL*8 :: p
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: mu
  ! Bad cases
  IF (N .LT. 0) THEN
    p = 1d0
    RETURN
  ELSE IF (mu .LT. 0d0) THEN
    p = -HUGE(1d0)
    RETURN
  END IF
  ! Catch special case mu = 0
  IF (mu .EQ. 0d0) THEN
    p = 0d0
    RETURN
  END IF
  ! General case
  ! Regularized gamma function P(N+1,mu)
  p = GAMMA_P(N+1d0,mu)
END FUNCTION Poisson1mCDF

!----------------------------------------
! Random number for Poisson distribution
!   mu          Average
FUNCTION RandPoisson(mu) RESULT(N)
  IMPLICIT NONE
  INTEGER :: N
  REAL*8, INTENT(IN) :: mu
  REAL*8 :: alpha,beta,k,c,u,u2,x,xsum,z
  REAL*8, PARAMETER :: PI = 3.1415926535897932d0
  ! Small mu case: simply iterate over terms.
  IF (mu .LE. 40d0) THEN
    u = rand()
    N = 0
    x = EXP(-mu)
    xsum = x
    DO WHILE (xsum .LT. u)
      N = N+1
      x = mu * x / N
      xsum = xsum + x
      ! Bad case: start over (can occur for u extremely
      ! close to 1).  Note 1-CDF(102|mu=40) = 7e-17
      ! => remaining terms too small.
      IF (N .GT. 102) THEN
        u = rand()
        N = 0
        x = EXP(-mu)
        xsum = x
      END IF
    END DO
    RETURN
  END IF
  ! Larger mu case: use rejection technique.
  ! See Atkinson, Appl. Statist. 28, 29 (1979).
  beta  = PI / SQRT(3*mu)
  alpha = beta*mu
  c = 0.767d0 - 3.36d0/mu
  k = LOG(c/beta) - mu
  DO WHILE (.TRUE.)
    u = rand()
    x = (alpha - LOG((1-u)/u)) / beta
    N = NINT(x)
    IF (N .LT. 0) CYCLE
    u2 = rand()
    z = alpha - beta*x
    IF (z+LOG(u2/(1+EXP(z))**2) .LE. k+N*LOG(mu)-LOG_GAMMA(N+1d0)) THEN
      RETURN
    END IF
  END DO
END FUNCTION RandPoisson


!----------BINOMIAL DISTRIBUTION----------------------------------------
! Binomial distribution has probability:
!   P(k|p,N) = [N!/(k! (N-k)!)] p^k (1-p)^(N-k)
!-----------------------------------------------------------------------

!----------------------------------------
! Probability for binomial distribution
!   k           Number of successful trials
!   p           Probability of success for each trial
!   N           Number of trials
ELEMENTAL FUNCTION BinomialPDF(k,p,N) RESULT(pdf)
  IMPLICIT NONE
  REAL*8 :: pdf
  INTEGER, INTENT(IN) :: k,N
  REAL*8, INTENT(IN) :: p
  INTEGER :: k0
  REAL*8 :: lnCnk
  ! Bad cases
  IF ((k .LT. 0) .OR. (k .GT. N)) THEN
    !pdf = -HUGE(pdf)
    pdf = 0d0
    RETURN
  END IF
  ! Special cases
  IF (p .LE. 0d0) THEN
    IF (k .EQ. 0) THEN
      pdf = 1d0
    ELSE
      pdf = 0d0
    END IF
    RETURN
  ELSE IF (p .GE. 1d0) THEN
    IF (k .EQ. N) THEN
      pdf = 1d0
    ELSE
      pdf = 0d0
    END IF
    RETURN
  END IF
  ! General case
  k0 = MIN(k,N-k)
  IF (N .LE. 30) THEN
    ! quick, exact for given condition
    pdf = BINOMIAL_COEFF(N,k) * p**k * (1-p)**(N-k)
  !ELSE IF (N .LE. 1000) THEN
  !  lnCnk = LOG_GAMMA(N+1d0) - LOG_GAMMA(k+1d0) - LOG_GAMMA((N-k)+1d0)
  !  pdf   = EXP(lnCnk) * p**k * (1-p)**(N-k)
  ! Overflow/underflow issues for N ~ O(10^3)
  ELSE
    lnCnk = LOG_GAMMA(N+1d0) - LOG_GAMMA(k+1d0) - LOG_GAMMA((N-k)+1d0)
    pdf   = EXP(lnCnk + k*LOG(p) + (N-k)*LOG(1-p))
  END IF
END FUNCTION BinomialPDF

!----------------------------------------
! Probability for binomial distribution
!   k           Number of successful trials
!   p           Probability of success for each trial
!   N           Number of trials
ELEMENTAL FUNCTION BinomialLogPDF(k,p,N) RESULT(lnpdf)
  IMPLICIT NONE
  REAL*8 :: lnpdf
  INTEGER, INTENT(IN) :: k,N
  REAL*8, INTENT(IN) :: p
  INTEGER :: k0
  REAL*8 :: lnCnk
  ! Bad cases
  IF ((k .LT. 0) .OR. (k .GT. N)) THEN
    lnpdf = -HUGE(1d0)
    RETURN
  END IF
  ! Special cases
  IF (p .LE. 0d0) THEN
    IF (k .EQ. 0) THEN
      lnpdf = 0d0
    ELSE
      lnpdf = -HUGE(1d0)
    END IF
    RETURN
  ELSE IF (p .GE. 1d0) THEN
    IF (k .EQ. N) THEN
      lnpdf = 0d0
    ELSE
      lnpdf = -HUGE(1d0)
    END IF
    RETURN
  END IF
  ! General case
  k0 = MIN(k,N-k)
  IF (N .LE. 30) THEN
    ! quick, exact for given condition
    lnCnk = LOG_BINOMIAL_COEFF(N,k)
  ELSE
    lnCnk = LOG_GAMMA(N+1d0) - LOG_GAMMA(k+1d0) - LOG_GAMMA((N-k)+1d0)
  END IF
  lnpdf = lnCnk + k*LOG(p) +(N-k)*LOG(1-p)
END FUNCTION BinomialLogPDF

!----------------------------------------
! CDF for binomial distribution
!   k           Number of successful trials
!   p           Probability of success for each trial
!   N           Number of trials
! NOTE: Not optimized for large N (scales as sqrt(N*p*q)).
ELEMENTAL FUNCTION BinomialCDF(k,p,N) RESULT(cdf)
  IMPLICIT NONE
  REAL*8 :: cdf
  INTEGER, INTENT(IN) :: k,N
  REAL*8, INTENT(IN) :: p
  REAL*8 :: ltsum,gtsum,pdf
  ! Note cdf = BetaRegularized[1-p,N-k,k+1]
  CALL BinomialSums(k,p,N,ltsum,gtsum,pdf)
  cdf = ltsum + pdf
END FUNCTION BinomialCDF

!----------------------------------------
! 1-CDF for binomial distribution
!   k           Number of successful trials
!   p           Probability of success for each trial
!   N           Number of trials
! NOTE: Not optimized for large N (scales as sqrt(N*p*q)).
ELEMENTAL FUNCTION Binomial1mCDF(k,p,N) RESULT(sf)
  IMPLICIT NONE
  REAL*8 :: sf
  INTEGER, INTENT(IN) :: k,N
  REAL*8, INTENT(IN) :: p
  REAL*8 :: ltsum,gtsum,pdf
  ! Note cdf = BetaRegularized[1-p,N-k,k+1]
  CALL BinomialSums(k,p,N,ltsum,gtsum,pdf)
  sf = gtsum
END FUNCTION Binomial1mCDF

!----------------------------------------
! Calculates various sums for binomial distribution.
! Utility routine for CDF and 1-CDF calculations.
!   k           Number of successful trials
!   p           Probability of success for each trial
!   N           Number of trials
! Output:
!   ltsum       Sum of PDFs for j=0,..,k-1
!   gtsum       Sum of PDFs for j=k+1,..,N
!   pdf         The PDF for k
! NOTE: Not optimized for large N (scales as sqrt(N*p*q)).
ELEMENTAL SUBROUTINE BinomialSums(k,p,N,ltsum,gtsum,pdf)
  IMPLICIT NONE
  REAL*8 :: cdf
  INTEGER, INTENT(IN) :: k,N
  REAL*8, INTENT(IN) :: p
  REAL*8, INTENT(OUT) :: ltsum,gtsum,pdf
  LOGICAL :: flipped
  INTEGER :: k0,I
  REAL*8 :: p0,q0,mu,xI,tmp
  
  ! Bad cases
  IF (k .LT. 0) THEN
    ltsum = 0d0
    gtsum = 1d0
    pdf   = 0d0
    RETURN
  ELSE IF (k .GT. N) THEN
    ltsum = 1d0
    gtsum = 0d0
    pdf   = 0d0
    RETURN
  END IF
  
  ! Can write CDF(k|p,N) = CDF(N-k|1-p,N) = CDF(k0,p0,N)
  ! Choose smaller of p or 1-p to work with
  IF (p .GT. 0.5d0) THEN
    p0 = 1-p
    q0 = p
    k0 = N-k
    flipped = .TRUE.
  ELSE
    p0 = p
    q0 = 1-p
    k0 = k
    flipped = .FALSE.
  END IF
  mu = N*p0
  
  ! Will find approximately smaller sum first.
  ! Special cases.
  IF (mu .EQ. 0d0) THEN
    IF (k0 .EQ. 0) THEN
      ltsum = 0d0
      gtsum = 0d0
      pdf   = 1d0
    ELSE
      ltsum = 1d0
      gtsum = 0d0
      pdf   = 0d0
    END IF
  ELSE IF (k0 .EQ. 0) THEN
    ltsum = 0d0
    pdf   = q0**N
    gtsum = 1d0 - pdf
  ELSE IF (k0 .EQ. N) THEN
    gtsum = 0d0
    pdf   = p0**N
    ltsum = 1d0 - pdf
  ! Work our way down from k0.
  ELSE IF (k0 .LE. mu) THEN
    I  = k0
    xI = EXP(LOG_BINOMIAL_COEFF(N,I) + I*LOG(p0) + (N-I)*LOG(q0))
    pdf   = xI
    ltsum = 0d0
    DO WHILE ((I .GT. 0) .AND. (xI .GT. EPSILON(1d0)*ltsum))
      I     = I-1
      xI    = xI * ((I+1)*q0) / ((N-I)*p0)
      ltsum = ltsum + xI
    END DO
    gtsum = 1d0 - pdf - ltsum
  ! Work our way up from k0.
  ELSE
    I  = k0
    xI = EXP(LOG_BINOMIAL_COEFF(N,I) + I*LOG(p0) + (N-I)*LOG(q0))
    pdf   = xI
    gtsum = 0d0
    DO WHILE ((I .LT. N) .AND. (xI .GT. EPSILON(1d0)*gtsum))
      I     = I+1
      xI    = xI * ((N-I+1)*p0) / (I*q0)
      gtsum = gtsum + xI
    END DO
    ltsum = 1d0 - pdf - gtsum
  END IF
  ! reverse sums if we reversed p & q, k & N-k
  IF (flipped) THEN
    tmp   = ltsum
    ltsum = gtsum
    gtsum = tmp
  END IF
END SUBROUTINE BinomialSums

!----------------------------------------
! Random number for binomial distribution
!   p           Probability of success for each trial
!   N           Number of trials
FUNCTION RandBinomial(p,N) RESULT(K)
  IMPLICIT NONE
  INTEGER :: K
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: p
  REAL*8, PARAMETER :: MUMAX = 20d0
  LOGICAL :: flipped
  REAL*8 :: p0,mu
  ! Bad/special cases
  IF ((N .LT. 0) .OR. (p .LT. 0d0) .OR. (p .GT. 1d0)) THEN
    K = -HUGE(K)
    RETURN
  ELSE IF (N .EQ. 0) THEN
    K = 0
    RETURN
  END IF
  ! Can write P(k|p,N) = P(N-k|1-p,N) = P(k0,p0,N)
  ! Choose smaller of p or 1-p to work with
  IF (p .GT. 0.5d0) THEN
    p0 = 1-p
    flipped = .TRUE.
  ELSE
    p0 = p
    flipped = .FALSE.
  END IF
  mu = N*p0
  ! Special case
  IF (mu .EQ. 0d0) THEN
    K = 0
  ! Small mu case: iterative calculation faster
  ELSE IF (mu .LE. MUMAX) THEN
    K = IteratedK(p0,N)
  ! Otherwise, use ~ constant time rejection algorithm
  ELSE
    K = RejectionK(p0,N)
  END IF
  IF (flipped) K = N - K
  
  CONTAINS
  !-----------------------------
  ! Use explicit iteration to find random K.
  ! Assumes 0 < p <= 0.5 and N > 0.
  FUNCTION IteratedK(p,N) RESULT(K)
    IMPLICIT NONE
    INTEGER :: K
    INTEGER, INTENT(IN) :: N
    REAL*8, INTENT(IN) :: p
    REAL*8 :: q,mu,xK,y,sum
    q  = 1-p
    mu = N*p
    y  = rand()
    ! Work our way up from K=0
    IF (y .LE. 0.9d0) THEN
      K = NINT(mu - 9*SQRT(mu*q))
      IF (K .LE. 0) THEN
        K  = 0
        xK = q**N
      ELSE
        xK = EXP(LOG_BINOMIAL_COEFF(N,K) + K*LOG(p) + (N-K)*LOG(q))
      END IF
      sum = xK
      DO WHILE (sum .LT. y)
        K   = K+1
        xK  = xK * ((N-K+1)*p) / (K*q)
        sum = sum + xK
      END DO
    ! Work backwards from tail
    ELSE
      y = 1-y
      K = NINT(mu + 9*SQRT(mu*q))
      IF (K .GE. N) THEN
        K  = N
        xK = p**N
      ELSE
        xK = EXP(LOG_BINOMIAL_COEFF(N,K) + K*LOG(p) + (N-K)*LOG(q))
      END IF
      sum = xk
      DO WHILE (sum .LT. y)
        K   = K-1
        xK  = xK * ((K+1)*q) / ((N-K)*p)
        sum = sum + xK
      END DO
    END IF
  END FUNCTION IteratedK
  
  !-----------------------------
  ! Use rejection technique to find random K, as
  ! described in Hormann, JSCS 46, 101 (1993).
  ! Assumes 0 < p <= 0.5 and N > 0.
  FUNCTION RejectionK(p,N) RESULT(K)
    IMPLICIT NONE
    INTEGER :: K
    INTEGER, INTENT(IN) :: N
    REAL*8, INTENT(IN) :: p
    INTEGER :: M
    REAL*8 :: q,mu,sd,r,a,b,c,alpha,ur,vr,urvr,u,v,lnf,lnf_accept
    q  = 1-p
    mu = N*p
    ! Rejection technique: binomial density is proportional to
    !   f(x) = Pb(floor(x)) / Pb(M) <= 1
    ! where Pb is the binomial probability and  M = floor((N+1)*p0)
    ! is the mode. Use dominating distribution with inverse
    !   G(u) = (2alpha/(0.5-|u|) + b)*u + c
    ! with u in [-0.5,0.5].  Appropriate choices of constants
    ! determined by Hormann, JSCS 46, 101 (1993).
    sd = SQRT(mu*q)
    M = INT((N+1)*p)
    r = p/q
    c = mu + 0.5d0
    b = 1.15d0 + 2.53d0*sd
    a = -0.0873d0 + 0.0248d0*b + 0.01d0*p
    alpha = (2.83d0 + 5.1d0/b) * sd
    !ur = 0.43d0
    vr = 0.92d0 - 4.2d0/b
    urvr = 0.86d0*vr
    ! Will generate point in u in [-0.5,0.5], v in [0,1].
    ! Use a little magic to reduce number of uniform random
    ! numbers needed (on average).
    DO WHILE (.TRUE.)
      v = rand()
      ! Box: u in [-ur,ur], v in [0,vr].
      ! Always accepted.  Uniform in [0,ur*vr] converted to
      ! uniform in [-ur,ur], avoiding second random number.
      IF (v .LE. urvr) THEN
        u = v/vr - 0.43d0
        K = FLOOR((2*a/(0.5d0-ABS(u))+b)*u+c)
        EXIT
      END IF
      ! Outside of above box.  Two regions to consider:
      ! 1) Point in v > vr.  Need random u (can keep v).
      IF (v .GE. vr) THEN
        u = rand() - 0.5d0
      ! 2) Point in v < vr, u in [-0.5,-ur] or [ur,0.5].
      ! Uniform in [2ur*vr,vr] converted to uniform in
      ! [-0.5,-ur]+[ur,0.5].  Need random v in [0,vr].
      ELSE
        u = v/vr - 0.93d0
        u = SIGN(0.5d0,u) - u
        v = vr*rand()
      END IF
      ! Transform U to K = floor(G(U)).
      K = FLOOR((2*a/(0.5d0-ABS(u))+b)*u+c)
      IF ((K .LT. 0) .OR. (K .GT. N)) CYCLE
      ! Acceptance condition: V < f(G(U)) G'(U) / alpha
      ! where f(G(U)) = f(K).  Rewrite as:
      !   alpha*V/G'(U) < f(G(U)) = f(K)
      lnf_accept = LOG(alpha*v / (a/(0.5d0-ABS(u))**2+b))
      ! Calculate f(K).
      ! This step not terribly well optimized: could use recursion
      ! for K ~ M, though LOG_GAMMA is fairly fast as an intrinsic.
      IF (N .LE. 30) THEN
        lnf = LOG_BINOMIAL_COEFF(N,K) - LOG_BINOMIAL_COEFF(N,M)         &
              + (K-M)*LOG(r)
      ELSE
        lnf = LOG_GAMMA(M+1d0) + LOG_GAMMA((N-M)+1d0)                   &
              - LOG_GAMMA(K+1d0) - LOG_GAMMA((N-K)+1d0)                 &
              + (K-M)*LOG(r)
      END IF
      IF (lnf .GE. lnf_accept) EXIT
    END DO
  END FUNCTION RejectionK
  
END FUNCTION RandBinomial



!=======================================================================
! ARRAY SORT/SEARCH
!=======================================================================

!-----------------------------------------------------------------------
! Sorts the given array using the heap sort method.
! Faster than quicksort when arrays are close to sorted already.
!   N           length of array
!   x           array of data to be sorted
! 
PURE SUBROUTINE HSort(N,x)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(INOUT) :: x(N)
  INTEGER :: I,Itmp
  REAL*8 :: xtmp
  DO I = N/2, 1, -1
    CALL HSort_siftdown(N,x,I,N)
  END DO
  DO I = N, 2, -1
    xtmp = x(1)
    x(1) = x(I)
    x(I) = xtmp
    CALL HSort_siftdown(N,x,1,I-1)
  END DO
  RETURN
    
  CONTAINS
  ! ------------------------------------------
  ! Utility function for heap sort
  PURE SUBROUTINE HSort_siftdown(N,x,Jlow,Jhigh)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N,Jlow,Jhigh
    REAL*8, INTENT(INOUT) :: x(N)
    INTEGER :: J,Jold
    REAL*8 :: x0
    x0 = x(Jlow)
    Jold = Jlow
    J = 2*Jlow
    DO WHILE (J .LE. Jhigh)
      IF ((J .LT. Jhigh) .AND. (x(J)) .LT. x(J+1)) THEN
        J = J+1
      END IF
      IF (x0 .GE. x(J)) EXIT
      x(Jold) = x(J)
      Jold = J
      J = 2*J
    END DO
    x(Jold) = x0
    RETURN
  END SUBROUTINE HSort_siftdown
END SUBROUTINE HSort


! ----------------------------------------------------------------------
! Searches the given array of size [1:N], assumed to be sorted in
! increasing order, for the given value using a binary search algorithm.
! Returns the index I such that x(I) <= x0 < x(I+1), 0 if x0 < x(1),
! or N if x(N) <= x0.
! 
! Input arguments:
!   N               Size of x array
!   x               Array of data to be searched (must be sorted)
!   x0              Value to search for
! Optional input arguments:
!   Istart          Index to start searching from
! 
PURE FUNCTION BSearch(N,x,x0,Istart) RESULT(index)
  IMPLICIT NONE
  INTEGER :: index
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: x(N),x0
  INTEGER, INTENT(IN), OPTIONAL :: Istart
  INTEGER :: Ilow,Ihigh,Imid,step,dir
  
  ! Check if in bounds
  IF (x0 .LT. x(1)) THEN
    index = 0
    RETURN
  ELSE IF (x0 .GE. x(N)) THEN
    index = N
    RETURN
  END IF
  
  ! If starting index given, bracket desired index
  IF (PRESENT(Istart)) THEN
    ! Find one bound
    IF (Istart .LE. 1) THEN
      Ilow = 1
      dir = +1
    ELSE IF (Istart .LT. N) THEN
      IF (x0 .GE. x(Istart)) THEN
        Ilow = Istart
        dir = +1
      ELSE
        Ihigh = Istart
        dir = -1
      END IF
    ELSE
      Ihigh = N
      dir = -1
    END IF
    ! Search up or down in increasing step sizes to find other bound
    step = 1
    IF (dir .GT. 0) THEN
      Ihigh = MIN(Ilow + step,N)
      DO WHILE (x0 .GT. x(Ihigh))
        Ilow  = Ihigh
        Ihigh = MIN(Ilow + step,N)
        step  = 2*step
      END DO
    ELSE IF (dir .LT. 0) THEN
      Ilow = MAX(Ihigh - step,1)
      DO WHILE (x0 .LT. x(Ilow))
        Ihigh = Ilow
        Ilow  = MAX(Ilow - step,1)
        step  = 2*step
      END DO
    END IF
  ! No starting index given, search entire array
  ELSE
    Ilow  = 1
    Ihigh = N
  END IF
  
  ! Binary search
  DO WHILE (Ihigh-Ilow .GT. 1)
    Imid = (Ilow+Ihigh)/2
    IF (x0 .GE. x(Imid)) THEN
      Ilow = Imid
    ELSE
      Ihigh = Imid
    END IF
  END DO
  index = Ilow
  
END FUNCTION BSearch



!=======================================================================
! STRING UTILITY FUNCTIONS
!=======================================================================

! ----------------------------------------------------------------------
! Changes string to uppercase.
! 
PURE SUBROUTINE ToUppercase(string)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(INOUT) :: string
  INTEGER, PARAMETER :: UA = ICHAR('A')
  INTEGER, PARAMETER :: la = ICHAR('a')
  INTEGER, PARAMETER :: UZ = ICHAR('Z')
  INTEGER, PARAMETER :: lz = ICHAR('z')
  INTEGER :: I,ival
  
  DO I = 1, LEN_TRIM(string)
    ival = ICHAR(string(I:I))
    IF ((ival .GE. la) .AND. (ival .LE. lz)) THEN
      string(I:I) = CHAR(ival+UA-la)
    END IF
  END DO
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Changes string to lowercase.
! 
PURE SUBROUTINE ToLowercase(string)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(INOUT) :: string
  INTEGER, PARAMETER :: UA = ICHAR('A')
  INTEGER, PARAMETER :: la = ICHAR('a')
  INTEGER, PARAMETER :: UZ = ICHAR('Z')
  INTEGER, PARAMETER :: lz = ICHAR('z')
  INTEGER :: I,ival

  DO I = 1, LEN_TRIM(string)
    ival = ICHAR(string(I:I))
    IF ((ival .GE. UA) .AND. (ival .LE. UZ)) THEN
      string(I:I) = CHAR(ival-UA+la)
    END IF
  END DO

END SUBROUTINE


! ----------------------------------------------------------------------
! Obtain the file suffix or an empty string if there is no suffix.
! 
PURE SUBROUTINE FileSuffix(file,suffix)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN) :: file
  CHARACTER*(*), INTENT(OUT) :: suffix
  INTEGER :: len,pos
  
  len = LEN_TRIM(file)
  pos = INDEX(file,'.',back=.TRUE.)
  
  ! No suffix
  IF (pos .EQ. 0) THEN
    suffix = ''
    RETURN
  END IF
  
  ! Get suffix
  suffix = file(pos+1:)
  
END SUBROUTINE


! ----------------------------------------------------------------------
! Obtain the file suffix or an empty string if there is no suffix.
! The file suffix is converted to lower case for easier parsing.
! 
PURE SUBROUTINE FileSuffixLC(file,suffix)
  IMPLICIT NONE
  CHARACTER*(*), INTENT(IN) :: file
  CHARACTER*(*), INTENT(OUT) :: suffix
  INTEGER :: len,pos
  
  len = LEN_TRIM(file)
  pos = INDEX(file,'.',back=.TRUE.)
  
  ! No suffix
  IF (pos .EQ. 0) THEN
    suffix = ''
    RETURN
  END IF
  
  ! Get suffix, convert to lower case for easier parsing
  suffix = file(pos+1:)
  CALL ToLowercase(suffix)
  
END SUBROUTINE



!=======================================================================
! TIME UTILITY FUNCTIONS
!=======================================================================

!-----------------------------------------------------------------------
! System time in seconds.
! 
FUNCTION GetTime() RESULT(t)
  IMPLICIT NONE
  REAL*8 :: t
  ! Resolution may vary by kind; ensure highest resolution
  INTEGER(KIND=8) :: count,count_rate
  CALL SYSTEM_CLOCK(count,count_rate)
  t = count*1d0 / count_rate
END FUNCTION


!-----------------------------------------------------------------------
! Time elapsed in seconds since TimeElapsedReset() was called.
! 
FUNCTION TimeElapsed() RESULT(t)
  IMPLICIT NONE
  REAL*8 :: t
  CALL TimeElapsedRoutine(t,.FALSE.)
END FUNCTION


!-----------------------------------------------------------------------
! Reset TimeElapsed timer.  Subroutine version of above.
! 
SUBROUTINE TimeElapsedReset()
  IMPLICIT NONE
  REAL*8 :: t
  CALL TimeElapsedRoutine(t,.TRUE.)
END SUBROUTINE


!-----------------------------------------------------------------------
! Time elapsed in seconds since the parameter 'reset' was set to true.
! Used by TimeElapsed() and TimeElapsedReset() routines, which should
! be used instead of direct calls to this routine.  Note that this is
! wall time, not CPU time.
! 
!   t           Time elpased
!   reset       Logical specifying if the counter should be
!               reset to the current cpu time.
! 
SUBROUTINE TimeElapsedRoutine(t,reset)
  IMPLICIT NONE
  REAL*8, INTENT(OUT) :: t
  LOGICAL,INTENT(IN) :: reset
  ! Resolution may vary by kind; ensure highest resolution
  INTEGER(KIND=8), SAVE :: icount = 0
  INTEGER(KIND=8) :: fcount,count_rate
  IF (reset) THEN
    CALL SYSTEM_CLOCK(icount)
    t = 0d0
    RETURN
  END IF
  CALL SYSTEM_CLOCK(fcount,count_rate)
  t = (fcount-icount)*1d0 / count_rate
END SUBROUTINE


!-----------------------------------------------------------------------
! CPU time elapsed in seconds since CPUTimeElapsedReset() was called.
! For multithreaded cases, this is the sum of time spent by all threads.
! 
FUNCTION CPUTimeElapsed() RESULT(t)
  IMPLICIT NONE
  REAL*8 :: t
  CALL CPUTimeElapsedRoutine(t,.FALSE.)
END FUNCTION


!-----------------------------------------------------------------------
! Reset CPUTimeElapsed timer.
! 
SUBROUTINE CPUTimeElapsedReset()
  IMPLICIT NONE
  REAL*8 :: t
  CALL CPUTimeElapsedRoutine(t,.TRUE.)
END SUBROUTINE


!-----------------------------------------------------------------------
! CPU time elapsed in seconds since the parameter 'reset' was set to
! true.  Used by CPUTimeElapsed() and CPUTimeElapsedReset() routines,
! which should be used instead of direct calls to this routine.
! 
!   t           Time elpased
!   reset       Logical specifying if the counter should be
!               reset to the current cpu time.
! 
SUBROUTINE CPUTimeElapsedRoutine(t,reset)
  IMPLICIT NONE
  REAL*8, INTENT(OUT) :: t
  LOGICAL,INTENT(IN) :: reset
  ! Resolution may vary by kind; ensure highest resolution
  REAL(KIND=8), SAVE :: t1
  REAL(KIND=8) :: t2
  IF (reset) THEN
    CALL CPU_TIME(t1)
    t = 0d0
  END IF
  CALL CPU_TIME(t2)
  t = t2 - t1
END SUBROUTINE



!#######################################################################
END MODULE
!#######################################################################


!#######################################################################
! END OF FILE
!#######################################################################

 
