!#######################################################################
! INTERFACE TESTING PROGRAM
! 
! Compile, assuming DDCalc0.o is already compiled:
!   gfortran -o InterfaceCheck InterfaceCheck.f90 DDCalc0.o
!   ifort -o InterfaceCheck InterfaceCheck.f90 DDCalc0.o
! ...or everything done in a single command (no DDCalc0.o necessary):
!   gfortran -fno-range-check DDCalc0.f90 -o InterfaceCheck InterfaceCheck.f90
!   ifort DDCalc0.f90 -o InterfaceCheck InterfaceCheck.f90
! 
! Run:
!   ./InterfaceCheck [--mfa|--mG|--msigma]
! The optional flags indicate which set of coupling parameters are
! to be used for input.
! 
! 
!   Created by Chris Savage
!   Nordita              (2014 -     )
! 
!#######################################################################

PROGRAM InterfaceCheck
  USE DDCALC0
  IMPLICIT NONE
  CHARACTER*32 :: arg
  INTEGER :: I,Narg,type
  REAL*8 :: M,xpSI,xnSI,xpSD,xnSD
  REAL*8 :: GpSI,GnSI,GpSD,GnSD,fp,fn,ap,an,sigmapSI,sigmanSI,sigmapSD,sigmanSD
  REAL*8 :: lnp
  TYPE(DetectorRateStruct) :: D
  ! These constants will be used to specify the type of input parameters
  INTEGER, PARAMETER :: TYPE_MG     = 1
  INTEGER, PARAMETER :: TYPE_MFA    = 2
  INTEGER, PARAMETER :: TYPE_MSIGMA = 3
  
  ! Parse command line options for this example program.
  ! Notably, determining how WIMP parameters will be specified.
  type = TYPE_MSIGMA
  Narg = IARGC()
  DO I = 1,Narg
    CALL GETARG(I,arg)
    IF (arg .EQ. '--mG') THEN
      type = TYPE_MG
    ELSE IF (arg .EQ. '--mfa') THEN
      type = TYPE_MFA
    ELSE IF (arg .EQ. '--msigma') THEN
      type = TYPE_MSIGMA
    END IF
  END DO
  
  ! NOTE: The module must be loaded via a 'USE DDCALC0' statement
  ! in any routine that is to make use of DDCalc0 routines (see
  ! the declaration area at the top of this routine).
  
  ! Initialize the module.
  CALL DDCalc0_Init()
  
  ! Initialize any experiments to be used.
  ! The argument indicates if extra sub-interval calculations should
  ! be performed.  Those calculations are required for maximum gap
  ! analyses, but are unnecessary for calculating total rates and
  ! likelihoods.  If .FALSE. is given, a no-background-subtraction
  ! p-value can still be calculated, but a Poisson is used instead
  ! of the maximum gap.  We show some maximum gap results below, so
  ! we must use .TRUE. here.
  CALL XENON100_2012_Init(.TRUE.)
  CALL LUX_2013_Init(.TRUE.)
  CALL DARWIN_Ar_2014_Init(.TRUE.)
  CALL DARWIN_Xe_2014_Init(.TRUE.)
  
  ! Can optionally specify a minimum recoil energy to be included in
  ! the rate calculations [keV].  Note the efficiency curves already
  ! account for detector and analysis thresholds regardless of this
  ! setting, so setting this to 0 keV (the default behavior when
  ! initialization is performed) does not imply that very low energy
  ! recoils actually contribute to the signal.
  ! EXAMPLE: Uncomment to set a minimum recoil energy of 3 keV.
  !CALL XENON100_2012_SetEmin(3d0)
  !CALL LUX_2013_SetEmin(3d0)
  !CALL DARWIN_Ar_2014_SetEmin(3d0)
  !CALL DARWIN_Xe_2014_SetEmin(3d0)
  
  ! Advanced usage:
  ! Explicit structure that is more configurable than above cases.
  ! Start by initializing the structure to the above LUX case.
  CALL LUX_2013_InitTo(D,.TRUE.)
  ! The DDCalc0_SetDetector() routine uses optional arguments that can
  ! be specified via keyword as shown before.  All the values used in
  ! the examples here are those used for the LUX case, so these
  ! specific calls are not actually necessary.
  ! 
  ! Set element to use.  Must be one of 14,18,32, or 54.  Spin-
  ! dependent interactions only implemented for 54 (xenon).
  CALL DDCalc0_SetDetector(D,Zelem=54)
  ! Change parameters from standard DARWIN argon case.
  CALL DDCalc0_SetDetector(D,mass=118d0,time=2*85.3d0,Nevents=1,background=0.64d0)
  ! Load efficiency curves from file.  First column is recoil energy
  ! [keV], the next column with values in [0,1] is the total detection
  ! efficiency.  Optionally (intervals=.TRUE.), additional columns are
  ! taken to be the detection efficiency for each interval between
  ! events.  If these columns are not available or should be ignored,
  ! set intervals=.FALSE.  There should be Nevents+1 intervals (not
  ! including total) if used.
  CALL DDCalc0_SetDetector(D,eff_file='example_efficiencies.dat',intervals=.TRUE.)
  ! Set the minimum recoil energy [keV] to include.
  CALL DDCalc0_SetDetector(D,Emin=0d0)
  
  ! TESTING:
  ! For comparison, set the D structure to represent the standard LUX
  ! case, but with a 3 keV minimum recoil energy.
  !CALL LUX_2013_InitTo(D,.TRUE.)
  !CALL DDCalc0_SetDetector(D,Emin=3d0)
  
  ! Write out directions for specifying input to this example program.
  ! WriteDescription is defined below.
  CALL WriteDescription(type)
  
  ! Loop over input to this example program.
  ! GetWIMPParams is defined below.
  DO WHILE (GetWIMPParams(type,M,xpSI,xnSI,xpSD,xnSD))
    
    WRITE(*,*)
    
    ! Set the WIMP parameters.
    ! There are three ways to specify the WIMP-nucleon couplings, with
    ! the WIMP mass [GeV] always the first argument:
    !   * DDCalc0_SetWIMP_mfa(m,fp,fn,ap,an)
    !     The standard couplings fp,fn [GeV^-2] & ap,an [unitless]
    !   * DDCalc0_SetWIMP_mG(m,GpSI,GnSI,GpSD,GnSD)
    !     The effective 4 fermion vertex couplings GpSI,GnSI,GpSD,GnSD
    !     [GeV^-2], related by:
    !         GpSI = 2 fp        GpSD = 2\sqrt{2} G_F ap
    !         GnSI = 2 fn        GnSD = 2\sqrt{2} G_F an
    !   * DDCalc0_SetWIMP_msigma(m,sigmapSI,sigmanSI,sigmapSD,sigmanSD)
    !     The WIMP-nucleon cross-sections [pb] (use a negative value
    !     to indicate the corresponding coupling should be negative).
    ! In the above, 'p' is for proton, 'n' is for neutron, 'SI' is for
    ! spin-independent, and 'SD' is for spin-dependent.
    SELECT CASE(type)
    CASE(TYPE_MG)
      CALL DDCalc0_SetWIMP_mG(M,xpSI,xnSI,xpSD,xnSD)
    CASE(TYPE_MFA)
      CALL DDCalc0_SetWIMP_mfa(M,xpSI,xnSI,xpSD,xnSD)
    CASE(TYPE_MSIGMA)
      CALL DDCalc0_SetWIMP_msigma(M,xpSI,xnSI,xpSD,xnSD)
    END SELECT
    
    ! Get the current WIMP parameters with the same signature and units
    ! as above.  The only difference is that WIMP-nucleon cross-sections
    ! are always positive (physical) values.
    CALL DDCalc0_GetWIMP_mG(M,GpSI,GnSI,GpSD,GnSD)
    CALL DDCalc0_GetWIMP_mfa(M,fp,fn,ap,an)
    CALL DDCalc0_GetWIMP_msigma(M,sigmapSI,sigmanSI,sigmapSD,sigmanSD)
    WRITE(*,'(A,1(1X,1PG12.4))') 'WIMP mass [GeV]     ',M
    WRITE(*,*)
    WRITE(*,'(A28,4(1X,A11))')     'WIMP-nucleon couplings          ',  &
        ' proton-SI ',' neutron-SI',' proton-SD ',' neutron-SD'
    WRITE(*,'(A28,4(1X,1PG11.4))') '  G [GeV^-2]                    ',  &
        GpSI,GnSI,GpSD,GnSD
    WRITE(*,'(A28,4(1X,1PG11.4))') '  f & a [GeV^-2,unitless]       ',  &
        fp,fn,ap,an
    WRITE(*,'(A28,4(1X,1PG11.4))') '  cross-section [pb]            ',  &
        sigmapSI,sigmanSI,sigmapSD,sigmanSD
    WRITE(*,*)
    
    ! Do rate calculations.
    ! After any change to the WIMP or halo parameters, perform the rate
    ! calculations necessary for predicted signals, likelihoods, and/or
    ! maximum gap statistics.
    CALL XENON100_2012_CalcRates()
    CALL LUX_2013_CalcRates()
    CALL DARWIN_Ar_2014_CalcRates()
    CALL DARWIN_Xe_2014_CalcRates()
    CALL DDCalc0_CalcRates(D)
    
    ! Header
    WRITE(*,'(A20,5(2X,A11))') '','XENON 2012 ',' LUX 2013  ',          &
        'DARWIN Ar  ','DARWIN Xe  ','(special)  '
    !WRITE(*,'(A20,5(1X,A12))') '','-----------','-----------',          &
    !    '-----------','-----------','-----------'
    
    ! Event quantities.
    ! The observed number of events (INTEGER).
    WRITE(*,'(A20,5(2X,1X,I6,4X))') 'Observed events                 ', &
        XENON100_2012_Events(),LUX_2013_Events(),                       &
        DARWIN_Ar_2014_Events(),DARWIN_Xe_2014_Events(),                &
        DDCalc0_Events(D)
    ! The average expected background.
    WRITE(*,'(A20,5(2X,1PG11.4))')  'Expected background             ', &
        XENON100_2012_Background(),LUX_2013_Background(),               &
        DARWIN_Ar_2014_Background(),DARWIN_Xe_2014_Background(),        &
        DDCalc0_Background(D)
    ! The average expected WIMP signal.
    WRITE(*,'(A20,5(2X,1PG11.4))')  'Expected signal                 ', &
        XENON100_2012_Signal(),LUX_2013_Signal(),                       &
        DARWIN_Ar_2014_Signal(),DARWIN_Xe_2014_Signal(),                &
        DDCalc0_Signal(D)
    ! The average expected WIMP spin-independent signal.
    WRITE(*,'(A20,5(2X,1PG11.4))')  '  spin-independent              ', &
        XENON100_2012_SignalSI(),LUX_2013_SignalSI(),                   &
        DARWIN_Ar_2014_SignalSI(),DARWIN_Xe_2014_SignalSI(),            &
        DDCalc0_SignalSI(D)
    ! The average expected WIMP spin-dependent signal.
    WRITE(*,'(A20,5(2X,1PG11.4))')  '  spin-dependent                ', &
        XENON100_2012_SignalSD(),LUX_2013_SignalSD(),                   &
        DARWIN_Ar_2014_SignalSD(),DARWIN_Xe_2014_SignalSD(),            &
        DDCalc0_SignalSD(D)
    
    ! The log-likelihoods for the current WIMP; note these are _not_
    ! multiplied by -2.  The likelihood is calculated using a Poisson
    ! given the observed number of events and expected signal+background.
    WRITE(*,'(A20,5(2X,1PG11.4))')  'Log-likelihood                  ', &
        XENON100_2012_LogLikelihood(),LUX_2013_LogLikelihood(),         &
        DARWIN_Ar_2014_LogLikelihood(),DARWIN_Xe_2014_LogLikelihood(),  &
        DDCalc0_LogLikelihood(D)
    
    ! The logarithm of the p-value, calculated without background
    ! subtraction, using either the maximum gap statistic or a Poisson
    ! statistic, depending on how the detector was initialized.  Note
    ! that this is actually a conservative upper _bound_ on the p-value
    ! in the event of an unknown background and is useful for excluding
    ! WIMP parameters.  However, since it is not a true p-value, it
    ! should not be interpreted as being related to any particular
    ! likelihood.
    WRITE(*,'(A20,5(2X,1PG11.4))')  'Max gap log(p-value)            ', &
        XENON100_2012_LogPValue(),LUX_2013_LogPValue(),                 &
        DARWIN_Ar_2014_LogPValue(),DARWIN_Xe_2014_LogPValue(),          &
        DDCalc0_LogPValue(D)
    
    ! The factor x by which the current WIMP cross- sections must be
    ! multiplied (sigma -> x*sigma, applied to all four WIMP-nucleon
    ! cross-sections) to achieve the given p-value (specified by its
    ! logarithm).  Useful for finding the no-background-subtraction
    ! exclusion limits.  For example, if setWIMP_msigma(100d0,10d0,10d0,
    ! 0d0,0d0) is called, then x*(10. pb) would be the SI cross-section
    ! at a WIMP mass of 100 GeV at which the experiment is excluded at
    ! the 90% CL (p=1-CL).
    lnp = LOG(0.1d0)
    WRITE(*,'(A20,5(2X,1PG11.4))')  'Max gap x for 90% CL ', &
        XENON100_2012_ScaleToPValue(lnp),LUX_2013_ScaleToPValue(lnp),   &
        DARWIN_Ar_2014_ScaleToPValue(lnp),DARWIN_Xe_2014_ScaleToPValue(lnp), &
        DDCalc0_ScaleToPValue(D,lnp)
    WRITE(*,'(A60)')  '  * Factor x such that sigma->x*sigma gives desired p-value'
    
    !WRITE(*,*)
    
  END DO
  
  
  CONTAINS

  ! --------------------------------------------
  ! Write a description of how input parameters should be specified
  SUBROUTINE WriteDescription(type)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: type 
    WRITE(*,'(A)') ''
    WRITE(*,'(A)') 'Enter WIMP parameters below.  Only the first two are necessary.'
    WRITE(*,'(A)') 'A blank line terminates input.  The parameters are:'
    WRITE(*,'(A)') ''
    SELECT CASE(type)
    CASE(TYPE_MG)
      WRITE(*,'(A)') '  M     WIMP mass [GeV]'
      WRITE(*,'(A)') '  GpSI  Spin-independent WIMP-proton effective coupling [GeV^-2]'
      WRITE(*,'(A)') '  GnSI  Spin-independent WIMP-neutron effective coupling [GeV^-2]'
      WRITE(*,'(A)') '  GpSD  Spin-dependent WIMP-proton effective coupling [GeV^-2]'
      WRITE(*,'(A)') '  GnSD  Spin-dependent WIMP-neutron effective coupling [GeV^-2]'
    CASE(TYPE_MFA)
      WRITE(*,'(A)') '  M     WIMP mass [GeV]'
      WRITE(*,'(A)') '  fp    Spin-independent WIMP-proton effective coupling [GeV^-2]'
      WRITE(*,'(A)') '  fn    Spin-independent WIMP-neutron effective coupling [GeV^-2]'
      WRITE(*,'(A)') '  ap    Spin-dependent WIMP-proton effective coupling [unitless]'
      WRITE(*,'(A)') '  an    Spin-dependent WIMP-neutron effective coupling [unitless]'
    CASE(TYPE_MSIGMA)
      WRITE(*,'(A)') '  M         WIMP mass [GeV]'
      WRITE(*,'(A)') '  sigmapSI  Spin-independent WIMP-proton cross-section [pb]'
      WRITE(*,'(A)') '  sigmanSI  Spin-independent WIMP-neutron cross-section [pb]'
      WRITE(*,'(A)') '  sigmapSD  Spin-dependent WIMP-proton cross-section [pb]'
      WRITE(*,'(A)') '  sigmanSD  Spin-dependent WIMP-neutron cross-section [pb]'
      WRITE(*,'(A)') ''
      WRITE(*,'(A)') 'Negative cross-section values can be given to indicate the'
      WRITE(*,'(A)') 'corresponding coupling should be taken to be negative.'
    END SELECT
    !WRITE(*,*) ''
  END SUBROUTINE


  ! --------------------------------------------
  ! Read WIMP parameters from standard input
  FUNCTION GetWIMPParams(type,M,xpSI,xnSI,xpSD,xnSD) RESULT(valid)
    IMPLICIT NONE
    LOGICAL :: valid
    INTEGER, INTENT(IN) :: type
    REAL*8, INTENT(OUT) :: M,xpSI,xnSI,xpSD,xnSD
    CHARACTER*256 :: line
    INTEGER :: I,ios
    REAL*8 :: x(5)
    
    valid = .FALSE.
    
    WRITE(*,'(A)') ''
    WRITE(*,'(A)') '------------------------------------------------------------'
    SELECT CASE(type)
    CASE(TYPE_MG)
      WRITE(*,'(A)') 'Enter values <M GpSI GnSI GpSD GnSD>:'
    CASE(TYPE_MFA)
      WRITE(*,'(A)') 'Enter values <M fp fn ap an>:'
    CASE(TYPE_MSIGMA)
      WRITE(*,'(A)') 'Enter values <M sigmapSI sigmanSI sigmapSD sigmanSD>:'
    END SELECT
    READ(*,'(A)') line
    IF (TRIM(line) .EQ. '') RETURN
    
    I = 6
    ios = -1
    DO WHILE ((I .GT. 2) .AND. (ios .NE. 0))
      I = I - 1
      READ(line,*,IOSTAT=ios) x(1:I)
    END DO
    IF (ios .NE. 0) RETURN
    
    valid = .TRUE.
    M    = x(1)
    xpSI = x(2)
    xnSI = x(3)
    xpSD = x(4)
    xnSD = x(5)
    IF (I .LT. 3) xnSI = xpSI
    IF (I .LT. 4) xpSD = 0d0
    IF (I .LT. 5) xnSD = xpSD
    
  END FUNCTION
  
END PROGRAM



!#######################################################################
! END OF FILE
!#######################################################################

