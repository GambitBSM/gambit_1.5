!#######################################################################
! DIRECT DETECTION RATES AND LIKELIHOODS (SIMPLE VERSION)
! Program to calculate the dark matter direct detection event rates
! and corresponding likelihoods/exclusion levels.  See DDCalc0run.use
! for instructions.
! 
! To compile (order is important!):
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
!#######################################################################

PROGRAM DDCalc0run
  USE DDCALC0
  IMPLICIT NONE
  CALL DDCalc0_Main()
END PROGRAM


!#######################################################################
! END OF FILE
!#######################################################################

 
