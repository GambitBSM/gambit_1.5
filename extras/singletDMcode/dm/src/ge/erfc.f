!=======================================================================
! Calculates erfc(x) to double precision.
! NOTE: This routine is provided for use with compilers that do not
!       provide an implementation themselves.  The intrinsic routine
!       should be used if available (it is probably much faster).
! 
! 
!    Created by Chris Savage (savage@fysik.su.se)
!    2011/05/23
! 
!=======================================================================
! 
      REAL*8 FUNCTION erfc(x)
      IMPLICIT NONE
      REAL*8 x,erf0,erfc0
      CALL erfsub(x,erf0,erfc0)
      erfc = erfc0
      END FUNCTION


