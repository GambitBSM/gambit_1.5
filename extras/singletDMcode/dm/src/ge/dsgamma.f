!=======================================================================
! Calculates gamma(x) to double precision.
!     gamma(x) = \int_0^\infty dt t^{x-1} e^{-t}
! Returns a large negative number where gamma is undefined
! (x=0,-1,-2,...).
! 
! 
!    Created by Chris Savage (savage@fysik.su.se)
!    2011/05/23
! 
!=======================================================================
! 
      REAL*8 FUNCTION dsgamma(x)
      IMPLICIT NONE
      REAL*8 x,x0,lngamma
      REAL*8 PI
      PARAMETER(PI=3.1415926535897932d0)
      ! We make use of the identity:
      !   Gamma(1-z) = pi*z / Gamma(1+z) / sin(pi*z)
      IF (x .GE. 1d0) THEN
        dsgamma = EXP(lngamma(x))
      ! Gamma is undefined when x=0,-1,-2,...
      ELSE IF (MOD(x,1d0) .EQ. 0d0) THEN
        dsgamma = -HUGE(x)
      ELSE
        x0 = 1d0 - x
        dsgamma = PI*x0 / (SIN(PI*x0) * EXP(lngamma(1d0+x0)))
      END IF
      END FUNCTION

