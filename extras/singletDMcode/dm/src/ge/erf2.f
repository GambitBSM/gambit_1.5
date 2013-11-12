!=======================================================================
! Error function between x1 and x2, i.e. erf(x2)-erf(x1).
! This routine accounts for the cases where loss of precision will
! result from canceling in an explicit subtraction of calls to ERF(x1)
! and ERF(x2).
! 
! 
!    Created by Chris Savage (savage@fysik.su.se)
!    2011/05/23
! 
!=======================================================================
! 
      REAL*8 FUNCTION erf2(x1,x2)
      IMPLICIT NONE
      REAL*8 x1,x2
      REAL*8 xc,delx
      REAL*8 SQRTPI
      PARAMETER(SQRTPI=1.7724538509055160d0)
      ! 1-eps is approximate level of cancelation at which to handle
      ! the difference more carefully
      REAL*8 EPS
      PARAMETER(EPS=0.03d0)
      ! Comment following line to use intrinsic erf and erfc
      ! Uncomment if intrinsic routines are not available
      !REAL*8 erf,erfc
      
      ! Opposite sign: no canceling to worry about here
      IF (x1*x2 .LE. 0d0) THEN
        erf2 = erf(x2) - erf(x1)
        RETURN
      END IF
      
      xc   = 0.5d0 * (x1+x2)
      delx = 0.5d0 * (x2-x1)
      
      ! Smaller arguments:
      !   |xc| < 1    --> |x1|,|x2| < 2
      ! Canceling is significant if |delx| < eps*|xc|
      IF ((ABS(xc) .LE. 1d0) .AND. (ABS(delx) .GT. EPS*ABS(xc))) THEN
        ! erf(x2) - erf(x1)
        erf2 = erf(x2) - erf(x1)
        RETURN
        
      ! At least one argument is "large": 
      !   |xc| > 1    --> |x1| > 1 or |x2| > 1
      ! Canceling is significant if |4*xc*delx| < eps
      ELSE IF ((ABS(xc) .GT. 1d0) .AND. (ABS(4*xc*delx) .GT. EPS)) THEN
        IF (xc .GT. 0d0) THEN
          ! Difference of complementary error function gives better
          ! precision here:
          ! erf(x2) - erf(x1) = erfc(x1) - erfc(x2)
          erf2 = -(erfc(x2) - erfc(x1))
        ELSE
          ! Difference of complementary error function gives better
          ! precision here (use symmetry of erf(x)):
          ! erf(x2) - erf(x1) = erf(-x1) - erf(-x2) = erfc(-x2) - erfc(-x1)
          erf2 = erfc(-x2) - erfc(-x1)
        END IF
        RETURN
      END IF
      
      ! If we reached this point, there is significant canceling.
      ! For these cases, the integrand in the error function does not
      ! change much over x1 < x < x2.  Taylor expand the integrand
      ! about xc = (x1+x2)/2 and integrate.  The following keeps terms
      ! up through tenth order.
      erf2 = 4 * delx * EXP(-xc**2) / SQRTPI
     &         * (1 + (2*xc**2 - 1)*delx**2 / 3
     &              + (4*xc**4 - 12*xc**2 + 3)*delx**4 / 30
     &              + (8*xc**6 - 60*xc**4 + 90*xc**2 - 15)*delx**6 / 630
     &              + (16*xc**8 - 224*xc**6 + 840*xc**4 - 840*xc**2
     &                          + 105)*delx**8 / 22680
     &           )
      
      END FUNCTION


