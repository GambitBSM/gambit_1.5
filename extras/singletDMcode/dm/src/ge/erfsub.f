!=======================================================================
! Calculates both erf(x) and erfc(x) to double precision.
! There is a negligible computational penalty to calculating both
! quantities in this routine compared to calculating either one alone.
! 
! 
!    Created by Chris Savage (savage@fysik.su.se)
!    2011/05/23
! 
!=======================================================================
! 
      SUBROUTINE erfsub(x,erf,erfc)
      IMPLICIT NONE
      REAL*8 x,erf,erfc
      REAL*8 x0
      ! Taylor expansion of erf(x)
      REAL*8 xsq,zk,series_sum
      ! Chebyshev polynomial for erfc(|x|)
      REAL*8 XMIN,XMAX,t,TA,TB,y,dk,dk1,dk2
      INTEGER K
      ! Misc
      REAL*8 SQRTPI
      PARAMETER(SQRTPI=1.7724538509055160d0)
      ! Ranges of x and t = 2/(2+x) for Chebyshev method:
      ! xmin = 1/2       ta = 2/(2+xmin) = 4/5
      ! xmax = \infty    tb = 2/(2+xmax) = 0
      !PARAMETER(XMIN=0.5d0,XMAX=HUGE(x),TA=2/(2+XMIN),TB=2/(2+XMAX))
      PARAMETER(TA=2/(2+0.5d0),TB=0d0)
      ! Chebyshev polynomial coefficients
      ! For single precision, use only C(0:10)
      INTEGER NC
      PARAMETER(NC=24)
      REAL*8 C(0:NC)
      DATA C / -1.5663741454481200d+00,
     &    -5.0589479331734064d-01, +2.0165332295700412d-02,
     &    +4.1507842704837142d-03, -6.9293930331318290d-04,
     &    -8.0613044480901668d-05, +2.5850999510766173d-05,
     &    +2.4532413619366540d-06, -1.0665628949056517d-06,
     &    -1.0956136201272647d-07, +4.6182837214229254d-08,
     &    +6.3050891373696668d-09, -1.9633309872353595d-09,
     &    -3.9793711163649333d-10, +7.3597120815386668d-11,
     &    +2.4619582503772174d-11, -1.7696477486824439d-12,
     &    -1.3810080809606214d-12, -4.6002066179630267d-14,
     &    +6.3804840465520990d-14, +9.9377184182455142d-15,
     &    -1.8569583472104992d-15, -7.7079563859333715d-16,
     &    -3.1621351296209296d-17, +3.4489452474809714d-17  /
      
      
      ! NOTE:
      ! Numerical Recipes uses a Chebyshev method to calculate erfc(x)
      ! over all x, then uses erf(x) = 1 - erfc(x).  The Chebyshev
      ! method works well for erfc(x), but there is a significant loss
      ! of precision when calculating erf(x) this way for small x
      ! (erfc is nearly 1 in this case, so the difference of erfc with
      ! 1 has a much smaller relative precision due to the canceling).
      ! Here, we ensure full precision over all x by calculating
      ! the (approximately) smaller of erf or erfc directly and using
      ! erf + erfc = 1 to find the other, thereby avoiding this
      ! loss of precision due to canceling.
      
      x0 = ABS(x)
      
      ! For |x| < 0.5, we calculate erf(x) using a Taylor expansion
      ! about x=0.
      !   single precision in 7 terms
      !   double precision in 12 terms
      IF (x0 .LE. 0.5d0) THEN
        xsq = x**2
        K = 0
        zk = 1
        series_sum = zk
        DO K = 1,12
          zk = -zk*xsq/K
          series_sum = series_sum + zk/(2*K+1)
        END DO
        erf = 2 * x * series_sum / SQRTPI
        erfc = 1 - erf
      
      ! For |x| > 0.5, we use a Chebyshev interpolation to calculate
      ! erfc(x).  We use this functional form:
      !     erfc(x) = t Exp[-x^2 + P(t)]
      ! where t = 2/(2+x) and P(t) is a polynomial obtained using the
      ! Chebyshev interpolation technique.  See Numerical Recipes
      ! (Press et al., 2007) for a description.  Note the range used
      ! here differs from Numerical Recipes, so the coefficients do
      ! not match.  Aside from the difference in range, the quantity
      ! 'y' (below) differs from the NR quantity 'ty' by a sign
      ! difference (y = -ty, if the limits were the same).
      ELSE
        ! Use t(x) = 2/(2+x), but we linearly rescale t into the
        ! quantity y such that -1 < y/2 < 1 over the range of t.
        ! The Chebyshev interpolating polynomial is actually
        ! constructed in the quantity y/2.  There is a sign
        ! ambiguity in the linear transformation from t to y that
        ! has been chosen differently between here and NR.  The
        ! two cases are related by y -> -y and C_k -> -C_k (k odd).
        t  = 2/(2+x0)
        !TA = 2/(2+XMIN)                   ! 0.8
        !TB = 2/(2+XMAX)                   ! 0.0
        y  = 2*(2*t - (TA+TB))/(TB-TA)     ! y = -(5t-2)
        ! Use recurrence relation to determine P(y)
        dk2 = 0
        dk1 = 0
        DO K = NC,1,-1
          dk = y*dk1 - dk2 + C(K)
          dk2 = dk1
          dk1 = dk
        END DO
        dk = (y*dk1 - 2*dk2 + C(0))/2
        ! This is erfc(|x|)
        erfc = t * EXP(-x**2 + dk)
        IF (x .LT. 0) THEN
          erf  = erfc - 1
          erfc = 1 - erf
        ELSE
          erf = 1 - erfc
        END IF
      END IF
      
      END SUBROUTINE


