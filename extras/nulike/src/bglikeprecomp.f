***********************************************************************
*** nulike_bglikeprecomp calls the calculation of the p-value for the
*** background in IceCube calculations based on Poissonian statistics.
***        
*** Author: Pat Scott (patscott@physics.mcgill.ca)
*** Date: May, July 2011
***********************************************************************

      subroutine nulike_bglikeprecomp

      implicit none
      include 'nulike.h'

      real*8 nulike_pval
      real*8 theta_sum
      real*8 events_sum
      theta_sum = sum(theta_BG)
      events_sum = sum(nEvents_inEAErrBins)

      BGpvalPoissonian = nulike_pval(events_sum,theta_sum,0.d0)
      pvalBGPoisComputed = .true.

      end subroutine nulike_bglikeprecomp
