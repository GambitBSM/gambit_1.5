**********************************************************************
*** subroutine dsnsigvgacont gives the number of photons above
*** the threshold egath (in GeV) times the annihilation cross section
*** into continuous gammas.
*** The result given is the dimensionless number
***     N_gamma * (sigma * v) / (10-29 cm^3 s^-1)
***
*** author: joakim edsjo, edsjo@physto.se
*** date: 00-09-03
**********************************************************************

      subroutine dsnsigvgacont(egath,nsigvgacont)
      implicit none
      real*8 nsigvgacont,dshaloyield,yieldga,egath,dssigmav
      integer istat
      include 'dsprep.h'

c-----------------------------------------------------------------------

      if (.not.dsprepcalled) then
        write(*,*) 'error in dsnsigvgacont: dsprep must be called',
     &    ' before any rate calculation'
        write(*,*) 'begins for every new model. program stopping.'
      endif

      yieldga=dshaloyield(egath,52,istat)
      nsigvgacont=yieldga*(dssigmav(0)/1.0d-29)

      end





