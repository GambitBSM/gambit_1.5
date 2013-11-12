**********************************************************************
*** subroutine dsnsigvgacdiff gives the number of contiuous gammas
*** per GeV at the energy ega (in GeV) times the annihilation cross
*** section.
*** The result given is the number
***     N_gamma * (sigma * v) / (10-29 cm^3 s^-1)
*** in units of GeV^-1.
***
*** author: joakim edsjo, edsjo@physto.se
*** date: 00-09-03
**********************************************************************

      subroutine dsnsigvgacdiff(ega,nsigvgacdiff)
      implicit none
      real*8 nsigvgacdiff,dshaloyield,yieldga,ega,dssigmav
      integer istat
      include 'dsprep.h'

c-----------------------------------------------------------------------

      if (.not.dsprepcalled) then
        write(*,*) 'error in dsnsigvgacontdiff: dsprep must be called',
     &    ' before any rate calculation'
        write(*,*) 'begins for every new model. program stopping.'
      endif

      yieldga=dshaloyield(ega,152,istat)
      nsigvgacdiff=yieldga*(dssigmav(0)/1.0d-29)

      end





