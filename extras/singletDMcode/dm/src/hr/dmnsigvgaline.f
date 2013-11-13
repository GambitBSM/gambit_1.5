**********************************************************************
*** subroutine dsnsigvgaline gives the number of photons times
*** the annihilation cross section into gamma gamma and Z gamma
*** respectively. The result given is the dimensionless number
***     N_gamma * (sigma * v) / (10-29 cm^3 s^-1)
***
*** author: joakim edsjo, edsjo@physto.se
*** date: 00-09-03
**********************************************************************

      subroutine dsnsigvgaline(nsigvgaga,nsigvgaz)
      implicit none
      real*8 nsigvgaga,nsigvgaz
      include 'dsprep.h'

c-----------------------------------------------------------------------

      if (.not.dsprepcalled) then
        write(*,*) 'error in dsnsigvgaline: dsprep must be called',
     &    ' before any rate calculation'
        write(*,*) 'begins for every new model. program stopping.'
      endif

c... gamma-gamma
      nsigvgaga=2.0d0*(sigv(28)/1.0d-29)
c... z gamma
      nsigvgaz=(sigv(29)/1.0d-29)
      end





