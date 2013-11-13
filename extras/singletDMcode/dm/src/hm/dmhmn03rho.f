****************************************************************
*** dark matter halo density profile in case of the          ***
*** navarro et al. (2003) model.                             ***
***                                                          ***  
*** radialdist = galactocentric distance in kpc              ***
*** an03  = length scale in kpc                              ***
*** rhoref = dark matter density in gev/cm**3 at the         ***
*** galactocentric distance Rref (in kpc)                    ***
***                                                          ***
*** Author: Piero Ullio (ullio@sissa.it)                     ***      
*** Date: 2004-01-12                                         ***
****************************************************************

      real*8 function dmhmn03rho(radialdist)

      implicit none
      real*8 radialdist,x

      include 'dshmcom.h'

      x=radialdist/an03

      dmhmn03rho=dexp(-2.d0/alphan03*(x**alphan03-1.d0))
      x=Rref/an03
      dmhmn03rho=dmhmn03rho
     &   *rhoref/dexp(-2.d0/alphan03*(x**alphan03-1.d0))

      return
      end









