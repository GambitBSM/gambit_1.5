****************************************************************
*** dark matter halo density profile in case of the          ***
*** burkert model.                                           ***
***                                                          ***  
*** radialdist = galactocentric distance in kpc              ***
*** ah  = length scale in kpc                                ***
*** rhoref = dark matter density in gev/cm**3 at the         ***
*** galactocentric distance Rref (in kpc)                    ***
***                                                          ***
*** Author: Piero Ullio (ullio@sissa.it)                     ***      
*** Date: 2004-01-12                                         ***
****************************************************************

      real*8 function dmhmburrho(radialdist)

      implicit none

      real*8 radialdist,x

      include 'dshmcom.h'
      x=radialdist/ah
      dmhmburrho=1.d0/(1.d0+x)/(1.d0+x**2)
      x=Rref/ah
      dmhmburrho=dmhmburrho*rhoref*(1.d0+x)*(1.d0+x**2)
      return
      end









