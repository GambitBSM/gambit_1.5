****************************************************************
*** dark matter halo density profile in case of the          ***
*** (alphah,beta,gamma) zhao model.                          ***
*** it is a double power law profile,  where -  gamma is     ***
***      the slope towards the galactic centre, - beta is    ***
***      the slope at large galactocentric distances and     ***
***      alphah determines the width of the transition zone. ***
***      e.g.: modified isothermal sphere profile = (2,2,0); ***
***            nfw profile = (1,3,1);                        ***
***            moore et al. profile = (1.5,3,1.5)            ***
***                                                          ***  
*** radialdist = galactocentric distance in kpc              ***
*** ah  = length scale in kpc                                ***
*** rhoref = dark matter density in gev/cm**3 at the         ***
*** galactocentric distance Rref (in kpc)                    ***
***                                                          ***
*** Author: Piero Ullio (ullio@sissa.it)                     ***      
*** Date: 2004-01-12                                         ***
****************************************************************

      real*8 function dmhmabgrho(radialdist)

      implicit none
      real*8 radialdist,fun,esp,rr

      include 'dshmcom.h'

      rr=radialdist
      fun=(rr/Rref)**(-gammah)
      esp=(betah-gammah)/alphah
      fun=fun/((1.d0+(rr/ah)**alphah)**esp)
      dmhmabgrho=rhoref*fun*((1.d0+(Rref/ah)**alphah)**esp)
 
      return
      end









