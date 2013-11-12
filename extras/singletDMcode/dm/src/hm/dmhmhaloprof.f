****************************************************************
***                                                          ***  
*** mod: 04-01-12 pu, this is obsolete and should not        *** 
*** be used anymore !!!!!!!!!!!!!                            ***
***                                                          ***  
*** dark matter halo density profile                         ***
*** it assumes that the halo                                 ***
***   1) is spherically symmetric                            ***
***   2) has a double power law profile: the                 ***
***      (alphah,beta,gamma) zhao model,  where - gammah is  ***
***      the slope towards the galactic centre, - beta is    ***
***      the slope at large galactocentric distances and     ***
***      alphah determines the width of the transition zone. ***
***      e.g.: modified isothermal sphere profile = (2,2,0); ***
***            nfw profile = (1,3,1);                        ***
***            moore et al. profile = (1.5,3,1.5)            ***
***                                                          ***  
*** rr = galactocentric distance in kpc                      ***
*** a  = length scale in kpc                                 ***
*** r_0 = galactocentric distance of the sun in kpc          ***
*** rho0 = local dark matter density (=val(r_0)) in gev/cm**3***
***                                                          ***
*** the profile is truncated at 10**-5 kpc assuming          ***
*** val(rr<10**-5 kpc) = val(rr=10**-5 kpc)                  ***
***                                                          ***
*** author: piero ullio (piero@tapir.caltech.edu)            ***
*** date: 00-07-13                                           ***
****************************************************************

      subroutine dshmhaloprof(rr,val)
      implicit none
      real*8 fun,esp,rr,val,rraux
      include 'dshmcom.h'
      if(rr.lt.1.d-5) then
        rraux=1.d-5
        fun=(rraux/r_0)**(-gammah)
        esp=(betah-gammah)/alphah
        fun=fun/((1.d0+(rraux/ah)**alphah)**esp)
        val=rho0*fun*((1.d0+(r_0/ah)**alphah)**esp)
      else
        fun=(rr/r_0)**(-gammah)
        esp=(betah-gammah)/alphah
        fun=fun/((1.d0+(rr/ah)**alphah)**esp)
        val=rho0*fun*((1.d0+(r_0/ah)**alphah)**esp)
      endif
      return
      end









