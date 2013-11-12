************************************************************************
*** function dshaib returns the yield from internal bremsstrahlung
*** for the current model.
***
*** The switch haib determines which calculation is used for
*** the internal bremsstrahlung.
*** Viable optons for the switch haib are:
***
*** haib       Meaning
*** ---        -------
*** susy       Include calculation for SUSY models
*** none       Include nothing (but keep the model-independent final
***            state radiation contribution that is automatically included 
***            in the tabulated PYTHIA results)
*** no_fsr     Subtract FSR contained in PYTHIA results 
*** user       Include calculation by user-supplied function
***
*** For the option user, replace the dummy function dsIByield_user.f with your
*** own one and link to that one instead of the default dummy function.
***
*** Inputs: egev - energy in GeV
***         yieldk - which kind of yield to calculate
***            (see dshaloyield for an explanation)
*** Output: yield (annihilation)**-1, differential yields are also GeV**-1
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2008-01-15
*** no_fsr option added by Torsten Bringmann, troms@physto.se
************************************************************************

      real*8 function dmhaib(egev,yieldk,istat)

      implicit none

      include 'dshacom.h'

      real*8 egev
      integer yieldk,istat
      real*8 dmibyield, dmibyield_fsr, dmibyield_user


      if (haib.eq.'none') then
         dmhaib=0.0d0
         return
      elseif (haib.eq.'no_fsr') then
         dmhaib= - dmibyield_fsr(egev,yieldk,istat)
      elseif (haib.eq.'susy') then
         dmhaib=dmibyield(egev,yieldk,istat)
      elseif (haib.eq.'user') then
         dmhaib=dmibyield_user(egev,yieldk,istat)
      else
         write(*,*) 'DS: Error in dmhaib: unknown option haib=',haib
      endif

      return
      end

         
