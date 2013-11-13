*****************************************************************************
***   function dsIByield_user is a dummy function to give the yield from 
***   internal bremsstrahlung (IB) for a user-defined model. 
***   
***   Input: egev - energy in GeV
***          yieldk - which yield to calculate
***             (see dshaloyield for an explanation)
***   Output: yield (annihilation)**1, differential also GeV**-1
***           istat is a flag for internal error handling 
***             (see dshaloyield for an explanation) 
***
*** Author: Torsten Bringmann, troms@physto.se
*** Date: 2008-02-10
*****************************************************************************

      real*8 function dmIByield_user(egev,yieldk,istat)
      implicit none
      include 'dssusy.h'
      include 'dshacom.h'
      include 'dshrcom.h'
      include 'dsidtag.h'
      include 'dsprep.h'
      include 'dsibcom.h'

      real*8 egev,yieldk
      integer istat


c...  see dsIByield for help on how to implement IB

      dmIByield_user=0.0d0
 

      end



















