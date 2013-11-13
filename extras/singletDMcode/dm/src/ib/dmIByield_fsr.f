*****************************************************************************
***   function dsIByield_fsr gives the photon yield from final state radiation 
***   (FSR) from the deacy of a hypothetical partical with mass 2*m0, such as
***   included in the Pythia runs. Just like in Pythia, only FSR from fermionic
***   final states is included
***   Input: egev - energy in GeV
***          yieldk - 52 for integrated or 152 for differential yield
***   Output: yield (annihilation)**1, differential also GeV**-1
***           istat is set as follows in case of errors
***   bit  decimal  reason
***     0        1  dsIBf_intdxdy failed
***     1        2  dsIBf_intdy failed
*** Author: Torsten Bringmann, troms@physto.se
*** Date: 2008-02-10
*****************************************************************************

      real*8 function dmIByield_fsr(egev,yieldk,istat)
      implicit none
c      include 'dssusy.h'
      include 'dshacom.h'
 

c------------------------ functions ------------------------------------

      real*8 dmIByieldone_fsr

c------------------------ variables ------------------------------------

      real*8 egev,yield,sv,dmsigmav
      real*8 FSRabr(4:12)               ! only fermions
      integer FSRch,istat,yieldk,itmp


      istat=0
      yield=0.0d0

      if (yieldk.eq.52.or.yieldk.eq.152) then 
        FSRabr(4)  = habr(15)     ! e+ e-
        FSRabr(5)  = habr(17)     ! mu+ mu-
        FSRabr(6)  = habr(19)     ! tau+ tau-
        FSRabr(7)  = habr(20)     ! u u-bar
        FSRabr(8)  = habr(21)     ! d d-bar
        FSRabr(9)  = habr(22)     ! c c-bar
        FSRabr(10) = habr(23)     ! s s-bar
        FSRabr(11) = habr(24)     ! t t-bar
        FSRabr(12) = habr(25)     ! b b-bar

c... at the moment, for some of these channels there are no Pythia runs.
c... For these, of course nothing should be subtracted
c... THIS IS A TEMPORARY FIX        
        FSRabr(4)  = 0d0
        FSRabr(7)  = 0d0
        FSRabr(8)  = 0d0
        FSRabr(10) = 0d0

c... now compute the contribution from each channel
        do 200 FSRch=4,12     
          if (FSRabr(FSRch).gt.0.0d0) then  
            yield=yield+FSRabr(FSRch)*
     &             dmIByieldone_fsr(hamwimp,egev,FSRch,yieldk,itmp)
           endif
  200   continue
      endif

      dmIByield_fsr=yield

      end
