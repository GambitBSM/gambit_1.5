*****************************************************************************
*** computes the envelope of the yield in channel IBch with the yield of
*** type yieldk due to the further decay of the particles in IBch
***
*** author: Torsten Bringmann (troms@physto.se)
*** date  : 2007-11-20
*** update: 2008-03-10 - better numerical convergence
*****************************************************************************

      real*8 function dmIBf_intdE(IBch,yieldk,x,mwimp,mp1,mp2)

      implicit none
      include 'dssusy.h'
      include 'dshacom.h'
      include 'dsprep.h'  
      include 'dsidtag.h'
      include 'dsibcom.h'


      real*8 a,b,x,mwimp,mp1,mp2,mptest
      real*8 dmIBf_intdy2
      integer yieldk,IBch

      real*8 dmf_int, dmIBintsel2
      external dmIBintsel2


      dmIBf_intdE=0.d0

c... set parameters needed by dmIBintsel
      if ((IBch.eq.1).or.
     &    (IBch.ge.5).and.(IBch.le.12)) then
          intch = -IBch
          intyield = yieldk
          ibcom_mx = mwimp
          ibcom_x = x
          ibcom_mp1 = mp1
          ibcom_mp2 = mp2
      else
          return
      endif


c... preliminary fix: dshayield uses different definitions of masses...
        if (IBch.eq.1) then
              mptest=map(chcomp(13))
           elseif (IBch.eq.5) then
              mptest=map(chcomp(17))
           elseif (IBch.eq.6) then
              mptest=map(chcomp(19))
           elseif ((IBch.ge.7).and.(IBch.le.10)) then
              mptest=map(chcomp(22))
           elseif (IBch.eq.11) then
              mptest=map(chcomp(24))
           elseif (IBch.eq.12) then
              mptest=map(chcomp(25))
         endif
         

c... determine integration limits;
c... (integration over xp = energy of IBch particle, in units of mwimp)
      a=1.03d0*x
      if (x.lt.(1.01d0*mptest/mwimp)) a=1.01d0*mptest/mwimp
c       if (x.lt.(1.1*mp1/mwimp)) a=1.1d0*mp1/mwimp
      b=0.999d0

      dmIBf_intdE=dmf_int(dmIBintsel2,a,b,IBacc)

      end



