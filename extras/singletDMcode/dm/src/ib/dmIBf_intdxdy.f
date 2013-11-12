*****************************************************************************
*** integrates a function f(x,y) over the kinematic variables
***
*** x = E_gamma/mx
*** y = (p+k)^2/(4 mx^2)
***
*** that appear for 3-body final states. (k is the photon momentum and p
*** that of one of the other final sates). The result is the
*** total photon yield per annihilation
***
*** called by dsIByield
***
*** based on paolo gondolos wxint.f routine.
*** author: Torsten Bringmann 2007-04-18
*****************************************************************************
      real*8 function dmIBf_intdxdy(IBch,x,mwimp,mp1,mp2)
      implicit none
      include 'dsprep.h'  
      include 'dsibcom.h'

      real*8 a,b,x,mwimp,mp1,mp2
      real*8 dmf_int, dmIBintsel2
      external dmIBintsel2
      integer IBch

      dmIBf_intdxdy=0D0

c... set parameters needed by dsIBintsel
      if ((IBch.ge.1.and.IBch.le.12).or.
     &    (IBch.ge.101.and.IBch.le.112)) then
          intch = IBch
      else
          return
      endif
      ibcom_mx=mwimp
      ibcom_mp1=mp1
      ibcom_mp2=mp2


c...determine integration limits
      a=x
      b=0.9999*(1d0-(mp1+mp2)**2/4./mwimp**2)
      if (a.ge.b) then
         return
      endif

      dmIBf_intdxdy=dmf_int(dmIBintsel2,a,b,IBacc)
  
      end
