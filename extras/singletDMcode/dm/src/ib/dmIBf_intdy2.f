*****************************************************************************
*** integrates a function f(z,y) over the kinematic variable 
***
*** y = (p+k)^2/(4 mx^2)
***
*** that appears for 3-body final states. (k is the photon momentum and p
*** that of one of the other final sates). The result is the
*** (dimensionless) differential p yield per annihilation
***
*** called by dsIByieldone
***
*** author: Torsten Bringmann (troms@physto.se)
*** date:   2007-10-20
*** update: 2008-03-10 - better numerical convergence
*****************************************************************************

      real*8 function dmIBf_intdy2(IBch,z,mwimp,mp1,mp2)

      implicit none
      include 'dsprep.h'  
      include 'dsibcom.h'

      real*8 z,y,mwimp,mp1,mp2,eta
      integer IBch

c... these parameters are needed by dqagse
      real*8 a,b,result,abserr
      integer neval,ier, limit
      parameter (limit=20)
      real*8 alist(limit),blist(limit),rlist(limit),elist(limit)
      integer iord(limit),last
   
      real*8 dmIBintsel
      external dmIBintsel


      dmIBf_intdy2=0.d0

c... set parameters needed by dmIBintsel
        if (IBch.eq.1.or.(IBch.ge.4.and.IBch.le.12)) then
          intch = -IBch
        else
          intch = 0
        endif
        ibcom_z = z
        ibcom_mx = mwimp
        ibcom_mp1 = mp1
        ibcom_mp2 = mp2

c...  if kinematical limits satisfied, determine integration limits
      if ((1d0.gt.z).and.(z.gt.(mp1/mwimp))) then 
         eta=(mp1/mwimp/2.)**2
         a=eta**2/((1 + eta - z)) + 
     -        ( (1. - Sqrt(1.-4.*eta/z**2))*(1. - z)*z )/
     -        (1 + eta - z)/2.
         b=eta**2/((1. + eta - z)) +
     -        ( (1. + Sqrt(1.-4.*eta/z**2))*(1. - z)*z )/
     -        (1. + eta - z)/2.
      else
         return
      endif

      call dqagse(dmIBintsel,a,b,1d-10,IBacc,20,result,abserr,
     &            neval,ier,alist,blist,rlist,elist,iord,last)

      if (ier.ne.0) then
        if (IBch.eq.4) iberr=ibset(iberr,2)
        if (IBch.ne.4) iberr=ibset(iberr,4)
      endif

      dmIBf_intdy2=result

 200  end
