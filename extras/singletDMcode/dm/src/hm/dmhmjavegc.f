****************************************************************
*** function dshmjavegc: average over the solid angle        ***
*** delta (sr) of the function dshmj(cospsi0) in case of the *** 
*** galactic center                                          ***
***                                                          ***
*** dshmj is the line of sight integral which enters in the  ***
*** computation of the gamma-ray and neutrino fluxes from    ***
*** pair annihilations of wimps in the halo.                 ***
***                                                          ***
*** see definition in e.g. bergstrom et al.,                 ***
***   phys. rev. d59 (1999) 043506                           ***
*** in case of the many unresolved clump scenario the term   ***
*** fdelta is factorized out                                 *** 
***                                                          ***  
*** it is valid for a spherical dark matter halo             ***
*** psi0, which must be 0.d0 is the angle between direction  ***
*** of observation and the direction of the galactic center; ***
*** cospsi0 is its cosine.                                   ***
***                                                          ***
*** author: piero ullio (piero@tapir.caltech.edu)            ***
*** date: 00-07-13                                           ***
****************************************************************


      real*8 function dshmjavegc(cospsi0in,deltain)
      implicit none
      include 'dshmcom.h'
      integer ier,iord,last,limit,neval
      real*8 abserr,alist,blist,elist,epsrel,result,rlist
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      real*8 cospsi0in,deltain
      real*8 cospsi0,pi,delta,par,rinf,rsup,rint(0:8)
      real*8 cosinf,cossup,checkres
      integer k
      real*8 cospsimin,cospsi,rr
      common/dshmjavegccom/cospsimin,cospsi,rr
      external dshmjavepar3,dshmjavepar5
      real*8 epsabs
      integer hmjavegcset,ndigits
      common/hmjavegcsetcom/epsabs,hmjavegcset,ndigits
      if(hmjavegcset.ne.123456) then
         ndigits=4
         epsabs=1.d-10          !numerical accuracy
         hmjavegcset=123456
      endif
      limit=5000
      epsrel=10.d0**(-ndigits)
      cospsi0=cospsi0in
      if(cospsi0.ne.1.d0) then
         write(*,*) 'dshmjavegc must be called for the galactic'
         write(*,*) 'center! it has been called instead for'
         write(*,*) 'cospsi =',cospsi0in
         stop
      endif   
      delta=deltain
      pi=4.d0*datan(1.d0)
      cospsimin=((2.d0*pi-delta)/(2.d0*pi))  
          !cosine of the cone's aperture angle
      if(cospsimin.le.0.d0) then
        write(*,*) 'dshmjavegc has been called for delta=',delta
        write(*,*) 'larger than 2*pi!'
        write(*,*) 'this routine does not work in this case'
        stop
      endif  
      par=0.d0
      cossup=1.d0
      cosinf=cospsimin
 11   call dqagseb(dshmjavepar3,cosinf,cossup,epsabs,epsrel,limit
     &  ,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      checkres=par+2.d0*pi*result
      if(dabs(checkres)/10.d0**ndigits.lt.epsabs) then
        epsabs=dabs(checkres)/10.d0**(ndigits+2) !too low accuracy
        goto 11
      elseif(dabs(checkres)/10.d0**(ndigits+2).gt.epsabs) then
        epsabs=dabs(checkres)/10.d0**(ndigits+2)  !too high accuracy 
      endif
      par=par+2.d0*pi*result
      rint(0)=0.d0
      rint(1)=1.d-5
      rint(5)=r_0*dsqrt(1.d0+cospsimin)*dsqrt(1.d0-cospsimin)
      do k=2,4
        rint(k)=10.d0**(dlog10(rint(1))+
     &                  (dlog10(rint(5))-dlog10(rint(1)))/4.d0*(k-1)) 
      enddo
      do k=0,4
        rinf=rint(k)
        rsup=rint(k+1)
 12     call dqagseb(dshmjavepar5,rinf,rsup,epsabs,epsrel,limit
     &   ,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        checkres=par+2.d0*pi*result
        if(dabs(checkres)/10.d0**ndigits.lt.epsabs) then
          epsabs=dabs(checkres)/10.d0**(ndigits+2) !too low accuracy
          goto 12
        elseif(dabs(checkres)/10.d0**(ndigits+2).gt.epsabs) then
          epsabs=dabs(checkres)/10.d0**(ndigits+2)  !too high accuracy 
        endif
        par=par+2.d0*pi*result
      enddo  
      dshmjavegc=par/8.5d0/delta*(rho0/0.3d0)**2
      return
      end
