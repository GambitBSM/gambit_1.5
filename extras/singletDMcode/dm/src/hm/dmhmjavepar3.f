****************************************************************
*** function integrated in dshmjavegc                        ***
***                                                          ***
*** author: piero ullio (piero@tapir.caltech.edu)            ***
*** date: 00-07-13                                           ***
****************************************************************

      real*8 function dshmjavepar3(cospsiin)
      implicit none
      include 'dshmcom.h'
      integer ier,iord,last,limit,neval
      real*8 abserr,alist,blist,elist,epsrel,result,rlist
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      real*8 cospsiin,par,rmin,rmax,rint(0:7),locinf,locsup,checkres
      integer k
      real*8 cospsimin,cospsi,rr
      common/dshmjavegccom/cospsimin,cospsi,rr
      external dshmjavepar4
      real*8 epsabs
      integer hmjavegc3set,ndigits
      common/hmjavegc3setcom/epsabs,hmjavegc3set,ndigits
      if(hmjavegc3set.ne.123456) then
         ndigits=4
         epsabs=1.d-10          !numerical accuracy
         hmjavegc3set=123456
      endif
      limit=5000
      epsrel=10.d0**(-ndigits)
      cospsi=cospsiin
      rmin=0.d0
      rmax=r_0*(cospsi-dsqrt(cospsi+cospsimin)*dsqrt(cospsi-cospsimin))
      rint(0)=rmin
      rint(5)=rmax*(1.d0-1.d-6)
      rint(6)=rmax
      do k=1,4
        rint(k)=rmax-10.d0**(dlog10(rmax)*(1.d0-1.d0/5.d0*k)
     &                       +dlog10(1.d-6*rmax)*(1.d0/5.d0*k))
      enddo
      par=0.d0
      do k=0,5
        locinf=rint(k)
        locsup=rint(k+1) 
 10     call dqagse(dshmjavepar4,locinf,locsup,epsabs,epsrel,limit,
     &    result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        checkres=par+result
        if(dabs(checkres)/10.d0**ndigits.lt.epsabs) then
          epsabs=dabs(checkres)/10.d0**(ndigits+2) !too low accuracy
          goto 10
        elseif(dabs(checkres)/10.d0**(ndigits+2).gt.epsabs) then
          epsabs=dabs(checkres)/10.d0**(ndigits+2)  !too high accuracy 
        endif
        par=par+result
      enddo
      rmin=r_0*(cospsi+dsqrt(cospsi+cospsimin)*dsqrt(cospsi-cospsimin))
      rmax=100.d0
      rint(0)=rmin
      rint(1)=rmin*(1.d0+1.d-6)
      rint(7)=rmax
      do k=2,6
        rint(k)=rmin+10.d0**(dlog10(rmin*1.d-6)*(1.d0-1.d0/6.d0*(k-1))
     &                       +dlog10(rmax-rmin)*(1.d0/6.d0*(k-1)))
      enddo
      do k=0,6
        locinf=rint(k)
        locsup=rint(k+1)
 11     call dqagse(dshmjavepar4,locinf,locsup,epsabs,epsrel,limit,
     &    result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        checkres=par+result
        if(dabs(checkres)/10.d0**ndigits.lt.epsabs) then
          epsabs=dabs(checkres)/10.d0**(ndigits+2) !too low accuracy
          goto 11
        elseif(dabs(checkres)/10.d0**(ndigits+2).gt.epsabs) then
          epsabs=dabs(checkres)/10.d0**(ndigits+2)  !too high accuracy 
        endif
        par=par+result
      enddo
      dshmjavepar3=par
      return
      end
