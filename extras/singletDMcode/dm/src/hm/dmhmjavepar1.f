****************************************************************
*** function integrated in dshmjave                          ***
***                                                          ***
*** author: piero ullio (piero@tapir.caltech.edu)            ***
*** date: 00-07-13                                           ***
****************************************************************

      real*8 function dmhmjavepar1(rrtmp)

      implicit none

      real*8 pi,rrtmp,wr,psi0r,linf
      real*8 cosinf,cossup,llsup
      integer ier,iord,last,limit,neval
      real*8 abserr,alist,blist,elist,epsrel,result,rlist
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      real*8 cospsi0,delta,rr
      common/dshmjavecom/cospsi0,delta,rr
      external dmhmjavepar2
      real*8 epsabs
      integer hmjaveset2,ndigits
      common/hmjaveset2com/epsabs,hmjaveset2,ndigits

      if(hmjaveset2.ne.123456) then
         ndigits=5
         epsabs=1.d-10          !numerical accuracy
         hmjaveset2=123456
      endif
      limit=5000
      epsrel=10.d0**(-ndigits)
      rr=rrtmp
      pi=4.d0*datan(1.d0)
      wr=dacos((2.d0*pi-delta)/(2.d0*pi))  !angle of the cone in rad
      psi0r=dacos(cospsi0)
      linf=psi0r-wr
      llsup=psi0r+wr
      if(linf.lt.0.d0) then
      linf=0.d0
      endif
      if(llsup.gt.pi) then
      llsup=pi
      endif
      cossup=dcos(linf)
      cosinf=dcos(llsup)
 10   call dqagseb(dmhmjavepar2,cosinf,cossup,epsabs,epsrel,limit
     &  ,result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      if(dabs(result)/10.d0**ndigits.lt.epsabs) then
        epsabs=dabs(result)/10.d0**(ndigits+2) !too low accuracy
        goto 10
      elseif(dabs(result)/10.d0**(ndigits+2).gt.epsabs) then
        epsabs=dabs(result)/10.d0**(ndigits+2)  !too high accuracy 
      endif
      dmhmjavepar1=result
      return
      end
