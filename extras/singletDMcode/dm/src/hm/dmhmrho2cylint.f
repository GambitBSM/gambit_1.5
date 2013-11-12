**********************************************************************
*** function which gives the integral of the square of the 
*** axisymmetric density profile dshmaxirho, normalized to the local 
*** halo density, over a cylindrical volume of radius rmax and height
*** 2*zmax, in the frame with the galactic center as its origin, 
***   i.e.:
***        2 \pi * int_{-zmax}^{+zmax} dz 
***              * int_0^{rmax} dr * r (dshmaxirho(r,zint)/rho0)^2
***        = 2 \pi * 2 int_0^{+zmax} dz 
***              * int_0^{rmax} dr * r (dshmaxirho(r,zint)/rho0)^2
***
*** zmax and rmax in kpc 
*** dshmrho2cylint in kpc^3
*** 
*** author: piero ullio (ullio@sissa.it)
*** date: 04-01-22
**********************************************************************

      real*8 function dshmrho2cylint(rmax,zmax)
      implicit none
      real*8 rmax,zmax
      real*8 zint,rmaxint
      common/dshmrho2intcom/zint,rmaxint
      real*8 up,low
      real*8 abserr,alist,blist,elist,epsrel,rlist,result
      integer ier,iord,last,limit,neval
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      external dshmrho2cylint1
      integer ndigits
      integer rho2cy1set
      real*8 epsabs
      common/rho2cy1com/epsabs,rho2cy1set,ndigits
      if(rho2cy1set.ne.123456) then
        ndigits=5 
        epsabs=1.d-14     !numerical accuracy
        rho2cy1set=123456
      endif
      epsrel=10.d0**(-ndigits)
      limit=5000
      rmaxint=rmax
      low=0.d0
      up=zmax
 20   call dqagse(dshmrho2cylint1,low,up,epsabs,epsrel,limit,
     &      result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      if(dabs(result)/10.d0**ndigits.lt.epsabs) then
        epsabs=dabs(result)/10.d0**(ndigits+2) !too low accuracy
        goto 20
      elseif(dabs(result)/10.d0**(ndigits+2).gt.epsabs) then
        epsabs=dabs(result)/10.d0**(ndigits+2)  !too high accuracy 
      endif
      dshmrho2cylint=result*16.d0*datan(1.d0)
      return
      end


      real*8 function dshmrho2cylint1(z)
      implicit none
      real*8 z
      real*8 up,low
      real*8 zint,rmaxint
      common/dshmrho2intcom/zint,rmaxint
      real*8 abserr,alist,blist,elist,epsrel,rlist,result
      integer ier,iord,last,limit,neval
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      external dshmrho2cylint2
      integer ndigits
      integer rho2cy2set
      real*8 epsabs
      common/rho2cy2com/epsabs,rho2cy2set,ndigits
      if(rho2cy2set.ne.123456) then
        ndigits=5 
        epsabs=1.d-14     !numerical accuracy
        rho2cy2set=123456
      endif
      epsrel=10.d0**(-ndigits)
      limit=5000
      zint=z
      low=0.d0
      up=rmaxint
 20   call dqagseb(dshmrho2cylint2,low,up,epsabs,epsrel,limit,
     &      result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      if(dabs(result)/10.d0**ndigits.lt.epsabs) then
        epsabs=dabs(result)/10.d0**(ndigits+2) !too low accuracy
        goto 20
      elseif(dabs(result)/10.d0**(ndigits+2).gt.epsabs) then
        epsabs=dabs(result)/10.d0**(ndigits+2)  !too high accuracy 
      endif
      dshmrho2cylint1=result
      return
      end



      real*8 function dshmrho2cylint2(r)
      implicit none
      include 'dshmcom.h'
      real*8 dshmaxirho,r
      real*8 zint,rmaxint
      common/dshmrho2intcom/zint,rmaxint
      dshmrho2cylint2=(dshmaxirho(r,zint)/rho0)**2*r
      return
      end

