****************************************************************
*** function integrated in dshmjavepar1                      ***
***                                                          ***
*** author: piero ullio (piero@tapir.caltech.edu)            ***
*** date: 00-07-13                                           ***
*** mod: 04-01-13 pu                                         ***
****************************************************************

      real*8 function dmhmjavepar2(cospsi)

      implicit none

      include 'dshmcom.h'

      real*8 y,val,pi,cospsi,sinpsi,dmhmsphrho,dmhmaxiprob
      real*8 sinpsi0,cosw,k1,par
      real*8 cospsi0,delta,rr
      common/dshmjavecom/cospsi0,delta,rr

      pi=4.d0*datan(1.d0)
      cosw=((2.d0*pi-delta)/(2.d0*pi))
      sinpsi=dsqrt(1.d0-cospsi**2)
      sinpsi0=dsqrt(1.d0-cospsi0**2)
      if(sinpsi.gt.1.d-16) then
        if(sinpsi0.gt.1.d-16) then
        k1=(cosw-cospsi*cospsi0)/(sinpsi*sinpsi0)
        else
        k1=(cosw-cospsi*cospsi0)/(sinpsi*1.d-16)
        endif
      else
        if(sinpsi0.gt.1.d-16) then
        k1=(cosw-cospsi*cospsi0)/(1.d-16*sinpsi0)
        else
        k1=(cosw-cospsi*cospsi0)/(1.d-16*1.d-16)
        endif
      endif
      if(k1.gt.1.d0) then
        par=0.d0
      else if(k1.le.-1.d0) then
        par=2.d0*pi
      else
        par=2*dacos(k1)
      endif
      if (hclumpy.eq.1) then ! smooth halo
        y=dsqrt(rr**2+r_0**2-2.d0*rr*r_0*cospsi)
        val=dmhmsphrho(y)/rho0 
        dmhmjavepar2=val**2*par
      elseif (hclumpy.eq.2) then ! clumpy halo
        y=dsqrt(rr**2+r_0**2-2.d0*rr*r_0*cospsi)
c just the spherical probability, you are not allowed to use here
c an axisymmetric one
        val=dmhmaxiprob(y,0.d0)/rho0 
        dmhmjavepar2=val*par
      else
        write(*,*) 'wrong value for hclumpy in dshmjavepar2'
        write(*,*) 'hclumpy =',hclumpy
        stop
      endif
      return
      end










