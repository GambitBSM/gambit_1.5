****************************************************************
*** function integrated in dshmjavepar3                      ***
***                                                          ***
*** author: piero ullio (piero@tapir.caltech.edu)            ***
*** date: 00-07-13                                           ***
****************************************************************

      real*8 function dshmjavepar4(rrin)
      implicit none
      include 'dshmcom.h'
      real*8 y,val,rrin,dshmsphrho,dshmaxiprob
      real*8 cospsimin,cospsi,rr
      common/dshmjavegccom/cospsimin,cospsi,rr
      rr=rrin
      y=dsqrt(rr**2+r_0**2-2.d0*rr*r_0*cospsi)
      if (hclumpy.eq.1) then ! smooth halo
        val=dshmsphrho(y)/rho0 
        dshmjavepar4=val**2
      elseif (hclumpy.eq.2) then ! clumpy halo
c just the spherical probability, you are not allowed to use here
c an axisymmetric one
        val=dshmaxiprob(y,0.d0)/rho0 
        dshmjavepar4=val
      else
        write(*,*) 'wrong value for hclumpy in dshmjavepar4'
        write(*,*) 'hclumpy =',hclumpy
        stop
      endif
      return
      end
