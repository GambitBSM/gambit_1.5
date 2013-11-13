****************************************************************
*** function integrated in dshmjavegc                        ***
*** integration in cos(phi) performed in here, the result    ***
*** still needs to be multiplied by 2 pi                     ***
***                                                          ***
*** author: piero ullio (ullio@sissa.it)                     ***
*** date: 01-10-10                                           ***
****************************************************************

      real*8 function dshmjavepar5(rrin)
      implicit none
      include 'dshmcom.h'
      real*8 rrin,val,rr,dshmsphrho,dshmaxiprob
      rr=rrin
      if (hclumpy.eq.1) then ! smooth halo
        val=dshmsphrho(rr)/rho0 
        dshmjavepar5=val**2
      elseif (hclumpy.eq.2) then ! clumpy halo
c just the spherical probability, you are not allowed to use here
c an axisymmetric one
        val=dshmaxiprob(rr,0.d0)/rho0 
        dshmjavepar5=val
      else
        write(*,*) 'wrong value for hclumpy in dshmjavepar5'
        write(*,*) 'hclumpy =',hclumpy
        stop
      endif
      if (rr.lt.r_0/2.d0) then
        dshmjavepar5=dshmjavepar5*rr/r_0*
     &    (dlog((1.d0+rr/r_0)/(1.d0-rr/r_0)))
      elseif (rr.ge.r_0/2.d0.and.rr.lt.r_0) then
        dshmjavepar5=dshmjavepar5*rr/r_0*
     &    (-dlog((1.d0-rr/r_0)/(1.d0+rr/r_0)))
      else
        write(*,*) 'dshmjavepar5 called with radius larger than r_0'
        write(*,*) 'radius, r_0 =',rr,r_0
        stop
      endif  
      return
      end
