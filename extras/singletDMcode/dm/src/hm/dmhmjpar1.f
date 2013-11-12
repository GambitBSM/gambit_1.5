****************************************************************
*** function integrated in dshmj                             ***
***                                                          ***
*** author: piero ullio (piero@tapir.caltech.edu)            ***
*** date: 00-07-13                                           ***
*** mod: 04-01-13 pu                                         ***
****************************************************************


      real*8 function dmhmjpar1(rr)

      implicit none
      include 'dshmcom.h'

      real*8 rr,y,val,dmhmsphrho,dmhmaxiprob
      real*8 cospsi0
      common/dshmjcom/cospsi0

      y=dsqrt(rr**2+r_0**2-2.d0*rr*r_0*cospsi0)
      if (hclumpy.eq.1) then ! smooth halo
        val=dmhmsphrho(y)/rho0 
        dmhmjpar1=val**2
      elseif (hclumpy.eq.2) then ! clumpy halo
c just the spherical probability, you are not allowed to use here
c an axisymmetric one
        val=dmhmaxiprob(y,0.d0)/rho0 
        dmhmjpar1=val
      else
        write(*,*) 'wrong value for hclumpy in dshmjavepar2'
        write(*,*) 'hclumpy =',hclumpy
        stop
      endif
      return
      end



