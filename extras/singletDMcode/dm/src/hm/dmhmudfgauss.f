***********************************************************************
*** The halo velocity profile in the Maxwell-Boltzmann (Gaussian)
*** approximation.
*** input: velocity relative to earth [ km s^-1 ]
*** output: f(u) / u [ (km/ s)^(-2) ]
*** date: april 6, 1999
*** Modified: 2004-01-29
***********************************************************************

      real*8 function dshmuDFgauss(u)
      implicit none
      include 'dshmcom.h'
      real*8 u

c...calculate f(u)/u
      dshmuDFgauss=sqrt(3./2.)/(vd_3d*v_sun)/sqrt(3.1415)*(
     &  exp(-1.5*(u-v_sun)**2/vd_3d**2)
     &  -exp(-1.5*(u+v_sun)**2/vd_3d**2))

      return
      end
