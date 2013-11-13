*****************************************************************************
*** auxiliary routine called by dsIBffdxdy
*** author: Torsten Bringmann, 2007-07-05
*****************************************************************************

      real*8 function dmIBffdxdy_1(x,y,m0,ml,msl1,msl2)
      real*8 x,y,m0,ml,msl1,msl2

      dmIBffdxdy_1 = 
     -      (-8*m0**2*(ml**4*(-1 + x) + 
     -      8*m0**4*(-1 + x)*(x**2 - 2*x*y + 2*y**2) + 
     -      4*m0**2*ml**2*(2*x**2 + 2*y - x*(1 + 2*y))))/
     -  ((3*ml**2 - 2*msl1**2 + m0**2*(-2 + 4*x - 4*y))*
     -    (3*ml**2 - 2*msl2**2 + m0**2*(-2 + 4*x - 4*y))*
     -    (ml**2 - 2*msl1**2 + m0**2*(-2 + 4*y))*
     -    (ml**2 - 2*msl2**2 + m0**2*(-2 + 4*y)))

      return
      end   ! dmIBffdxdy_1
