*****************************************************************************
*** auxiliary routine called by dsIBwhdxdy
*** author: Torsten Bringmann, 2007-04-18
*****************************************************************************

      real*8 function dmIBwhdxdy_9(x,y,m0,mw,mhc,mc1,mc2)
      real*8 x,y,m0,mw,mhc,mc1,mc2

      dmIBwhdxdy_9 = 
     -   (8*m0**3*mc2*(mw**4*(1 - 4*x) + 
     -      4*m0**2*(mhc**2*x + 4*m0**2*(-1 + x)*(x - y))*y - 
     -      mw**2*(mhc**2*x + 4*m0**2*(-x + 2*x**2 + 2*y - 5*x*y))))/
     -  (mw**2*(-2*mc1**2 + mhc**2 + 2*mw**2 + m0**2*(-2 + 4*x-4*y))*
     -    (-2*mc2**2 + mhc**2 + 2*mw**2 + m0**2*(-2 + 4*x - 4*y))*
     -    (-2*mc1**2 + mhc**2 + m0**2*(-2 + 4*y))*
     -    (-2*mc2**2 + mhc**2 + m0**2*(-2 + 4*y)))
      return
      end   ! dmIBwhdxdy_9
