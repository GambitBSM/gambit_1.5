*****************************************************************************
*** auxiliary routine called by dsIBwhdxdy
*** author: Torsten Bringmann, 2007-04-18
*****************************************************************************

      real*8 function dmIBwhdxdy_7(x,y,m0,mw,mhc,mc1,mc2)
      real*8 x,y,m0,mw,mhc,mc1,mc2

      dmIBwhdxdy_7 = 
     -   (8*(2*m0**2*mw**6 + m0**4*
     -       (mhc**2*mw**2*x + mw**4*(-1 + 8*x - 16*y)) - 
     -      16*m0**8*(-1 + x)*(x - y)*y + 
     -      4*m0**6*(-(mhc**2*x*y) + 
     -         mw**2*(-x + 2*x**2 + 2*y - 9*x*y + 8*y**2))))/
     -  (mw**2*(-2*mc1**2 + mhc**2 + 2*mw**2 + m0**2*(-2 + 4*x-4*y))*
     -    (-2*mc2**2 + mhc**2 + 2*mw**2 + m0**2*(-2 + 4*x - 4*y))*
     -    (-2*mc1**2 + mhc**2 + m0**2*(-2 + 4*y))*
     -    (-2*mc2**2 + mhc**2 + m0**2*(-2 + 4*y)))
      return
      end   ! dmIBwhdxdy_7
