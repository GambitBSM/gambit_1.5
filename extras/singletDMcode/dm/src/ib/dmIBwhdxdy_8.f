*****************************************************************************
*** auxiliary routine called by dsIBwhdxdy
*** author: Torsten Bringmann, 2007-04-18
*****************************************************************************

      real*8 function dmIBwhdxdy_8(x,y,m0,mw,mhc,mc1,mc2)
      real*8 x,y,m0,mw,mhc,mc1,mc2

      dmIBwhdxdy_8 = 
     -   (-8*m0**2*mc1*mc2*(mw**4 + 
     -      4*m0**2*(mhc**2*x + 4*m0**2*(-1 + x)*(x - y))*y - 
     -      mw**2*(mhc**2*x + m0**2*(8*x**2 + 8*y - 4*x*(1 + y)))))/
     -  (mw**2*(-2*mc1**2 + mhc**2 + 2*mw**2 + m0**2*(-2 + 4*x - 4*y))*
     -    (-2*mc2**2 + mhc**2 + 2*mw**2 + m0**2*(-2 + 4*x - 4*y))*
     -    (-2*mc1**2 + mhc**2 + m0**2*(-2 + 4*y))*
     -    (-2*mc2**2 + mhc**2 + m0**2*(-2 + 4*y)))
      return
      end   ! dmIBwhdxdy_8
