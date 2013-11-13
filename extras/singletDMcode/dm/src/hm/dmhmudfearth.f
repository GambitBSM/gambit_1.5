****************************************************************
*** Dark matter halo velocity profile as seen from the       ***
*** Earth. Compared to dshmuDF, this routine also includes   ***
*** the possibility to use distribution functions where      ***
*** solar system diffusion is included.                      ***
***                                                          ***
*** This routine gives back u*DF(u) in units of (km/s)^(-2)  ***
*** DF(u) = int dOmega DF(abs(v)) where DF(abs(v-vector)) is ***
*** the three-dimensional distribution function in the halo  ***
*** and v = v_us + u with u being the velocity relative us   ***
*** (Earth/Sun).                                             ***
*** Note: u*DF(u) is the same as f(u)/u, where f(u) is the   ***
*** one-dimensional distribution function as defined in e.g. ***
*** Gould, ApJ 321 (1987) 571.                               ***
***                                                          ***
*** dshmuDF is normalized such that                          ***
*** int_0^\infty u*dshmuDF du = int_0^\infty u^2 DF(u) du =  ***
*** int_0^\infty f(u) du = 1                                 ***
***                                                          ***
*** Input:  u = Speed in km/s                                ***
*** Output: u*DF(u) in (km/s)^(-2)                           ***
***                                                          ***
*** Calls other routines depending of choice of velocity     ***
*** distribution function (as set by veldfearth in the   ***
*** common blocks in dshmcom.h)                              ***
*** Author: Joakim Edsjo                                     ***
*** Date: 2004-01-29                                         ***
****************************************************************

      real*8 function dshmuDFearth(u)
      implicit none
      include 'dshmcom.h'
      real*8 u,dshmuDFgauss,dshmuDFearthtab,dsntdkfoveru

      if (veldfearth.eq.'gauss'.or.veldfearth.eq.'iso') then
        dshmudfearth=dshmudfgauss(u)
      elseif(veldfearth.eq.'sdbest') then
        dshmudfearth=dshmudfearthtab(u,1)
      elseif(veldfearth.eq.'sdcons') then
        dshmudfearth=dshmudfearthtab(u,2)
      elseif(veldfearth.eq.'sducons') then
        dshmudfearth=dshmudfearthtab(u,3)
      elseif(veldfearth.eq.'sdgauss') then
        dshmudfearth=dshmudfearthtab(u,4)
      elseif(veldfearth.eq.'user') then
        dshmudfearth=dshmudfearthtab(u,5)
      elseif(veldfearth.eq.'dk') then
        dshmudfearth=dsntdkfoveru(u)
      else
        write(*,*) 'ERROR in dshmuDFearth: Unknown velocity type,'
        write(*,*) 'veldfearth = ',veldfearth
        stop
      endif
      
      return
      end









