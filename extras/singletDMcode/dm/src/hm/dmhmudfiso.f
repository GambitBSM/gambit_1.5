****************************************************************
***                                                          ***
*** function which gives u*DF(u) where:                      ***
***                                                          ***
***    u is the modulus of \vec{u} = \vec{v} - \vec{v}_{ob}  ***
***    with \vec{v} the 3-d velocity of a WIMP in the        ***
***    galactic frame, and \vec{v}_{ob} the projection on    ***
***    the frame you are considering                         ***
***                                                          ***
***    DF(u) = int dOmega DF(\vec{v}), where                 ***
***    DF(\vec{v}) is the halo local velocity distribution   ***
***    function in the galactic frame                        ***
***                                                          ***
*** the function implemented here is valid for:              ***
***    a) an isothermal sphere profile                       ***
***    b) an isotropic profile, i.e.                         ***
***       DF(\vec{v}) = DF(|\vec{v}|)                        ***
*** condition b) implies that the integral is performed by   ***
*** setting |\vec{v}|^2 = u^2 + |\vec{v}_ob|^2               ***
***                      + 2*cos(alpha)*|\vec{v}_ob|*u       ***
*** and then integrating in d(cos(alpha))                    ***
***                                                          ***
*** u in km s**-1                                            ***
*** u*DF(u) in km**-2 s**2                                   ***
***                                                          ***
*** Author: Piero Ullio                                      ***
*** Date: 2004-01-29                                         ***
****************************************************************

      real*8 function dshmuDFiso(u)
      implicit none
      include 'dshmcom.h'
      real*8 pi,v0,x,y,u,exp1,exp2,par
ccc
      pi=4.d0*datan(1.d0)
      v0=dsqrt(2.d0/3.d0)*vd_3d
      x=u/v0
      y=v_obs/v0
      par=1.d0/v0**3/pi**1.5d0
      if(x*y.lt.1.d-8) then
        par=par*0.5d0*dexp(-x**2-y**2)*
     &   (1.d0+(2.d0*x*y)**2/6.d0+(2.d0*x*y)**4/120.d0)
      else
        exp1=-(x-y)**2
        exp2=-(x+y)**2
        par=par*(dexp(exp1)-dexp(exp2))/(4.d0*x*y)
      endif  
      par=par*4.d0*pi*u ! km^-2 s^2
      dshmuDFiso=par
      return
      end

