*****************************************************************************
*** The function dsIBhhdxdy gives the full analytical expressions for
*** the differential IB photon yield (dxdy) from H+H- final states, 
*** normalized to the  *total* neutralino annihilation rate
***
*** The kinematic variables x,y are
***
*** x = E_gamma/mx
*** y = (p+k)^2/(4 mx^2),
***
*** where p denotes the H+ momentum and k the photon momentum.
*** (note that the expressions above and below only apply to the v->0 limit;
***  in this limit the lowest order result vanishes)
***
*** author: Torsten Bringmann (bringman@sissa.it)
*** date: 2007-05-01
*****************************************************************************

      real*8 function dmIBhhdxdy(IBch,x,y)

      implicit none
      include 'dssusy.h'
      include 'dshacom.h'
      include 'dsprep.h'
      include 'dshmcom.h'
      include 'dsidtag.h'

      integer IBch
      real*8 x,y,tmpresult,dmsigmav,sv

      character*12 memory    ! to suppress multiple error messages
      save memory
      data memory /'____________'/


 
c------------- couplings and masses relevant for IB -------------------

      real*8 CHH21, CHH22
      real*8 mc1, mc2, m0, mhc

      dmIBhhdxdy=0.0d0
      tmpresult=0.0d0

      if (IBch.ne.3) return

c...  this is the total annihilation cross section; the call also makes
c...  sure that all couplings were computed

      sv=dmsigmav(0)

c...set up masses

      m0   = mass(kn(1))
      mhc  = mass(khc)
      mc1  = mass(kcha(1))
      mc2  = mass(kcha(2))
      
c...set up couplings
 
      CHH21 = (abs(gr(khc,kn(1),kcha(1)))**2 
     &         + abs(gl(khc,kn(1),kcha(1)))**2)
      CHH22 = (abs(gr(khc,kn(1),kcha(2)))**2 
     &         + abs(gl(khc,kn(1),kcha(2)))**2)


c...import full IB expression for |M|**2 from form/mathematica:

      tmpresult=
     -  -16*m0**2*(-(mhc**4*(-1 + x)) + 16*m0**4*(-1 + x)*(x - y)*y + 
     -    4*m0**2*mhc**2*(x - 2*y + 2*x*y))*
     -  (CHH21/
     -      ((-2*mc1**2 + 3*mhc**2 + m0**2*(-2 + 4*x - 4*y))*
     -        (-2*mc1**2 + mhc**2 + m0**2*(-2 + 4*y))) + 
     -     CHH22/
     -      ((-2*mc2**2 + 3*mhc**2 + m0**2*(-2 + 4*x - 4*y))*
     -        (-2*mc2**2 + mhc**2 + m0**2*(-2 + 4*y))))**2


c...check that the result is positive
      if (0.gt.tmpresult) then
        if (m0**2*tmpresult.lt.(-1D-15)
     -     .and.(idtag.ne.memory)) then
          write(*,*) '*****'
          write (*,*) 'Error in dmIBhhdxdy: ',
     -                'negative |M|^2 for model ',
     -                 idtag,'.'
          write (*,*) ' Setting corresponding contributions to zero...'
          write(*,*) '*****'
        endif
        return
      endif

c... Since the 2-body annihilation amplitude vanishes for v->0,
c... we now calculate the 3-body cross section and divide it by the 
c... *total* neutralino annihilation cross section:

      tmpresult = gev2cm3s*(alphem/pi)*tmpresult/(32*pi)/sv

      dmIBhhdxdy=tmpresult

      return
      end

