**********************************************************************
*** function dssigmav returns the annihilation cross section
*** sigma v at p=0 for neutralino-neutralino annihilation.
*** if partch=0, the full sigma v is obtained and if partch>0, the
*** cross section to channel partch is obtained.
*** units: cm^-3 s^-1
**********************************************************************

      function dmsigmav(partch)
      implicit none
      real*8 dmsigmav
      integer partch
      include 'dssusy.h'
      include 'dsandwcom.h'
      include 'dshacom.h'
      include 'dsidtag.h'
      include 'dsprep.h'
      real*8 dsandwdcosnn
      integer i


      if (newmodelsigmav) then
c         wtot=dsandwdcosnn(0.0d0,0.0d0,kn(1),kn(1))
c         abr(1) = prtial(22)/wtot ! cc-bar
c         abr(2) = prtial(25)/wtot ! bb-bar
c         abr(3) = prtial(24)/wtot ! tt-bar
c         abr(4) = prtial(19)/wtot ! tau+ tau-
c         abr(5) = prtial(13)/wtot ! w+ w-
c         abr(6) = prtial(12)/wtot ! z0 z0
c         abr(7) = prtial(5)/wtot ! h10 h30
c         abr(8) = prtial(8)/wtot ! z0 h10
c         abr(9) = prtial(9)/wtot ! z0 h20
c         abr(10)= prtial(11)/wtot ! w+ h- / w- h+
c         abr(11)= prtial(6)/wtot ! h20 h30
c         abr(12)= prtial(26)/wtot ! gluon gluon
c         abr(13)= prtial(17)/wtot ! mu+ mu-
c         abr(14)= prtial(29)/wtot ! z gamma
         mx=mass(kn(1))
c... sigma v = w / (4*E_1^2) = wtot / (2*mx**2) since integration
c... over cos theta gives factor of 2.
         sigmav = 0.38937966d-27*3.d10*wtot/(2.d0*mx**2) ! in cm^3/s
         do i=1,30
            sigv(i) = prtial(i)/wtot*sigmav
         enddo
         newmodelsigmav=.false.
      endif
      
      dmsigmav=sigmav
      if (partch.gt.0) dmsigmav=sigv(partch)
      return

      end
