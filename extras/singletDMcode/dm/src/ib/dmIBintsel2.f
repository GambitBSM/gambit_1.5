*****************************************************************************
*** auxiliary function that selects integrand for integration routines
*** author: Torsten Bringmann 2008-03-12
*****************************************************************************

      real*8 function dmIBintsel2(xint) 

      implicit none
      include 'dsibcom.h'

      real*8 xint,xintcut
      real*8 dmIBf_intdy,dmIBf_intdy2,dmhayield
      integer istat,ch

c... intch = 1..12 for IB photons
c... intch = 101..112 for FSR photons
c... intch = -1...-12 for IB positrons

      if ((intch.ge.1.and.intch.le.12).or.
     &    (intch.ge.101.and.intch.le.112)) then
         dmIBintsel2 = 
     &      dmIBf_intdy(intch,xint,ibcom_mx,ibcom_mp1,ibcom_mp2)
      elseif (intch.eq.-4) then
         dmIBintsel2 =
     &      dmIBf_intdy2(4,xint,ibcom_mx,ibcom_mp1,ibcom_mp2)
      elseif ((intch.eq.-1).or.
     &        (intch.le.-5.and.intch.ge.-12)) then
c... set channel for dshayield ! quark final states preliminary!!!
         if (intch.eq.-1) then
              ch=13
           elseif (intch.eq.-5) then
              ch=17
           elseif (intch.eq.-6) then
              ch=19
           elseif ((intch.le.-7).and.(intch.ge.-10)) then
              ch=22
           elseif (intch.eq.-11) then
              ch=24
           elseif (intch.eq.-12) then
              ch=25
           else
              dmIBintsel2=0.d0
              return
         endif
         xintcut=xint                  ! for higher energies, take the IB 
         if (xint.ge.0.8) xintcut=0.8  ! contribution at xp=IRcut as a
                                       ! *conservative* estimate for the yield;
                                       ! a full treatment would have to include 
                                       ! virtual photons
         dmIBintsel2 =
     &        dmIBf_intdy2(-intch,xintcut,ibcom_mx,ibcom_mp1,ibcom_mp2)*
     &          dmhayield(xint*ibcom_mx,ibcom_x*ibcom_mx,
     &                    ch,intyield,istat)
      else
         dmIBintsel2=0.d0
      endif
 
      end



