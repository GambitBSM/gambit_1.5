***********************************************************************
*** Routine dsIBselect selects which channels to include based on the
*** mass degeneracies. Used when ibhow=2.
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2007-10-19
***********************************************************************

      subroutine dmibselect

      implicit none
      include 'dssusy.h'
      include 'dsprep.h'
      include 'dsibcom.h'

      real*8 mn0,zg1,dmabsq
      integer i

      mn0=mx
!      zg1=dmabsq(neunmx(1,1))+dmabsq(neunmx(1,2))

c...Reset all channels
      do i=1,12
            IBflag(i)=0
      enddo 

c...Hard to make 100% safe criterion, uncomment with care
c      if (zg1.gt.0.9d0.and.m0.gt.100) then
c        IBflag(1)=1 ! W+W-
c      endif
      IBflag(1)=1 ! W+W- ! always on for safety

c...Hard to make 100% safe criterion, uncomment with care
c      if (zg1.gt.0.0001d0.and.mn0.gt.800.0d0) then
c        IBflag(2)=1 ! W+H- and W-H+ ! always on for safety
c      endif
c     IBflag(2)=1 ! W+H- and W-H+ ! always on for safety

c...Never important, don't include
      IBflag(3)=0 ! H+H-

c      if (min(mass(kse1),mass(kse2)).lt.ibmfr*mn0) then
c         IBflag(4)=1 ! e+ e-
c      endif     

c      if (min(mass(ksmu1),mass(ksmu2)).lt.ibmfr*mn0) then
c         IBflag(5)=1 ! mu+ mu-
c      endif     

c      if (min(mass(kstau1),mass(kstau2)).lt.ibmfr*mn0) then
c         IBflag(6)=1 ! tau+ tau-
c      endif     

c      if (min(mass(ksu1),mass(ksu2)).lt.ibmfr*mn0) then
c         IBflag(7)=1 ! u u-bar
c      endif     

c      if (min(mass(ksd1),mass(ksd2)).lt.ibmfr*mn0) then
c         IBflag(8)=1 ! d d-bar
c      endif     

c      if (min(mass(ksc1),mass(ksc2)).lt.ibmfr*mn0) then
c         IBflag(9)=1 ! c c-bar
c      endif     

c      if (min(mass(kss1),mass(kss2)).lt.ibmfr*mn0) then
c         IBflag(10)=1 ! s s-bar
c      endif     

c      if (min(mass(kst1),mass(kst2)).lt.ibmfr*mn0) then
c         IBflag(11)=1 ! t t-bar
c      endif     

c      if (min(mass(ksb1),mass(ksb2)).lt.ibmfr*mn0) then
c         IBflag(12)=1 ! b b-bar
c      endif     


      return

      end





        

