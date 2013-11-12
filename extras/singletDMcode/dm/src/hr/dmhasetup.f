*****************************************************************************
***   subroutine dshasetup prepares the common blocks with annihilation
***   channel branching ratios and Higgs decay widths for halo yield
***   calculations.
***   Note: This routine needs to be called for each each model before
***   dshaloyield is called.
***   This routine is the interface between SUSY and the halo annihilation
***   routines.
***
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 08-01-15
*****************************************************************************

      subroutine dmhasetup
      implicit none
      include 'dssusy.h'
      include 'dshacom.h'
      include 'dshrcom.h'
      include 'dsidtag.h'
      include 'dsprep.h'
      include 'dsibcom.h'

      integer i

c----------------------------------------------- set-up common variables


c...Transfer to Halo annihilation common blocks
c...For a description of the channels, see dshayield.f
      do i=1,30
         habr(i)=sigv(i)/sigmav
      enddo

      
c...Internal Bremsstrahlung
c...See dshaib for more info on viable options.
      haib='none'


c...Mass of WIMP
      hamwimp = mx

      end



















