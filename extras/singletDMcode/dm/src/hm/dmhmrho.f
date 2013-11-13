****************************************************************
*** dark matter halo density profile                         ***
***                                                          ***
*** Input:  r = galactocentric distance in kpc               ***
*** Output: dshmrho = density in GeV/cm^3                    ***
***                                                          ***
*** links to dshmsphrho where the profile is calculated      ***
*** Author: Joakim Edsjo                                     ***
*** Date: 2000-09-02                                         ***
*** mod: 04-01-13 pu                                         ***
****************************************************************

      real*8 function dshmrho(r)
      implicit none
      real*8 r,dshmsphrho
      dshmrho=dshmsphrho(r)
      return
      end









