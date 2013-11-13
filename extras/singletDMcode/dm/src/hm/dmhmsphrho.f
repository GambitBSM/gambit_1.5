****************************************************************
*** spherically symmetric dark matter halo density profile   ***
***                                                          ***
*** Input:  radialdist = radial coordinate in kpc            ***
***         in a spherically symmetric coordinate system     ***
***         centered in the Galactic Center                  ***
***         if radialdist lower than the cut radius rhcut,   ***
***         radialdist is shifted to rhcut                   ***
*** Output: dshmsphrho = density in GeV/cm^3                 ***
***   e.g.: local halo density = rho0 = dshmsphrho(r_0)      ***
***                                                          ***
*** Author: Piero Ullio (ullio@sissa.it)                     ***      
*** Date: 2004-01-12                                         ***
****************************************************************

      real*8 function dmhmsphrho(radialdist)

      implicit none
      include 'dshmcom.h'

      real*8 radialdist,rr,dmhmabgrho,dmhmburrho,dmhmn03rho
     &  ,dmhmnumrho
c
      if(radialdist.lt.rhcut) then
        rr=rhcut
      else
        rr=radialdist
      endif  
c
      if(halotype.eq.'albega') then
c alpha-beta-gamma parametrization
        dmhmsphrho=dmhmabgrho(rr)
      elseif(halotype.eq.'burkert') then
c burkert parametrization
        dmhmsphrho=dmhmburrho(rr) 
      elseif(halotype.eq.'n03') then
c navarro et al.(2003) parametrization
        dmhmsphrho=dmhmn03rho(rr) 
      elseif(halotype.eq.'numerical') then
c profile read from disc
        dmhmsphrho=dmhmnumrho(rr) 
      else
c no other option set yet  
        write(*,*) 'dmhmsphrho called wrong option halotype'
        write(*,*) 'halotype = ',halotype
        write(*,*) 'program stopped'
        stop
      endif
      return
      end
