****************************************************************
*** axisymmetric dark matter halo density profile            ***
***                                                          ***
*** Input:  radcoord = radial coordinate in kpc              ***
***         vertcoord = vertical coordinate in kpc           ***
***         in a cylindrical coordinate system centered in   *** 
***         the Galactic Center                              ***
*** Output: dshmaxirho = density in GeV/cm^3                 ***
***   e.g.: local halo density = rho0 = dshmaxirho(r_0,0.d0) ***
***                                                          ***
*** Author: Piero Ullio (ullio@sissa.it)                     ***      
*** Date: 2004-01-12                                         ***
****************************************************************

      real*8 function dshmaxirho(radcoord,vertcoord)
      implicit none
      include 'dshmcom.h'
      real*8 radcoord,vertcoord,radialdist,dshmsphrho,
     &   dshmboerrhoaxi
c
      if(haloshape.eq.'spherical') then
c spherical halo density profile
        radialdist=dsqrt(radcoord**2+vertcoord**2)
        dshmaxirho=dshmsphrho(radialdist)
        return
      elseif (haloshape.eq.'boer') then
        dshmaxirho=dshmboerrhoaxi(radcoord,vertcoord,2)
      else
c no other option set yet  
        write(*,*) 'dshmaxirho called wrong option haloshape'
        write(*,*) 'haloshape = ',haloshape
        write(*,*) 'program stopped'
        stop
      endif
      return
      end









