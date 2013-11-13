****************************************************************
*** axisymmetric probability distribution function for       ***
*** equal and small (i.e. unresolved) dark matter clumps     ***
***                                                          ***
*** Input:  radcoord = radial coordinate in kpc              ***
***         vertcoord = vertical coordinate in kpc           ***
***         in a cylindrical coordinate system centered in   *** 
***         the Galactic Center                              ***
*** Output: dshmaxiprob in GeV/cm^3, i.e. for the moment     ***
***         the normalization has to be such that            ***
***         dshmaxiprob(r_0,0.d0)=local halo density=rho0    ***
***                                                          ***
*** Author: Piero Ullio (ullio@sissa.it)                     ***      
*** Date: 2004-01-12                                         ***
****************************************************************

      real*8 function dmhmaxiprob(radcoord,vertcoord)

      implicit none
      include 'dshmcom.h'
      real*8 radcoord,vertcoord,radialdist,dmhmsphrho
c
      if(probshape.eq.'spherhalo') then
c spherical probability following the halo density profile
        radialdist=dsqrt(radcoord**2+vertcoord**2)
        dmhmaxiprob=dmhmsphrho(radialdist)
        return
      else
c no other option set yet  
        write(*,*) 'dmhmaxirho called wrong option probshape'
        write(*,*) 'probshape = ',probshape
        write(*,*) 'program stopped'
        stop
      endif
      return
      end









