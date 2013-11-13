****************************************************************
*** A symmetrized (avergaed over phi) version of de Boers profile
*** as reported in astro-ph/0508617.
*** The profile has a triaxial smooth halo and two
*** rings of dark matter.
***
*** Input: r - cylindrical radius from galactic center (kpc)
***        z - height above galactic plane (kpc)
***        how - 1: returns the average <rho> over phi
***              2: returns the average sqrt(<rho^2>) suitable
***                 for symmetrization for annihilation rates (like pbar)
*** The paramters of the profile are set in dshmset.
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date:   2005-12-08
****************************************************************

      real*8 function dshmboerrhoaxi(r,z,how)
      implicit none
      real*8 r,z,x,y,rhotmp,phi
      real*8 dshmboerrho
      integer i,n,how


      include 'dshmcom.h'
      include 'dssusy.h' 

c...Now take the average over phi
c...Since both the triaxial halo and rings have elliptic symmetry, it is
c...enough to take the average over only one half of a full turn
c...(in principle only a quarter is enough, but we would then need
c...to do it from max to min of each elliptic contribution, which is
c...not possible here).
      n=10
      rhotmp=0.0d0
      do i=1,n
        phi=pi*(dble(i)-0.5d0)/dble(n)
        x=r*cos(phi)
        y=r*sin(phi)
        if (how.eq.1) then
          rhotmp=rhotmp+dshmboerrho(x,y,z)
        else
          rhotmp=rhotmp+dshmboerrho(x,y,z)**2
        endif
      enddo

      if (how.eq.1) then
        dshmboerrhoaxi=rhotmp/dble(n)
      else
        dshmboerrhoaxi=sqrt(rhotmp/dble(n))
      endif

      return
      end









