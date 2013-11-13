**********************************************************************
*** function dshrgacsusy gives the susy dependent term in the
*** flux of gamma-rays with continuum energy spectrum above the 
*** threshold egath (gev) from neutralino annihilation in the halo.
***
*** dshrgacsusy is dimensionless
*** 
*** the flux in a solid angle delta in the direction psi0 is given by:
***   cm^-2 s^-1 sr^-1 * dshrgacsusy(egath,istat)
***   * dshmjave(cospsi0,delta) * delta (in sr)
***
*** the flux per solid angle in the direction psi0 is given by:
***   cm^-2 s^-1 sr^-1 * dshrgacsusy(egath,istat)
***   * dshmj(cospsi0)
*** 
*** in case of a clumpy halo the factor fdelta has to be added 
***
*** author: joakim edsjo, edsjo@physto.se
*** modified: piero ullio (piero@tapir.caltech.edu) 00-07-13
***           Joakim Edsjo (edsjo@physto.se) 03-01-21, factor of 1/2
***           in annihilation rate added
**********************************************************************

      real*8 function dmhrgacsusy(egath,istat)

      implicit none

      include 'dssusy.h'
      include 'dshmcom.h'
      include 'dshacom.h'
      include 'dsprep.h'

      real*8 phiga,dmhaloyield,egath,dnsigma,dmsigmav
      integer istat
c...pg dshaloyield calls hrsetup which assigns mx in dshacom.h
      phiga=dmhaloyield(egath,52,istat)
      dnsigma=phiga*(dmsigmav(0)/1.0d-29)*0.5d0 ! JE Corr 03-01-21
      dmhrgacsusy=1.879d-11*dnsigma*(10.0d0/mx)**2
      end
