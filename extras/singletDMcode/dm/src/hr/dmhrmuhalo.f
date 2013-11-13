**********************************************************************
*** function dshrmuhalo calculates the flux of diffuse neutrino-
*** induced muons from neutralino annihilation in the halo.
*** the flux given, is the total flux above a given threshold.
*** there are some approximations going on for the higgses assuming
*** de/dx for muons are constant to simplify the integration. the
*** errors for this shouldn't be too big.
*** units: km^-2 yr^-1 sr^-1 (if jpsi is given as 1st arg.)
***        km^-2 yr^-1       (if jpsi*delta is given as 1st arg.)
*** dnsigma is also returned, which is the dimensionless number
***   N_mu * (sigma v) / (10^-29 cm^3 s^-1)
*** author: joakim edsjo, edsjo@physto.se
*** date: 98-05-07
*** modified: 00-09-03
***           Joakim Edsjo (edsjo@physto.se) 03-01-21, factor of 1/2
***           in annihilation rate added
**********************************************************************

      real*8 function dshrmuhalo(jpsi,eth,dnsigma,istat)
      implicit none

      include 'dssusy.h'
      include 'dshmcom.h'
      include 'dsidtag.h'
      include 'dshacom.h'
      include 'dsprep.h'

c---------------------------------------------------------------------

      real*8 phimu,phimutot,dshaloyield,eth,
     &  jpsi,dnsigma,dssigmav
      integer istat

c---------------------------------------------------------------------

      phimu=dshaloyield(eth,73,istat)

      dnsigma=phimu*(dssigmav(0)/1.0d-29)*0.5d0 ! JE Corr 03-01-21
      phimutot=dnsigma*jpsi*1.0d0*1.879d-11*
     &  (10.0d0/mx)**2
c...this give the flux per sr. ! convert to km^-2 yr^-1
      phimutot=phimutot*3.15d17       ! km^-2 yr^-1 sr^-1

      dshrmuhalo=phimutot*(rhox/rho0)**2

      end

