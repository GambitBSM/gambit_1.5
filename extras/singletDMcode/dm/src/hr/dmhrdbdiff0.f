      real*8 function dshrdbdiff0(td,solarmod,how)

**********************************************************************
*** function dshrdbdiff0 calculates the differential flux of
*** antideuterons for the kinetic energy per nucleon td as a result of
*** neutralino annihilation in the halo.
*** dshrdbdiff0 uses the unrescaled local density
*** see dshrdbardiff for rescaling the local density
***   input:
***     td - antideuteron kinetic energy per nucleon in gev
***     solarmod - 0 no solar modulation
***     how - 1 calculate t_diff only for requested momentum
***           2 tabulate t_diff for first call and use table for
***             subsequent calls
***           3 as 2, but also write the table to disk as pbtd.dat
***             at the first call
***           4 read table from disk on first call, and use that for
***             subsequent calls
*** units: gev^-1 cm^-2 s^-1 sr^-1
*** author: joakim edsjo
*** date: 98-02-10
*** modified: joakim edsjo, edsjo@physto.se
*** modified: 98-07-13, 00-07-19 paolo gondolo
**********************************************************************

      implicit none

      include 'dshmcom.h'
      include 'dshacom.h'
      include 'dsprep.h'
      include 'dspbcom.h'

c----------------------------------------------------------------------
      real*8 td,tdb,dqde,mp,pi,vtdiff,md,nucleon,pcoal,
     &  difftime,ed,edb,pd,pdb,dshaloyielddb,phiap,dsdbtd15,dssigmav
      integer istat,solarmod,how
      parameter (mp=0.93827231d0,pi=3.141592653589793238d0)
      parameter (md=1.875612762d0)

      pcoal=0.058d0
      istat=0
c----------------------------------------------------------------------

c...solar modulation: calculate kinetic energy at heliosphere boundary
c...antiproton flux before modulation at kinetic energy tdb
      if (solarmod.eq.1) then  ! solar modulation
        nucleon=2.d0
        ed=nucleon*td+md
        edb=ed+smod
        tdb=td+smod/nucleon
        pd=dsqrt(dabs(ed**2-md**2))
        pdb=dsqrt(dabs(edb**2-md**2))
        phiap=(4.d0/3.d0*pcoal**3/pd)*md/mp**2
     &          *dshaloyielddb(tdb,159,istat)
     &          *pd**2/pdb**2
      else               ! no solar modulation
        nucleon=2.d0
        tdb=td
        edb=nucleon*td+md
        pdb=dsqrt(dabs(edb**2-md**2))
        phiap=(4.d0/3.d0*pcoal**3/pdb)*md/mp**2
     &          *dshaloyielddb(td,159,istat)  ! JE corr 030423
      endif


c...antideuteron flux, factor of 1/2 inserted
      dqde = 0.5d0*(rho0/mx)**2*dssigmav(0)*phiap/(4.0d0*pi)
      difftime = dsdbtd15(tdb,how)
      vtdiff = (pdb/edb*2.99792458d10)*difftime*1.d15
      dshrdbdiff0 = vtdiff*dqde

      return
      end
