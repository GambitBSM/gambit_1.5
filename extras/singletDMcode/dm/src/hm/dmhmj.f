****************************************************************
*** function dshmj: line of sight integral which enters in the ***
*** computation of the gamma-ray and neutrino fluxes from    ***
*** pair annihilations of wimps in the halo.                 ***
***                                                          ***
*** see definition in e.g. bergstrom et al.,                 ***
***   phys. rev d59 (1999) 043506                            ***
*** in case of the many unresolved clump scenario the term   ***
*** fdelta is factorized out                                 *** 
***                                                          ***  
*** it is valid for a spherical dark matter halo             ***
*** psi0 is the angle between direction of observation       ***
*** and the direction of the galactic center; cospsi0 is     ***
*** its cosine.                                              ***
***                                                          ***
*** author: piero ullio (piero@tapir.caltech.edu)            ***
*** date: 00-07-13                                           ***
****************************************************************

      real*8 function dmhmj(cospsi0in)

      implicit none

      include 'dshmcom.h'

      integer ier,iord,last,limit,neval
      real*8 abserr,alist,blist,elist,epsabs,epsrel,result,rlist
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      real*8 cospsi0in,par,sd,logrmin,logrmax,inter,locinf,
     & locsup
      integer k
      real*8 cospsi0
      common/dshmjcom/cospsi0
      external dmhmjpar1

      cospsi0=cospsi0in
      epsabs=1.d-10     !numerical accuracy
      epsrel=1.d-10
      limit=5000
      sd=1.d-5
      par=0.d0
      logrmin=dlog10(sd)
      logrmax=dlog10(r_0)
      inter=(logrmax-logrmin)/6.d0
      par=0.d0
      do k=0,5
        locinf=r_0-10.d0**(logrmax-inter*k)
        locsup=r_0-10.d0**(logrmax-inter*(k+1)) 
        call dqagse(dmhmjpar1,locinf,locsup,epsabs,epsrel,limit,
     &    result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        par=par+result
      enddo
      locinf=r_0-sd
      locsup=r_0+sd 
      call dqagse(dmhmjpar1,locinf,locsup,epsabs,epsrel,limit,
     &  result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      logrmin=dlog10(sd)
      logrmax=dlog10(100.d0)
      par=par+result
      inter=(logrmax-logrmin)/7.d0
      do k=0,6
        locinf=r_0+10.d0**(logrmin+inter*k)
        locsup=r_0+10.d0**(logrmin+inter*(k+1)) 
        call dqagse(dmhmjpar1,locinf,locsup,epsabs,epsrel,limit,
     &    result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        par=par+result
      enddo
      dmhmj=par/8.5d0*(rho0/0.3d0)**2
      return
      end


