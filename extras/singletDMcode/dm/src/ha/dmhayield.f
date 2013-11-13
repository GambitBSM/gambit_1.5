*****************************************************************************
*** function dshayield calculates the yield above threshold
*** or the differential flux, for the
*** fluxtype given by yieldk, according to the following table.
***
*** particle       integrated yield     differential yield
*** --------       ----------------     ------------------
*** positron                     51                    151
*** cont. gamma                  52                    152
*** nu_mu and nu_mu-bar          53                    153
*** antiproton                   54                    154
*** cont. gamma w/o pi0          55                    155
*** nu_e and nu_e-bar            56                    156
*** nu_tau and nu_tau-bar        57                    157
*** pi0                          58                    158
*** nu_mu and nu_mu-bar          71                    171 (same as 53/153)
*** muons from nu at creation    72                    172
*** muons from nu at detector    73                    173
***
*** The channel numbers should be the new channel numbers as listed below
*** S# is neutral scalar number #
*** S+ is charged scalar
***
*** Ch No  Particles                 Old Ch No   Old chi   New chcomp
*** -----  ---------                 ---------   -------   ----------
***  1     S1 S1                      -          -         Not done yet
***  2     S1 S2                      -          -
***  3     S2 S2                      -          -
***  4     S3 S3                      -          -
***  5     S1 S3                      7          -
***  6     S2 S3                      11         -
***  7     S- S+                      -          -
***  8     S1 Z                       8          -
***  9     S2 Z                       9          -
*** 10     S3 Z	                      -          -
*** 11     W- S+ and W+ S-            10         -
*** 12     Z0 Z0 	              6          6
*** 13     W+ W-                      5          5
*** 14     nu_e nu_e-bar              -          -
*** 15     e+ e-                      -          -
*** 16     nu_mu nu_mu-bar            -          -

*** 17     mu+ mu-                    13         7
*** 18     nu_tau nu_tau-bar	      -          -
*** 19     tau+ tau-	              4          4
*** 20     u u-bar                    -          -
*** 21     d d-bar                    -          -
*** 22     c c-bar                    1          1
*** 23     s s-bar                    -          -
*** 24     t t-bar                    3          3
*** 25     b b-bar                    2          2
*** 26     gluon gluon                12         8
*** 27     q q gluon (not implemented yet, put to zero)
*** 28     gamma gamma (1-loop)
*** 29     Z0 gamma (1-loop)
*** 30     hh
***
*** the units are (annihilation)**-1
*** for the differential yields, the units are the same plus gev**-1.
***
*** Note 1. The correct data files need to be loaded. This is handled by
*** a call to dshainit. It is done automatically here upon first call.
***
*** Note 2. These routines do not contain internal bremsstrahlung (IB)
*** contributions (except those final state radiations (FSR) that are
*** included in the Pythia runs). The full IB contributions are added
*** in dshaloyield.f with calls to dshaib.f. The function that
*** (in the case of neutralinos) corresponds to dshayield is dshaIByieldone.f.
***
*** author: joakim edsjo (edsjo@physto.se)
*** date: 98-01-26
*** modified: 08-01-15
*****************************************************************************

      real*8 function dmhayield(mwimp,emuthr,ch,yieldk,istat,mtop)

      implicit none

      include 'dshacom.h'
      include 'dsidtag.h'
      include 'dsio.h'

c------------------------ variables ------------------------------------

      real*8 mwimp,emuthr,mp1,mp2,e1,e2,yield,mtop
      integer ch,chold,istat,yieldk,fltyp,i,j,fi,chi

      logical chok

      logical first
      data first /.true./
      save first

c------------------------ functions ------------------------------------

      real*8 dmhayieldf,dmhayields,dmhayieldee,dmhayieldhh

c-----------------------------------------------------------------------



      chi=chcomp(ch) ! convert from new to compact channel numbers


      istat=0
      haerr=0
      mp1=0.d0
      mp2=0.d0

      dmhayield=0.0

      call dmhadec(yieldk,fltyp,fi)

c--------------------------------------- if first call, load tables

c      write(*,*) 'dshayield called with:'
c      write(*,*) '  mwimp = ',mwimp
c      write(*,*) '  emuthr = ',emuthr
c      write(*,*) '  ch = ',ch

      if (first) then
        do i=1,ntype
          do j=1,2
            yieldtype(j,i)=0
          enddo
        enddo

        first=.false.
        call dmhainit(yieldk)
      endif

      if (yieldtype(fltyp,fi).eq.0) then
        call dmhainit(yieldk)
      endif

!trick by rruiz
      map(3) = mtop

c-----------------------------------------------------------------------
      
      if (chi.ge.1.and.chi.le.8) then                ! "fundamental" channel
        chok = .false.
        if (mwimp.ge.map(chi)) chok = .true.
        if ((chi.eq.3.or.chi.eq.5.or.chi.eq.6).and.
     &     mwimp.gt.(0.99*map(chi))) chok = .true.
        if (chok) then
          dmhayield=dmhayieldf(mwimp,emuthr,chi,yieldk,istat)
        else
          dmhayield=0.0
c          istat=mod(istat,2)+2
c            write(6,5000) 'error in dmhayield: channel ',ch,' is not',
c     &        ' energetically allowed.'
c            write(6,*) 'Mass of WIMP: ',mwimp
c            write(6,*) 'Mass of each final state particle: ',map(chi)
c            write(6,*) 'Energy: ',emuthr
c            write(6,*) 'model: ',idtag
        endif
        
      elseif (ch .eq. 30 .and. yieldk .eq. 152) then ! hh channel differential yield

        dmhayield = dmhayieldhh(emuthr/mwimp,mwimp)

      elseif (ch .eq. 15 .and. yieldk .eq. 152) then ! e+e- channel differential yield

        dmhayield = dmhayieldee(emuthr/mwimp,mwimp,0.511d-3)
      
      else                                           ! "complex" channel

c...determine masses of the annihilation particles
        if (ch.eq.5) then
c          mp1=has0m(1)         ! S10 mass
c          mp2=has0m(3)         ! S30 mass
c        elseif (ch.eq.8) then
c          mp1=map(6)        ! z0 mass
c          mp2=has0m(1)         ! S10 mass
c        elseif (ch.eq.9) then
c          mp1=map(6)        ! z0 mass
c          mp2=has0m(2)         ! S20 mass
c        elseif (ch.eq.11) then
c          mp1=map(5)        ! w+- mass
c          mp2=hascm         ! S+- mass
c        elseif (ch.eq.6) then
c          mp1=has0m(2)         ! S20 mass
c          mp2=has0m(3)         ! S30 mass
         elseif (ch.eq.29) then
           mp1=map(6)        ! z0 mass
           mp2=0.0d0         ! gamma mass
         endif

c...if energetically allowed channel, go on...
        if (mwimp.ge.0.995d0*((mp1+mp2)/2.0d0)) then

c...calculate the energy of the annihilation particles
          e1=((2.0d0*mwimp)**2-mp2**2+mp1**2)/(4.0d0*mwimp)
          e2=2.0d0*mwimp-e1
          e1=max(e1,mp1+0.001d0)
          e2=max(e2,mp2+0.001d0)

c...check different annihilation channels
          yield=0.0d0

c---------- h10 h30 channel ----------
c          if (ch.eq.5) then
c            yield=yield+dmhayields(e1,emuthr,1,yieldk,istat)
c            yield=yield+dmhayields(e2,emuthr,3,yieldk,istat)
c          endif

c---------- z0 h10 channel ----------
c          if (ch.eq.8) then
c            yield=yield+0.5d0*dmhayieldf(e1,emuthr,6,yieldk,istat)
c            yield=yield+dmhayields(e2,emuthr,1,yieldk,istat)
c          endif

c---------- z0 h20 channel ----------
c          if (ch.eq.9) then
c            yield=yield+0.5d0*dmhayieldf(e1,emuthr,6,yieldk,istat)
c            yield=yield+dmhayields(e2,emuthr,2,yieldk,istat)
c          endif

c---------- w+h- w-h+ channel ----------
c...this calculation gives a mean of the two channels w+h- & w-h+
c          if (ch.eq.11) then
c            yield=yield+0.5d0*dmhayieldf(e1,emuthr,5,yieldk,istat)
c            yield=yield+dmhayields(e2,emuthr,4,yieldk,istat)
c          endif

c---------- h20 h30 channel ----------
c          if (ch.eq.6) then
c            yield=yield+dmhayields(e1,emuthr,2,yieldk,istat)
c            yield=yield+dmhayields(e2,emuthr,3,yieldk,istat)
c          endif

c---------- z0 gamma channel ----------
           if (ch.eq.29) then
             yield=yield+0.5d0*dmhayieldf(e1,emuthr,6,yieldk,istat)
           endif

           dmhayield=yield

        else   ! not energetically allowed channel
          dmhayield=0.0d0
          istat=mod(istat,2)+2
c          if (prtlevel.gt.0) then
c            write(6,5000) 'error in dshayield: channel ',ch,' is not',
c     +        ' energetically allowed.'
c            write(6,*) 'Mass of WIMP: ',mwimp
c            write(6,*) 'Mass of final state particles: ',mp1,mp2
c            write(6,*) 'Energy: ',emuthr
c            write(6,*) 'model: ',idtag
c          endif
        endif

      endif

      if (haerr.gt.0.and.prtlevel.gt.0) then
        write(*,*) 'warning in dshayield for model ',idtag,
     &    ', yield ',yieldk,' and channel ',ch
c        write(*,*) '  the integration over higgs decay angles ran',
c     &    'into numerical problems.'
        write(*,*) '  the results can only be trusted as a lower',
     &    'bound.'
        istat=ibset(istat,2)
      endif

 5000 format(' ',a,i2,a,a)

      end




