*****************************************************************************
*** function dshayieldfth integrates dshadydth over cos theta.
*** it is the yield from particle 1 (which decays from m0)
*** that is calculated. particle one corresponds to channel ch.
*** units: (annihilation)**-1
*****************************************************************************

      real*8 function dmhayieldfth(e0,m0,mp1,mp2,emuthr,ch,yieldk,
     &  istat)
      implicit none
      include 'dshacom.h'
      include 'dsidtag.h'
      include 'dsio.h'

c------------------------ variables ------------------------------------

      real*8 e0,m0,mp1,mp2,emuthr,m00,e00
      integer ch,istat,yieldk

      real*8 sum,pi,dth,elab
      parameter (pi=3.141592653589793238d0,dth=0.001d0)
      logical wb,chok

c------------------------ functions ------------------------------------

      external dmhadydth
      real*8 dmhadydth,dmhayield_int

c-----------------------------------------------------------------------

      dmhayieldfth=0.d0
      wb=.true.
      m00=m0
      e00=e0
c...take care of slight numerical inaccuracy problems
      if (m00.lt.(mp1+mp2)) then
        if (m00.gt.0.99*(mp1+mp2)) then
          m00=(mp1+mp2)*1.0001
          if(e0.lt.m00) e00=m00
        else
          write(*,*) 'error in dmhayieldfth: a particle with mass ',m0
          write(*,*) 'is let to decay to two particles with mass ',
     &      mp1,' and ',mp2
          write(*,*) 'which is not energetically allowed.'
          write(*,*) 'particle 1 has code ',ch
          write(*,*) 'the yield from this decay is put to 0.0'
          write(*,*) 'model: ',idtag
          dmhayieldfth=0.0
          return
        endif
      endif

c...calculate boost parameters of parent for common block
      phibep=sqrt(e00**2-m00**2)/e00
      phigap=e00/m00

c...set common block variables
      phim0=m00
      phie0=e00
      phim1=mp1
      phim2=mp1
      phie1=(m00**2+mp1**2-mp2**2)/(2.0d0*m00) ! cm energy of 1
      phie2=(m00**2+mp2**2-mp1**2)/(2.0d0*m00) ! cm energy of 2
      phieth=emuthr
      phich=ch
      phifk=yieldk

      if (yieldk.eq.54.or.yieldk.eq.154) then ! antiproton
        phimp=0.93827231
        elab=phieth+phimp

        if (elab.le.phigap*phimp) then
          phitype=1
          phicthmax=-sqrt(max((phigap**2*phimp**2-elab**2)/
     &      (phimp**2*phigap**2*phibep**2),0.0d0))
          phicthmax=max(-1.0d0,phicthmax-dth)
          phicthmin=max((phieth-phie1*phigap)/
     &      (phigap*phibep*sqrt(phie1**2-phimp**2)),-1.0d0)
        else
          phitype=2
          phicthmax=1.0d0
          phicthmin=max((phieth-phie1*phigap)/
     &      (phigap*phibep*sqrt(phie1**2-phimp**2)),-1.0d0)
        endif
      else  ! all others where mass can be neglected
        phimp=0.0d0
        phitype=0
        phicthmax=1.0d0
        phicthmin=max((phieth-phie1*phigap)/
     &    (phie1*phigap*phibep),-1.0d0)
      endif

      if (phicthmin.ne.-1.0d0) then
        phicthmin=min(phicthmin+dth,1.0d0)
      endif

      if (yieldk.le.100) then ! integrated yields ok to go all the way
        phicthmax=1.0d0
      endif

      if (phicthmin.ge.phicthmax) then
        dmhayieldfth=0.0d0
        return
      endif


c...check if mneu within correct bounds
      chok = .true.
      if ((ch.eq.1.or.ch.eq.2.or.ch.eq.4).and.phie1.lt.lb(ch))
     &  chok=.false.
      if ((ch.eq.3.or.ch.eq.5.or.ch.eq.6).and.phie1.lt.(0.99*lb(ch)))
     &  chok=.false.

      if (.not.chok) then
        if (prtlevel.gt.0) then
          write(6,*)
          write(6,5000)
     +     'warning in dmhayieldfth: a cm energy of ',
     +      2.0d0*phie1,' gev wants to be used,'
          write(6,5010) 'while the lower bound for channel ',ch,
     +      ' is ',2.0d0*lb(ch),' gev.'
          write(6,*) 'the yield is put to 0 for these too low energies.'
          write(6,*) 'the results can thus only be trusted as a',
     +      ' lower bound.'
        endif
        wb=.false.
        dmhayieldfth=0.0d0
        istat=(istat/2)*2+1
      endif

      if (phie1.gt.ub(ch)) then
        if (prtlevel.gt.0) then
          write(6,*)
          write(6,5000) 'warning in dmhayieldfth: a cm energy of ',
     +      2.0d0*phie1,' gev wants to be used,'
          write(6,5010) 'while the upper bound for channel ',ch,
     +      ' is ',2.0d0*ub(ch),' gev.'
          write(6,5020) 'a cm energy of ',2.0d0*ub(ch),' gev is used',
     +      ' instead for these too high energies.'
          write(6,*) 'the results can thus only be trusted as a',
     +      ' lower bound.'
          endif
        istat=(istat/2)*2+1
      endif

      if (wb) then
        if (phicthmin.lt.0.99d0.and.phicthmax.gt.0.99d0) then
          if (dmhadydth(0.99d0).eq.0.0d0) then
            sum=dmhayield_int(dmhadydth,0.99d0,phicthmax)+
     &        dmhayield_int(dmhadydth,phicthmin,0.99d0)
          else
            sum=dmhayield_int(dmhadydth,phicthmin,phicthmax)
          endif
        else
          sum=dmhayield_int(dmhadydth,phicthmin,phicthmax)
        endif

c        eps=0.01
c        sum=0.0
c        call gadap(0.0,thu,dmhadydth,eps,sum)
c
c        if (sum.lt.1.0.and.sum.gt.0.0) then
c          write(6,*) 'sum less then 1.0 in gadap (phiith)'
c          eps=0.01*sum
c          eps0=eps
c  100     call gadap(0.0,thu,dmhadydth,eps,sum)
c          if ((eps0/sum).gt.0.01) then
c            eps0=eps0/2.0
c            eps=eps0
c          endif
c          if ((eps0*2/sum).gt.0.01) goto 100
c        endif

        dmhayieldfth=sum  ! JE Corr 080114: 1.d-15 removed

      endif

 5000 format(' ',a,f8.2,a)
 5010 format(' ',a,i2,a,f8.2,a)
 5020 format(' ',a,f8.2,a,a)

      end
