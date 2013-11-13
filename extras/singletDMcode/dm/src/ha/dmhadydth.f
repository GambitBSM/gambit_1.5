*****************************************************************************
*** function dshadydth is the differential yield dyield/dcostheta in the
*** cm system boosted to the lab system (including proper jacobians if
*** we are dealing with a differential yield.
*** the function should be integrated from -1 to 1, by e.g.
*** the routine gadap.
*** units: (annihilation)**-1
*****************************************************************************

      real*8 function dmhadydth(cth)
      implicit none
      include 'dshacom.h'

c------------------------ variables ------------------------------------

      real*8 cth,yield,jac,elab,yield1,yield2
      real*8 e1cm,ethcm,dtldt0,tcm(2),dtdt(2)
      integer istat

c------------------------ functions ------------------------------------
      real*8 dmhayieldf
c-----------------------------------------------------------------------

c...calculate the energy of particle 1 in the cm system
      e1cm=phie1  ! is this reasonable even for t b-bar ?

c... calculate the cm energy for this angle that corresponds to requested
c... lab energy phieth (m_mu, m_ep etc neglected, but not m_p)
      dtldt0=0.d0
      if (phitype.eq.0) then ! mass of particle neglected
        dtldt0=abs(phigap*(1.0d0+phibep*cth))
        ethcm=phieth/dtldt0
      elseif (phitype.eq.1) then ! mass not neglected and elab<gamma m
        if (cth.le.phicthmax) then
          elab=phieth+phimp
          tcm(1)=(elab-
     &        phibep*cth*sqrt(max(elab**2-phigap**2*phimp**2+phigap**2
     &        *phibep**2*cth**2*phimp**2,0.0d0)))/
     &        (phigap*(1.0d0-phibep*cth)*(1.0d0+phibep*cth))
     &        -phimp    ! kinetic energy in parent cm system
          dtdt(1)=abs(phigap*(1.0d0+phibep*cth*(tcm(1)+phimp)/
     &      sqrt(tcm(1)**2+2.0d0*tcm(1)*phimp)))

          tcm(2)=(elab+
     &        phibep*cth*sqrt(max(elab**2-phigap**2*phimp**2+phigap**2
     &        *phibep**2*cth**2*phimp**2,0.0d0)))/
     &        (phigap*(1.0d0-phibep*cth)*(1.0d0+phibep*cth))
     &        -phimp    ! kinetic energy in parent cm system
          dtdt(2)=abs(phigap*(1.0d0+phibep*cth*(tcm(2)+phimp)/
     &      sqrt(tcm(2)**2+2.0d0*tcm(2)*phimp)))
c          write(*,*) 'tcm(1) = ',tcm(1),'  dtdt(1) = ',dtdt(1)
c          write(*,*) 'tcm(2) = ',tcm(2),'  dtdt(2) = ',dtdt(2)
        else
          dmhadydth=0.0d0
c          write(*,*) 'cth = ',cth,' > phicthmax'
          return
        endif
      else  ! phitype.eq.2  ! mass not neglected, elab>gamma m
        elab=phieth+phimp
        tcm(1)=(elab-
     &        phibep*cth*sqrt(max(elab**2-phigap**2*phimp**2+phigap**2
     &        *phibep**2*cth**2*phimp**2,0.0d0)))/
     &        (phigap*(1.0d0-phibep*cth)*(1.0d0+phibep*cth))
     &        -phimp    ! kinetic energy in parent cm system
        dtdt(1)=abs(phigap*(1.0d0+phibep*cth*(tcm(1)+phimp)/
     &    sqrt(tcm(1)**2+2.0d0*tcm(1)*phimp)))
      endif

c... get the yield at this energy
      if (phitype.eq.0) then
        yield=dmhayieldf(e1cm,ethcm,phich,phifk,istat)
        if (phifk.ge.100) then  ! differential yields
          jac=1.0d0/dtldt0
          yield=yield*jac
        endif
      elseif (phitype.eq.1) then
        yield1=dmhayieldf(e1cm,tcm(1),phich,phifk,istat)
        yield2=dmhayieldf(e1cm,tcm(2),phich,phifk,istat)
        if (phifk.ge.100) then  ! differential yields
          jac=1.0d0/dtdt(1)
          yield1=yield1*jac
          jac=1.0d0/dtdt(2)
          yield2=yield2*jac
        endif
        yield=yield1+yield2
      else  ! phytype.eq.2
        yield=dmhayieldf(e1cm,tcm(1),phich,phifk,istat)
        if (phifk.ge.100) then  ! differential yields
          jac=1.0d0/dtdt(1)
          yield=yield*jac
        endif
      endif

c... take care of factors depending on energy in neutrino-nucleus cross
c... section and muon range
      if (phifk.eq.72.or.phifk.eq.172) then
        yield=yield*dtldt0   ! nu-nucleon cross section
      elseif (phifk.eq.73.or.phifk.eq.173) then
        yield=yield*(dtldt0**2)  ! nu-nucleon cross section and mu range
      endif

c      if (phifk.eq.72.or.phifk.eq.73) then
c        write(*,*) 'warning in dmhadydth: integrated yields for',
c     &    ' channel ',phifk,' has to be checked!'
c      endif

c...  include a factor of 1/2 to get average of yield,
c...  (1/(int_-1^1 1 dcth)=1/2)
c...  another factor of 1/2 from that we only consider
c...  one particle.

      dmhadydth=0.25d0*yield
c      write(*,*) 'cth = ',cth,'  dmhadydth = ',0.25d0*yield

      end











