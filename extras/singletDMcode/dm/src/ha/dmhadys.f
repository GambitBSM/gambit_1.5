*****************************************************************************
*** function dshadys is the differential yield dyield/dcostheta in the
*** cm system boosted to the lab system (including proper jacobians if
*** we are dealing with a differential yield. all decay channels of the
*** scalar boson in question are summed.
*** the function should be integrated from -1 to 1, by e.g.
*** the routine gadap.
*** units: (annihilation)**-1
*** author: joakim edsjo (edsjo@physto.se)
*** date: 1998
*** modified: 98-04-15
*****************************************************************************

      real*8 function dmhadys(cth)
      implicit none
      include 'dshacom.h'

c------------------------ variables ------------------------------------

      real*8 cth,yield,jac,elab,yield1,yield2,mp,ep
      real*8 e1cm,ethcm,dtldt0,tcm(2),dtdt(2)
      integer istat,hdi(7),ch

      data hdi/22,25,24,19,13,12,17/ ! fundamental channels

c------------------------ constants ------------------------------------

      real*8 tmin,pmin
      parameter(tmin=1.0d-5)   ! to avoid numerical roundoff problems
      parameter(pmin=1.0d-5)   ! to avoid numerical roundoff problems

c------------------------ functions ------------------------------------
      real*8 dmhayieldf
c-----------------------------------------------------------------------

c...calculate the energy of particle 1 in the cm system
      e1cm=has0m(phihno)/2.0d0

c... calculate the cm energy for this angle that corresponds to requested
c... lab energy phieth (m_mu, m_ep etc neglected, but not m_p)
      dtldt0=0.d0
      if (phitype.eq.0) then ! mass of particle neglected
        dtldt0=abs(phigap*(1.0d0+phibep*cth))
        ethcm=phieth/dtldt0
      elseif (phitype.eq.1) then ! mass not neglected and elab<gamma m
        if (cth.le.phicthmax) then
          elab=phieth+phimp
          tcm(1)=max((elab-
     &        phibep*cth*sqrt(max(elab**2-phigap**2*phimp**2+phigap**2
     &        *phibep**2*cth**2*phimp**2,0.0d0)))/
     &        (phigap*(1.0d0-phibep*cth)*(1.0d0+phibep*cth))
     &        -phimp,tmin)    ! kinetic energy in parent cm system
          tcm(2)=max((elab+
     &        phibep*cth*sqrt(max(elab**2-phigap**2*phimp**2+phigap**2
     &        *phibep**2*cth**2*phimp**2,0.0d0)))/
     &        (phigap*(1.0d0-phibep*cth)*(1.0d0+phibep*cth))
     &        -phimp,tmin)    ! kinetic energy in parent cm system

          if (phifk.ge.100) then ! differential yields
            dtdt(1)=abs(phigap*(1.0d0+phibep*cth*(tcm(1)+phimp)/
     &        sqrt(max(tcm(1)**2+2.0d0*tcm(1)*phimp,pmin**2))))
            dtdt(2)=abs(phigap*(1.0d0+phibep*cth*(tcm(2)+phimp)/
     &        sqrt(max(tcm(2)**2+2.0d0*tcm(2)*phimp,pmin**2))))
          else
            tcm(1)=max(tmin,tcm(1))
            tcm(2)=max(tmin,tcm(2))
          endif
c          write(*,*) 'tcm(1) = ',tcm(1),'  dtdt(1) = ',dtdt(1)
c          write(*,*) 'tcm(2) = ',tcm(2),'  dtdt(2) = ',dtdt(2)
        else
          dmhadys=0.0d0
c          write(*,*) 'cth = ',cth,' > phicthmax'
          return
        endif
      else  ! phitype.eq.2  ! mass not neglected, elab>gamma m
        elab=phieth+phimp
        tcm(1)=max((elab-
     &        phibep*cth*sqrt(max(elab**2-phigap**2*phimp**2+phigap**2
     &        *phibep**2*cth**2*phimp**2,0.0d0)))/
     &        (phigap*(1.0d0-phibep*cth)*(1.0d0+phibep*cth))
     &        -phimp,tmin)    ! kinetic energy in parent cm system
        dtdt(1)=abs(phigap*(1.0d0+phibep*cth*(tcm(1)+phimp)/
     &    sqrt(max(tcm(1)**2+2.0d0*tcm(1)*phimp,pmin**2))))
      endif

c... get the yield at this energy
      yield=0.0d0

      do ch=1,7
        mp=map(ch)
        ep=e1cm
        if (has0br(hdi(ch),phihno).gt.0.0d0) then
          if(ep.lt.0.99d0*mp) then
            write(*,*) 'error in dmhadys: ep < 0.99 mp'
            write(*,*) '  channel = ',ch
            write(*,*) '  mp = ',mp
            write(*,*) '  ep = ',ep
          else
            ep=max(ep,mp*1.0001d0)
          endif


          if (ep.gt.ub(ch)) then
            write(*,*) 'error in dmhadys: too high ep = ',ep
          endif

          if (ep.lt.lb(ch)) goto 150

          if (phitype.eq.0) then
            yield1=has0br(hdi(ch),phihno)*
     &        dmhayieldf(e1cm,ethcm,ch,phifk,istat)
            if (phifk.ge.100) then  ! differential yields
              jac=1.0d0/dtldt0
              yield1=yield1*jac
            endif
            yield=yield+yield1
          elseif (phitype.eq.1) then
            yield1=dmhayieldf(e1cm,tcm(1),ch,phifk,istat)
            yield2=dmhayieldf(e1cm,tcm(2),ch,phifk,istat)
            if (phifk.ge.100) then  ! differential yields
              jac=1.0d0/dtdt(1)
              yield1=yield1*jac
              jac=1.0d0/dtdt(2)
              yield2=yield2*jac
            endif
            yield=yield+has0br(hdi(ch),phihno)*(yield1+yield2)
          else  ! phitype.eq.2
            yield1=has0br(hdi(ch),phihno)
     &        *dmhayieldf(e1cm,tcm(1),ch,phifk,
     &      istat)
            if (phifk.ge.100) then  ! differential yields
              jac=1.0d0/dtdt(1)
              yield1=yield1*jac
            endif
            yield=yield+yield1
          endif
  150   continue
        endif
      enddo

c... take care of factors depending on energy in neutrino-nucleus cross
c... section and muon range
      if (phifk.eq.72.or.phifk.eq.172) then
        yield=yield*dtldt0   ! nu-nucleon cross section
      elseif (phifk.eq.73.or.phifk.eq.173) then
        yield=yield*(dtldt0**2)  ! nu-nucleon cross section and mu range
      endif

c      if (phifk.eq.72.or.phifk.eq.73) then
c        write(*,*) 'warning in dmhadys: integrated yields for',
c     &    ' channel ',phifk,' has to be checked!'
c      endif

c...  include a factor of 1/2 to get average of yield,
c...  (1/(int_-1^1 1 dcth)=1/2).

      dmhadys=0.5d0*yield
c      write(*,*) 'cth = ',cth,'  dshadys = ',0.25d0*yield

      end











