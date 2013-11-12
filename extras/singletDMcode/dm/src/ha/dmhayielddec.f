*****************************************************************************
*** function dshayielddec integrates dshadys over cos theta.
*** the scalar decay channels are summed in dshadys
*** This routine is for the fundamental channels from S1, S2 and S3 decay.
*** units: (annihilation)**-1
*****************************************************************************

      real*8 function dmhayielddec(eh,hno,emuth,yieldk,istat)
      implicit none
      include 'dshacom.h'
      include 'dsidtag.h'

c------------------------ variables ------------------------------------

      real*8 eh,emuth
      real*8 elab,e1
      integer istat,ch,yieldk,hdi(7),hno

      real*8 sum,pi,dth
      parameter (pi=3.141592653589793238d0)
      parameter (dth=0.00001d0)  ! safety when integrating in theta

      data hdi/22,25,24,19,13,12,17/
c------------------------ functions ------------------------------------

      external dmhadys
      real*8 dmhadys,dmhayield_int

c-----------------------------------------------------------------------

c... check that there are no problems
      do ch=1,7
c...take care of slight numerical inaccuracy problems
        if (has0br(hdi(ch),hno).gt.0.0d0) then
          if (0.99d0*2.0d0*map(ch).gt.has0m(hno)) then
            write(*,*) 'error in dmhayielddec: a particle with mass ',
     &        has0m(hno)
            write(*,*) 'is let to decay to two particles with mass ',
     &        map(ch)
            write(*,*) 'which is not energetically allowed.'
            write(*,*) 'particle 1 has code ',ch
            write(*,*) 'the yield from this decay is put to 0.0'
            write(*,*) 'model: ',idtag
          endif
        endif
      enddo

c...set common block variables: boost parameters etc.
      phibep=sqrt(eh**2-has0m(hno)**2)/eh
      phigap=eh/has0m(hno)
      phihno=hno
      phifk=yieldk
      phieth=emuth
      e1=has0m(hno)/2.0d0

      if (yieldk.eq.54.or.yieldk.eq.154) then ! antiproton
        phimp=0.93827231
        elab=phieth+phimp

        if (elab.le.phigap*phimp) then
          phitype=1
          phicthmax=-sqrt(max((phigap**2*phimp**2-elab**2)/
     &      (phimp**2*phigap**2*phibep**2),0.0d0))
          phicthmax=max(-1.0d0,phicthmax-dth)
          phicthmin=max((phieth-e1*phigap)/
     &      (phigap*phibep*sqrt(e1**2-phimp**2)),-1.0d0)
        else
          phitype=2
          phicthmax=1.0d0
          phicthmin=max((phieth-e1*phigap)/
     &      (phigap*phibep*sqrt(e1**2-phimp**2)),-1.0d0)
        endif
      else  ! all others where mass can be neglected
        phimp=0.0d0
        phitype=0
        phicthmax=1.0d0
        phicthmin=max((phieth-e1*phigap)/
     &    (e1*phigap*phibep),-1.0d0)
      endif

      if (phicthmin.ne.-1.0d0) then
        phicthmin=min(phicthmin+dth,1.0d0)
      endif

      if (yieldk.le.100) then ! integrated yields ok to go all the way
        phicthmax=1.0d0
      endif

      if (phicthmin.ge.phicthmax) then
        dmhayielddec=0.0d0
        return
      endif

c      write(*,*) 'phitype = ',phitype
c      write(*,*) 'phicthmin = ',phicthmin
c      write(*,*) 'phicthmax = ',phicthmax

        sum=dmhayield_int(dmhadys,phicthmin,phicthmax)

c        eps=1.0d-3
c        call dgadap(phicthmin,phicthmax,dmhadys,eps,sum)

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

c        dmhayielddec=sum

        dmhayielddec=sum
c        write(*,*) 'dmhayielddec called with'
c        write(*,*) '    e0 = ',e0
c        write(*,*) '    m0 = ',m0
c        write(*,*) '    emuthr = ',emuthr
c        write(*,*) '    mp1 = ',mp1
c        write(*,*) '    mp2 = ',mp2
c        write(*,*) '    ch = ',ch
c        write(*,*) '  result is ',sum

 5000 format(' ',a,f8.2,a)
 5010 format(' ',a,i2,a,f8.2,a)
 5020 format(' ',a,f8.2,a,a)

      end




