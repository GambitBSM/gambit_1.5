*****************************************************************************
*** The function dsIBffdxdy gives the full analytical expressions for
*** the differential IB photon yield (dxdy) from fermion final states, 
*** normalized to the annihilation rate into fermion pairs f fbar
***
*** The kinematic variables x,y are
***
*** x = E_gamma/mx
*** y = (p+k)^2/(4 mx^2),
***
*** where p denotes the fermion momentum and k the photon momentum.
*** (note that the expressions above and below only apply to the v->0 limit)
***
*** author: Torsten Bringmann (bringman@sissa.it)
*** date: 2007-07-05
*****************************************************************************

      real*8 function dmIBffdxdy(IBch,x,y)

      implicit none
      include 'dssusy.h'
      include 'dshacom.h'
      include 'dsprep.h'
      include 'dshmcom.h'
      include 'dsidtag.h'


      real*8 x,y,tmpresult,tmpdecay,dmsigmav,sv
      real*8 msq2bodyds
      integer IBch, kf, ksf1, ksf2                ! IB channel and particle codes
      real*8 dmIBffdxdy_1,dmIBffdxdy_2,dmIBffdxdy_3
      real*8 dmIBffdxdy_4,dmIBffdxdy_5,dmIBffdxdy_6
      real*8 dmIBffdxdy_7,dmIBffdxdy_8,dmIBfsrdxdy

      character*12 memory    ! to suppress multiple error messages
      save memory
      data memory /'____________'/

c      real*8 msq2body        ! to check the lowest order annihilation rate


c------------- couplings and masses relevant for IB -------------------

      real*8 C11, C12, C151, C152, C21, C22, C251, C252, CZ5, CH
      real*8 m0,msf2, msf1, mz, mf, mh03, Gh03, GZ, zf, hf

    
      dmIBffdxdy=0d0
      tmpresult=0d0
      tmpdecay=0d0

      if ((IBch.lt.4).or.(IBch.gt.12)) return

c...make sure we computed the branching ratios
      sv=dmsigmav(0)

c...determine fermion and sfermion particle codes
     
      if (IBch.eq.4) then
        kf   = ke
        ksf1 = kse(1)
        ksf2 = kse(2)
      endif
      if (IBch.eq.5) then
        kf   = kmu
        ksf1 = ksmu(1)
        ksf2 = ksmu(2)
      endif
      if (IBch.eq.6) then
        kf   = ktau
        ksf1 = kstau(1)
        ksf2 = kstau(2)
      endif
      if (IBch.eq.7) then
        kf   = ku
        ksf1 = ksu(1)
        ksf2 = ksu(2)
      endif
      if (IBch.eq.8) then
        kf   = kd
        ksf1 = ksd(1)
        ksf2 = ksd(2)
      endif
      if (IBch.eq.9) then
        kf   = kc
        ksf1 = ksc(1)
        ksf2 = ksc(2)
      endif
      if (IBch.eq.10) then
        kf   = ks
        ksf1 = kss(1)
        ksf2 = kss(2)
      endif
      if (IBch.eq.11) then
        kf   = kt
        ksf1 = kst(1)
        ksf2 = kst(2)
      endif
      if (IBch.eq.12) then
        kf   = kb
        ksf1 = ksb(1)
        ksf2 = ksb(2)
      endif

c...set up masses and widths

      m0   = mass(kn(1))
      mz   = mass(kz)
      mh03 = mass(kh3)
      Gh03 = width(kh3)
      GZ   = width(kz)

      mf   = mass(kf)
      msf1 = mass(ksf1)
      msf2 = mass(ksf2)
      Gh03 = width(kh3)
  
c...set up couplings
 
      C11   = abs(gl(ksf1,kf,kn(1)))**2 - 
     -        abs(gr(ksf1,kf,kn(1)))**2
      C12   = abs(gl(ksf2,kf,kn(1)))**2 - 
     -        abs(gr(ksf2,kf,kn(1)))**2
      C151  = abs(gl(ksf1,kf,kn(1)))**2 + 
     -        abs(gr(ksf1,kf,kn(1)))**2
      C152  = abs(gl(ksf2,kf,kn(1)))**2 + 
     -        abs(gr(ksf2,kf,kn(1)))**2
      C21   = imag(gr(ksf1,kf,kn(1))*conjg(gl(ksf1,kf,kn(1))))
      C22   = imag(gr(ksf2,kf,kn(1))*conjg(gl(ksf2,kf,kn(1))))
      C251  = dble(gr(ksf1,kf,kn(1))*conjg(gl(ksf1,kf,kn(1))))
      C252  = dble(gr(ksf2,kf,kn(1))*conjg(gl(ksf2,kf,kn(1))))
      CZ5   = dble(gl(kz,kn(1),kn(1)))*
     -        dble(gr(kz,kf,kf)-gl(kz,kf,kf))
      CH    = imag(gr(kh3,kn(1),kn(1)))*imag(gr(kh3,kf,kf))


c...import IB expressions for |M|**2 from form/mathematica; 
c...for light leptons and quarks, take the simplified expression with mf=0:

      if ((IBch.eq.4).or.(IBch.eq.5).or.(IBch.eq.7).or.
     -    (IBch.eq.8).or.(IBch.eq.10)) then
            tmpresult=
     -        -2*m0**6*(-1 + x)*(x**2 - 2*x*y + 2*y**2)*
     -       ((C11/((msf1**2 + m0**2*(1 - 2*y))*
     -             (msf1**2 + m0**2*(1 - 2*x + 2*y))) + 
     -         C12/
     -          ((msf2**2 + m0**2*(1 - 2*y))*
     -          (msf2**2 + m0**2*(1 - 2*x + 2*y))))**2 + 
     -        (C151/
     -          ((msf1**2 + m0**2*(1 - 2*y))*
     -            (msf1**2 + m0**2*(1 - 2*x + 2*y))) + 
     -          C152/
     -          ((msf2**2 + m0**2*(1 - 2*y))*
     -            (msf2**2 + m0**2*(1 - 2*x + 2*y))))**2)
           endif

c...for heavy leptons and quarks, add
c...contributions from all linear independent combinations of coupling
c...constants separately, utilizing various symmetries

      if ((IBch.eq.6).or.(IBch.eq.9).or.(IBch.eq.11)
     -    .or.(IBch.eq.12)) then

          zf=(mz**2*(16*m0**4+GZ**2*mz**2-8*m0**2*mz**2
     -       +mz**4))/(-4*m0**2 + mz**2)**2
          hf=((-16*m0**4-Gh03**2*mh03**2+8*m0**2*mh03**2
     -       -mh03**4)*mf)/(2.*m0*(4*m0**2 - mh03**2))

      tmpresult=
     -     C11*C12*dmIBffdxdy_1(x,y,m0,mf,msf1,msf2)
     -  +  C11**2*dmIBffdxdy_1(x,y,m0,mf,msf1,msf1)/2.
     -  +  C12**2*dmIBffdxdy_1(x,y,m0,mf,msf2,msf2)/2.
     -
     -  +  C151*C152*dmIBffdxdy_2(x,y,m0,mf,msf1,msf2)
     -  +  C151**2*dmIBffdxdy_2(x,y,m0,mf,msf1,msf1)/2.
     -  +  C152**2*dmIBffdxdy_2(x,y,m0,mf,msf2,msf2)/2
     -
     -  +  C21*C22*dmIBffdxdy_3(x,y,m0,mf,msf1,msf2)
     -  +  C21**2*dmIBffdxdy_3(x,y,m0,mf,msf1,msf1)/2.
     -  +  C22**2*dmIBffdxdy_3(x,y,m0,mf,msf2,msf2)/2.
     -
     -  +  C251*C252*dmIBffdxdy_4(x,y,m0,mf,msf1,msf2)
     -  +  C251**2*dmIBffdxdy_4(x,y,m0,mf,msf1,msf1)/2.
     -  +  C252**2*dmIBffdxdy_4(x,y,m0,mf,msf2,msf2)/2.
     -
     -  +  C151*C252*dmIBffdxdy_5(x,y,m0,mf,msf1,msf2)
     -  +  C151*C251*dmIBffdxdy_5(x,y,m0,mf,msf1,msf1)
     -  +  C152*C252*dmIBffdxdy_5(x,y,m0,mf,msf2,msf2)
     -  +  C152*C251*dmIBffdxdy_5(x,y,m0,mf,msf2,msf1)
     -
     -  +  C151*CZ5*dmIBffdxdy_6(x,y,m0,mf,msf1)/zf
     -  +  C152*CZ5*dmIBffdxdy_6(x,y,m0,mf,msf2)/zf
     -  +  C151*CH*dmIBffdxdy_6(x,y,m0,mf,msf1)/hf
     -  +  C152*CH*dmIBffdxdy_6(x,y,m0,mf,msf2)/hf
     -
     -  +  C251*CZ5*dmIBffdxdy_7(x,y,m0,mf,msf1)/zf
     -  +  C252*CZ5*dmIBffdxdy_7(x,y,m0,mf,msf2)/zf
     -  +  C251*CH*dmIBffdxdy_7(x,y,m0,mf,msf1)/hf
     -  +  C252*CH*dmIBffdxdy_7(x,y,m0,mf,msf2)/hf
     -
     -  +  dmIBffdxdy_8(x,y,CZ5,CH,m0,mf,mz,mh03,GZ,Gh03)

      endif   ! contribution from heavy fermions


c... Since the branching ratios were calculated with dmsigmav, do
c... the same for the photon multiplicity for consistency.

      if (IBch.eq.4) msq2bodyds=dmsigmav(15)
      if (IBch.eq.5) msq2bodyds=dmsigmav(17)
      if (IBch.eq.6) msq2bodyds=dmsigmav(19)
      if (IBch.eq.7) msq2bodyds=dmsigmav(20)/3.   ! color factor for quarks
      if (IBch.eq.8) msq2bodyds=dmsigmav(21)/3.
      if (IBch.eq.9) msq2bodyds=dmsigmav(22)/3.
      if (IBch.eq.10) msq2bodyds=dmsigmav(23)/3.
      if (IBch.eq.11) msq2bodyds=dmsigmav(24)/3.
      if (IBch.eq.12) msq2bodyds=dmsigmav(25)/3.
      msq2bodyds=msq2bodyds*32*pi*m0**2/sqrt(1-mf**2/m0**2)/
     -           gev2cm3s


c... for comparison, this is the analytically obtained 2-body 
c... annihilation amplitude squared (in the same v->0 limit, 
c... and up to color factors, as considered for the 3-body case) :
c 
c      msq2body=  
c     -  4*(m0**4 - m0**2*mf**2)*
c     -   (C21/(m0**2 - mf**2 + msf1**2) + 
c     -      C22/(m0**2 - mf**2 + msf2**2))**2 + 
c     -  m0**2*((2*C251*m0 + C151*mf)/(m0**2 - mf**2 + msf1**2) + 
c     -      (2*C252*m0 + C152*mf)/(m0**2 - mf**2 + msf2**2))**2 + 
c     -  (4*CZ5*m0**2*mf*((2*C251*m0 + C151*mf)/
c     -        (m0**2 - mf**2 + msf1**2) + 
c     -       (2*C252*m0 + C152*mf)/(m0**2 - mf**2 + msf2**2) + 
c     -       (CZ5*mf)/mz**2)*(-4*m0**2 + mz**2)**2)/
c     -   (mz**2*(GZ**2*mz**2 + (-4*m0**2 + mz**2)**2)) + 
c     -  (8*CH*m0**3*(2*CH*m0 + 
c     -       (-4*m0**2 + mh03**2)*
c     -        ((2*C251*m0 + C151*mf)/(m0**2 - mf**2 + msf1**2) + 
c     -          (2*C252*m0 + C152*mf)/(m0**2 - mf**2 + msf2**2)) + 
c     -       (2*CZ5*mf*(-(Gh03*GZ*mh03*mz*(4*m0**2 - mz**2)) + 
c     -            (-4*m0**2 + mh03**2)*(-4*m0**2 + mz**2)**2))/
c     -        (mz**2*(GZ**2*mz**2 + (-4*m0**2 + mz**2)**2))))/
c     -   (Gh03**2*mh03**2 + (-4*m0**2 + mh03**2)**2)
c 
c... compare msq2body to msq2bodyds:
c      if (((1d0-msq2body/msq2bodyds)**2.ge.0.01)
c     -     .and.(memory.ne.idtag)) then
c        memory=idtag
c        write(*,*) '*****'
c        write(*,*) 'model ',idtag,' - ', 'warning from ',
c     -             'dmIBffdxdy:'
c        write(*,*) 'IB contribution from channel',IBch,
c     -             ' could be up to a factor',msq2bodyds/msq2body,
c     -             ' higher!'
c        write(*,*) '*****'
c       endif

c...check that result is positive
      if (0.gt.tmpresult) then
        if (m0**2*tmpresult.lt.(-1D-12)
     -     .and.(idtag.ne.memory)) then
          write(*,*) '*****'
          write (*,*) 'Error in dmIBffdxdy (channel:',IBch,
     -                '): negative |M|^2 for model ',idtag,'.'
          write (*,*) 'Setting corresponding contributions to zero...'
          write(*,*) '*****'
        endif
        return
      endif

c...take into account different electric charges for quarks

      if ((IBch.eq.7).or.(IBch.eq.9).or.(IBch.eq.11))
     -    tmpresult=tmpresult*4/9.
      if ((IBch.eq.8).or.(IBch.eq.10).or.(IBch.eq.12))
     -    tmpresult=tmpresult/9.


c... The photon multiplicity is given by the ratio of the squared
c... amplitudes, times a phase space factor:

      tmpresult = (alphem/pi)*m0**2/sqrt(1-mf**2/m0**2)
     -            *(tmpresult/msq2bodyds)


c... Finally, subtract the model-independent part which should already have
c... been taken into account by Pythia

      if (IBch.eq.5.or.IBch.eq.6.or.                ! do this only for those 
     -    IBch.eq.9.or.IBch.eq.11.or.IBch.eq.12)    ! channels for which Pythia
     -                                              ! runs exist.
     -                                              ! THIS IS A TEMPORARY FIX
     -    tmpresult = tmpresult - dmIBfsrdxdy(IBch,x,y)


      dmIBffdxdy=tmpresult
      return

      end

