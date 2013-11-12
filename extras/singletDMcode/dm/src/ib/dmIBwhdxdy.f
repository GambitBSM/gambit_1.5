*****************************************************************************
*** The function dsIBwhdxdy gives the full analytical expressions for
*** the differential IB photon yield (dxdy) from W+H- and W-H+ final states, 
*** normalized to the  annihilation rate into W+H- and W-H+
***
*** The kinematic variables x,y are
***
*** x = E_gamma/mx
*** y = (p+k)^2/(4 mx^2),
***
*** where p denotes the W momentum and k the photon momentum.
*** (note that the expressions above and below only apply to the v->0 limit)
***
*** author: Torsten Bringmann (bringman@sissa.it)
*** date: 2007-05-01
*****************************************************************************


      real*8 function dmIBwhdxdy(IBch,x,y)

      implicit none
      include 'dssusy.h'
      include 'dshacom.h'
      include 'dsprep.h'
      include 'dshmcom.h'
      include 'dsidtag.h'

      integer IBch
      real*8 x,y,tmpresult,dmsigmav,sv,msq2bodyds
      real*8 dmIBwhdxdy_1,dmIBwhdxdy_2,dmIBwhdxdy_3
      real*8 dmIBwhdxdy_4,dmIBwhdxdy_5,dmIBwhdxdy_6
      real*8 dmIBwhdxdy_7,dmIBwhdxdy_8,dmIBwhdxdy_9

      character*12 memory    ! to suppress multiple error messages
      save memory
      data memory /'____________'/


c      real*8 msq2body        ! to check the lowest order annihilation rate

 
c------------- couplings and masses relevant for IB -------------------

      complex*8 CHW11, CHW12, CHW21, CHW22, CHW3, 
     -          CHWA11, CHWA12, CHWA21, CHWA22
      real*8 mc1, mc2, mh03, Gh03, mhc, m0, mw

      dmIBwhdxdy=0d0
      tmpresult=0d0

      if (IBch.ne.2) return

c...make sure we computed the branching ratios
      sv=dmsigmav(0)
      if ((1D-6).gt.(dmsigmav(11)/sv)) return  ! continue only for non-zero br

c...set up masses

      m0   = mass(kn(1))
      mw   = mass(kw)
      mhc  = mass(khc)
      mc1  = mass(kcha(1))
      mc2  = mass(kcha(2))
      mh03 = mass(kh3)
      Gh03 = width(kh3)
      
c...set up couplings
  
      CHW11  = gr(khc,kn(1),kcha(1))*conjg(gr(kw,kn(1),kcha(1))) -
     -         gl(khc,kn(1),kcha(1))*conjg(gl(kw,kn(1),kcha(1)))
      CHW12  = gr(khc,kn(1),kcha(2))*conjg(gr(kw,kn(1),kcha(2))) -
     -         gl(khc,kn(1),kcha(2))*conjg(gl(kw,kn(1),kcha(2)))
      CHW21  = gr(khc,kn(1),kcha(1))*conjg(gl(kw,kn(1),kcha(1))) -
     -         gl(khc,kn(1),kcha(1))*conjg(gr(kw,kn(1),kcha(1)))
      CHW22  = gr(khc,kn(1),kcha(2))*conjg(gl(kw,kn(1),kcha(2))) -
     -         gl(khc,kn(1),kcha(2))*conjg(gr(kw,kn(1),kcha(2)))
      CHWA11 = gr(khc,kn(1),kcha(1))*conjg(gr(kw,kn(1),kcha(1))) +
     -         gl(khc,kn(1),kcha(1))*conjg(gl(kw,kn(1),kcha(1)))
      CHWA12 = gr(khc,kn(1),kcha(2))*conjg(gr(kw,kn(1),kcha(2))) +
     -         gl(khc,kn(1),kcha(2))*conjg(gl(kw,kn(1),kcha(2)))
      CHWA21 = gr(khc,kn(1),kcha(1))*conjg(gl(kw,kn(1),kcha(1))) +
     -         gl(khc,kn(1),kcha(1))*conjg(gr(kw,kn(1),kcha(1)))
      CHWA22 = gr(khc,kn(1),kcha(2))*conjg(gl(kw,kn(1),kcha(2))) +
     -         gl(khc,kn(1),kcha(2))*conjg(gr(kw,kn(1),kcha(2)))
      CHW3 = gr(kw,khc,kh3)*(gr(kh3,kn(1),kn(1))-gl(kh3,kn(1),kn(1)))/
     -         cmplx(4*m0**2-mh03**2,mh03*Gh03)


c...import IB expressions for |M|**2 from form/mathematica; add  
c...contributions from all linear independent combinations of coupling
c...constants separately, utilizing various symmetries

      tmpresult=
     -   + 2*dble(CHW11*conjg(CHW12))
     -       *dmIBwhdxdy_1(x,y,m0,mw,mhc,mc1,mc2)
     -   + (abs(CHW11))**2
     -       *dmIBwhdxdy_1(x,y,m0,mw,mhc,mc1,mc1)
     -   + (abs(CHW12))**2
     -       *dmIBwhdxdy_1(x,y,m0,mw,mhc,mc2,mc2)
     -
     -   + 2*dble(CHW11*conjg(CHW22))
     -       *dmIBwhdxdy_2(x,y,m0,mw,mhc,mc1,mc2)
     -   + 2*dble(CHW11*conjg(CHW21))
     -       *dmIBwhdxdy_2(x,y,m0,mw,mhc,mc1,mc1)
     -   + 2*dble(CHW12*conjg(CHW22))
     -       *dmIBwhdxdy_2(x,y,m0,mw,mhc,mc2,mc2)
     -   + 2*dble(CHW12*conjg(CHW21))
     -       *dmIBwhdxdy_2(x,y,m0,mw,mhc,mc2,mc1)
     - 
     -   + 2*dble(CHW11*conjg(CHW3))
     -       *dmIBwhdxdy_3(x,y,m0,mw,mhc,mc1)
     -   + 2*dble(CHW12*conjg(CHW3))
     -       *dmIBwhdxdy_3(x,y,m0,mw,mhc,mc2)
     -
     -   + 2*dble(CHW21*conjg(CHW22))
     -       *dmIBwhdxdy_4(x,y,m0,mw,mhc,mc1,mc2)
     -   + (abs(CHW21))**2
     -       *dmIBwhdxdy_4(x,y,m0,mw,mhc,mc1,mc1)
     -   + (abs(CHW22))**2
     -       *dmIBwhdxdy_4(x,y,m0,mw,mhc,mc2,mc2)
     -
     -   + 2*dble(CHW21*conjg(CHW3))
     -       *dmIBwhdxdy_5(x,y,m0,mw,mhc,mc1)
     -   + 2*dble(CHW22*conjg(CHW3))
     -       *dmIBwhdxdy_5(x,y,m0,mw,mhc,mc2)
     -
     -   + (abs(CHW3))**2*dmIBwhdxdy_6(x,y,m0,mw,mhc)
     -
     -   + 2*dble(CHWA11*conjg(CHWA12))
     -       *dmIBwhdxdy_7(x,y,m0,mw,mhc,mc1,mc2)
     -   + (abs(CHWA11))**2
     -       *dmIBwhdxdy_7(x,y,m0,mw,mhc,mc1,mc1)
     -   + (abs(CHWA12))**2
     -       *dmIBwhdxdy_7(x,y,m0,mw,mhc,mc2,mc2)
     -
     -   + 2*dble(CHWA21*conjg(CHWA22))
     -       *dmIBwhdxdy_8(x,y,m0,mw,mhc,mc1,mc2)
     -   + (abs(CHWA21))**2
     -       *dmIBwhdxdy_8(x,y,m0,mw,mhc,mc1,mc1)
     -   + (abs(CHWA22))**2
     -       *dmIBwhdxdy_8(x,y,m0,mw,mhc,mc2,mc2)
     -
     -   + 2*dble(CHWA11*conjg(CHWA22))
     -       *dmIBwhdxdy_9(x,y,m0,mw,mhc,mc1,mc2)
     -   + 2*dble(CHWA11*conjg(CHWA21))
     -       *dmIBwhdxdy_9(x,y,m0,mw,mhc,mc1,mc1)
     -   + 2*dble(CHWA12*conjg(CHWA22))
     -       *dmIBwhdxdy_9(x,y,m0,mw,mhc,mc2,mc2)
     -   + 2*dble(CHWA12*conjg(CHWA21))
     -       *dmIBwhdxdy_9(x,y,m0,mw,mhc,mc2,mc1)


c... Since the branching ratios were calculated with dmsigmav, do
c... the same for the photon multiplicity for consistency.

      msq2bodyds=dmsigmav(11)*32*pi*m0**2/gev2cm3s/
     -           sqrt(1-(mw**2+mhc**2)/2./m0**2
     -                +(mw**2-mhc**2)**2/16./m0**4)


c... for comparison, this is the analytically obtained 2-body
c... annihilation amplitude squared (in the same v->0 limit
c... as considered for the 3-body case) :
c 
c      msq2body=  
c     -   ((16*m0**4 + (mhc**2 - mw**2)**2 - 8*m0**2*(mhc**2 + mw**2))*
c     -    ((m0**2*dble(CHW11*conjg(CHW11)) + 
c     -      2*m0*mc1*dble(CHW11*conjg(CHW21)) + 
c     -         mc1**2*dble(CHW21*conjg(CHW21)))/
c     -       (-2*m0**2 - 2*mc1**2 + mhc**2 + mw**2)**2 + 
c     -      (2*m0**2*dble(CHW11*conjg(CHW12)) + 
c     -        2*m0*mc1*dble(CHW12*conjg(CHW21)) + 
c     -         2*m0*mc2*dble(CHW11*conjg(CHW22)) + 
c     -         2*mc1*mc2*dble(CHW21*conjg(CHW22)))
c     -        /((-2*m0**2 - 2*mc1**2 + mhc**2 + mw**2)*
c     -         (-2*m0**2 - 2*mc2**2 + mhc**2 + mw**2)) + 
c     -      (m0**2*dble(CHW12*conjg(CHW12)) +
c     -          2*m0*mc2*dble(CHW12*conjg(CHW22)) + 
c     -         mc2**2*dble(CHW22*conjg(CHW22)))/
c     -       (-2*m0**2 - 2*mc2**2 + mhc**2 + mw**2)**2 + 
c     -      (2*m0**2*dble(CHW11*conjg(CHW3)) +
c     -         2*m0*mc1*dble(CHW21*conjg(CHW3)))/
c     -       (-2*m0**2 - 2*mc1**2 + mhc**2 + mw**2) + 
c     -      (2*m0**2*dble(CHW12*conjg(CHW3)) +
c     -         2*m0*mc2*dble(CHW22*conjg(CHW3)))/
c     -       (-2*m0**2 - 2*mc2**2 + mhc**2 + mw**2) + 
c     -      m0**2*dble(CHW3*conjg(CHW3))))/(2.*mw**2)
c 
c      msq2body = 2*msq2body        !W+H- & W-H+ final states
c 
c... check  lowest order result
c      if (((1d0-msq2body/msq2bodyds)**2.ge.0.001)
c     -     .and.(memory.ne.idtag)) then
c        memory=idtag
c        write(*,*) '*****'
c        write(*,*) 'model ',idtag,' - ', 'warning from ',
c     -             'dmIBwhdxdy:'
c        write(*,*) 'IB contribution from channel WH',
c     -             ' could be up to a factor',msq2bodyds/msq2body,
c     -             ' higher!'
c        write(*,*) '*****'
c       endif


c...check that result is positive
      if (0.gt.tmpresult) then
        if (m0**2*tmpresult.lt.(-1D-10)
     -     .and.(idtag.ne.memory)) then
          write(*,*) '*****'
          write (*,*) 'Error in dmIBwhdxdy: ',
     -                'negative |M|^2 for model ',idtag,'.'
          write (*,*) 'Setting corresponding contributions to zero...'
          write(*,*) '*****'
        endif
        return
      endif


c... The photon multiplicity is given by the ratio of the squared
c... amplitudes, times a phase space factor:

      tmpresult = (alphem/pi)*m0**2*(tmpresult/msq2bodyds)/
     -     sqrt(1-(mw**2+mhc**2)/2./m0**2+(mw**2-mhc**2)**2/16./m0**4)


      dmIBwhdxdy=tmpresult

      return
      end

