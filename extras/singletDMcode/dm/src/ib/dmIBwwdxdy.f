*****************************************************************************
*** The function dsIBwwdxdy gives the full analytical expressions for
*** the differential IB photon yield (dxdy) from W+W- final states, 
*** normalized to the  annihilation rate into W+W-
***
*** The kinematic variables x,y are
***
*** x = E_gamma/mx
*** y = (p+k)^2/(4 mx^2),
***
*** where p denotes the W+ momentum and k the photon momentum.
*** (note that the expressions above and below only apply to the v->0 limit)
***
*** author: Torsten Bringmann (bringman@sissa.it)
*** date: 2007-05-01
*****************************************************************************

      real*8 function dmIBwwdxdy(IBch,x,y)

      implicit none
      include 'dssusy.h'
      include 'dshacom.h'
      include 'dsprep.h'
      include 'dshmcom.h'
      include 'dsidtag.h'

 
      integer IBch
      real*8 x,y,tmpresult,dmsigmav,sv,msq2bodyds
      real*8 dmIBwwdxdy_1,dmIBwwdxdy_2,dmIBwwdxdy_3
      real*8 dmIBwwdxdy_4,dmIBwwdxdy_5

      character*12 memory    ! to suppress multiple error messages
      save memory
      data memory /'____________'/

c      real*8 msq2body        ! to check the lowest order annihilation rate

 
c------------- couplings and masses relevant for IB -------------------

      real*8 Cllww1, Crrww1, Clrww1, Clrww51, 
     &       Cllww2, Crrww2, Clrww2, Clrww52
  
      real*8 mc1, mc2, m0, mw

      dmIBwwdxdy=0d0
      tmpresult=0d0

      if (IBch.ne.1) return

c...set up masses

      m0   = x
      mw   = mass(kw)
      mc1  = mass(kcha(1))
      mc2  = mass(kcha(2))
      
c...set up couplings
 
      Cllww1  = abs(gl(kw,kn(1),kcha(1)))**2
      Crrww1  = abs(gr(kw,kn(1),kcha(1)))**2
      Clrww1  = imag(gl(kw,kn(1),kcha(1))*conjg(gr(kw,kn(1),kcha(1))))
      Clrww51 = dble(gl(kw,kn(1),kcha(1))*conjg(gr(kw,kn(1),kcha(1))))
      Cllww2  = abs(gl(kw,kn(1),kcha(2)))**2
      Crrww2  = abs(gr(kw,kn(1),kcha(2)))**2
      Clrww2  = imag(gl(kw,kn(1),kcha(2))*conjg(gr(kw,kn(1),kcha(2))))
      Clrww52 = dble(gl(kw,kn(1),kcha(2))*conjg(gr(kw,kn(1),kcha(2))))
     

c...import IB expressions for |M|**2 from form/mathematica; add  
c...contributions from all linear independent combinations of coupling
c...constants separately, utilizing various symmetries

      tmpresult=
     -     (Cllww1+Crrww1)*(Cllww2+Crrww2)
     -       *dmIBwwdxdy_1(x,y,m0,mw,mc1,mc2)
     -   + (Cllww1+Crrww1)**2/2.
     -       *dmIBwwdxdy_1(x,y,m0,mw,mc1,mc1)
     -   + (Cllww2+Crrww2)**2/2.
     -       *dmIBwwdxdy_1(x,y,m0,mw,mc2,mc2)
     - 
     -   + (Cllww1*Crrww2+Crrww1*Cllww2)
     -       *dmIBwwdxdy_2(x,y,m0,mw,mc1,mc2)
     -   + Cllww1*Crrww1*dmIBwwdxdy_2(x,y,m0,mw,mc1,mc1)
     -   + Cllww2*Crrww2*dmIBwwdxdy_2(x,y,m0,mw,mc2,mc2)
     -
     -   + (Cllww1 + Crrww1)*Clrww52
     -       *dmIBwwdxdy_3(x,y,m0,mw,mc1,mc2)
     -   + (Cllww2 + Crrww2)*Clrww51
     -       *dmIBwwdxdy_3(x,y,m0,mw,mc2,mc1)
     -   + (Cllww1 + Crrww1)*Clrww51
     -       *dmIBwwdxdy_3(x,y,m0,mw,mc1,mc1)
     -   + (Cllww2 + Crrww2)*Clrww52
     -       *dmIBwwdxdy_3(x,y,m0,mw,mc2,mc2)
     -
     -   + Clrww1*Clrww2*dmIBwwdxdy_4(x,y,m0,mw,mc1,mc2)
     -   + Clrww1*Clrww1*dmIBwwdxdy_4(x,y,m0,mw,mc1,mc1)/2.
     -   + Clrww2*Clrww2*dmIBwwdxdy_4(x,y,m0,mw,mc2,mc2)/2.
     -
     -   + Clrww51*Clrww52*dmIBwwdxdy_5(x,y,m0,mw,mc1,mc2)
     -   + Clrww51*Clrww51*dmIBwwdxdy_5(x,y,m0,mw,mc1,mc1)/2.
     -   + Clrww52*Clrww52*dmIBwwdxdy_5(x,y,m0,mw,mc2,mc2)/2.

c... Since the branching ratios were calculated with dmsigmav, do
c... the same for the photon multiplicity for consistency.

      msq2bodyds=dmsigmav(13)*32*pi*m0**2/sqrt(1-mw**2/m0**2)/
     -           gev2cm3s

c... for comparison, this is the analytically obtained 2-body
c... annihilation amplitude squared (in the same v->0 limit,
c... and up to color factors, as considered for the 3-body case) :
c
c      msq2body=  
c     -   4*m0**2*((m0**2 - mw**2)*
c     -     ((Cllww1 + Crrww1)/(m0**2 + mc1**2 - mw**2) + 
c     -        (Cllww2 + Crrww2)/(m0**2 + mc2**2 - mw**2))**2 + 
c     -    (2*(4*m0**4 - 4*m0**2*mw**2 + 3*mw**4)*
c     -       ((Clrww1*mc1)/(m0**2 + mc1**2 - mw**2) + 
c     -          (Clrww2*mc2)/(m0**2 + mc2**2 - mw**2))**2)/mw**4)
c 
c... compare msq2body to msq2bodyds:
c      if (((1d0-msq2body/msq2bodyds)**2.ge.0.001)
c     -     .and.(memory.ne.idtag)) then
c        memory=idtag
c        write(*,*) '*****'
c        write(*,*) 'model ',idtag,' - ', 'warning from ',
c     -             'dmIBwwdxdy:'
c        write(*,*) 'IB contribution from WW',
c     -             ' could be up to a factor',msq2bodyds/msq2body,
c     -             ' higher!'
c        write(*,*) '*****'
c       endif


c...check that result is positive
      if (0.gt.tmpresult) then
        if (m0**2*tmpresult.lt.(-1D-10)
     -     .and.(idtag.ne.memory)) then
          write(*,*) '*****'
          write (*,*) 'Error in dmIBwwdxdy:',
     -                ' negative |M|^2 for model ',idtag,'.'
          write (*,*) 'Setting corresponding contributions to zero...'
          write(*,*) '*****'
        endif
        return
      endif


c... The photon multiplicity is given by the ratio of the squared
c... amplitudes, times a phase space factor:

      tmpresult = (alphem/pi)*m0**2/sqrt(1-mw**2/m0**2)
     -            *(tmpresult/msq2bodyds)

      dmIBwwdxdy=tmpresult
      return
             
      end
