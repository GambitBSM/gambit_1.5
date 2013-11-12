!         -*- mode: fortran -*-
!######################################################################*
!                       i n c l u d e     f i l e                      *
!######################################################################*

************************************************************************
***                             susy.h                               ***
***         this piece of code is needed as a separate file          ***
***               the rest of the code 'includes' susy.h             ***
c----------------------------------------------------------------------c
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994,1995
c    10-nov-95 complex vertex constants
c  modified by joakim edsjo (edsjo@teorfys.uu.se) 97-02-11
** general variables
*   particle codes
      integer knue,ke,knumu,kmu,knutau,ktau,ku,kd,kc,ks,kt,kb,kgamma,
     &     kw,kz,kgluon,kh1,kh2,kh3,khc,ksnue,kse1,kse2,ksnumu,ksmu1,
     &     ksmu2,ksnutau,kstau1,kstau2,ksu1,ksu2,ksd1,ksd2,ksc1,ksc2,
     &     kss1,kss2,kst1,kst2,ksb1,ksb2,kn1,kn2,kn3,kn4,kcha1,kcha2,
     &     kgluin,kgold0,kgoldc
      parameter (knue=1,ke=2,knumu=3,kmu=4,knutau=5,ktau=6)
      parameter (ku=7,kd=8,kc=9,ks=10,kt=11,kb=12)
      parameter (kgamma=13,kw=14,kz=15,kgluon=16)
      parameter (kh1=17,kh2=18,kh3=19,khc=20)
      parameter (ksnue=21,kse1=22,kse2=23,ksnumu=24,ksmu1=25,
     &     ksmu2=26,ksnutau=27,kstau1=28,kstau2=29)
      parameter (ksu1=30,ksu2=31,ksd1=32,ksd2=33,ksc1=34,ksc2=35,
     &     kss1=36,kss2=37,kst1=38,kst2=39,ksb1=40,ksb2=41)
      parameter (kn1=42,kn2=43,kn3=44,kn4=45,kcha1=46,kcha2=47,
     &     kgluin=48)
      parameter (kgold0=49,kgoldc=50)
c  knu=(nue,numu,nutau)   kl=(e,mu,tau)    kqu=(u,c,t)    kqd=(d,s,b)
c  ksqu=(~u1,~c1,~t1,~u2,~c2,~t2)    ksqd=(~d1,~s1,~b1,~d2,~s2,~b2)
      integer kse(2),ksmu(2),kstau(2),ksu(2),ksd(2),ksc(2),kss(2),
     &     kst(2),ksb(2),kn(4),kcha(2),knu(3),kl(3),kqu(3),kqd(3),
     &     ksnu(6),ksl(6),ksqu(6),ksqd(6)
      character*8 pname(0:50)
      common /pacodes/ kse,ksmu,kstau,ksu,ksd,ksc,kss,kst,ksb,kn,kcha,
     &     knu,kl,kqu,kqd,ksnu,ksl,ksqu,ksqd,pname
* program switches
      integer higloop,neuloop,bsgqcd
      real*8 msquarks,msleptons
      common /switch/ msquarks,msleptons,
     &     higloop,neuloop,bsgqcd
* set sfermion masses by hand
      real*8 massup1(3),massup2(3),thetamixu(3),
     & massdn1(3),massdn2(3),thetamixd(3),
     & masssn(3),masssl1(3),masssl2(3),thetamixsl(3)
      common/sfermionmass/ massup1,massup2,thetamixu,
     & massdn1,massdn2,thetamixd,
     & masssn,masssl1,masssl2,thetamixsl
* indirect rates variables
*      real*8 emuth,thmumax
*      common /indrates/emuth,thmumax
* model parameters
      real*8 ckms12,ckms23,ckms13,ckmdelta
      real*8 tanbe,mu,m2,m1,m3,ma,mass2u(3),mass2q(3),
     & mass2d(3),mass2l(3),mass2e(3),asoftu(3),asoftd(3),asofte(3)
      common /sckm/ ckms12,ckms23,ckms13,ckmdelta
      common /spar/ tanbe,mu,m2,m1,m3,ma,mass2u,mass2q,mass2d,
     & mass2l,mass2e,asoftu,asoftd,asofte
* useful global variables + roption
      integer lsp,kln
      real*8 pi,s2thw,sinthw,costhw,cosbe,sinbe,
     & cos2be,sin2be,zg,lgzg,delrho,alph3mz,gev2cm3s
      character*5 roption
      common /iuseful/ lsp,kln
      common /ruseful/ pi,s2thw,sinthw,costhw,
     & cosbe,sinbe,cos2be,sin2be,zg,lgzg,delrho,alph3mz,gev2cm3s
      common /cuseful/ roption
* coupling constants
      real*8 g2weak,gyweak,g3stro,alphem,alph3,yukawa(12),
     &     lam1,lam2,lam3,lam4,lam5,lam6,lam7
      common /couplingconstants/ g2weak,gyweak,g3stro,alphem,
     & alph3,yukawa,lam1,lam2,lam3,lam4,lam5,lam6,lam7
* mass spectrum
      real*8 mass(0:50),runmass(0:50)
      real*8 mcmc,mbmb,mtmt
      common /mspctm/ mass,runmass,mcmc,mbmb,mtmt
c* masses at mb scale
c      real*8 mbb,mcb
c      common /massmb/mbb,mcb
* decay widths
      real*8 width(0:50),hdwidth(32,4)
      common /widths/ width,hdwidth
* partial decay widths
      real*8 prtial(54)
      common /partials/ prtial
* vertices
      complex*16  gl(50,50,50),gr(50,50,50)
      common /vrtxs/ gl,gr
* quantum numbers
      real*8 ncolor(12),wiso3(12),echarg(12)
      complex*16 ckm(3,3)
      common /qnum/ ckm,ncolor,wiso3,echarg
      include 'dsio.h'
* msugra variables
      real*8 m0var,mhfvar,a0var,tgbetavar,sgnmuvar
      common/sugrainput/m0var,mhfvar,a0var,tgbetavar,sgnmuvar
* DarkSUSY version and root directory
      character*50 dmversion
      character*128 dmroot
      common/dmver/dmversion,dmroot
* rr
      real*8 ylsp_tau,ylsp_b,ylsp_t
      common/yklsp/ylsp_tau,ylsp_b,ylsp_t
      logical debug
      common/debugging/debug
      integer suwar
      common/suwarning/suwar
      integer aswar
      common/aswarning/aswar
      integer hdecwar
      common/hdecaywarning/hdecwar
      real*8 ftp_s, ftn_s
      common/hadronic/ftp_s, ftn_s
c...Init models
      integer IDmod 
      logical first_idmod
      common /idetec/IDmod, first_idmod
c save common blocks
      save /pacodes/,/switch/,/sckm/,/spar/,/iuseful/,/ruseful/,
     &  /cuseful/,/couplingconstants/,/mspctm/,/widths/,/partials/,
     &  /qnum/,/sugrainput/,
     &  /dmver/,/yklsp/, /debugging/,/suwarning/,/aswarning/,
     &  /hdecaywarning/,/hadronic/,/idetec/
***                                                                 ***
************************* end of susy.h *******************************




