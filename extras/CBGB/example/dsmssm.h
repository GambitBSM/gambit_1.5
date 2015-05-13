*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dsmssm.h                              ***
***         this piece of code is needed as a separate file          ***
***              the rest of the code 'includes' dsmssm.h            ***
c----------------------------------------------------------------------c
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994,1995
c    10-nov-95 complex vertex constants
c  modified by joakim edsjo (edsjo@teorfys.uu.se) 97-02-11
c  derived from dssusy.h 2008-01-22 by paolo gondolo
c  modified by paolo gondolo (paolo@physics.utah.edu) 08-01-30
c  modified by pat scott 08-12-03, 09-10-20, 11-04-27
c   for including running masses

c
c    variables common to both standard and mssm models
c

*  number of particle species
      integer numpartspecies
      parameter (numpartspecies=50)
*   particle codes
* IMPORTANT: the first 12 particle species must be the SM leptons and
* quarks and must be listed in the order below (because they are 
* sometimes referenced directly)
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
      character*8 pname(0:numpartspecies)
      common /pacodes/ kse,ksmu,kstau,ksu,ksd,ksc,kss,kst,ksb,kn,kcha,
     &     knu,kl,kqu,kqd,ksnu,ksl,ksqu,ksqd,pname
* mass spectrum
      real*8 mass(0:numpartspecies)
      real*8 runmass(0:numpartspecies)
      real*8 mu2gev,md2gev,ms2gev,mcmc,mbmb,mtmt
      common /mspctm/ mass,runmass,mu2gev,md2gev,ms2gev,
     & mcmc,mbmb,mtmt
* decay widths
      real*8 width(0:numpartspecies)
      common /widths/ width
* number of spin+color states
      integer kdof(0:numpartspecies)
      common /intdof/ kdof
* vertices
      complex*16 gl(numpartspecies,numpartspecies,numpartspecies),
     &     gr(numpartspecies,numpartspecies,numpartspecies)
      common /vrtxs/ gl,gr
c save common blocks
      save /pacodes/,/mspctm/,/widths/

c
c    standard model variables
c

* useful global variables + roption
      real*8 s2thw,sinthw,costhw,delrho,alph3mz,GFermi
      real*8 s2wmz,swmz,cwmz
      character*5 roption
      common /smruseful/ s2thw,sinthw,costhw,delrho,alph3mz,GFermi,
     &   s2wmz,swmz,cwmz
      common /smcuseful/ roption
* coupling constants
      real*8 g2weak,gyweak,g3stro,alphem,alph3,yukawa(12)
      real*8 g2wmz,gywmz
      common /couplingconstants/ g2weak,gyweak,g3stro,alphem,
     & alph3,yukawa,g2wmz,gywmz
* quantum numbers
      real*8 ncolor(12),wiso3(12),echarg(12)
      common /qnum/ ncolor,wiso3,echarg
* mixings
      real*8 ckms12,ckms23,ckms13,ckmdelta
      complex*16 ckm(3,3)
      common /sckm/ ckms12,ckms23,ckms13,ckmdelta
      common /mixing/ ckm
c save common blocks
      save /sckm/,/smruseful/,/smcuseful/,
     &  /couplingconstants/,
     &  /qnum/,/mixing/

c
c    minimal supersymmetric model variables
c

* type of model (MSSM-7, mSUGRA, etc.) (same meaning as ModSel_Model in SLHA)
* 0=full MSSM, 1=mSUGRA
      integer modeltype
      common /mssmtype/ modeltype
* model parameters
      real*8 tanbe,mu,m2,m1,m3,ma,mass2u(3),mass2q(3),
     & mass2d(3),mass2l(3),mass2e(3),asoftu(3),asoftd(3),asofte(3)
      common /mssmpar/ tanbe,mu,m2,m1,m3,ma,mass2u,mass2q,mass2d,
     & mass2l,mass2e,asoftu,asoftd,asofte
* program switches
      integer higloop,neuloop,bsgqcd,higwid
      real*8 msquarks,msleptons
      common /mssmswitch/ msquarks,msleptons,
     &     higloop,neuloop,bsgqcd,higwid
* set sfermion masses by hand
      real*8 massup1(3),massup2(3),thetamixu(3),
     & massdn1(3),massdn2(3),thetamixd(3),
     & masssn(3),masssl1(3),masssl2(3),thetamixsl(3)
      common/sfermionmass/ massup1,massup2,thetamixu,
     & massdn1,massdn2,thetamixd,
     & masssn,masssl1,masssl2,thetamixsl
* useful global variables + roption
      integer lsp,kln
      real*8 cosbe,sinbe,cos2be,sin2be,zg,lgzg
      common /mssmiuseful/ lsp,kln
      common /mssmruseful/ cosbe,sinbe,cos2be,sin2be,zg,lgzg
* coupling constants
      real*8 lam1,lam2,lam3,lam4,lam5,lam6,lam7
      common /mssmcouplingconstants/ lam1,lam2,lam3,lam4,lam5,lam6,lam7
* decay widths
      real*8 hdwidth(32,4)
      common /mssmwidths/ hdwidth
* mixings
      real*8 alpha,mix_stop,mix_sbot,mix_stau
      complex*16 neunmx(4,4),chaumx(2,2),chavmx(2,2),
     & slulmx(3,3),sldlmx(6,3),sldrmx(6,3),
     & squlmx(6,3),squrmx(6,3),sqdlmx(6,3),sqdrmx(6,3)
      common /mssmmixing/ neunmx,chaumx,chavmx,
     & slulmx,sldlmx,sldrmx,
     & squlmx,squrmx,sqdlmx,sqdrmx,alpha,mix_stop,mix_sbot,mix_stau
* msugra variables
      real*8 m0var,mhfvar,a0var,tgbetavar,sgnmuvar
      common/sugrainput/m0var,mhfvar,a0var,tgbetavar,sgnmuvar
* flag indicating if first or later call of dsralph34loop
      logical first_dsralph34loop
      common/sufirsts/first_dsralph34loop
* cross section ratios for HiggsBounds
c no longer needed (PG 20110606)
c      real*8 xsec_lep_hjZ_ratio(3),xsec_lep_hjhi_ratio(3,3)
c      common /xseclepratio/ xsec_lep_hjZ_ratio,
c     &                      xsec_lep_hjhi_ratio
c save common blocks
      save /mssmpar/,/mssmswitch/,/mssmiuseful/,/mssmruseful/,
     &  /mssmcouplingconstants/,
     &  /mssmmixing/,/sugrainput/,/sufirsts/

***                                                                 ***
************************ end of dsmssm.h ******************************
