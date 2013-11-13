*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dswacom.h                               ***
***         this piece of code is needed as a separate file          ***
***             the rest of the code 'includes' dswacom.h              ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@physto.se) 96-03-22
c  modified: 00-08-16, 08-04-01
c....wasim - simulation result tables
c...wamax is the number of tables to maximally have in memory at the
c...same time. For most users, the default below is OK, but if you
c...want to use all tables simultaneously, wamax has to be increased
c...to 52.
      integer wamax,walast(2)
      parameter(wamax=6) ! number of tables to load simultaneously
      real*8 lb,ub,mi,thindex,zindex,dth,dz
      real phiint,phidiff
      integer milow,thn,zn,yload
      character*128 wadir
      character waftype
      character*40 wabase
      common/wasim/lb(14),ub(14),mi(19),thindex(-1:90,2),
     &  zindex(-1:50,2),dth(-1:90),dz(-1:50),
     &  yload(2,26),walast,
     &  milow(14),thn,zn,
     &  wadir,waftype,wabase
      common/wasim2/ phiint(-1:90,0:50,19,13,2,wamax),
     &  phidiff(-1:90,-1:50,19,13,2,wamax)

c...wagen - general stuff
      integer ch2chi(29),chi2chii(14),chii2chi(13),chi2ch(14)
      common /wagen/ch2chi,chi2chii,chii2chi,chi2ch

c....wapar - variable passing to some routines
      real*8 phim0,phim1,phim2,phie0,phieth,phithm
      integer phichi,phiwh,phifk,phifv   ! phiwh=1 - sun, 2 - earth
      common/wapar/phim0,phim1,phim2,phie0,phieth,phithm,phichi,
     &  phiwh,phifk,phifv
c...wainfo - tag etc.
      integer waerr,waistat
      common/wainfo/waerr,waistat

c...wasim3 - Information about simulation channels
      real*8 map(14)
      common /wasim3/map

c...wabranch - annihilation branching rates and Scalar decay rates
      real*8 wabr(29),was0br(29,3),wascbr(15),was0m(3),wascm,wamwimp
      logical dswasetupcalled
      common /wabranch/wabr,was0br,wascbr,was0m,wascm,wamwimp,
     &  dswasetupcalled

c save common block
      save /wasim/,/wapar/,/wainfo/,/wasim2/,/wasim3/,/wagen/,
     &  /wabranch/
***                                                                 ***
************************** end of muoncom.h ***************************






