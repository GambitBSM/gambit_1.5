*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dshacom.h                               ***
***         this piece of code is needed as a separate file          ***
***             the rest of the code 'includes' dshacom.h              ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@physics.berkeley.edu) 98-01-26
c  modified: 00-08-16
c....simres - simulation result tables
      real*8 lb,ub,mi,zindex,dz,ndec
      real phiint,phidiff,phidiffPPPC
      integer yieldtype,milow,thn,zn,ntype
      character*128 hadir
      character haftype
      parameter(zn=250,       ! Number of bins in z
     &          ndec=10.0d0,  ! Number of decades tabulated
     &          ntype=23)     ! Code for highest yield type
      common/hasim/lb(8),ub(8),mi(18),
     &  zindex(-1:zn,2),dz(-1:zn),
     &  yieldtype(2,ntype),milow(8),thn,
     &  hadir,haftype
      common/hasim2/phiint(0:zn,18,8,ntype),
     &  phidiff(-1:zn,18,8,ntype),
     &  phidiffPPPC(11160, 14)

c....phi2par - variable passing to some routines
      real*8 phim0,phim1,phim2,phie0,phieth,phigap,phibep,
     &  phie1,phie2,phicthmin,phicthmax,phimp
      integer phich,phifk,phitype,phihno
      common/hapar/phim0,phim1,phim2,phie0,phieth,
     &  phigap,phibep,phie1,phie2,phicthmin,phicthmax,phimp,
     &  phitype,phich,phifk,phihno

c....hainfo - option info switches for the program
      integer haexhi,hasmooth,haerr,haistat
      common /hainfo/haexhi,hasmooth,haerr,haistat

c...hach - channel conversion blocks
c...chcomp converts from new channel numbers to compressed channel numbers
c...Currently these are the old channel numbers, before Pythia runs
c...have been updated
      integer chcomp(30) 
      common /hach/chcomp

c...habranch - annihilation branching rates and Scalar decay rates
c...Also, switches for IB (internal bremsstrahlung) are set here.
      real*8 habr(30),has0br(30,3),hascbr(15),has0m(3),hascm,hamwimp
      character*8 haib
      logical dshasetupcalled
      common /habranch/habr,has0br,hascbr,has0m,hascm,hamwimp,
     &  haib,dshasetupcalled

c...hasim3 - Information about simulation channels
      real*8 map(8),mqu,mqd,mqs
      common /hasim3/map,mqu,mqd,mqs

c save common block
      save /hasim/,/hapar/,/hainfo/,/hasim2/,/hach/,/habranch/,
     & /hasim3/

***                                                                 ***
************************** end of dshacom.h *****************************





