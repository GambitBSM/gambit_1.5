*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dsmucom.h                               ***
***         this piece of code is needed as a separate file          ***
***             the rest of the code 'includes' dsmucom.h              ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@teorfys.uu.se) 96-03-22
c  modified: 00-08-16
c....mfsim - simulation result tables
      real*8 lb,ub,mi,thindex,zindex,dth,dz
      real phiint,phidiff
      integer flxtype,milow,thn,zn
      character*128 mfdir
      character mfftype
      character*4 mfmed
      common/musim/lb(6),ub(6),mi(18),thindex(-1:90,2),
     &  zindex(-1:50,2),dth(-1:90),dz(-1:50),
     &  flxtype(2,3),milow(6),thn,zn,
     &  mfdir,mfftype,mfmed
      common/musim2/ phiint(-1:90,0:50,18,6,2,3),
     &  phidiff(-1:90,-1:50,18,6,2,3)

c....mfpar - variable passing to some routines
      real*8 phim0,phim1,phim2,phie0,phieth,phithm
      integer phich,phiwh,phifk,phifv   ! phiwh=1 - sun, 2 - earth
      common/mupar/phim0,phim1,phim2,phie0,phieth,phithm,phich,
     &  phiwh,phifk,phifv
c...mfinfo - tag etc.
      integer mferr,mfistat
      common/muinfo/mferr,mfistat

c save common block
      save /musim/,/mupar/,/muinfo/,/musim2/
***                                                                 ***
************************** end of muoncom.h ***************************






