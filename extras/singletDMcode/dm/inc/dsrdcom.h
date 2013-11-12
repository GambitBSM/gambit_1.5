*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dsrdcom.h                               ***
***         this piece of code is needed as a separate file          ***
***            the rest of the code 'includes' dsrdcom.h               ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@cfpa.berkeley.edu) 96-03-22
c  modified: 98-03-03
c....rdmgev - relic and coannihilating masses in gev's
      integer tharsi
c      parameter (tharsi=50)
      parameter (tharsi=1000)  ! New slepton coannihilation value
      real*8 mco(tharsi),mdof(tharsi),rgev(tharsi),rwid(tharsi)
      integer nco,nres
      common /rdmgev/ mco,mdof,rgev,rwid,nco,nres
c....rdpth - momentum threshold to take into account when integrating
      integer nth,incth(0:tharsi+1)
      real*8 pth(0:tharsi+1)
      common /rdpth/ pth,incth,nth
c....rddof - degrees of freedom (3*300 real, 3 integer)
      real*8 tgev(300),fh(300),fg(300)
      integer nf,khi,klo
      common /rddof/ tgev,fh,fg,nf,khi,klo
c....rderrors - error flags (2 integer)
      integer rderr,rdwar,rdinit
      common /rderrors/ rderr,rdwar,rdinit
c....rdlun - logical units (2 integer)
      integer rdlulog,rdluerr
      common /rdlun/ rdlulog,rdluerr
c....rdrate - tabulated invariant rate (3*nrmax real,3*nrmax+3)
      integer nrmax,ispl
      parameter (nrmax=2000) 
      real*8 pp(nrmax),yy(nrmax),yy2(nrmax)
      integer indx(nrmax),nlo,nhi,nr
      common /rdrate/ pp,yy,yy2,indx,ispl(0:tharsi+1),nlo,nhi,nr
c...rdinfo - info
      character*12 rdtag
      common /rdinfo/rdtag
c...rdpars - modifiable parameters
      real*8 cosmin,waccd,dpminr,dpthr,wdiffr,wdifft,
     &  hstep,hmin,compeps,xinit,xfinal,umax,cfr,pmax
      common /rdpars/cosmin,waccd,dpminr,dpthr,wdiffr,wdifft,
     &  hstep,hmin,compeps,xinit,xfinal,umax,cfr,pmax
c...rdpadd - how points are added
      real*8 pdivr,dpres
      integer nlow,nhigh,npres,nthup,cthtest,spltest
      common /rdpadd/pdivr,dpres,nlow,nhigh,npres,nthup,
     &  cthtest,spltest
c...rdlims
      real*8 plow(1:2*tharsi),phigh(1:2*tharsi)
      integer nlim
      common /rdlims/plow,phigh,nlim
c...rdswitch - switches
      integer thavint,rdprt
      common /rdswitch/thavint,rdprt
c....save common blocks
      save /rdmgev/,/rddof/,/rderrors/,/rdlun/,/rdrate/,/rdpth/,
     &  /rdinfo/,/rdpars/,/rdpadd/,/rdlims/,/rdswitch/
***                                                                 ***
************************** end of dsrdcom.h *****************************







