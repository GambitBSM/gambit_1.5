*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                           diacom.h                               ***
***         this piece of code is needed as a separate file          ***
***            the rest of the code 'includes' diacom.h              ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@teorfys.uu.se) 95-10-20
c  modified: 96-03-23
      complex*16 aa(0:1,-1:1,-1:1,-1:1)
      real*8 mp1,mp2,mk,mp3,mp4,ppl,pmi,epl,emi,mw,mz,pp,kpl,
     &  kmi,efpl,efmi,kk,e1,e2,e3,e4,s,t,u,wd(2,-2:2,-2:2)
      common /diacom/ aa,
     &  mp1,mp2,mk,mp3,mp4,ppl,pmi,epl,emi,mw,mz,pp,kpl,kmi,
     &  efpl,efmi,kk,e1,e2,e3,e4,s,t,u,wd
c save common block
      save /diacom/
***                                                                 ***
************************* end of diacom.h *****************************


