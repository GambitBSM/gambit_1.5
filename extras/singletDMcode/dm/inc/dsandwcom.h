*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                          dsandwcom.h                               ***
***         this piece of code is needed as a separate file          ***
***          the rest of the code 'includes' dsandwcom.h               ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@physto.se) 97-09-09
c  modified: 01-09-12
c....dwcom - common block needed for the dwdcos optimization
      real*8 parts(6,6,54)
      logical tur(6,6,54)
      common /dwcom1/ parts,tur

      real*8 mcofr
      integer coann,coproc,slcode,copart
      logical incglue,incgaga,incgaz
      common /dwcom2/ mcofr,coann,incglue,incgaga,incgaz,
     &  coproc,slcode,copart
      save /dwcom1/,/dwcom2/
***                                                                 ***
************************** end of dwcom.h *****************************







