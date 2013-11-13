*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                             dsge.h                               ***
***         this piece of code is needed as a separate file          ***
***               the rest of the code 'includes' susy.h             ***
c----------------------------------------------------------------------c
c  derived from dssusy.h 2008-01-22, paolo gondolo

* useful global variables
      real*8 pi
      real*8 gev2cm3s     ! added by /tb
      common /ruseful/ pi, gev2cm3s
c save common blocks
      save /ruseful/
***                                                                 ***
************************* end of dsge.h *******************************
