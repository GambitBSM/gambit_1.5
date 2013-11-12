*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dshrcom.h                               ***
***         this piece of code is needed as a separate file          ***
***             the rest of the code 'includes' dshrcom.h              ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@physics.berkeley.edu) 98-02-10
c  modified: 98-03-08
      real*8 ehere
      common /hrint/ehere
	  integer hristat
	  common /hrstat/hristat
c save common block
      save /hrint/,/hrstat/
***                                                                 ***
************************** end of dshrcom.h *****************************




