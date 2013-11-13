*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                             dsio.h                               ***
***         this piece of code is needed as a separate file          ***
***               the rest of the code 'includes' dsio.h             ***
c----------------------------------------------------------------------c
c  author: paolo gondolo 1999
* program i/o units
      integer prtlevel,lulog,luerr,luout
      common /dsio/ prtlevel,lulog,luerr,luout
c save common blocks
      save /dsio/
***                                                                 ***
************************* end of dsio.h *******************************
