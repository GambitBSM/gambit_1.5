*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            idtag.h                               ***
***         this piece of code is needed as a separate file          ***
***               the rest of the code 'includes' idtag.h            ***
***            where reference to a model id tag is needed           ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@teorfys.uu.se) 96-12-12

* model id tag
      character*12 idtag
      common /tag/ idtag

***************************** end **************************************
