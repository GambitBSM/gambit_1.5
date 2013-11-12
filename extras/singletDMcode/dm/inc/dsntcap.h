*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dsntcap.h                             ***
***         this piece of code is needed as a separate file          ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@physto.se), 2003-11-27


c...common block with available tables for capture rates
      integer nc,ntea,ntsu
      parameter(nc=2500,  ! number of entries in the tables
     &          ntea=6,   ! number of earth and sun tables that can be
     &          ntsu=6)   ! loaded simultaneously
      real*8 ctabea(0:nc,ntea),ctabsusi(0:nc,ntsu),ctabsusd(0:nc,ntsu)
      common /dsntcap/ctabea,ctabsusi,ctabsusd
      save /dsntcap/

c...common block with information about what is loaded and with what
      integer nealoaded,nsuloaded
      character*200 fileea(ntea),filesu(ntsu)
      common /dsntcap2/fileea,filesu,nealoaded,nsuloaded
      save /dsntcap2/

************************* end of dsntcap.h *****************************






