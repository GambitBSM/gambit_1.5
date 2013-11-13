c     -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                       dsgalpropcom.h                             ***
***       this piece of code is needed as a separate file            ***
***        the rest of the code 'includes' dsgalpropcom.h            ***
c----------------------------------------------------------------------c
c---  author: e.a.baltz 2/21/2006
c---  common blocks included by galprop routines ---
      real*8 kegpgf(1000,1000),kepgpgf(1000,1000),epgpgf(1000,1000),
     &     pbgpgf(1000,1000),dbgpgf(1000,1000)
      integer gpgfunit,gpnumdecade,gpnumin,gpnumout
      character*200 gpgffile,gpgaldeffile,gpgaldeftmpfile
      character*4 gpmodtag
      logical gpgfread
      common /galpropgreen/kegpgf,kepgpgf,epgpgf,pbgpgf,dbgpgf,
     +     gpgfunit,gpnumdecade,gpnumin,gpnumout,
     +     gpgffile,gpgaldeffile,gpgaldeftmpfile,gpgfread,gpmodtag
      
      save /galpropgreen/
