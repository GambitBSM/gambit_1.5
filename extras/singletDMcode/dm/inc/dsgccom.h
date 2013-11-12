c     -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dsgccom.h                             ***
***         this piece of code is needed as a separate file          ***
***             the rest of the code 'includes' dsgccom.h            ***
c----------------------------------------------------------------------c
c---  common blocks included by other galactic-center routines ---
      real*8 mbh,d,vhd,lambda_diff
      common /dsgcpar/ mbh,d,vhd,lambda_diff
c
      real*8 lphibar,RS,lD,lmbh_md,lrhoc,lrhoDbar
      common /dsgcuseful/ lphibar,RS,lD,lmbh_md,lrhoc,lrhoDbar
c
      real*8 lRin,lRspike,lRplateau
      common /dsgcspikepar/lRin,lRspike,lRplateau
c     
      integer npmax
      parameter(npmax=500)
      real*8 ye(0:npmax),emin,lemin,emax,lemax,dle
      integer npt
      common /dsyetable/ ye,emin,lemin,emax,lemax,dle,npt

      save /dsgcpar/,/dsgcuseful/,/dsgcspikepar/,/dsyetable/
