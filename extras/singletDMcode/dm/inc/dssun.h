*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dssun.h                               ***
***         this piece of code is needed as a separate file          ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@physto.se), 2003-11-24, 2005-11-24

      integer sdmax,sdn,sdnne
      parameter(sdmax=3000)
      real*8 sdr(sdmax),sdm(sdmax),sdrho(sdmax),sdphi(sdmax),
     &  sdmfr(16,sdmax),sdcdens(sdmax,0:2)
      real*8 sdrne(sdmax),sdne(sdmax)
      real*8 sdabund(16),sdaa(16),sdma(16),r_sun,m_sun,cd_sun(0:2)
      common /dssun/sdr,sdm,sdrho,sdphi,sdmfr,sdcdens,
     &  sdabund,sdaa,sdma,
     &  r_sun,m_sun,cd_sun,sdrne,sdne,sdn,sdnne
      save /dssun/

      character cdt
      common /dssuncd/cdt
      save /dssuncd/


*************************** end of dssun.h *****************************






