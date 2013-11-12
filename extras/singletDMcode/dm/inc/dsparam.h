*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dsparam.h                             ***
***         this piece of code is needed as a separate file          ***
*** included here are parameters needed for several other codes      ***
c----------------------------------------------------------------------c
c  author: joakim edsjo (edsjo@physto.se), 2006-03-21

      real*8 m_p,m_n,N_a,c_0
      parameter(
     &  m_p=0.938271998d0,  ! proton mass [GeV]
     &  m_n=0.9396d0,       ! neutron mass [GeV]
     &  N_a=6.022d23,       ! Avogradros number [number / mol ]
     &  c_0=299792.458d0)   ! speed of light [km/s]

*************************** end of dsparam.h *****************************






