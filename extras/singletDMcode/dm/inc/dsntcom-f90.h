!         -*- mode: fortran -*-
!######################################################################*
!                       i n c l u d e     f i l e                      *
!######################################################################*

!************************************************************************
!***                            dsntcom.h                             ***
!***         this piece of code is needed as a separate file          ***
!----------------------------------------------------------------------c
!  author: joakim edsjo (edsjo@physto.se), 2000-08-16

real*8 :: tausu,csu,tauea,cea,ntarateea,ntaratesu,ceadk,ntmx,gtot10
common /ntres/tausu,csu,tauea,cea,ntarateea,ntaratesu,ceadk, &
 ntmx,gtot10

integer :: ntcalcmet,ntlambda,nttab
common /ntpara/ntcalcmet,ntlambda,nttab


!*************************** end of dsntcom.h *****************************






