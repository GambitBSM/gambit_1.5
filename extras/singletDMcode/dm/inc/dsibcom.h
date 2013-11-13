*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                           dsibcom.h                              ***
***         this piece of code is needed as a separate file          ***
***             the rest of the code 'includes' dsibcom.h            ***
c----------------------------------------------------------------------c
* IB switches
      integer IBflag(1:12),IBerr,IBhow,IBhow_pos
      real*8 IBacc,IBmfr
      common /IBswitches/ IBacc,IBmfr,IBflag,IBerr,IBhow,IBhow_pos

* global variables needed for better performance of integration routines
      integer intch, intyield
      real*8  ibcom_x, ibcom_z, ibcom_mx, ibcom_mp1, ibcom_mp2
      common /IBintvars/ ibcom_x, ibcom_z, ibcom_mx, 
     &                   ibcom_mp1, ibcom_mp2, intch, intyield
