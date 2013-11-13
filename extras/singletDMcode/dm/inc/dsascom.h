*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

**************************************************************
*** include file                                           ***
*** this file contains common blocks for dsas files        ***
*** AUTHOR: Piero Ullio, piero@physto.se                   ***
*** Date: 99-06-22                                         ***
**************************************************************

c family indices (set in dsinit)
      integer ivfam(0:50)
      common /asfamily/ ivfam
c internal family indices
      integer iifam(4)
      common /asfamilyint/ iifam
c mass parameters
      real*8 mass1,mass2,mass3,mass4
      common /asparmass/ mass1,mass2,mass3,mass4
c degrees of freedom
      real*8 gg1,gg2
      common /asdegfree/ gg1,gg2
c symmetry factor of the final state
      integer s34
c and type of fermions in the final state
      integer  q3,q4
      common /aspartype/ q3,q4,s34
c kinematical input values
      real*8 p12,costheta
      common /askin/ p12,costheta
c kinematical values defined in first step
      real*8 Svar,Tvar,Uvar,k34,ep1,ep2,ek3,ek4,scal(4,4)
      common /askinder/ Svar,Tvar,Uvar,k34,ep1,ep2,ek3,ek4,
     & scal
c terms to compute a fermion line
      complex*16 ASxpl(4),ASxpr(4),ASyl,ASyr
      common /asampli/ ASxpl,ASxpr,ASyl,ASyr
c terms to compute a fermion line
      complex*16 ASxplc(6,4),ASxprc(6,4),ASylc(6),ASyrc(6),
     & colfactor(6,6)
      common /asamplic/ ASxplc,ASxprc,ASylc,ASyrc,colfactor
c some roundoff terms
      real*8 thstep,fertoll
      common /astoll/ thstep,fertoll 
