c     -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                            dsepcom.h                               ***
***         this piece of code is needed as a separate file          ***
***             the rest of the code 'includes' dshacom.h              ***
c----------------------------------------------------------------------c
c---  common blocks included by other positron routines ---
      real*8 l_h,r_e,k27,tau16,n_c,alphexp
      integer dsepdiffmethod
      logical epload,epcl_load
      common /eppar/ l_h,r_e,k27,tau16,n_c,alphexp,dsepdiffmethod,
     &  epload,epcl_load
c
      real*8 k0tau,sumtol,rtol,vtol,rwid,av,ametric,
     +     tableep,table,etable,vtable,wtable
      character*10 epc
      common /dsepcoms/ k0tau,sumtol,rtol,vtol,rwid,av,ametric,
     +     tableep(10,22001),table(22001),etable(13001),vtable(13001),
     +     wtable(13001),epc
c     
      integer npmax
      parameter(npmax=500)
      real*8 dnde(0:npmax),emin,lemin,lemax,dle
      integer npt
      common /dseptable/ dnde,emin,lemin,lemax,dle,npt

      real*8 r_gc
      common /ephaloint/ r_gc
	  
      real*8 mselo(10),msehi(3),
     +     msatab(10),msbtab(10),msctab(10),
     +     mswtab(3),msxtab(3),msytab(3)
      common /epmstab/ mselo,msehi,msatab,msbtab,msctab,
     +     mswtab,msxtab,msytab
      
      save /eppar/,/dsepcoms/,/dseptable/,/ephaloint/,/epmstab/
