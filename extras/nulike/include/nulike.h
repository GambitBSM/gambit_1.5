!*         -*- mode: fortran -*-
!************************************************************************
!***                           nucommon.h                             ***
!----------------------------------------------------------------------c
!  author: Pat Scott (patscott@physics.mcgill.ca), June 12, 2014

      include "nuconst.h"

      interface
        subroutine nulike_bounds(analysis_name_in, mwimp, annrate, 
     &   nuyield, Nsignal_predicted, NBG_expected, Ntotal_observed, 
     &   lnlike, pvalue, liketype, pvalFromRef, referenceLike, dof,
     &   context)  BIND(C)
          use iso_c_binding, only: c_ptr, c_char
          implicit none
          include "nuconst.h"
          integer Ntotal_observed, liketype
          real*8 Nsignal_predicted, NBG_expected, lnlike, pvalue, referenceLike,
     &           dof, mwimp, annrate, nuyield
          logical pvalFromRef
          character(kind=c_char), dimension(nulike_clen) :: analysis_name_in
          type(c_ptr) :: context
          external nuyield
        end subroutine
      end interface

!*************************** end of nulike.h *****************************

