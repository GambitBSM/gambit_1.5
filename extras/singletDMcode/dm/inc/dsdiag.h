*         -*- mode: fortran -*-

* diag.h
* global declarations for the Diag routines
* this file is part of Diag
* last modified 23 Nov 07 th
* adapted to darksusy 04 June 08 pg

* The maximum dimension of a matrix, needed for allocating internal
* memory, i.e. the routines handle at most MAXDIM-by-MAXDIM matrices.

      integer MAXDIM
      parameter (MAXDIM=16)


* A matrix is considered diagonal if the sum of the squares
* of the off-diagonal elements is less than EPS.  SYM_EPS is
* half of EPS since only the upper triangle is counted for
* symmetric matrices.
* 52 bits is the mantissa length for IEEE double precision.

      real*8 EPS,SYM_EPS,DBL_EPS
      parameter (EPS=2D0**(-102))
      parameter (SYM_EPS=2D0**(-103))
      parameter (DBL_EPS=2D0**(-52))


