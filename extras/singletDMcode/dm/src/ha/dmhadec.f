*****************************************************************
*** suboutine dmhadec decomposes yieldcode yieldk to flyxtype
*** fltype and fi
*****************************************************************

      subroutine dmhadec(yieldk,fltyp,fi)
      implicit none

      integer yieldk,fltyp,fi

      if (yieldk.lt.100) then
        fltyp=1      ! integrated yields
        fi=yieldk-50
      else
        fltyp=2      ! differential yields
        fi=yieldk-150
      endif

      if (fltyp.ne.1.and.fltyp.ne.2) then
        write(*,*) 'inconsistent yield code: ',yieldk
        stop
      endif

      if (fi.lt.1.or.(fi.gt.8.and.fi.lt.21).or.fi.gt.23) then
        write(*,*) 'inconsistent yield code: ',yieldk
        stop
      endif
      end
