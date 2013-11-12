      function dsrndsgn(idum)
      implicit none
      real*8 dsrndsgn,dsrnd1,r
      integer idum
      r = dsrnd1(idum)
      if (r.le.0.5d0) then
         dsrndsgn = -1.d0
      else
         dsrndsgn = +1.d0
      endif
      return
      end
