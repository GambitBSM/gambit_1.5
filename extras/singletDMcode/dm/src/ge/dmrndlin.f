      function dsrndlin(idum,a,b)
      implicit none
      real*8 dsrndlin,a,b,dsrnd1
      integer idum
      dsrndlin = a + (b-a) * dsrnd1(idum)
      return
      end
