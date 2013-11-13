      function dsrndlog(idum,a,b)
      implicit none
      real*8 dsrndlog,a,b,dsrnd1
      integer idum
      if (a*b.le.0.d0) then
         call dswrite(0,0,'dsrndlog: interval includes zero')
         stop
      endif
      dsrndlog = a * exp( dsrnd1(idum) * log(b/a) )
      return
      end
