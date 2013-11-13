      function dsfval(a)
      implicit none
      real*8 dsfval,fval
      character*(*) a
      character*80 message
      read (a,*,err=10) fval
      dsfval=fval
      return
 10   write (message,1000) a
 1000 format('error while reading real number ',a)
      call dswrite(0,0,message)
      stop
      end
