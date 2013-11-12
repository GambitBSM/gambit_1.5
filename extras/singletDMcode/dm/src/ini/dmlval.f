      function dslval(a)
      implicit none
      logical dslval,lval
      character*(*) a
      character*80 message
      read (a,*,err=10) lval
      dslval=lval
      return
 10   write (message,1000) a
 1000 format('error while reading logical ',a)
      call dswrite(0,0,message)
      stop
      end
