      function dsival(a)
      implicit none
      integer dsival,ival
      character*(*) a
      character*80 message
      read (a,*,err=10) ival
      dsival=ival
      return
 10   write (message,1000) a
 1000 format ('error while reading integer number ',a)
      call dswrite(0,0,message)
      stop
      end
