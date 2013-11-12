      character*12 function dsf2s(x)
      implicit none
      real*8 x
      character*12 a

      write(a,'(e12.6)') x
      dsf2s=a
      return
      end
