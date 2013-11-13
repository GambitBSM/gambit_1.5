      character*8 function dsi2s(x)
      implicit none
      integer x
      character*8 a

      write(a,'(i8)') x
      dsi2s=a
      return
      end
