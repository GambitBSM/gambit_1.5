      subroutine dslowcase(a)
      implicit none
      character*(*) a
      character*26 low,upp
      integer i,j,l
      data upp/'abcdefghijklmnopqrstuvwxyz'/
      data low/'abcdefghijklmnopqrstuvwxyz'/
      l = len(a)
      do i=1,l
         j=index(upp,a(i:i))
         if (j.ne.0) a(i:i) = low(j:j)
      enddo
      return
      end
