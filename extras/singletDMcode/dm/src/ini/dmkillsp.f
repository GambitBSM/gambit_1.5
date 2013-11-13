      function dskillsp(a1,a2)
      implicit none
      character*(*) a1,a2
      integer dskillsp,i,j,l
      l = len(a1)
      j = 1
      do i=1,l
         if (a1(i:i).ne.' ') then
            a2(j:j) = a1(i:i)
            j = j+1
         endif
      enddo
      dskillsp = j-1
      do i=j,len(a2)
         a2(i:i) = '\0'
      enddo
      return
      end
