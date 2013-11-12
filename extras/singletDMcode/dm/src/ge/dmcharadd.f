***********************************************************************
*** dscharadd takes a string str, adds a string add to it
*** and returns the concatenated string with spaces removed.
*** Author: Joakim Edsjo, edsjo@physto.se
*** Date: 2004-01-19
***********************************************************************

      subroutine dscharadd(str,add)

      character*(*) str
      character*(*) add
 
      character*500 tmp
      integer i,j,tl

      tmp=str//add

      tl=500

      do i=1,tl
 40     if (tmp(i:i).eq.' ') then
          tl=tl-1
          do j=i,tl
            tmp(j:j)=tmp(j+1:j+1)
          enddo
          if (tl.eq.i) goto 50
          goto 40
        endif
      enddo
 50   continue

      str=tmp

      return
      end

