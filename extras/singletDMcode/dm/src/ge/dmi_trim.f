      function dsi_trim(s)
      implicit none
      integer dsi_trim,n
      character*(*) s
      n=len(s)
      if (n.eq.0) goto 20
 10   if (s(n:n).eq.' ') then
         n=n-1
         if (n.gt.0) goto 10
      endif
 20   dsi_trim=n
      return
      end
