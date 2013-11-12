      function dscval(a)
      implicit none
      character*80 dscval
      character*(*) a
      integer i,j
      i = index(a,'''')
      j = index(a(i+1:),'''')
      dscval=a(i+1:j)
      return
      end
