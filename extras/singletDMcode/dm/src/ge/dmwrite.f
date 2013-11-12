      subroutine dswrite(level,opt,message)
c_______________________________________________________________________
c  handle writing onto standard output in darksusy
c  input:
c    level - print message if prtlevel is >= level
c    opt   - print (1) the model tag or not (0)
c    message - string containing the text to print
c  common:
c    'dsio.h' - i/o units numbers and prtlevel
c  author: paolo gondolo 1999
c=======================================================================
      implicit none
      include 'dsio.h'
      character*(*) message
      character*12 dsidtag
      integer level,i,dsi_trim,opt
      if (prtlevel.lt.level) return
      i=dsi_trim(message)
      if (opt.eq.1) then
        if (i.gt.0) write (luout,1000) '('//dsidtag()//') ',
     &    message(1:i)
      else
        if (i.gt.0) write (luout,1000) message(1:i)
      endif
 1000 format (1x,a,a)
      return
      end
