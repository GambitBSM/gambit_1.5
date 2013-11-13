c_______________________________________________________________________
c  abs squared of a complex*16 number.
c  called by dwdcos.
c  author: paolo gondolo (gondolo@lpthe.jussieu.fr) 1994
c=======================================================================
      function dsabsq(z)
      implicit none
      real*8 dsabsq
      complex*16 z
      dsabsq=dreal(z)**2+dimag(z)**2
      end
