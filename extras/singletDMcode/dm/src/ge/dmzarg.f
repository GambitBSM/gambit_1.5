************************************************************************
*** dszarg gives the argument (or phase) of a complex number [radians].
************************************************************************
      real*8 function dszarg(z)
      implicit none
      include 'dssusy.h'

      complex*16 z
      real*8 a,b

      a=dreal(z)
      b=dimag(z)

      if (a.eq.0.0d0) then
        if (b.eq.0.0d0) then
          dszarg=0.0d0
        else
          dszarg=pi/2.0d0*sign(1.0d0,b)
        endif
        return
      endif
      
      dszarg=atan(b/a)

      return
      end

