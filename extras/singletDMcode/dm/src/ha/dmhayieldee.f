c********************************************************************
c astro-ph/0410359, see their eq.(2)

      function dmhayieldee(x,mx,ml)

      real*8 dmhayieldee,x,mx,ml
      real*8 alpha,pi

      pi=3.141592653589793238d0
      alpha=1.d0/128.d0

      if(x.lt.1.d0-1.d-16) then
      dmhayieldee=alpha/pi*(x*x-2.d0*x+2.d0)/x*dlog((1.d0-x)*(mx/ml)**2)
      else
       dmhayieldee=0.d0
      endif

      dmhayieldee=dmhayieldee/mx

      return
      end
c**********************
