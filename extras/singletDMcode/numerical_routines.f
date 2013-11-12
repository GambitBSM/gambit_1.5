      FUNCTION factln(n)
      INTEGER n
      REAL factln
CU    USES gammln
      REAL a(100),gammln
      SAVE a
      DATA a/100*-1./
      if (n.lt.0) pause 'negative factorial in factln'
      if (n.le.99) then
        if (a(n+1).lt.0.) a(n+1)=gammln(n+1.)
        factln=a(n+1)
      else
        factln=gammln(n+1.)
      endif
      return
      END

      function dirichlet(x1,k)
      real*8 :: x1
      integer :: k
      dirichlet = 1.
      end function

      FUNCTION factrl(n)
      INTEGER n
      REAL factrl
CU    USES gammln
      INTEGER j,ntop
      REAL a(33),gammln
      SAVE ntop,a
      DATA ntop,a(1)/0,1./
      if (n.lt.0) then
        pause 'negative factorial in factrl'
      else if (n.le.ntop) then
        factrl=a(n+1)
      else if (n.le.32) then
        do j=ntop+1,n
          a(j+1)=j*a(j)
       end do
        ntop=n
        factrl=a(n+1)
      else
        factrl=exp(gammln(n+1.))
      endif
      return
      END


      FUNCTION dfridr(func,x,h,err)
      INTEGER NTAB
      REAL*8 dfridr,err,h,x,func,CON,CON2,BIG,SAFE
      PARAMETER (CON=1.4,CON2=CON*CON,BIG=1.E30,NTAB=10,SAFE=2.)
      EXTERNAL func
!CU    USES func
      INTEGER i,j
      REAL*8 errt,fac,hh,a(NTAB,NTAB)
      if(h.eq.0.) pause 'h must be nonzero in dfridr'
      hh=h
      a(1,1)=(func(x+hh)-func(x-hh))/(2.0*hh)
      err=BIG
      do 12 i=2,NTAB
        hh=hh/CON
        a(1,i)=(func(x+hh)-func(x-hh))/(2.0*hh)
        fac=CON2
        do 11 j=2,i
          a(j,i)=(a(j-1,i)*fac-a(j-1,i-1))/(fac-1.)
          fac=CON2*fac
          errt=max(abs(a(j,i)-a(j-1,i)),abs(a(j,i)-a(j-1,i-1)))
          if (errt.le.err) then
            err=errt
            dfridr=a(j,i)
          endif
11      continue
        if(abs(a(i,i)-a(i-1,i-1)).ge.SAFE*err)return
12    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software V,3.

c---------------------------

        function rombint(f,a,b,tol)
c  Rombint returns the integral from a to b of using Romberg integration.
c  The method converges provided that f(x) is continuous in (a,b).
c  f must be double precision and must be declared external in the calling
c  routine.  tol indicates the desired relative accuracy in the integral.
c
        parameter (MAXITER=40,MAXJ=5)
        implicit double precision (a-h,o-z)
        dimension g(MAXJ+1)
        double precision f
        external f
c        
        h=0.5d0*(b-a)
        gmax=h*(f(a)+f(b))
        g(1)=gmax
        nint=1
        error=1.0d20
        i=0
10        i=i+1
          if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol))
     2      go to 40
c  Calculate next trapezoidal rule approximation to integral.
          g0=0.0d0
            do 20 k=1,nint
            g0=g0+f(a+(k+k-1)*h)
20        continue
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1.0d0
            do 30 j=1,jmax
c  Use Richardson extrapolation.
            fourj=4.0d0*fourj
            g1=g0+(g0-g(j))/(fourj-1.0d0)
            g(j)=g0
            g0=g1
30        continue
          if (abs(g0).gt.tol) then
            error=1.0d0-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        go to 10
40      rombint=g0
        if (i.gt.MAXITER.and.abs(error).gt.tol)
     2    write(*,*) 'Rombint failed to converge; integral, error=',
     3    rombint,error
        return
        end

c---------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Fifth order Runge-Kutta Method with Adaptive Stepsize.
C     Integrate 'func' with parameters array fp(np) which 
C     contains any extra parameters other than the integration variable
C     from a to b, with initial step size dxinit and fractional accuracy
C     eps.
C
C     In other words,
C          _b
C         /
C        |
C        |  FUNC(x,fp)dx
C        |
C        /
C       _ a
C
C     fp of length np, (i.e. real*8 fp(np)) contains all variables other
C     than the integration variable, say, x.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      real*8 function INTEGRATE(FUNC,fp,np,a,b,dxinit,eps)
      
      implicit none
      integer np 
      real*8 FUNC
      real*8 a, b, eps, dxinit, fp(np)
      external FUNC

      integer maxsteps
      parameter(maxsteps=1.d8)

      real*8  x, dx, dxnext, y, dydx, yscale, ds2dm
      integer  Nstep
      external ds2dm

      x     = a
      dx    = dxinit
      y     = 0.d0
      Nstep = 0

      do while ((x-b)*(b-a).lt.0.d0.and.Nstep.lt.maxsteps)
        Nstep = Nstep + 1
        dydx = FUNC(x,fp,np)

        yscale = max(abs(y) + abs(dx*dydx), 1.d-12)
        if ((x+dx-b)*(x+dx-a).gt.0.d0)  ! If stepsize overshoots, decrease it.
     1    dx = b - x

        call RUNGE5VAR(y,dydx,x,dx,eps,yscale,dxnext,FUNC,fp,np)

        dx = dxnext
      enddo

      if (Nstep.ge.maxsteps)
     1  write (*,*) 'WARNING: failed to converge in INTEGRATE.'

      INTEGRATE = y

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Individual Step: From J. S. Bullock
c
c     Fifth-order Runge-Kutta step with monitoring of local truncation error
c     to ensure accuracy and adjust stepsize.  Input are the dependent
c     variable y and its derivative dydx at the starting value of the
c     independent variable x.  Also input are the stepsize to be attempted
c     htry, the required accuracy eps, and the value yscale, against which the
c     error is scaled.  On output, y and x are replaced by their new values.
c     hdid is the stepsize that was actually accomplished, and hnext is the
c     estimated next stepsize.  DERIVS is the user-supplied routine that
c     computes right-hand-side derivatives.  The argument fparm is for an
c     optional second argument to DERIVS (NOT integrated over).
c
c
c     by M.A.K. Gross
c
      SUBROUTINE RUNGE5VAR(y,dydx,x,htry,eps,yscale,hnext,DERIVS,
     1                     fp,np)
      implicit none

      integer np 
      real*8 fp(np)
      real*8 eps,hnext,htry,x,dydx,y,yscale,DERIVS
      external DERIVS

      real*8 errmax,h,hold,htemp,xnew,yerr,ytemp

      real*8 safety,pgrow,pshrink,errcon
      parameter (safety  =  0.9d0)
      parameter (pgrow   = -0.2d0)
      parameter (pshrink = -0.25d0)
      parameter (errcon  =  1.89d-4)

      h = htry                         ! Set stepsize to initial accuracy.
      errmax = 10.d0
      do while (errmax.gt.1.d0)
        call RUNGE(y,dydx,x,h,ytemp,yerr,DERIVS,fp,np)

        errmax = abs(yerr/yscale)/eps   ! Scale relative to required accuracy.
        if (errmax.gt.1.d0) then        ! Truncation error too large; reduce h
          htemp = safety*h*(errmax**pshrink)
          hold = h
          h = sign(max(abs(htemp),0.1d0*abs(h)),h)  ! No more than factor of 10
          xnew = x + h
          if (xnew.eq.x) then
            write (*,*) 'WARNING: ',
     1                  'Stepsize underflow in RUNGE5VAR().'
            h = hold
            errmax = 0.d0
          endif
        endif
      enddo
c
c     Step succeeded.  Compute estimated size of next step.
c
      if (errmax.gt.errcon) then
        hnext = safety*h*(errmax**pgrow)
      else
        hnext = 5.d0 * h                ! No more than factor of 5 increase.
      endif
      x = x + h

      y = ytemp

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Integrator Step: From J.S.Bullock.
c
c     Given values for a variable y and its derivative dydx known at x, use
c     the fifth-order Cash-Karp Runge-Kutta method to advance the solution
c     over an interval h and return the incremented variables as yout.  Also
c     return an estimate of the local truncation error in yout using the
c     embedded fourth order method.  The user supplies the routine
c     DERIVS(x,y,dydx), which returns derivatives dydx at x.
c
c     by M.A.K. Gross
c
      SUBROUTINE RUNGE(y,dydx,x,h,yout,yerr,DERIVS,fp,np)
      implicit none

      integer np 
      real*8 h,x,dydx,y,yerr,yout,DERIVS,fp(np)

      external DERIVS

      real*8 ak3, ak4, ak5 ,ak6

      real*8 a2,a3,a4,a5,a6
      real*8 c1,c3,c4,c6,dc1,dc3,dc4,dc5,dc6
      parameter(a2  =    0.2d0)
      parameter(a3  =    0.3d0)
      parameter(a4  =    0.6d0)
      parameter(a5  =    1.d0)
      parameter(a6  =    0.875d0)
      parameter(c1  =   37.d0/378.d0)
      parameter(c3  =  250.d0/621.d0)
      parameter(c4  =  125.d0/594.d0)
      parameter(c6  =  512.d0/1771.d0)
      parameter(dc1 = c1 -  2825.d0/27648.d0)
      parameter(dc3 = c3 - 18575.d0/48384.d0)
      parameter(dc4 = c4 - 13525.d0/55296.d0)
      parameter(dc5 = -277.d0/14336.d0)
      parameter(dc6 = c6 -     0.25d0)

      ak3 = DERIVS(x+a3*h,fp,np)
      ak4 = DERIVS(x+a4*h,fp,np)
      ak5 = DERIVS(x+a5*h,fp,np)
      ak6 = DERIVS(x+a6*h,fp,np)
c
c     Estimate the fifth order value.
c
      yout = y + h*(c1*dydx + c3*ak3 + c4*ak4  + c6*ak6)
c
c     Estimate error as difference between fourth and fifth order
c
      yerr = h*(dc1*dydx + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6)

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

***********************
c-----------------------------------------------------
c-----------------------------------------------------
cc Spline fit subroutines

      SUBROUTINE spline_v(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      real*8 yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=100010)
      INTEGER i,k
      real*8 p,qn,sig,un,u(NMAX)
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
C  (C) Copr. 1986-92 Numerical Recipes Software =$j*m,).
      END
c--------------------------------------------------
      Subroutine splint_v(xa,ya,y2a,n,x,y)
      Integer n
      real*8 x,y,xa(n),y2a(n),ya(n)
      Integer k,khi,klo
      real*8 a,b,h

      klo=1
      khi=n
 1    if(khi-klo .gt. 1)then
       k=(khi+klo)/2
       if(xa(k) .gt. x)then
        khi=k
       else
        klo=k
       endif
       goto 1
      endif
      h=xa(khi)-xa(klo)
      if(h .eq. 0.)pause'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     &       ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end
c---------------------------------------------------
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
c-----------------------------------------------------------
      FUNCTION gasdev(idum)
      INTEGER idum
      REAL*8 gasdev
CU    USES ran1
      INTEGER iset
      REAL*8 fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.
c---------------------------------------------------
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.
c--------------------------------------------------------
      SUBROUTINE choldc(a,n,np,p)
      INTEGER n,np
      REAL*8 a(np,np),p(n)
      INTEGER i,j,k
      REAL*8 sum
      do 13 i=1,n
        do 12 j=i,n
          sum=a(i,j)
          do 11 k=i-1,1,-1
            sum=sum-a(i,k)*a(j,k)
11        continue
          if(i.eq.j)then
            if(sum.le.0.)pause 'choldc failed'
            p(i)=sqrt(sum)
          else
            a(j,i)=sum/p(i)
          endif
12      continue
13    continue
      return
      END
c---------------------------------------------------

      SUBROUTINE jacobi(a,n,np,d,v,nrot)
      INTEGER n,np,nrot,NMAX
      REAL*8 a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j
      REAL*8 c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+
     *g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      pause 'too many iterations in jacobi'
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.
c---------------------
	SUBROUTINE moment(data1,n,ave,adev,sdev,var,skew,curt)
	INTEGER n
	REAL*8 adev,ave,curt,sdev,skew,var,data1(n)
	INTEGER j
	REAL*8 p,s,ep
	if(n.le.1)pause 'n must be at least 2 in moment'
	s=0.
	do j=1,n
	   s=s+data1(j)
	enddo
	ave=s/n
	adev=0.
	var=0.
	skew=0.
	curt=0.
	ep=0.
	do j=1,n
	   s=data1(j)-ave
	   ep=ep+s
	   adev=adev+abs(s)
	   p=s*s
	   var=var+p
	   p=p*s
	   skew=skew+p
	   p=p*s
	   curt=curt+p
	enddo
	adev=adev/n
	var=(var-ep**2/n)/(n-1)
	sdev=sqrt(var)
	if(var.ne.0.)then
	   skew=skew/(n*sdev**3)
	   curt=curt/(n*var**2)-3.
	else
	   pause 'no skew or kurtosis when zero variance in moment'
	endif
	return
	END
C       (C) Copr. 1986-92 Numerical Recipes Software V,3.
	
c------------------------------------------
      FUNCTION gammp(a,x)
      REAL*8 a,gammp,x
CU    USES gcf,gser
      REAL*8 gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammp'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.
c------------------------------------------------
      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
CU    USES gammln
      INTEGER i
      REAL*8 an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.
c-------------------------------------------------
      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
CU    USES gammln
      INTEGER n
      REAL*8 ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.
c-----------------------------------------------
      FUNCTION gammln(xx)
      REAL*8 gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
      end do 
      gammln=tmp+log(stp*ser/x)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.
c--------------------------------------------------


      FUNCTION rtbis(func,x1,x2,xacc)
      INTEGER JMAX
      REAL*8 rtbis,x1,x2,xacc,func
      EXTERNAL func
      PARAMETER (JMAX=200)
      INTEGER j
      REAL*8 dx,f,fmid,xmid
      fmid=func(x2)
      f=func(x1)
      if(f*fmid.ge.0.) pause 'root must be bracketed in rtbis'
      if(f.lt.0.)then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbis+dx
        fmid=func(xmid)
        if(fmid.le.0.)rtbis=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
11    continue
      pause 'too many bisections in rtbis'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.

c------------------------------------------

      SUBROUTINE zbrac(func,x1,x2,succes)
      INTEGER NTRY
      REAL*8 x1,x2,func,FACTOR
      EXTERNAL func
      PARAMETER (FACTOR=1.6,NTRY=1000)
      INTEGER j
      REAL*8 f1,f2
      LOGICAL succes
      if(x1.eq.x2)pause 'you have to guess an initial range in zbrac'
      f1=func(x1)
      f2=func(x2)
      succes=.true.
      do 11 j=1,NTRY
        if(f1*f2.lt.0.)return
        if(abs(f1).lt.abs(f2))then
          x1=x1+FACTOR*(x1-x2)
          f1=func(x1)
        else
          x2=x2+FACTOR*(x2-x1)
          f2=func(x2)
        endif
11    continue
      succes=.false.
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.

c--------------------------------------

      SUBROUTINE zbrak(fx,x1,x2,n,xb1,xb2,nb)
      INTEGER n,nb
      REAL*8 x1,x2,xb1(nb),xb2(nb),fx
      EXTERNAL fx
      INTEGER i,nbb
      REAL*8 dx,fc,fp,x
      nbb=0
      x=x1
      dx=(x2-x1)/n
      fp=fx(x)
      do 11 i=1,n
        x=x+dx
        fc=fx(x)
        if(fc*fp.lt.0.) then
          nbb=nbb+1
          xb1(nbb)=x-dx
          xb2(nbb)=x
          if(nbb.eq.nb)goto 1
        endif
        fp=fc
11    continue
1     continue
      nb=nbb
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.

c--------------------------------------------

      FUNCTION zbrent(func,x1,x2,tol)
      INTEGER ITMAX
      REAL*8 zbrent,tol,x1,x2,func,EPS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      REAL*8 a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))pause
     *'root must be bracketed for zbrent'
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b)
11    continue
      pause 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.

c------------------------------
      SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      real*8 arr(n)
      PARAMETER (M=7,NSTACK=50)

      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      real*8 a,temp
      jstack=0
      l=1
      ir=n
          
 1    if(ir-l.lt.M)then
         do j=l+1,ir
            a=arr(j)
            do i=j-1,l,-1
               if(arr(i).le.a)goto 2
               arr(i+1)=arr(i)
            enddo 
            i=l-1
 2          arr(i+1)=a
         enddo

         if(jstack.eq.0)return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         temp=arr(k)
         arr(k)=arr(l+1)
         arr(l+1)=temp
         if(arr(l).gt.arr(ir))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
         endif
         if(arr(l+1).gt.arr(ir))then
            temp=arr(l+1)
            arr(l+1)=arr(ir)
            arr(ir)=temp
         endif
         if(arr(l).gt.arr(l+1))then
            temp=arr(l)
            arr(l)=arr(l+1)
            arr(l+1)=temp
         endif
         i=l+1
         j=ir
         a=arr(l+1)
 3       continue
         i=i+1
         if(arr(i).lt.a)goto 3
 4       continue
         j=j-1
         if(arr(j).gt.a)goto 4
         if(j.lt.i)goto 5
         temp=arr(i)
         arr(i)=arr(j)
         arr(j)=temp
         goto 3
 5       arr(l+1)=arr(j)
         arr(j)=a
         jstack=jstack+2
         if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
         if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
      goto 1
      END

c------------------------------------------

      SUBROUTINE sort3(n,ra,rb,rc,wksp,iwksp)
      INTEGER n,iwksp(n)
      real*8 ra(n),rb(n),rc(n),wksp(n)
CU    USES indexx
      INTEGER j
      call indexx(n,ra,iwksp)
      do 11 j=1,n
        wksp(j)=ra(j)
11    continue
      do 12 j=1,n
        ra(j)=wksp(iwksp(j))
12    continue
      do 13 j=1,n
        wksp(j)=rb(j)
13    continue
      do 14 j=1,n
        rb(j)=wksp(iwksp(j))
14    continue
      do 15 j=1,n
        wksp(j)=rc(j)
15    continue
      do 16 j=1,n
        rc(j)=wksp(iwksp(j))
16    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.

c-------------------------------------------

      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      real*8 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      real*8 a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.

c-----------------------------------------------------

      SUBROUTINE qtrap(func,a,b,s)
      INTEGER JMAX
      REAL a,b,func,s,EPS
      EXTERNAL func
      PARAMETER (EPS=1.e-6, JMAX=20)
CU    USES trapzd
      INTEGER j
      REAL olds
      olds=0.
      do 11 j=1,JMAX
        call trapzd(func,a,b,s,j)
        if (j.gt.5) then
          if (abs(s-olds).lt.EPS*abs(olds).or.(s.eq.0..and.olds.eq.0.)) 
     *return
        endif
        olds=s
11    continue
      pause 'too many steps in qtrap'
      END

      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END

      FUNCTION erf_v(x)
      REAL*8 erf_v,x
CU    USES gammp
      REAL*8 gammp
      if(x.lt.0.)then
        erf_v=-gammp(.5d0,x**2)
      else
        erf_v=gammp(.5d0,x**2)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.

cc--------------------------------------------------------------

      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      REAL*8 a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
c     Linear equation solution by Gauss-Jordan elimination, equation 
c     (2.1.1) above. a(1:n,1:n) is an input matrix stored in an array 
c     of physical dimensions np by np. b(1:n,1:m) is an input matrix 
c     containing the m right-hand side vectors, stored in an array of 
c     physical dimensions np by mp. On output, a(1:n,1:n) is replaced 
c     by its matrix inverse, and b(1:n,1:m) is replaced by the 
c     corresponding set of solution vectors.
c     Parameter: NMAX is the largest anticipated value of n.

      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)    
c     The integer arrays ipiv, indxr, and indxc are used
      REAL*8 big,dum,pivinv    !for bookkeeping on the pivoting.

		do j=1,n
         ipiv(j)=0
		enddo 
         do i=1,n   
            big=0.   
            do j=1,n    
               if(ipiv(j).ne.1)then    
                  do k=1,n
                     if (ipiv(k).eq.0) then
                        if (abs(a(j,k)).ge.big)then
                           big=abs(a(j,k))
                           irow=j
                           icol=k
                        endif
                     endif
                  enddo 
               endif
            enddo 
            ipiv(icol)=ipiv(icol)+1
c     We now have the pivot element, so we interchange rows, if needed, 
c     to put the pivot element on the diagonal. The columns are not 
c     physically interchanged, only relabeled:
c     indxc(i), the column of the ith pivot element, is the ith column 
c     that is reduced, while indxr(i) is the row in which that pivot 
c     element was originally located. If indxr(i) =indxc(i) there is an 
c     implied column interchange. With this form of bookkeeping, the
c     solution b's will end up in the correct order, and the inverse 
c     matrix will be scrambled by columns.
c
            if (irow.ne.icol) then
            do l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
         enddo 
         do l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
         enddo 
      endif
      indxr(i)=irow            
      indxc(i)=icol    
      if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
      pivinv=1./a(icol,icol)
      a(icol,icol)=1.
      do l=1,n
         a(icol,l)=a(icol,l)*pivinv
      enddo 
      do l=1,m
         b(icol,l)=b(icol,l)*pivinv
      enddo 
      do  ll=1,n    !Next, we reduce the rows...
         if(ll.ne.icol)then     !...except for the pivot one, of course.
            dum=a(ll,icol)
            a(ll,icol)=0.
            do l=1,n
               a(ll,l)=a(ll,l)-a(icol,l)*dum
            enddo 
            do l=1,m
               b(ll,l)=b(ll,l)-b(icol,l)*dum
            enddo 
         endif
      enddo 
      enddo    
      do l=n,1,-1    
         if(indxr(l).ne.indxc(l))then    
            do k=1,n    
               dum=a(k,indxr(l))   
               a(k,indxr(l))=a(k,indxc(l))
               a(k,indxc(l))=dum
            enddo 
         endif
      enddo 
      return    !And we are done.
      END
	
c----------------------------------------------

      SUBROUTINE vegas(region,ndim,fxn,init,ncall,itmx,nprn,tgral,sd,
     *chi2a)
      INTEGER init,itmx,ncall,ndim,nprn,NDMX,MXDIM
      real*8 tgral,chi2a,sd,region(2*ndim),fxn,ALPH,TINY
      PARAMETER (ALPH=1.5,NDMX=50,MXDIM=10,TINY=1.e-30)
      EXTERNAL fxn
CU    USES fxn,ran2,rebin
      INTEGER i,idum,it,j,k,mds,nd,ndo,ng,npg,ia(MXDIM),kg(MXDIM)
      real*8 calls,dv2g,dxg,f,f2,f2b,fb,rc,ti,tsi,wgt,xjac,xn,xnd,xo,
     *d(NDMX,MXDIM),di(NDMX,MXDIM),dt(MXDIM),dx(MXDIM),r(NDMX),x(MXDIM),
     *xi(NDMX,MXDIM),xin(NDMX),ran2
      DOUBLE PRECISION schi,si,swgt
      COMMON /ranno/ idum
!$omp threadprivate(/ranno/)
      SAVE
      if(init.le.0)then
        mds=1
        ndo=1
        do 11 j=1,ndim
          xi(1,j)=1.
11      continue
      endif
      if (init.le.1)then
        si=0.d0
        swgt=0.d0
        schi=0.d0
      endif
      if (init.le.2)then
        nd=NDMX
        ng=1
        if(mds.ne.0)then
          ng=(ncall/2.+0.25)**(1./ndim)
          mds=1
          if((2*ng-NDMX).ge.0)then
            mds=-1
            npg=ng/NDMX+1
            nd=ng/npg
            ng=npg*nd
          endif
        endif
        k=ng**ndim
        npg=max(ncall/k,2)
        calls=float(npg)*float(k)
        dxg=1./ng
        dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-1.)
        xnd=nd
        dxg=dxg*xnd
        xjac=1./calls
        do 12 j=1,ndim
          dx(j)=region(j+ndim)-region(j)
          xjac=xjac*dx(j)
12      continue
        if(nd.ne.ndo)then
          do 13 i=1,max(nd,ndo)
            r(i)=1.
13        continue
          do 14 j=1,ndim
            call rebin(ndo/xnd,nd,r,xin,xi(1,j))
14        continue
          ndo=nd
        endif
        if(nprn.ge.0) write(*,200) ndim,calls,it,itmx,nprn,ALPH,mds,nd,
     *(j,region(j),j,region(j+ndim),j=1,ndim)
      endif
      do 28 it=1,itmx
        ti=0.
        tsi=0.
        do 16 j=1,ndim
          kg(j)=1
          do 15 i=1,nd
            d(i,j)=0.
            di(i,j)=0.
15        continue
16      continue
10      continue
          fb=0.
          f2b=0.
          do 19 k=1,npg
            wgt=xjac
            do 17 j=1,ndim
              xn=(kg(j)-ran2(idum))*dxg+1.
              ia(j)=max(min(int(xn),NDMX),1)
              if(ia(j).gt.1)then
                xo=xi(ia(j),j)-xi(ia(j)-1,j)
                rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
              else
                xo=xi(ia(j),j)
                rc=(xn-ia(j))*xo
              endif
              x(j)=region(j)+rc*dx(j)
              wgt=wgt*xo*xnd
17          continue
            f=wgt*fxn(x,wgt)
            f2=f*f
            fb=fb+f
            f2b=f2b+f2
            do 18 j=1,ndim
              di(ia(j),j)=di(ia(j),j)+f
              if(mds.ge.0) d(ia(j),j)=d(ia(j),j)+f2
18          continue
19        continue
          f2b=sqrt(f2b*npg)
          f2b=(f2b-fb)*(f2b+fb)
          if (f2b.le.0.) f2b=TINY
          ti=ti+fb
          tsi=tsi+f2b
          if(mds.lt.0)then
            do 21 j=1,ndim
              d(ia(j),j)=d(ia(j),j)+f2b
21          continue
          endif
        do 22 k=ndim,1,-1
          kg(k)=mod(kg(k),ng)+1
          if(kg(k).ne.1) goto 10
22      continue
        tsi=tsi*dv2g
        wgt=1./tsi
        si=si+dble(wgt)*dble(ti)
        schi=schi+dble(wgt)*dble(ti)**2
        swgt=swgt+dble(wgt)
        tgral=si/swgt
        chi2a=max((schi-si*tgral)/(it-.99d0),0.d0)
        sd=sqrt(1./swgt)
        tsi=sqrt(tsi)
        if(nprn.ge.0)then
          write(*,201) it,ti,tsi,tgral,sd,chi2a
          if(nprn.ne.0)then
            do 23 j=1,ndim
              write(*,202) j,(xi(i,j),di(i,j),i=1+nprn/2,nd,nprn)
23          continue
          endif
        endif
        do 25 j=1,ndim
          xo=d(1,j)
          xn=d(2,j)
          d(1,j)=(xo+xn)/2.
          dt(j)=d(1,j)
          do 24 i=2,nd-1
            rc=xo+xn
            xo=xn
            xn=d(i+1,j)
            d(i,j)=(rc+xn)/3.
            dt(j)=dt(j)+d(i,j)
24        continue
          d(nd,j)=(xo+xn)/2.
          dt(j)=dt(j)+d(nd,j)
25      continue
        do 27 j=1,ndim
          rc=0.
          do 26 i=1,nd
            if(d(i,j).lt.TINY) d(i,j)=TINY
            r(i)=((1.-d(i,j)/dt(j))/(log(dt(j))-log(d(i,j))))**ALPH
            rc=rc+r(i)
26        continue
          call rebin(rc/xnd,nd,r,xin,xi(1,j))
27      continue
28    continue
      return
200   FORMAT(/' input parameters for vegas:  ndim=',i3,'  ncall=',
     *f8.0/28x,'  it=',i5,'  itmx=',i5/28x,'  nprn=',i3,'  alph=',
     *f5.2/28x,'  mds=',i3,'   nd=',i4/(30x,'xl(',i2,')= ',g11.4,' xu(',
     *i2,')= ',g11.4))
201   FORMAT(/' iteration no.',I3,': ','integral =',g14.7,'+/- ',g9.2/
     *' all iterations:   integral =',g14.7,'+/- ',g9.2,
     *' chi**2/it''n =',g9.2)
202   FORMAT(/' data for axis ',I2/'    X       delta i       ',
     *'   x       delta i       ','    x       delta i       ',/(1x,
     *f7.5,1x,g11.4,5x,f7.5,1x,g11.4,5x,f7.5,1x,g11.4))
      END

c-------------------------------------------------

      SUBROUTINE rebin(rc,nd,r,xin,xi)
      INTEGER nd
      real*8 rc,r(*),xi(*),xin(*)
      INTEGER i,k
      real*8 dr,xn,xo
      k=0
      xo=0.
      dr=0.
      do 11 i=1,nd-1
1       if(rc.gt.dr)then
          k=k+1
          dr=dr+r(k)
        goto 1
        endif
        if(k.gt.1) xo=xi(k-1)
        xn=xi(k)
        dr=dr-rc
        xin(i)=xn-(xn-xo)*dr/r(k)
11    continue
      do 12 i=1,nd-1
        xi(i)=xin(i)
12    continue
      xi(nd)=1.
      return
      END
c---------------------------------------                                                                                                                  
      FUNCTION rtbis_1(func,x1,x2,xacc)
      INTEGER JMAX
      REAL*8 rtbis_1,x1,x2,xacc,func
      EXTERNAL func
      PARAMETER (JMAX=200)
      INTEGER j
      REAL*8 dx,f,fmid,xmid
      fmid=func(x2)
      f=func(x1)
      if(f*fmid.ge.0.) pause 'root must be bracketed in rtbis'
      if(f.lt.0.)then
         rtbis_1=x1
        dx=x2-x1
      else
        rtbis_1=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
         dx=dx*.5
         xmid=rtbis_1+dx
         fmid=func(xmid)
        if(fmid.le.0.)rtbis_1=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
 11         continue
      pause 'too many bisections in rtbis'
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.                                                                                                       
c----------------------
      SUBROUTINE zbrac_1(func,x1,x2,succes)
      INTEGER NTRY
      REAL*8 x1,x2,func,FACTOR
      EXTERNAL func
      PARAMETER (FACTOR=1.6,NTRY=1000)
      INTEGER j
      REAL*8 f1,f2
      LOGICAL succes

      if(x1.eq.x2)pause 'you have to guess an initial range in zbrac'
      f1=func(x1)
      f2=func(x2)
      succes=.true.
      do 11 j=1,NTRY
        if(f1*f2.lt.0.)return
        if(abs(f1).lt.abs(f2))then
          x1=x1+FACTOR*(x1-x2)
          f1=func(x1)
        else
          x2=x2+FACTOR*(x2-x1)
          f2=func(x2)
        endif
 11         continue
      succes=.false.
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software V,3.                                                                                                                     



