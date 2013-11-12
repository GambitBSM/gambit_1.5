
c#######################################################################
c
c   one- and two-dimensional adaptive gaussian integration routines.
c
c **********************************************************************

      subroutine dgadap(a0,b0,f,eps0,sum)
c
c   purpose           - integrate a function f(x)
c   method            - adaptive gaussian
c   usage             - call gadap(a0,b0,f,eps,sum)
c   parameters  a0    - lower limit (input,real)
c               b0    - upper limit (input,real)
c               f     - function f(x) to be integrated. must be
c                       supplied by the user. (input,real function)
c               eps0  - desired relative accuracy. if sum is small eps
c                       will be absolute accuracy instead. (input,real)
c               sum   - calculated value for the integral (output,real)
c   precision         - single (see below)
c   req'd prog's      - f
c   author            - t. johansson, lund univ. computer center, 1973
c   reference(s)      - the australian computer journal,3 p.126 aug. -71
c
c  made real*8 by j. edsjo 97-01-17
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      include 'dsidtag.h'
      common/gadap1/ num,ifu
      external f
      dimension a(300),b(300),f1(300),f2(300),f3(300),s(300),n(300)
    1 format(16h gadap:i too big)
      dsum(f1f,f2f,f3f,aa,bb)=5./18.*(bb-aa)*(f1f+1.6*f2f+f3f)
      eps=eps0
      if(eps.lt.1.0d-8) eps=1.0d-8
      red=1.3
      l=1
      i=1
      sum=0.
      c=sqrt(15.)/5.
      a(1)=a0
      b(1)=b0
      f1(1)=f(0.5*(1+c)*a0+0.5*(1-c)*b0)
      f2(1)=f(0.5*(a0+b0))
      f3(1)=f(0.5*(1-c)*a0+0.5*(1+c)*b0)
      ifu=3
      s(1)=  dsum(f1(1),f2(1),f3(1),a0,b0)
  100 continue
      l=l+1
      n(l)=3
      eps=eps*red
      a(i+1)=a(i)+c*(b(i)-a(i))
      b(i+1)=b(i)
      a(i+2)=a(i)+b(i)-a(i+1)
      b(i+2)=a(i+1)
      a(i+3)=a(i)
      b(i+3)=a(i+2)
      w1=a(i)+(b(i)-a(i))/5.
      u2=2.*w1-(a(i)+a(i+2))/2.
      f1(i+1)=f(a(i)+b(i)-w1)
      f2(i+1)=f3(i)
      f3(i+1)=f(b(i)-a(i+2)+w1)
      f1(i+2)=f(u2)
      f2(i+2)=f2(i)
      f3(i+2)=f(b(i+2)+a(i+2)-u2)
      f1(i+3)=f(a(i)+a(i+2)-w1)
      f2(i+3)=f1(i)
      f3(i+3)=f(w1)
      ifu=ifu+6
      if(ifu.gt.5000) goto 125
      s(i+1)=  dsum(f1(i+1),f2(i+1),f3(i+1),a(i+1),b(i+1))
      s(i+2)=  dsum(f1(i+2),f2(i+2),f3(i+2),a(i+2),b(i+2))
      s(i+3)=  dsum(f1(i+3),f2(i+3),f3(i+3),a(i+3),b(i+3))
      ss=s(i+1)+s(i+2)+s(i+3)
      i=i+3
      if(i.gt.300)goto 120
      sold=s(i-3)
      if(abs(sold-ss).gt.eps*(1.+abs(ss))/2.) goto 100
      sum=sum+ss
      i=i-4
      n(l)=0
      l=l-1
  110 continue
      if(l.eq.1) goto 130
      n(l)=n(l)-1
      eps=eps/red
      if(n(l).ne.0) goto 100
      i=i-1
      l=l-1
      goto 110
  120 write(*,*) 'warning in dgadap: too high i (>300)'
      write(*,*) '  for model ',idtag
      goto 130
 125  write(*,*) 'warning in dgadap: too many function calls (>5000)'
      write(*,*) '  for model ',idtag
  130 return
      end






