***********************************************************************
*** function dsbessek1 returns the value of the modified bessel       ***
*** function of the second kind of order 1 times exp(x).            ***
*** works for positive real x                                       ***
*** coefficients from abramovitz and stegun.                        ***
*** e-mail: edsjo@physto.se                                 ***
*** date: 98-04-29                                                  ***
***********************************************************************

      real*8 function dsbessek1(x)
      implicit none
      real*8 x,y
      real*8 p1,p2,p3,p4,p5,p6,p7,
     &  q1,q2,q3,q4,q5,q6,q7
      save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
      real*8 dsbessei1
      data p1,p2,p3,p4,p5,p6,p7/
     &  1.0d0,0.15443144d0,-0.67278579d0,-0.18156897d0,
     &  -0.1919402d-1,-0.110404d-2,-0.4686d-4/
      data q1,q2,q3,q4,q5,q6,q7/
     &  1.25331414d0,0.23498619d0,-0.3655620d-1,
     &  0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/
      
      if (x.le.2.0d0) then
        y=x*x/4.0d0
        dsbessek1=exp(x)*((log(x/2.0d0)*exp(abs(x))*dsbessei1(x))+
     &       (1.0d0/x)*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))))
      else
        y=2.0d0/x
        dsbessek1=(1.0d0/sqrt(x))*(q1+y*(q2+y*(q3+
     &    y*(q4+y*(q5+y*(q6+y*q7))))))
      endif

      return
      end
