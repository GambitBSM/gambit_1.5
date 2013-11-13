***********************************************************************
*** function bessk2 returns the value of the modified bessel        ***
*** function of the second kind of order 2 times exp(x)             ***
*** works for positive real x                                       ***
*** recurrence relation                                             ***
*** e-mail: gondolo@mppmu.mpg.de                                    ***
*** date: 00-07-07                                                  ***
***********************************************************************

      real*8 function dsbessek2(x)
      implicit none
      real*8 x
      real*8 dsbessek1,dsbessek0
      dsbessek2 = 2.d0/x*dsbessek1(x)+dsbessek0(x)
      return
      end
