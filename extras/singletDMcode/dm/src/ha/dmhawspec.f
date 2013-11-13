**********************************************************************
*** subroutine dshawspec dumps the spectrum for which dshayield_int failed
*** to the file haspec.dat
**********************************************************************

      subroutine dmhawspec(f,a,b,n)
      implicit none

      real*8 f,a,b,y,x
      integer i,n
c      parameter(n=5000)
      external f

      open(unit=15,file='haspec.dat',status='unknown',
     &  form='formatted')
      do i=0,n
        x=a+(b-a)*dble(i)/dble(n)
        y=f(x)
        write(15,*) x,y
      enddo

      close(15)

      return
      end
