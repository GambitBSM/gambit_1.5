!Example driver program for CTA Galactic Centre likelihood code
!Pat Scott patscott@physics.mcgill.ca
!May 24 2013

module example_helper

use types
use cta_gc

implicit none
double precision :: ms, deltalnlike
Type(BFSet) :: BFs   

contains

  double precision function cta_likelihood_difference(logsigv)

    double precision :: logsigv, sigv
    sigv = 10.**logsigv 
    cta_likelihood_difference = ctalike(BFs, sigv, ms) + deltalnlike

  end function cta_likelihood_difference


end module example_helper


program cta_example

use example_helper

implicit none

  double precision, parameter :: tol=1.d-6, msmin = 32.d0, msmax = 5000.d0
  integer, parameter :: nms = 200

  double precision :: zbrent, CL, erfinv, lhs, t
  integer :: i

  feedback = 0

  call init_ctalike

  call clear(BFs)
  BFs%b=1.d0

  CL = 90.0d0
  deltalnlike = erfinv(CL*0.01d0, (100.d0-CL)*0.01d0)**2
  do i = 1, nms
    ms = msmin * 10.d0**(dble(i-1)/dble(nms-1) * (log10(msmax) - log10(msmin)))
    lhs = zbrent(cta_likelihood_difference,-35.d0,-10.d0,tol)
    write(*,*) ms, lhs 
  enddo

end program

