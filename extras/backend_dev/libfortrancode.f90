!
! Dummy-code used as an example of an external Fortran library
! that can be accessed in GAMBIT via a C++ backend
!
! Author: Anders Kvellestad
! Date:   2013-03-26
!


module fortrancode
  
  implicit none

  ! Some variables
  integer :: testInt = 10 

  ! A (useless) common block, just to test the C++ backend
  integer, dimension(3) :: dummy_array = (/1,2,3/)
  double precision :: dummy_double = 1.2345
  common /commonBlock/ dummy_double, dummy_array

contains

  subroutine printMe(arg)
    implicit none
    double precision, dimension(3) :: arg
    write(*,*) "fortrancode: This is function 'printMe'."
    write(*,*) "fortrancode: Got array:", arg(1), arg(2), arg(3)
  end subroutine printMe


  integer function total(x, y)
    integer x, y
    write (*,*) "fortrancode: This is function 'total'."
    write (*,*) "fortrancode: Received arguments:",x,y
    total = x + y
    write (*,*) "fortrancode: Will return result:", total
    return
  end function total

end module fortrancode
