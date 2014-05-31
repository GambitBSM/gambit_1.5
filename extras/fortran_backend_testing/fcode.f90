!
! Dummy-code used as an example of an external Fortran code
! that can be accessed in GAMBIT via a C++ backend
!
! Author: Anders Kvellestad
! Date:   2013-06-25
!


module fcode
  implicit none

contains

  !
  ! Dummy subroutine 
  !
  subroutine addOneSubr(x)
    implicit none
    integer x 
    write(*,*) "-- Hello, this is the 'addOneSubr' subroutine."
    x = x + 1
  end subroutine addOneSubr

  !
  ! Dummy function
  !
  integer function addOneFunc(x)
    implicit none
    integer x
    write(*,*) "-- Hello, this is the 'addOneFunc' function."
    addOneFunc = x + 1
    return
  end function addOneFunc

end module fcode
