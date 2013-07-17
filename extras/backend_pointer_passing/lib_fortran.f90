!
! Dummy-code used as an example of an external Fortran library
! that can be accessed in GAMBIT via a C++ backend
!
! Author: Anders Kvellestad
!

subroutine testMe()
  write(*,*) "This is testMe. Not much happens here." 
end subroutine testMe


subroutine runMe(f, i)
  external f
  real f,f_res
  integer i

  write(*,*) "This is runMe. Calling externalRoutine with arguments:",i
  f_res = f(i)
  write(*,*) "This is runMe. Got result:",f_res
  
  !call f(i)
end subroutine runMe


subroutine externalRoutine(i)
  write(*,*) "This is externalRoutine called with arguments:",i
end subroutine externalRoutine


real function externalFunction(i)
  write(*,*) "This is externalFunction called with arguments:",i
  externalFunction = i*3.14
  return
end function



