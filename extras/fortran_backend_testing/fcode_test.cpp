
#include <iostream>
#include "backend_fcode.hpp"

int main()
{
  std::cout << std::endl;
 
  /* Test of backend for Fortran libraray 'libfcode.so' */

  int x = 1;

  /* First test the Fortran function 'addOneFunc'. 
   * Can either use the functor (in Functown) or the function pointer directly. */
  std::cout << "x = " << x << std::endl;
  x = GAMBIT::Backends::LibFcode::Functown::addOneFunc(x);
  //x = GAMBIT::Backends::LibFcode::addOneFunc(x);
  std::cout << "x = " << x << std::endl;

  /* Then test the Fortran subroutine 'addOneSubr': */
  GAMBIT::Backends::LibFcode::Functown::addOneSubr(x);
  //GAMBIT::Backends::LibFcode::addOneSubr(x);
  std::cout << "x = " << x << std::endl;
  std::cout << std::endl;

  return 0;
}
