#include "backend_lib_fortran.hpp"
#include <iostream>



int main()
{


  int arg = 10;

  /* Use the 'runMe' functor */
  GAMBIT::Backends::LibFortran::Functown::runMe.calculate( GAMBIT::Backends::LibFortran::externalFunction, 10 );

  std::cout << std::endl;
  
  /* Use the 'runMe' function pointer */
  GAMBIT::Backends::LibFortran::runMe( GAMBIT::Backends::LibFortran::externalFunction, GAMBIT::byVal(arg) );
  //GAMBIT::Backends::LibFortran::runMe( GAMBIT::Backends::LibFortran::externalFunction, 10 );
  
  return 0;
}

