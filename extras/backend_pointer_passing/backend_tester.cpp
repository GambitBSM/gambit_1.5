#include "backend_lib_fortran.hpp"
#include <iostream>


double some_other_function(int &input)
{
  std::cout << "  This is some_other_function, invoked with argument " << input << std::endl;
  return input * 2.0;
}


int main()
{


  int arg = 10;

  /* Use the 'runMe' functor */
  GAMBIT::Backends::LibFortran::Functown::runMe.calculate( GAMBIT::Backends::LibFortran::externalFunction, 20 );

  std::cout << std::endl;
  
  /* Use the 'runMe' function pointer, passing in a pointer to externalFunction */
  GAMBIT::Backends::LibFortran::runMe( GAMBIT::Backends::LibFortran::externalFunction, GAMBIT::byVal(arg) );
  GAMBIT::Backends::LibFortran::runMe( GAMBIT::Backends::LibFortran::externalFunction, 10 );

  std::cout << std::endl;

  /* Use the 'runMe' function pointer, passing in a pointer to some_other_function */
  GAMBIT::Backends::LibFortran::runMe( &some_other_function, GAMBIT::byVal(arg) );
  GAMBIT::Backends::LibFortran::runMe( &some_other_function, 10 );
  
  return 0;
}

