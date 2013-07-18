/*
 * Testing backends
 * 
 * \author Anders Kvellestad
 * \date 2013-03-26
 * 
 * Modified: 2013-07-18
 */

#include <iostream>
#include "backend_rollcall.hpp"

int main()
{
  std::cout << std::endl;

  
  /* Test of backends to C++ libraray 'libfirst.so' */
  
  double myResult;

  /* Initialize and perform calculation with the first backend 
   * (in GAMBIT::Backend::LibFirst), and then report on the result. */

  
  std::cout << "BACKEND: LibFirst" << std::endl;
  std::cout << "-----------------" << std::endl;

  
  GAMBIT::Backend::LibFirst::initialize(2);
  GAMBIT::Backend::LibFirst::someFunction();

  std::cout << std::endl;
  myResult = GAMBIT::Backend::LibFirst::returnResult();
  std::cout << "Function 'returnResult' returned: " << myResult << std::endl; 
  

  /* Check the corresponding result using the second backend 
   * (in GAMBIT::Backend::LibFirstCopy).
   * This should be 0 as we have not calculated anything using this backend.
   * But if the two backends infact point to the same 'instance' of the library, 
   * this result will be identical to the one reported by the first backend. */
  
 
  std::cout << std::endl;
  std::cout << "BACKEND: LibFirstCopy" << std::endl;
  std::cout << "---------------------" << std::endl;
  
  // GAMBIT::Backend::LibFirstCopy::someFunction();

  std::cout << std::endl;
  myResult = GAMBIT::Backend::LibFirstCopy::returnResult();
  std::cout << "Function 'returnResult' returned: " << myResult << std::endl; 
  
 
  std::cout << std::endl;
  return 0;
}
