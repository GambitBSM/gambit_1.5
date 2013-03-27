/*
 * Testing backends
 * 
 * \author Anders Kvellestad
 * \date 2013-03-26
 * 
 */

#include <iostream>
#include "backend_rollcall.hpp"

int main()
{
  std::cout << std::endl;

  double myResult;
  
  
  /* Initialize and perform calculation with the first backend 
     (in GAMBIT::Backend::LibFirst), and then report on the result. */

  std::cout << "FIRST BACKEND:" << std::endl;

  GAMBIT::Backend::LibFirst::initialize(1);
  GAMBIT::Backend::LibFirst::someFunction();

  std::cout << std::endl;
  myResult = GAMBIT::Backend::LibFirst::returnResult();
  std::cout << "Function 'returnResult' returned: " << myResult << std::endl; 


  /* Check the corresponding result using the second backend 
     (in GAMBIT::Backend::LibFirstCopy).
     This should be 0 as we have not calculated anything using this backend.
     But if the two backends infact point to the same 'instance' of the library, 
     this result will be identical to the one reported by the first backend. 
     
     To test this, in backend_general.hpp change the line
     pHandle = dlmopen(LM_ID_NEWLM, LIBPATH, RTLD_LAZY);
     to 
     pHandle = dlopen(LIBPATH, RTLD_LAZY);  */
  
  std::cout << std::endl;
  std::cout << "SECOND BACKEND:" << std::endl;

  std::cout << std::endl;
  myResult = GAMBIT::Backend::LibFirstCopy::returnResult();
  std::cout << "Function 'returnResult' returned: " << myResult << std::endl; 


  std::cout << std::endl;
  return 0;
}
