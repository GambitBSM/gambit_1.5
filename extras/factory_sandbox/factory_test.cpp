/*
 * Testing class loading
 * 
 * \author Pat Scott
 * \date 2013-06-23
 * 
 * Adapted from Anders Kvellestad's use_backends for testing backend macros
 */

#include <iostream>
#include "backend_libfirst.hpp"

int main()
{
  std::cout << std::endl;

  
  /* Test of backends to C++ libraray 'libfirst.so' */
  
  double myResult;
  typedef GAMBIT::Classes::LibFirst::v1_0::cattleClass Cow;

  /* Initialize and perform calculation with LibFirst backend 
   * (in GAMBIT::Backend::LibFirst), and then report on the result. */

  
  std::cout << "BACKEND: LibFirst" << std::endl;
  std::cout << "-----------------" << std::endl;
  
  GAMBIT::Backend::LibFirst::initialize(2);
  std::cout << "Hmmm OK, moving right along..." << std::endl;
  GAMBIT::Backend::LibFirst::someFunction();

  std::cout << std::endl;
  myResult = GAMBIT::Backend::LibFirst::returnResult();
  std::cout << "Function 'returnResult' returned: " << myResult << std::endl << std::endl; 
  
  std::cout << "Creating a proprietary cattleClass object using libFirst..." << std::endl;
  Cow* Crusher = GAMBIT::Backend::LibFirst::cattleClassFactory(1000.0);

  std::cout << "Checking how much I can get for it:" << std::endl;
  std::cout << Crusher->valueMe(4.45) << " dollars.  Not bad. ";
  std::cout << "OK, time to cash it in!" << std::endl;
  delete Crusher;

  std::cout << std::endl;
  std::cout << "No cattle were harmed in the making of this example.  Some coffee beans might have been though." << std::endl;
  std::cout << std::endl;
  return 0;
}
