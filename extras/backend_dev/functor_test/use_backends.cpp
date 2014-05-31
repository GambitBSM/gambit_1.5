/*
 * Testing backends
 * 
 * \author Anders Kvellestad
 * \date 2013-03-26
 * 
 * Modified: 2013-04-05
 */

#include <iostream>
#include "backend_rollcall.hpp"
#include "backendfunctors.hpp"

int main()
{
  std::cout << std::endl;

  /* Simulate core/module/backend interaction: */

  /* The module creates local pointers to backend functors: 
   * 
   * Problem:
   * How can it be determined what template arguments to give?
   * (The backend_functor2 class is templated both on the return type
   * and on the argument list.) The return type is known from 
   * START_BACKEND_REQ, but the argument list is not. Can this list 
   * be communicated from the backend functor object in any way, so that
   * START_BACKEND_REQ can create a fully specialized functor pointer? 
   */
  GAMBIT::backend_functor2<void, int> *LibFirst_initialize;
  GAMBIT::backend_functor2<void> *LibFirst_someFunction;
  GAMBIT::backend_functor2<double> *LibFirst_returnResult;
  GAMBIT::backend_functor2<double,int> *LibFirst_doAll;

  GAMBIT::backend_functor2<double> *LibFirst_getSomeDouble;
  GAMBIT::backend_functor2<void,double> *LibFirst_setSomeDouble;

  /* The dependency resolver connects these pointers to the 
   * corresponding functors within the backend: */
  LibFirst_initialize   = &(GAMBIT::Backend::LibFirst::Functown::initialize);
  LibFirst_someFunction = &(GAMBIT::Backend::LibFirst::Functown::someFunction);
  LibFirst_returnResult = &(GAMBIT::Backend::LibFirst::Functown::returnResult);
  LibFirst_doAll        = &(GAMBIT::Backend::LibFirst::Functown::doAll);

  LibFirst_getSomeDouble = &(GAMBIT::Backend::LibFirst::Functown::getSomeDouble);
  LibFirst_setSomeDouble = &(GAMBIT::Backend::LibFirst::Functown::setSomeDouble);


  /* Test of backends to C++ libraray 'libfirst.so' */
  
  double myResult;

  std::cout << "BACKEND: LibFirst" << std::endl;
  std::cout << "-----------------" << std::endl;

  /* Test identification methods of a functor object */
  std::cout << std::endl;
  std::cout << "==== Testing identification methods of the 'initialize' functor: ====" << std::endl;
  std::cout << std::endl;
  std::cout << "Name           : " << LibFirst_initialize->name()       << std::endl 
            << "Capability     : " << LibFirst_initialize->capability() << std::endl 
            << "Return type    : " << LibFirst_initialize->type()       << std::endl 
            << "Origin         : " << LibFirst_initialize->origin()     << std::endl
            << "Origin version : " << LibFirst_initialize->version()    << std::endl;


  /* Use backend functors */
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "==== Testing the basic backend functors: ====" << std::endl;  
  (*LibFirst_initialize)(2);
  (*LibFirst_someFunction)();

  std::cout << std::endl;
  myResult = (*LibFirst_returnResult)();
  std::cout << "Function 'returnResult' returned: " << myResult << std::endl; 
  

  /* Test get/set functors */
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "==== Testing the 'get/set' functors for library variables: ====" << std::endl;  
  double testdouble;
  
  (*LibFirst_setSomeDouble)(1.23456);
  testdouble = (*LibFirst_getSomeDouble)();
  std::cout << std::endl;
  std::cout << "testdouble: " << testdouble << std::endl; 
  std::cout << std::endl;

  myResult = (*LibFirst_returnResult)();
  std::cout << "Function 'returnResult' returned: " << myResult << std::endl; 


  /* Test functor for convenience function */
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "==== Testing functors for convenience functions: ====" << std::endl;  
  myResult = (*LibFirst_doAll)(3);
  std::cout << "Function 'returnResult' returned: " << myResult << std::endl; 


  std::cout << std::endl;
  return 0;
}
