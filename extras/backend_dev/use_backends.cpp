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

  
  /* Test of backends to C++ libraray 'libfirst.so' */
  
  double myResult;

  /* Initialize and perform calculation with the first backend 
     (in GAMBIT::Backend::LibFirst), and then report on the result. */

  std::cout << "BACKEND: LibFirst" << std::endl;
  std::cout << "-----------------" << std::endl;

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
  std::cout << "BACKEND: LibFirstCopy" << std::endl;
  std::cout << "---------------------" << std::endl;
  
  std::cout << std::endl;
  myResult = GAMBIT::Backend::LibFirstCopy::returnResult();
  std::cout << "Function 'returnResult' returned: " << myResult << std::endl; 




  /* Test of backend to Fortran libraray 'libfortrancode.so' */

  std::cout << std::endl << std::endl;
  std::cout << "BACKEND: LibFortranCode" << std::endl;
  std::cout << "-----------------------" << std::endl;
  std::cout << std::endl;
  

  std::cout << "Testing get/set of int variable 'testInt':" << std::endl;
  std::cout << std::endl;

  int testInt;

  testInt = GAMBIT::Backend::LibFortranCode::getTestInt();
  std::cout << "testInt: " << testInt << std::endl;

  GAMBIT::Backend::LibFortranCode::setTestInt(20);
  testInt = GAMBIT::Backend::LibFortranCode::getTestInt();
  std::cout << "testInt: " << testInt << std::endl;
  

  std::cout << std::endl;
  std::cout << "Testing subroutine 'printMe':" << std::endl;
  std::cout << std::endl;
  
  double array[3] = {1.0,2.0,3.0};
  GAMBIT::Backend::LibFortranCode::printMe(array);


  std::cout << std::endl;
  std::cout << "Testing function 'total':" << std::endl;
  std::cout << std::endl;

  double myTotal = GAMBIT::Backend::LibFortranCode::total(2, 3);
  std::cout << "Function 'total' returned: " << myTotal << std::endl;

  
  std::cout << std::endl;
  std::cout << "Testing get/set of common block 'commonBlock':" << std::endl;
  std::cout << std::endl;

  commonBlock_type myCommonBlock;
  myCommonBlock = GAMBIT::Backend::LibFortranCode::getCommonBlock();
  std::cout << "commonBlock.dummyDouble  : " << myCommonBlock.dummyDouble << std::endl;
  std::cout << "commonBlock.dummyArray[0]: " << myCommonBlock.dummyArray[0] << std::endl;
  std::cout << "commonBlock.dummyArray[1]: " << myCommonBlock.dummyArray[1] << std::endl;
  std::cout << "commonBlock.dummyArray[2]: " << myCommonBlock.dummyArray[2] << std::endl;

  myCommonBlock.dummyArray[1] = 200;
  
  GAMBIT::Backend::LibFortranCode::setCommonBlock(myCommonBlock);

  myCommonBlock = GAMBIT::Backend::LibFortranCode::getCommonBlock();
  std::cout << std::endl;
  std::cout << "commonBlock.dummyDouble  : " << myCommonBlock.dummyDouble << std::endl;
  std::cout << "commonBlock.dummyArray[0]: " << myCommonBlock.dummyArray[0] << std::endl;
  std::cout << "commonBlock.dummyArray[1]: " << myCommonBlock.dummyArray[1] << std::endl;
  std::cout << "commonBlock.dummyArray[2]: " << myCommonBlock.dummyArray[2] << std::endl;


  std::cout << std::endl;
  return 0;
}
