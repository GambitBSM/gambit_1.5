#include <iostream>
#include "dlfcn.h"
#include <string>
#include <stdlib.h>


#include "GAMBIT_wrapper_U.hpp"
#include "GAMBIT_wrapper_T.hpp"
#include "GAMBIT_wrapper_X.hpp"
#include "GAMBIT_wrapper_functions.hpp"

using std::string;
using std::cout;
using std::endl;




// --------------------------------------------

int main(int argc, char * argv[])
{
  

  cout << endl;

  //
  // Library path
  //
  
  string libName = "../after_BOSS/lib.so";
  
  
  //
  // Load library
  //

  cout << "Trying to load library from:" << endl;
  cout << libName.c_str() << endl;
  void * pHandle = dlopen(libName.c_str(), RTLD_LAZY);
  if(pHandle)
  {
    cout << "Succeeded in loading library!" << endl;
  }
  else
  {
    cout << dlerror() << endl;
    cout << "Failed loading library!" << endl;
    exit(1);
  }
  cout << endl;
  

  //
  // Get pointers to library symbols and initialize the 
  // function pointers declared in GAMBIT_wrapper_*.hpp
  //

  void * temp_ptr;

  // U factory 1
  temp_ptr = dlsym(pHandle, "_Z9Factory_Uv");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  Factory_U_1 = reinterpret_cast<Abstract_U* (*)()>( temp_ptr );

  // T factory 1
  temp_ptr = dlsym(pHandle, "_Z9Factory_Tv");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  Factory_T_1 = reinterpret_cast<Abstract_T* (*)()>( temp_ptr );

  // T factory 2
  temp_ptr = dlsym(pHandle, "_Z9Factory_Tid");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  Factory_T_2 = reinterpret_cast<Abstract_T* (*)(int,double)>( temp_ptr );

  // X factory 1
  temp_ptr = dlsym(pHandle, "_Z9Factory_Xv");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  Factory_X_1 = reinterpret_cast<Abstract_X* (*)()>( temp_ptr );

  // X factory 2
  // temp_ptr = dlsym(pHandle, "_Z9Factory_XR10Abstract_T");
  temp_ptr = dlsym(pHandle, "_Z9Factory_X10Abstract_T");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  Factory_X_2 = reinterpret_cast<Abstract_X* (*)(Abstract_T&)>( temp_ptr );

  // // Container factories
  // temp_ptr = dlsym(pHandle, "_Z21Factory_Container_intv");
  // if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  // Factory_Container_int = reinterpret_cast<Abstract_Container<int>* (*)()>( temp_ptr );

  // temp_ptr = dlsym(pHandle, "_Z19Factory_Container_Tv");
  // if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  // Factory_Container_T   = reinterpret_cast<Abstract_Container<Abstract_T>* (*)()>( temp_ptr );

  // temp_ptr = dlsym(pHandle, "_Z19Factory_Container_Xv");
  // if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  // Factory_Container_X   = reinterpret_cast<Abstract_Container<Abstract_X>* (*)()>( temp_ptr );

  
  // printT function
  temp_ptr = dlsym(pHandle, "_Z13printT_GAMBITR10Abstract_T");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  printT_GAMBIT = reinterpret_cast<void (*)(Abstract_T&)>( temp_ptr );

  // doubleT function
  temp_ptr = dlsym(pHandle, "_Z14doubleT_GAMBITR10Abstract_T");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  doubleT_GAMBIT = reinterpret_cast<void (*)(Abstract_T&)>( temp_ptr );
  
  // doubleX_byVal function
  temp_ptr = dlsym(pHandle, "_Z20doubleX_byVal_GAMBITR10Abstract_X");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  doubleX_byVal_GAMBIT = reinterpret_cast<Abstract_X* (*)(Abstract_X&)>( temp_ptr );



  // ------------
  // Test library
  // ------------


  cout << endl;
  cout << "=======================" << endl;
  cout << endl;

  T_gambit t;
  X_gambit x;

  // // Testing various ways of assigning values to class variables

  // t.printMe();

  // t.i = 99;
  // t.d = 1.11;
  // t.printMe();
  
  // x.t.printMe();

  // x.setT(t);
  // x.t.printMe();

  // t.i = 1000;
  // x.t = t;
  // x.t.printMe();

  // cout << endl;
  // cout << "=======================" << endl;
  // cout << endl;


  // // Testing global library functions 
  // printT(t);
  // printT(x.t);

  // doubleT(x.t);
  // x.t.printMe();

  // X_gambit x2 = doubleX_byVal(x);
  // x2.t.printMe();

  // cout << endl;
  // cout << "=======================" << endl;
  // cout << endl;


  // // Testing template class 'Container'

  // Container_gambit<int> ci;
  // ci.printMsg();
  
  // // cout << "\nHI!\n" << endl;
  // Container_gambit<X_gambit> cx;
  // cout << "\nHI2!\n" << endl;
  // cx.printMsg();
  
  // cx.var.t.printMe();
  // cx.var.t = t;
  // cx.var.t.printMe();
  
  // cout << endl;
  // cout << "=======================" << endl;
  // cout << endl;


  // // Testing member_variable system


  // cx.var._print_member();
  // cx._print_member();
  // x.t._print_member();
  // x._print_member();


  // cout << endl;
  // cout << "=======================" << endl;
  // cout << endl;


  // Testing refTest function:
  int a = 555;

  cout << " a = " << a << endl;
  t.printMe();
  cout << endl;
  
  x.refTest(t, a);

  cout << " a = " << a << endl;
  t.printMe();
  
  cout << endl;
  cout << "=======================" << endl;
  cout << endl;


  // Testing new constructors
  T_gambit t2(99,9.9);
  t2.printMe();
  cout << endl;

  X_gambit x3(t2);

  x3.t.printMe();
  
  cout << endl;
  cout << "=======================" << endl;
  cout << endl;


  // Testing handling of input and return types
  int i = 10;
  int* ip = &i;
  int** ipp = &ip;
  double d = 8.8;
  cout << " i = " << **ipp << endl;
  x3.testFunc(&t, t2, ipp, d);
  cout << " i = " << **ipp << endl;

  cout << endl;
  cout << "=======================" << endl;
  cout << endl;


  // Testing inheritance
  cout << endl;
  (*t.BEptr).memberFunc();
  t.memberFunc();
  cout << endl;

  cout << endl;
  cout << "=======================" << endl;
  cout << endl;

  //
  // Done
  // 

  // dlclose(pHandle);

  return 0;
}
