#include <iostream>
#include "dlfcn.h"
#include <string>
#include <stdlib.h>


#include "GAMBIT_wrapper_classes.hpp"
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

  // T factory
  Factory_T = reinterpret_cast<Abstract_T* (*)()>( dlsym(pHandle, "_Z9Factory_Tv") );

  // X factory
  Factory_X = reinterpret_cast<Abstract_X* (*)()>( dlsym(pHandle, "_Z9Factory_Xv") );

  // Container factories
  Factory_Container_int = reinterpret_cast<Abstract_Container<int>* (*)()>( dlsym(pHandle, "_Z21Factory_Container_intv") );
  Factory_Container_T   = reinterpret_cast<Abstract_Container<Abstract_T>* (*)()>( dlsym(pHandle, "_Z19Factory_Container_Tv") );
  Factory_Container_X   = reinterpret_cast<Abstract_Container<Abstract_X>* (*)()>( dlsym(pHandle, "_Z19Factory_Container_Xv") );

  
  // printT function
  printT_WRAPPER = reinterpret_cast<void (*)(Abstract_T&)>( dlsym(pHandle, "_Z14printT_WRAPPERR10Abstract_T") );

  // doubleT function
  doubleT_WRAPPER = reinterpret_cast<void (*)(Abstract_T&)>( dlsym(pHandle, "_Z15doubleT_WRAPPERR10Abstract_T") );
  
  // doubleX_byVal function
  doubleX_byVal_WRAPPER = reinterpret_cast<Abstract_X* (*)(Abstract_X&)>( dlsym(pHandle, "_Z21doubleX_byVal_WRAPPERR10Abstract_X") );


  //
  // Test library
  //

  cout << endl;
  cout << "=======================" << endl;
  cout << endl;
  
  // Testing various ways of assigning values to class variables

  T_gambit test_t;
  test_t.printMe();

  test_t.i = 99;
  test_t.d = 1.11;
  test_t.printMe();
  
  X_gambit test_x;
  test_x.t.printMe();

  test_x.setT(test_t);
  test_x.t.printMe();

  test_t.i = 1000;
  test_x.t = test_t;
  test_x.t.printMe();

  cout << endl;
  cout << "=======================" << endl;
  cout << endl;


  // Testing global library functions 

  printT(test_t);
  printT(test_x.t);

  doubleT(test_x.t);
  test_x.t.printMe();

  X_gambit test_x2 = doubleX_byVal(test_x);
  test_x2.t.printMe();

  cout << endl;
  cout << "=======================" << endl;
  cout << endl;


  // Testing template class 'Container'

  Container_gambit<int> ci;
  ci.printMsg();
  
  Container_gambit<X_gambit> cx;
  cx.printMsg();
  
  cx.var.t.printMe();
  cx.var.t = test_t;
  cx.var.t.printMe();
  
  cout << endl;
  cout << "=======================" << endl;
  cout << endl;


  // Testing member_variable system


  cx.var._print_member();
  cx._print_member();
  test_x.t._print_member();
  test_x._print_member();


  cout << endl;
  cout << "=======================" << endl;
  cout << endl;

  
  // Testing convertion of multiple pointers

  T_gambit t1;
  T_gambit t2;

  t1.i = 10;
  t1.d = 1.11;
  t2.i = 20;
  t2.d = 2.22;

  T_gambit *t1_ptr;
  T_gambit **t1_ptr2;

  t1_ptr  = &t1;
  t1_ptr2 = &t1_ptr;

  T_gambit *t2_ptr;
  T_gambit **t2_ptr2;

  t2_ptr  = &t2;
  t2_ptr2 = &t2_ptr;


  (*t1_ptr).printMe();
  (*t2_ptr).printMe();
  cout << "passing in pointers: " << t1_ptr << "  -  " << t2_ptr << endl;
  // test_x.t.printMe();
  test_x.setT2(t1_ptr2, t2_ptr2);
  // test_x.t.printMe();
  cout << "returned pointers  : " << t1_ptr << "  -  " << t2_ptr << endl;
  (*t1_ptr).printMe();
  (*t2_ptr).printMe();


  //
  // Done
  // 

  cout << endl;
  return 0;
}
