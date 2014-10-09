#include <iostream>
#include "dlfcn.h"
#include <string>
#include <stdlib.h>


#include "loadedtypes_minimal_0_1.hpp"

using namespace Gambit::Backends::Minimal_0_1;

// #include "boss_loaded_classes.hpp"

// #include "GAMBIT_wrapper_X_def.hpp"
// #include "GAMBIT_wrapper_Y_def.hpp"

// #include "GAMBIT_wrapper_typedefs.hpp"




using std::string;
using std::cout;
using std::endl;




// --------------------------------------------

int main(int argc, char * argv[])
{
  

  cout << endl;

  
  //
  // Load libraries
  //

  string libName = "../libminimal.so";

  cout << "Trying to load library from: " << libName.c_str() << endl;
  
  void * pHandle;
  if( (pHandle = dlopen(libName.c_str(), RTLD_LAZY | RTLD_GLOBAL)) == NULL)
  {
    cout << "Fail!" << endl; 
  }
  
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
  // Get pointers to library symbols
  //

  void * temp_ptr;

  // class X
  X (*Factory_X_1)();
  temp_ptr = dlsym(pHandle, "_Z9Factory_Xv");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  Factory_X_1 = reinterpret_cast<X (*)()> (temp_ptr);

  X (*Factory_X_2)(int);
  temp_ptr = dlsym(pHandle, "_Z9Factory_Xi");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  Factory_X_2 = reinterpret_cast<X (*)(int)> (temp_ptr);


  // class Y
  Y (*Factory_Y_1)();
  temp_ptr = dlsym(pHandle, "_Z9Factory_Yv");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  Factory_Y_1 = reinterpret_cast<Y (*)()> (temp_ptr);

  Y (*Factory_Y_2)(X);
  temp_ptr = dlsym(pHandle, "_Z9Factory_YR8X_GAMBIT");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  Factory_Y_2 = reinterpret_cast<Y (*)(X)> (temp_ptr);

  // Abstract_Y* (*Factory_Y_2_abs)(Abstract_X*);
  // temp_ptr = dlsym(pHandle, "_Z13Factory_Y_absP10Abstract_X");
  // if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  // Factory_Y_2_abs = reinterpret_cast<Abstract_Y* (*)(Abstract_X*)> (temp_ptr);



  // ------------
  // Test library
  // ------------


  cout << endl;
  cout << "=======================" << endl;
  cout << endl;


  X x1 = Factory_X_1();
  x1.i = 10;
  cout << "  x1.i = " << x1.i << endl;

  X& x2 = x1.return_ref_this();
  x2.i = 20;

  cout << "  x1.i = " << x1.i << endl;
  cout << "  x2.i = " << x2.i << endl;

  // X x3(x1);
  X x3 = x2;

  x3.i = 30;
  cout << "  x3.i = " << x3.i << endl;
  cout << "  x1.i = " << x1.i << endl;


  X x4 = x3 + x1;
  cout << "  x4.i = " << x4.i << endl;


  cout << endl;
  cout << endl;
  cout << endl;


  // X x3 = Factory_X_2(100);


  // X x2 = Factory_X_2(10);

  Y y1 = Factory_Y_1();
  Y y2 = Factory_Y_2(x2);
  // Y y2( Factory_Y_2_abs(x2.BEptr) );


  cout << "  x1.i = " << x1.i << endl;
  cout << "  x2.i = " << x2.i << endl;
  cout << "  y1.x.i = " << y1.x.i << endl;
  cout << "  y2.x.i = " << y2.x.i << endl;

  cout << "  ... " << endl;

  // X x3 = y2.get_x();
  x3.i += 5;
  y2.set_x(x3);

  cout << "  x3.i = " << x3.i << endl;
  cout << "  y2.x.i = " << y2.x.i << endl;


  cout << "  ... " << endl;


  X* xptr = x3.return_ptr_this();
  xptr->i += 5;

  cout << "  x3.i = " << x3.i << endl;

  y2.set_x_ptr(xptr);
  
  cout << "  y2.x.i = " << y2.x.i << endl;


  cout << "  ... " << endl;


  x3.set_yptr(&y2);

  Y* yptr = x3.get_yptr();

  yptr->x.i = 999;

  cout << "  y2.x.i = " << y2.x.i << endl;





  //
  // Done
  // 

  cout << endl;
  cout << "=======================" << endl;
  cout << endl;

  // dlclose(pHandle);

  return 0;
}
