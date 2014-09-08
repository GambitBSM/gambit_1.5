#include <iostream>
#include "dlfcn.h"
#include <string>
#include <stdlib.h>

#include "GAMBIT_wrapper_X.hpp"
#include "GAMBIT_wrapper_Y.hpp"
#include "GAMBIT_wrapper_typedefs.hpp"


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



  // ------------
  // Test library
  // ------------


  cout << endl;
  cout << "=======================" << endl;
  cout << endl;


  X x1 = Factory_X_1();
  X x2 = Factory_X_2(10);

  Y y1 = Factory_Y_1();
  Y y2 = Factory_Y_2(x2);


  cout << "  x1.i = " << x1.i << endl;
  cout << "  x2.i = " << x2.i << endl;
  cout << "  y1.x.i = " << y1.x.i << endl;
  cout << "  y2.x.i = " << y2.x.i << endl;

  cout << "  ... " << endl;

  X x3 = y2.get_x();
  x3.i += 5;
  y2.set_x(x3);

  cout << "  x3.i = " << x3.i << endl;
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
