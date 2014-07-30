#include <iostream>
#include "dlfcn.h"
#include <string>
#include <stdlib.h>

// #include "GAMBIT_wrapper_Event.h"
// #include "GAMBIT_wrapper_Particle.h"
#include "GAMBIT_wrapper_Hist.h"
#include "GAMBIT_wrapper_Pythia.h"
#include "GAMBIT_wrapper_typedefs.h"


using std::string;
using std::vector;
using std::cout;
using std::endl;




// --------------------------------------------

int main(int argc, char * argv[])
{
  

  cout << endl;

  
  //
  // Load libraries
  //

  // string libName2 = "../../pythia8186/lib/liblhapdfdummy.so";

  // cout << "Trying to load library from: " << libName2.c_str() << endl;
  
  // void * pHandle2;
  // if( (pHandle2 = dlopen(libName2.c_str(), RTLD_LAZY | RTLD_GLOBAL)) == NULL)
  // {
  //   cout << "Fail!" << endl; 
  // }
  // if(pHandle2)
  // {
  //   cout << "Succeeded in loading library!" << endl;
  // }
  // else
  // {
  //   cout << dlerror() << endl;
  //   cout << "Failed loading library!" << endl;
  //   exit(1);
  // }
  // cout << endl;


  char* pythiaPath;
  pythiaPath = getenv("PYTHIA8DIR");
  string libName;

  if(pythiaPath != NULL)
  {
    libName = pythiaPath;
    libName += "/lib";
  }
  else
    libName = "..";
  libName += "/libpythia8.so";

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





  // exit(1);  

  //
  // Get pointers to library symbols and initialize the 
  // function pointers declared in GAMBIT_wrapper_*.hpp
  //

  void * temp_ptr;

  // class Pythia
  temp_ptr = dlsym(pHandle, "_Z14Factory_PythiaSsb");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  Factory_Pythia_0 = reinterpret_cast<Pythia8::Abstract_Pythia* (*)(std::string, bool)> (temp_ptr);



  // class Hist
  temp_ptr = dlsym(pHandle, "_Z12Factory_Histv");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  Factory_Hist_0 = reinterpret_cast<Pythia8::Abstract_Hist* (*)()> (temp_ptr);

  temp_ptr = dlsym(pHandle, "_Z12Factory_HistSsidd");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  Factory_Hist_1 = reinterpret_cast<Pythia8::Abstract_Hist* (*)(string, int, double, double)> (temp_ptr);



  // // class Particle
  // temp_ptr = dlsym(pHandle, "_Z16Factory_Particlev");
  // if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  // Factory_Particle_0 = reinterpret_cast<Pythia8::Abstract_Particle* (*)()> (temp_ptr);



  // // class Event
  // temp_ptr = dlsym(pHandle, "_Z13Factory_Eventi");
  // if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  // Factory_Event_0 = reinterpret_cast<Pythia8::Abstract_Event* (*)(int)> (temp_ptr);


  // ------------
  // Test library
  // ------------


  cout << endl;
  cout << "=======================" << endl;
  cout << endl;

  // Generator. Process selection. Tevatron initialization. Histogram.
  Pythia pythia("../../pythia8186/xmldoc", true);

  pythia.readString("Beams:idB = -2212", true);
  pythia.readString("Beams:eCM = 1960.", true);
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on", true);
  pythia.readString("PhaseSpace:mHatMin = 80.", true);
  pythia.readString("PhaseSpace:mHatMax = 120.", true);
  pythia.init();
  Hist pTZ("dN/dpTZ", 100, 0., 100.);
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 1000; ++iEvent) {
    if (!pythia.next()) continue;
    // Loop over particles in event. Find last Z0 copy. Fill its pT.
    int iZ = 0;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].id() == 23) iZ = i;
    pTZ.fill( pythia.event[iZ].pT(), 1.0 );
  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  // cout << mult;  <-- This won't work as we haven't defined the << operator for Hist


  //
  // Done
  // 

  cout << endl;
  cout << "=======================" << endl;
  cout << endl;

  // dlclose(pHandle);

  return 0;
}
