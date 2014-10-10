#include <iostream>
#include "dlfcn.h"
#include <string>
#include <stdlib.h>

#include "GAMBIT_wrapper_Particle.h"
#include "GAMBIT_wrapper_Event.h"
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



  // class Vec4
  temp_ptr = dlsym(pHandle, "_Z12Factory_Vec4dddd");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  Factory_Vec4_0 = reinterpret_cast<Pythia8::Abstract_Vec4* (*)(double, double, double, double)> (temp_ptr);





  // class Particle
  temp_ptr = dlsym(pHandle, "_Z16Factory_Particlev");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  Factory_Particle_0 = reinterpret_cast<Pythia8::Abstract_Particle* (*)()> (temp_ptr);



  // class Event
  temp_ptr = dlsym(pHandle, "_Z13Factory_Eventi");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  Factory_Event_0 = reinterpret_cast<Pythia8::Abstract_Event* (*)(int)> (temp_ptr);

  temp_ptr = dlsym(pHandle, "_Z13Factory_Eventv");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  Factory_Event_0_overload_1 = reinterpret_cast<Pythia8::Abstract_Event* (*)()> (temp_ptr);


  // ------------
  // Test library
  // ------------


  cout << endl;
  cout << "=======================" << endl;
  cout << endl;


  // {
  //   cout << endl;
  //   cout << "START SCOPE 1" << endl;
  //   cout << endl;


  //   cout << "Constructing Vec4 'p1'" << endl;
  //   Vec4 p1(0.0, 0.0, 0.0, 0.0);  

  //   {
  //     cout << endl;
  //     cout << "START SCOPE 2" << endl;
  //     cout << endl;

  //     cout << endl;
  //     cout << "Constructing Vec4 'p2'" << endl;
  //     Vec4 p2(1.0, 2.0, 3.0, 4.0);
      
  //     cout << endl;
  //     cout << "Exectuing 'Vec4 p3 = p1 += p2'" << endl;
  //     p1 += p2;
  //     Vec4& p1_ref = p1;

  //     p1_ref = Vec4(2.0, 2.0, 2.0, 2.0);

  //     p1_ref = 9.9;

  //     cout << endl;
  //     cout << "END SCOPE 2" << endl;
  //     cout << endl;
  //   }
  
  //   cout << "  p1 : " << p1.px() << " , " << p1.py() << " , " << p1.pz() << " , " << p1.e() << endl;


  //   cout << endl;
  //   cout << "END SCOPE 1" << endl;
  //   cout << endl;
  // }

  // PROBLEM: p1 is never deleted, while p2 is deleted twice.




  Pythia pythia("../../pythia8186/xmldoc", false);

  // Shorthand for the event record in pythia.
  Event& event = pythia.event;

  // Read in commands from external file.
  pythia.readFile("./main03.cmnd");

  // Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Initialize.
  pythia.init();

  // Book histograms.
  Hist pThard("process pT scale", 100, 0., 500.);
  Hist nFinal("final particle multiplicity", 100, -0.5, 1599.5);
  Hist nCharged("charged particle multiplicity", 100, -0.5, 799.5);
  Hist dndy("dn/dy for charged particles", 100, -10., 10.);
  Hist dndeta("dn/d(eta) for charged particles", 100, -10., 10.);
  Hist dndpT("dn/dpT for charged particles", 100, 0., 10.);
  Hist epCons("deviation from energy-momentum conservation", 100, 0., 1e-6);

  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }

    // Fill hard scale of event.
    pThard.fill( pythia.info.pTHat() );

    // Loop over final particles in the event.
    int  nFin = 0;
    int  nChg = 0;
    Vec4 pSum(0.0, 0.0, 0.0, 0.0);


    for (int i = 0; i < event.size(); ++i) if (event[i].isFinal()) {

      // Analyze all particles.
      nFin++;

      pSum += event[i].p();  // SEGFAULT!
    
      cout << endl;
      cout << "TEST:  Start" << endl;
      cout << "TEST:  statement: Particle p1" << endl;
      Particle p1;
      cout << "TEST:  statement: p1.id(123)" << endl;
      p1.id(123);
      cout << "TEST:  p1.id() = " << p1.id() << endl;
      cout << "TEST:  statement: event[i] = p1" << endl;
      event[i] = p1;
      cout << "TEST:  event[i].id() = " << event[i].id() << endl;
      cout << "TEST:  End" << endl;



      // Analyze charged particles and fill histograms.
      if (event[i].isCharged()) {
        ++nChg;
        dndy.fill( event[i].y() );
        dndeta.fill( event[i].eta() );
        dndpT.fill( event[i].pT() );
      }

    // End of particle loop. Fill global properties.
    }

    nFinal.fill( nFin );
    nCharged.fill( nChg );
    pSum /= event[0].e();
    double epDev = abs(pSum.e() - 1.) + abs(pSum.px()) + abs(pSum.py())
      + abs(pSum.pz());
    epCons.fill(epDev);

  // End of event loop.
  }
   

  // Final statistics. Normalize and output histograms.
  pythia.stat();
  dndy   *=  5. / nEvent;
  dndeta *=  5. / nEvent;
  dndpT  *= 10. / nEvent;
  // cout << pThard << nFinal << nCharged << dndy << dndeta << dndpT << epCons;  //  <-- Haven't loaded the global << operator yet...


  //
  // Done
  // 

  cout << endl;
  cout << "=======================" << endl;
  cout << endl;

  // dlclose(pHandle);

  return 0;
}
