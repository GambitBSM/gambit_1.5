#include <iostream>
#include "dlfcn.h"
#include <string>
#include <stdlib.h>

#include <cmath>

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



  // class Vec4
  temp_ptr = dlsym(pHandle, "_Z12Factory_Vec4dddd");
  if(!temp_ptr) { cout << dlerror() << endl; exit(1); }
  Factory_Vec4_0 = reinterpret_cast<Pythia8::Abstract_Vec4* (*)(double, double, double, double)> (temp_ptr);



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


  Pythia pythia("../../pythia8186/xmldoc", false);

  // Shorthand for the event record in pythia.
  Event& event = pythia.event;

  // Read in commands from external file.
  pythia.readFile("./main04.cmnd");

  // Extract settings to be used in the main program.
  int    nEvent    = pythia.mode("Main:numberOfEvents");
  int    nAbort    = pythia.mode("Main:timesAllowErrors");
 
  // Initialize.
  pythia.init();

  // Book histograms: multiplicities and mean transverse momenta.
  Hist yChg("rapidity of charged particles; all",      100, -10., 10.);
  Hist nChg("number of charged particles; all",        100, -0.5, 799.5);
  Hist nChgSD("number of charged particles; single diffraction",
                                                       100, -0.5, 799.5);
  Hist nChgDD("number of charged particles, double diffractive",
                                                       100, -0.5, 799.5);
  Hist nChgCD("number of charged particles, central diffractive",
                                                       100, -0.5, 799.5);
  Hist nChgND("number of charged particles, non-diffractive",
                                                       100, -0.5, 799.5);
  Hist pTnChg("<pt>(n_charged) all",                   100, -0.5, 799.5);
  Hist pTnChgSD("<pt>(n_charged) single diffraction",  100, -0.5, 799.5);
  Hist pTnChgDD("<pt>(n_charged) double diffraction",  100, -0.5, 799.5);
  Hist pTnChgCD("<pt>(n_charged) central diffraction", 100, -0.5, 799.5);
  Hist pTnChgND("<pt>(n_charged) non-diffractive   ",  100, -0.5, 799.5);

  // Book histograms: ditto as function of separate subsystem mass.
  Hist mLogInel("log10(mass), by diffractive system",  100, 0., 5.);
  Hist nChgmLog("<n_charged>(log10(mass))",            100, 0., 5.);
  Hist pTmLog("<pT>_charged>(log10(mass))",            100, 0., 5.);

  // Book histograms: elastic/diffractive.
  Hist tSpecEl("elastic |t| spectrum",              100, 0., 1.);
  Hist tSpecElLog("elastic log10(|t|) spectrum",    100, -5., 0.);
  Hist tSpecSD("single diffractive |t| spectrum",   100, 0., 2.);
  Hist tSpecDD("double diffractive |t| spectrum",   100, 0., 5.);
  Hist tSpecCD("central diffractive |t| spectrum",  100, 0., 5.);
  Hist mSpec("diffractive mass spectrum",           100, 0., 100.);
  Hist mLogSpec("log10(diffractive mass spectrum)", 100, 0., 4.);

  // Book histograms: inelastic nondiffractive.
  double pTmax = 20.;
  double bMax  = 4.;
  Hist pTspec("total pT_hard spectrum",             100, 0., pTmax);
  Hist pTspecND("nondiffractive pT_hard spectrum",  100, 0., pTmax);
  Hist bSpec("b impact parameter spectrum",         100, 0., bMax);
  Hist enhanceSpec("b enhancement spectrum",        100, 0., 10.);
  Hist number("number of interactions",             100, -0.5, 99.5);
  Hist pTb1("pT spectrum for b < 0.5",              100, 0., pTmax);
  Hist pTb2("pT spectrum for 0.5 < b < 1",          100, 0., pTmax);
  Hist pTb3("pT spectrum for 1 < b < 1.5",          100, 0., pTmax);
  Hist pTb4("pT spectrum for 1.5 < b",              100, 0., pTmax);
  Hist bpT1("b spectrum for pT < 2",                100, 0., bMax);
  Hist bpT2("b spectrum for 2 < pT < 5",            100, 0., bMax);
  Hist bpT3("b spectrum for 5 < pT < 15",           100, 0., bMax);
  Hist bpT4("b spectrum for 15 < pT",               100, 0., bMax);


  // Begin event loop.
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if too many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n";
      break;
    }


    // Extract event classification.
    int code = pythia.info.code();
    
    // Charged multiplicity and mean pT: all and by event class.
    int nch = 0;
    double pTsum = 0.;
    for (int i = 1; i < event.size(); ++i)
    if (event[i].isFinal() && event[i].isCharged()) {
      yChg.fill( event[i].y() );
      ++nch;
      pTsum += event[i].pT();
    }
    nChg.fill( nch );
    if (nch > 0) pTnChg.fill( nch, pTsum/nch);
    if (code == 103 || code == 104) {
      nChgSD.fill( nch );
      if (nch > 0) pTnChgSD.fill( nch, pTsum/nch);
    } else if (code == 105) {
      nChgDD.fill( nch );
      if (nch > 0) pTnChgDD.fill( nch, pTsum/nch);
    } else if (code == 106) {
      nChgCD.fill( nch );
      if (nch > 0) pTnChgCD.fill( nch, pTsum/nch);
    } else if (code == 101) {
      nChgND.fill( nch );
      if (nch > 0) pTnChgND.fill( nch, pTsum/nch);
      double mLog = log10( event[0].m() );
      mLogInel.fill( mLog );
      nChgmLog.fill( mLog, nch );
      if (nch > 0) pTmLog.fill( mLog, pTsum / nch );
    }

    // Charged multiplicity and mean pT: per diffractive system.
    for (int iDiff = 0; iDiff < 3; ++iDiff)
    if ( (iDiff == 0 && pythia.info.isDiffractiveA())
      || (iDiff == 1 && pythia.info.isDiffractiveB())
      || (iDiff == 2 && pythia.info.isDiffractiveC()) ) {
      int ndiff = 0;
      double pTdiff = 0.;
      int nDoc = (iDiff < 2) ? 4 : 5;
      for (int i = nDoc + 1; i < event.size(); ++i)
      if (event[i].isFinal() && event[i].isCharged()) {
        // Trace back final particle to see which system it comes from.
        int k = i;
        do k = event[k].mother1();
        while (k > nDoc);
        if (k == iDiff + 3) {
          ++ndiff;
          pTdiff += event[i].pT();
        }
      }
      // Study diffractive mass spectrum.
      double mDiff = event[iDiff+3].m();
      double mLog  = log10( mDiff);
      mLogInel.fill( mLog );
      nChgmLog.fill( mLog, ndiff );
      if (ndiff > 0) pTmLog.fill( mLog, pTdiff / ndiff );
      mSpec.fill( mDiff );
      mLogSpec.fill( mLog );
    }

    // Study pT spectrum of all hard collisions, no distinction.
    double pT = pythia.info.pTHat();
    pTspec.fill( pT );

    // Study t distribution of elastic/diffractive events.
    if (code > 101) {
      double tAbs = abs(pythia.info.tHat());
      if (code == 102) {
        tSpecEl.fill(tAbs);
        tSpecElLog.fill(log10(tAbs));
      }
      else if (code == 103 || code == 104) tSpecSD.fill(tAbs);
      else if (code == 105) tSpecDD.fill(tAbs);
      else if (code == 106) {

        // Anders: HACK! (To avoid operator error - result is nonsense)
        // double t1Abs = abs( (event[3].p() - event[1].p()).m2Calc() );
        // double t2Abs = abs( (event[4].p() - event[2].p()).m2Calc() );
        double t1Abs = abs( (event[3].p()).m2Calc() );
        double t2Abs = abs( (event[4].p()).m2Calc() );

        tSpecCD.fill(t1Abs);
        tSpecCD.fill(t2Abs);
      }

    // Study nondiffractive inelastic events in (pT, b) space.
    } else {
      double b = pythia.info.bMPI();
      double enhance = pythia.info.enhanceMPI();
      int nMPI = pythia.info.nMPI();
      pTspecND.fill( pT );
      bSpec.fill( b );
      enhanceSpec.fill( enhance );
      number.fill( nMPI );
      if (b < 0.5) pTb1.fill( pT );
      else if (b < 1.0) pTb2.fill( pT );
      else if (b < 1.5) pTb3.fill( pT );
      else pTb4.fill( pT );
      if (pT < 2.) bpT1.fill( b );
      else if (pT < 5.) bpT2.fill( b );
      else if (pT < 15.) bpT3.fill( b );
      else bpT4.fill( b );
    }

  // End of event loop.
  }

  // Final statistics and histograms.
  pythia.stat();
  pTnChg   /= nChg;
  pTnChgSD /= nChgSD;
  pTnChgDD /= nChgDD;
  pTnChgCD /= nChgCD;
  pTnChgND /= nChgND;
  nChgmLog /= mLogInel;
  pTmLog   /= mLogInel;


  //
  // Done
  // 

  cout << endl;
  cout << "=======================" << endl;
  cout << endl;

  // dlclose(pHandle);

  return 0;
}



