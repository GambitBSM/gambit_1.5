//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for DarkSUSY_generic_wimp_6.1.1 backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Torsten Bringmann
///          (torsten.bringmann@fys.uio.no)
///  \date 2020
///
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/DarkSUSY_generic_wimp_6_1_1.hpp"
#include "gambit/Utils/file_lock.hpp"
#include "gambit/Utils/mpiwrapper.hpp"

//#define DARKSUSY_DEBUG

// Initialisation function (definition)
BE_INI_FUNCTION
{
  // Initialize DarkSUSY if run for the first time
  bool static scan_level = true;

  if (scan_level)
  {

    // Do the call to dsinit one-by-one for each MPI process, as DarkSUSY loads up
    // HiggsBounds, which writes files at init then reads them back in later.
    Utils::FileLock mylock("DarkSUSY_" STRINGIFY(SAFE_VERSION) "_init");
    mylock.get_lock();
    dsinit();
    mylock.release_lock();

    //// Initialize yield tables for use in cascade decays (initialize more if needed)
    // This makes sure that different processes later don't read the yield tables
    // from disk simultaneously
    int istat=0;
    double mdm=100.0, egev=10.0;
    int pdg=5, yieldpdg, diff=1;
    char*hel =  (char *)"0";
    
    yieldpdg = 22;// gamma rays
    dsanyield_sim(mdm,egev,pdg,hel,yieldpdg,diff,istat);
    yieldpdg = -2212; //antiprotons
    dsanyield_sim(mdm,egev,pdg,hel,yieldpdg,diff,istat);
    yieldpdg = -11; //positrons
    dsanyield_sim(mdm,egev,pdg,hel,yieldpdg,diff,istat);

    scan_level = false;

  }

// unlike in the MSSM case, we don't actually initialize a generic WIMP
// model "point" here. We thus just use basic, "model-independent"
// DarkSUSY functinalities.
    
    

}
END_BE_INI_FUNCTION

// Convenience functions (definitions)
BE_NAMESPACE
{
   // No specific convenience functions added so far, as this backend is
   // currently only used in DarkBit_standalone. Eventually, as more
   // functionality is requested, such convenienc functions may have to be added.

}
END_BE_NAMESPACE
