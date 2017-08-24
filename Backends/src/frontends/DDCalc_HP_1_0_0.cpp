//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend source for the DDCalc_HP backend.
///
///  Actual implementation of DDCalc_HP ini function.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Ankit Beniwal
///          (ankit.beniwal@adelaide.edu.au)
///  \date 2017 Aug
///
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/DDCalc_HP_1_0_0.hpp"
#include "gambit/Utils/standalone_error_handlers.hpp"

#include <map>
#include <stdexcept>

// File-local globals
BE_NAMESPACE
{
  // Map returning detector index given an analysis name.
  std::map<str,int> ex_map;
  // DM model and halo model singleton indices.
  int WIMP, Halo;
}
END_BE_NAMESPACE


// Initialisation function
BE_INI_FUNCTION
{

  // Halo model parameters and pointers to their entries in the Params map.
  static double rho0_eff = 0.4, vrot = 235, v0 = 235, vesc = 550;
  static safe_ptr<LocalMaxwellianHalo> LocalHaloParameters_ptr;

  // Fraction of DM
  double fraction = *Dep::RD_fraction;

  // Scan-level initialization -----------------------------
  static bool scan_level = true;
  if (scan_level)
  {
    // Initialize halo and WIMP models
    WIMP = DDCalc_InitWIMP();
    Halo = DDCalc_InitHalo();

    // Initialize experiments
    if (*InUse::DDCalc_Experiment)
    {
      ex_map["XENON100_2012"] = XENON100_2012_Init(false);
      ex_map["LUX_2013"] = LUX_2013_Init(false);
      ex_map["SuperCDMS_2014"] = SuperCDMS_2014_Init(false);
      ex_map["SIMPLE_2014"] = SIMPLE_2014_Init(false);
      ex_map["LUX_2016"] = LUX_2016_Init(false);
      ex_map["PandaX_2016"] = PandaX_2016_Init(false);
      ex_map["LUX_2015"] = LUX_2015_Init(false);
      ex_map["PICO_2L"] = PICO_2L_Init(false);
      ex_map["PICO_60_F"] = PICO_60_F_Init(false);
      ex_map["PICO_60_I"] = PICO_60_I_Init(false);
      //ex_map["DARWIN_Ar"] = DARWIN_Ar_Init(false);
      //ex_map["DARWIN_Xe"] = DARWIN_Xe_Init(false);
    }

    // Save safe pointers to local halo parameters.
    LocalHaloParameters_ptr = Dep::LocalHalo.safe_pointer();
  }
  scan_level = false;

  // Point-level initialization ----------------------------

  // Change DM parameters
  // Expected input is: f = G/2 where G is the effective 4 fermion DM-nucleon coupling
  DDCalc_SetWIMP_higgsportal(WIMP,*Dep::mwimp,(Dep::DD_couplings->gps)/2,(Dep::DD_couplings->gns)/2,
                    (Dep::DD_couplings->gpa)/2,(Dep::DD_couplings->gna)/2);

  // Change halo parameters.
    bool halo_changed = false;

    if (LocalHaloParameters_ptr->rho0 * fraction != rho0_eff) {rho0_eff = LocalHaloParameters_ptr->rho0 * fraction; halo_changed = true;}
    if (LocalHaloParameters_ptr->vrot != vrot)                {vrot     = LocalHaloParameters_ptr->vrot;            halo_changed = true;}
    if (LocalHaloParameters_ptr->v0   != v0)                  {v0       = LocalHaloParameters_ptr->v0;              halo_changed = true;}
    if (LocalHaloParameters_ptr->vesc != vesc)                {vesc     = LocalHaloParameters_ptr->vesc;            halo_changed = true;}

    if (halo_changed)
    {
      DDCalc_SetSHM(Halo,rho0_eff,vrot,v0,vesc);

      // Log stuff if in debug mode
      #ifdef DDCALC_DEBUG
        logger() << "Updated DDCalc_HP halo parameters:" << EOM;
        logger() << "    rho0 [GeV/cm^3]     = " << LocalHaloParameters_ptr->rho0 << EOM;
        logger() << "    rho0_eff [GeV/cm^3] = " << rho0_eff << EOM;
        logger() << "    vrot [km/s]         = " << vrot << EOM;
        logger() << "    v0   [km/s]         = " << v0   << EOM;
        logger() << "    vesc [km/s]         = " << vesc << EOM;
      #endif
    }

}
END_BE_INI_FUNCTION

// Convenience functions
BE_NAMESPACE
{
  // Convenience function for returning detector index given an analysis name.
  int DDCalc_Experiment(const str& ex)
  {
    int result = -1;
    try { result = ex_map.at(ex); }
    catch(std::out_of_range) { backend_error().raise(LOCAL_INFO, "Unknown experiment requested from DDCalc_HP."); }
    return result;
  }

  // Convenience function for calling CalcRates with internally-initialised WIMP and halo objects.
  void DDCalc_CalcRates_simple(const int& D) { DDCalc_CalcRates(D, WIMP, Halo); }
}
END_BE_NAMESPACE
