//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for the plc 3.0 backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2019 Aug
//
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/plc_3_0.hpp"

// Set the path to the plc_2.0 directory:
// This can either done at compile-time or during a scan via the rules-section
#ifndef PLC2_PATH
#define PLC2_PATH "EMPTY"
#endif

// Set the path to the plc_3.0 directory:
// This can either done at compile-time or during a scan via the rules-section
#ifndef PLC3_PATH
#define PLC3_PATH "EMPTY"
#endif

// Before we define the function, registered in the frontend-header let's define some macros,
// since the structure of these functions are redundant

// Check if we will calll "plc_loglike_NAME()" in the scan. If so, then run "initialize_NAME()"
// First check if the respective path (plc2_location or plc3_location) is already set.
#define CHECK_AND_ACTIVATE(NAME,VERSION)                                                 \
                                                                                         \
if(*InUse::CAT(plc_loglike_,NAME))                                                       \
{                                                                                        \
  if ( CAT_3(plc,VERSION,_location).empty() )                                            \
  {                                                                                      \
    CAT(set_planck_path_,VERSION) ( CAT_3(planck,VERSION,_path) );                       \
  }                                                                                      \
  CAT(initialize_,NAME) ();                                                              \
}                                                                                        \

// Definition of "void initialize_NAME()"  [Helper functions - Not registered in the header]
#define PLC_INITIALIZE(NAME,TYPE,VERSION,LOC,PRINT_NAME)                                 \
                                                                                         \
  void CAT(initialize_,NAME) ()                                                          \
  {                                                                                      \
    std::cout << "Loading the " << PRINT_NAME << " likelihood " << std::endl;            \
    std::string full_path = CAT_3(plc,VERSION,_location);                                \
    full_path += LOC;                                                                    \
    if (not Utils::file_exists(full_path))                                               \
    {                                                                                    \
      std::string ErrMssg = "Error while loading likelihood in plc_3.0:\n\nThe file \'"; \
      ErrMssg += full_path;                                                              \
      ErrMssg += "\' does not exist. Hence we will fail to load it.";                    \
      backend_error().raise(LOCAL_INFO,ErrMssg.c_str());                                 \
    }                                                                                    \
    char* clik_path = (char*)full_path.c_str();                                          \
    CAT(TYPE,_map)[STRINGIFY(NAME)] = CAT(TYPE,_init) (clik_path, NULL);                 \
  }

// Definition of "double plc_loglike_NAME(double* cl_and_pars)"
#define PLC_CLIK_LOGLIKE(NAME,TYPE)                                                      \
                                                                                         \
  double CAT(plc_loglike_,NAME)(double* cl_and_pars)                                     \
  {                                                                                      \
    const std::string name = STRINGIFY(NAME);                                            \
    logger() << "Calling \""<< STRINGIFY(CAT(plc_loglike_,NAME)) << "\":" << EOM;        \
    if (CAT(TYPE,_map).count(name) == 0)                                                 \
    {                                                                                    \
      std::string mssg = "Could not find the Likelihood object \"";                      \
      mssg += name;                                                                      \
      mssg += "\" in the map of activated likelihoods.";                                 \
      backend_error().raise(LOCAL_INFO,mssg.c_str());                                    \
    }                                                                                    \
    double res;                                                                          \
    clik_error* plc_Error = initError();                                                 \
    res =  CAT(TYPE,_compute) (CAT(TYPE,_map)[name], cl_and_pars, &plc_Error);           \
    if (isError(plc_Error))                                                              \
    {                                                                                    \
      std::string ErrMssg = "Calling \"";                                                \
      ErrMssg +=  std::string(STRINGIFY(CAT(plc_loglike_,NAME))) + "\" ";                \
      ErrMssg += "was not successful.\n\n";                                              \
      char forwardedErr[4096];                                                           \
      stringError(byVal(forwardedErr), plc_Error);                                       \
      ErrMssg += forwardedErr;                                                           \
      logger() << ErrMssg << "\n\nPoint will be invalidated.\n" << EOM;                  \
      invalid_point().raise(ErrMssg.c_str());                                            \
    }                                                                                    \
    else                                                                                 \
    {                                                                                    \
      logger() << "Calling \""<< STRINGIFY(CAT(plc_loglike_,NAME)) << "\"";              \
      logger() << " was successfull. Got " << res << EOM;                                \
    }                                                                                    \
    cleanupError(&plc_Error);                                                            \
    return res;                                                                          \
  }                                                                                      \

// Macros end here.

BE_NAMESPACE
{
  // Static varibales to be kept alive in the frontend for the scan.
  //  -> Prefix to the plc_2.0 data folder
  static std::string plc2_location;
  //  -> Prefix to the plc_3.0 data folder
  static std::string plc3_location;
  //  -> Container for all activated likelihood objects
  static std::map<std::string, clik_object* > clik_map;
  static std::map<std::string, clik_lensing_object* > clik_lensing_map;

  void set_planck_path_2(std::string& planck_path)
  {
    // Check if the path is set (e.g. it is different to "EMPTY")
    if (planck_path.compare("EMPTY") == 0)
    {
      std::string message = "The location of the planck data does not seem to be set properly!\n\n";
      message += "Either change the PLC2_PATH macro in Backends/include/... .../plc_3_0.hpp or include\n\n";
      message += "  - capability:   plc_3_0_init\n";
      message += "    options:\n";
      message += "      plc_2.0_path: /YOUR/PATH/TO/plc_2.0/\n\n";
      message += "in the rules section of your input.";
      backend_error().raise(LOCAL_INFO,message.c_str());
    }

    // If the path does not already end with "/", then add it.
    if (planck_path.back() != '/')
      planck_path.push_back('/');

    plc2_location = planck_path;
  }

  void set_planck_path_3(std::string& planck_path)
  {
    // Check if the path is set (e.g. it is different to "EMPTY")
    if (planck_path.compare("EMPTY") == 0)
    {
      std::string message = "The location of the planck data does not seem to be set properly!\n\n";
      message += "Either change the PLC3_PATH macro in Backends/include/... .../plc_3_0.hpp or include\n\n";
      message += "  - capability:   plc_3_0_init\n";
      message += "    options:\n";
      message += "      plc_3.0_path: /YOUR/PATH/TO/plc_3.0/\n\n";
      message += "in the rules section of your input.";
      backend_error().raise(LOCAL_INFO,message.c_str());
    }

    // If the path does not already end with "/", then add it.
    if (planck_path.back() != '/')
      planck_path.push_back('/');

    plc3_location = planck_path;
  }

  // "Write" the code for all "void initialize_NAME()"
  PLC_INITIALIZE(highl_TT_2015,clik,2,"hi_l/plik/plik_dx11dr2_HM_v18_TT.clik","high-l TT (PR2 - 2015)")
  PLC_INITIALIZE(highl_TTTEEE_2015,clik,2,"hi_l/plik/plik_dx11dr2_HM_v18_TTTEEE.clik","high-l TTTEEE (PR2 - 2015)")
  PLC_INITIALIZE(highl_TT_lite_2015,clik,2,"hi_l/plik_lite/plik_lite_v18_TT.clik","high-l TT-lite (PR2 - 2015)")
  PLC_INITIALIZE(highl_TTTEEE_lite_2015,clik,2,"hi_l/plik_lite/plik_lite_v18_TTTEEE.clik","high-l TTTEEE-lite (PR2 - 2015)")
  PLC_INITIALIZE(lowl_TEB_2015,clik,2,"low_l/bflike/lowl_SMW_70_dx11d_2014_10_03_v5c_Ap.clik","low-l TEB (PR2 - 2015)")
  PLC_INITIALIZE(lowl_TT_2015,clik,2,"low_l/commander/commander_rc2_v1.1_l2_29_B.clik","low-l TT (PR2 - 2015)")
  PLC_INITIALIZE(lensing_2015,clik_lensing,2,"lensing/smica_g30_ftl_full_pp.clik_lensing","lensing (PR2 - 2015)")

  PLC_INITIALIZE(highl_TT_2018,clik,3,"hi_l/plik/plik_rd12_HM_v22_TT.clik","high-l TT (PR3 - 2018)")
  PLC_INITIALIZE(highl_TTTEEE_2018,clik,3,"hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik","high-l TTTEEE (PR3 - 2018)")
  PLC_INITIALIZE(highl_TT_lite_2018,clik,3,"hi_l/plik_lite/plik_lite_v22_TT.clik","high-l TT-lite (PR3 - 2018)")
  PLC_INITIALIZE(highl_TTTEEE_lite_2018,clik,3,"hi_l/plik_lite/plik_lite_v22_TTTEEE.clik","high-l TTTEEE-lite (PR3 - 2018)")
  PLC_INITIALIZE(lowl_TT_2018,clik,3,"low_l/commander/commander_dx12_v3_2_29.clik","low-l TT (PR3 - 2018)")
  PLC_INITIALIZE(lowl_EE_2018,clik,3,"low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik","low-l EE (PR3 - 2018)")
  PLC_INITIALIZE(lensing_2018,clik_lensing,3,"lensing/smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8.clik_lensing","lensing (PR3 - 2018)")
  PLC_INITIALIZE(lensing_marged_2018,clik_lensing,3,"lensing/smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8_CMBmarged.clik_lensing","marginalized lensing (PR3 - 2018)")

  // "Write" the code for all "double plc_loglike_NAME(double* cl_and_pars)"
  PLC_CLIK_LOGLIKE(highl_TTTEEE_2015,clik)
  PLC_CLIK_LOGLIKE(highl_TT_2015,clik)
  PLC_CLIK_LOGLIKE(highl_TTTEEE_lite_2015,clik)
  PLC_CLIK_LOGLIKE(highl_TT_lite_2015,clik)
  PLC_CLIK_LOGLIKE(lowl_TEB_2015,clik)
  PLC_CLIK_LOGLIKE(lowl_TT_2015,clik)
  PLC_CLIK_LOGLIKE(lensing_2015,clik_lensing)

  PLC_CLIK_LOGLIKE(highl_TTTEEE_2018,clik)
  PLC_CLIK_LOGLIKE(highl_TT_2018,clik)
  PLC_CLIK_LOGLIKE(highl_TTTEEE_lite_2018,clik)
  PLC_CLIK_LOGLIKE(highl_TT_lite_2018,clik)
  PLC_CLIK_LOGLIKE(lowl_TT_2018,clik)
  PLC_CLIK_LOGLIKE(lowl_EE_2018,clik)
  PLC_CLIK_LOGLIKE(lensing_2018,clik_lensing)
  PLC_CLIK_LOGLIKE(lensing_marged_2018,clik_lensing)
}
END_BE_NAMESPACE

BE_INI_FUNCTION
{
  // Most of the Backend initialisation is only relevant for the first parameter point
  static bool scan_level = true;
  if (scan_level)
  {
    std::string planck2_path;
    std::string planck3_path;

    // Check if "plc_prefix" is in the runOptions. If so, then we assume that
    // the folders "plc_2.0" and "plc_3.0" are contained in it.
    // If this is not the case we will check for "plc_2.0_path" and "plc_3.0_path"
    // in case they have completely different locations.
    if (runOptions->hasKey("plc_prefix"))
    {
      std::string plc_pre = runOptions->getValue<std::string>("plc_prefix");
      if (plc_pre.back() != '/')
        plc_pre.push_back('/');
      planck2_path = plc_pre + "plc_2.0/";
      planck3_path = plc_pre + "plc_3.0/";
    }
    else
    {
      planck2_path = runOptions->getValueOrDef<std::string>(PLC2_PATH,"plc_2.0_path");
      planck3_path = runOptions->getValueOrDef<std::string>(PLC3_PATH,"plc_3.0_path");
    }

    // Check whcih lilelihood is used and initialise it when needed.
    CHECK_AND_ACTIVATE(highl_TTTEEE_2015,2)
    CHECK_AND_ACTIVATE(highl_TT_2015,2)
    CHECK_AND_ACTIVATE(highl_TTTEEE_lite_2015,2)
    CHECK_AND_ACTIVATE(highl_TT_lite_2015,2)
    CHECK_AND_ACTIVATE(lowl_TEB_2015,2)
    CHECK_AND_ACTIVATE(lowl_TT_2015,2)
    CHECK_AND_ACTIVATE(lensing_2015,2)

    CHECK_AND_ACTIVATE(highl_TTTEEE_2018,3)
    CHECK_AND_ACTIVATE(highl_TT_2018,3)
    CHECK_AND_ACTIVATE(highl_TTTEEE_lite_2018,3)
    CHECK_AND_ACTIVATE(highl_TT_lite_2018,3)
    CHECK_AND_ACTIVATE(lowl_TT_2018,3)
    CHECK_AND_ACTIVATE(lowl_EE_2018,3)
    CHECK_AND_ACTIVATE(lensing_2018,3)
    CHECK_AND_ACTIVATE(lensing_marged_2018,3)
  }
  scan_level = false;
}
END_BE_INI_FUNCTION

// undefine all macros
#undef CHECK_AND_ACTIVATE
#undef PLC_CLIK_LOGLIKE
#undef PLC_INITIALIZE

//#undef PLC2_PATH
//#undef PLC3_PATH
