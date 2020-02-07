//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for the plc 2.0 backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2019 July, Aug, Nov
///  \date 2020 Feb
//
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/plc_2_0.hpp"

// Set the path to the plc_2.0 directory:
// This can either done at compile-time or during a scan via the rules-section
#ifndef PLC2_PATH
#define PLC2_PATH "EMPTY"
#endif

// Before we define the function, registered in the frontend-header let's define some macros,
// since the structure of these functions are redundant

// Check if we will calll "plc_loglike_NAME()" in the scan. If so, then run "initialize_NAME()"
#define CHECK_AND_ACTIVATE(NAME)                                                                          \
                                                                                                          \
if(*InUse::CAT(plc_loglike_,NAME))                                                                        \
{                                                                                                         \
  CAT(initialize_,NAME) ();                                                                               \
}                                                                                                         \

// Definition of "void initialize_NAME()"  [Helper functions - Not registered in the header]
// The wildcard "VERSION" is not relevant for this version of the frontend
#define PLC_INITIALIZE(NAME,TYPE,VERSION,LOC,PRINT_NAME)                                                  \
                                                                                                          \
  void CAT(initialize_,NAME) ()                                                                           \
  {                                                                                                       \
    std::cout << "Loading the " << PRINT_NAME << " likelihood " << std::endl;                             \
    std::string full_path = plc_location;                                                                 \
    full_path += LOC;                                                                                     \
    if (not Utils::file_exists(full_path))                                                                \
    {                                                                                                     \
      std::string ErrMssg = "Error while loading likelihood in plc_2.0:\n\nThe file \'";                  \
      ErrMssg += full_path;                                                                               \
      ErrMssg += "\' does not exist. Hence we will fail to load it.\n(Hint: Is PLC2_PATH set properly?)"; \
      backend_error().raise(LOCAL_INFO,ErrMssg.c_str());                                                  \
    }                                                                                                     \
    char* clik_path = (char*)full_path.c_str();                                                           \
    clik_error* plc_Error = initError();                                                                  \
    auto destructor = []( CAT(TYPE,_object)* clikid){ CAT(TYPE,_cleanup)(&(clikid)); };                   \
    CAT(TYPE,_map)[STRINGIFY(NAME)] =                                                                     \
      std::shared_ptr<CAT(TYPE,_object)>(CAT(TYPE,_init)(clik_path, &plc_Error), destructor);             \
    cleanupError(&plc_Error);                                                                             \
  }

// Definition of "double plc_loglike_NAME(double* cl_and_pars)"
#define PLC_CLIK_LOGLIKE(NAME,TYPE)                                                                       \
                                                                                                          \
  double CAT(plc_loglike_,NAME)(double* cl_and_pars)                                                      \
  {                                                                                                       \
    const std::string name = STRINGIFY(NAME);                                                             \
    logger() << "Calling \""<< STRINGIFY(CAT(plc_loglike_,NAME)) << "\":" << EOM;                         \
    if (CAT(TYPE,_map).count(name) == 0)                                                                  \
    {                                                                                                     \
      std::string mssg = "Could not find the Likelihood object \"";                                       \
      mssg += name;                                                                                       \
      mssg += "\" in the map of activated likelihoods.";                                                  \
      backend_error().raise(LOCAL_INFO,mssg.c_str());                                                     \
    }                                                                                                     \
    double res;                                                                                           \
    clik_error* plc_Error = initError();                                                                  \
    res =  CAT(TYPE,_compute) (CAT(TYPE,_map)[name].get(), cl_and_pars, &plc_Error);                      \
    if (isError(plc_Error))                                                                               \
    {                                                                                                     \
      std::string ErrMssg = "Calling \"";                                                                 \
      ErrMssg +=  std::string(STRINGIFY(CAT(plc_loglike_,NAME))) + "\" ";                                 \
      ErrMssg += "was not successful.\n\n";                                                               \
      char forwardedErr[4096];                                                                            \
      stringError(byVal(forwardedErr), plc_Error);                                                        \
      ErrMssg += forwardedErr;                                                                            \
      logger() << ErrMssg << "\n\nPoint will be invalidated." << EOM;                                     \
      invalid_point().raise(ErrMssg.c_str());                                                             \
    }                                                                                                     \
    else                                                                                                  \
    {                                                                                                     \
      logger() << "Calling \""<< STRINGIFY(CAT(plc_loglike_,NAME)) << "\"";                               \
      logger() << " was successfull. Got " << res << EOM;                                                 \
    }                                                                                                     \
    cleanupError(&plc_Error);                                                                             \
    return res;                                                                                           \
  }                                                                                                       \

// Macros end here.

BE_NAMESPACE
{
  // Static varibales to be kept alive in the frontend for the scan.
  //  -> Prefix to the plc_2.0 data folder
  static std::string plc_location;
  //  -> Container for all activated likelihood objects
  static std::map<std::string, std::shared_ptr<clik_object> > clik_map;
  static std::map<std::string, std::shared_ptr<clik_lensing_object> > clik_lensing_map;
  //  -> Array containing the required lmax.
  //     The order is [phiphi TT EE BB TE TB EB]
  static std::array<int,7> lmax_array{-1,-1,-1,-1,-1,-1,-1};

  void set_planck_path(std::string& planck_path)
  {
    // Check if the path is set (e.g. it is different to "EMPTY")
    if (planck_path.compare("EMPTY") == 0)
    {
      std::string message = "The location of the planck data does not seem to be set properly!\n\n";
      message += "Either change the PLC2_PATH macro in Backends/include/... .../plc_2_0.hpp or include\n\n";
      message += "  - capability:   plc_2_0_init\n";
      message += "    options:\n";
      message += "      plc_2.0_path: /YOUR/PATH/TO/plc_2.0/\n\n";
      message += "in the rules section of your input.";
      backend_error().raise(LOCAL_INFO,message.c_str());
    }

    // If the path does not already end with "/", then add it.
    if (planck_path.back() != '/')
      planck_path.push_back('/');

    plc_location = planck_path;
  }

  // "Write" the code for all "void initialize_NAME()"
  PLC_INITIALIZE(highl_TT_2015,clik,2,"hi_l/plik/plik_dx11dr2_HM_v18_TT.clik","high-l TT (PR2 - 2015)")
  PLC_INITIALIZE(highl_TTTEEE_2015,clik,2,"hi_l/plik/plik_dx11dr2_HM_v18_TTTEEE.clik","high-l TTTEEE (PR2 - 2015)")
  PLC_INITIALIZE(highl_TT_lite_2015,clik,2,"hi_l/plik_lite/plik_lite_v18_TT.clik","high-l TT-lite (PR2 - 2015)")
  PLC_INITIALIZE(highl_TTTEEE_lite_2015,clik,2,"hi_l/plik_lite/plik_lite_v18_TTTEEE.clik","high-l TTTEEE-lite (PR2 - 2015)")
  PLC_INITIALIZE(lowl_TEB_2015,clik,2,"low_l/bflike/lowl_SMW_70_dx11d_2014_10_03_v5c_Ap.clik","low-l TEB (PR2 - 2015)")
  PLC_INITIALIZE(lowl_TT_2015,clik,2,"low_l/commander/commander_rc2_v1.1_l2_29_B.clik","low-l TT (PR2 - 2015)")
  PLC_INITIALIZE(lensing_2015,clik_lensing,2,"lensing/smica_g30_ftl_full_pp.clik_lensing","lensing (PR2 - 2015)")

  // "Write" the code for all "double plc_loglike_NAME(double* cl_and_pars)"
  PLC_CLIK_LOGLIKE(highl_TTTEEE_2015,clik)
  PLC_CLIK_LOGLIKE(highl_TT_2015,clik)
  PLC_CLIK_LOGLIKE(highl_TTTEEE_lite_2015,clik)
  PLC_CLIK_LOGLIKE(highl_TT_lite_2015,clik)
  PLC_CLIK_LOGLIKE(lowl_TEB_2015,clik)
  PLC_CLIK_LOGLIKE(lowl_TT_2015,clik)
  PLC_CLIK_LOGLIKE(lensing_2015,clik_lensing)

  void plc_required_Cl(int& lmax, bool& needs_tCl, bool& needs_pCl)
  {
    lmax = -1;
    for (const auto& it: lmax_array)
    {
      lmax = std::max(lmax, it);
    }

    needs_tCl = (lmax_array[1] > -1);

    auto begin = lmax_array.begin() + 2;
    auto end = lmax_array.end();
    needs_pCl = std::any_of(begin,end,[](int& i){return i > -1;});
  }
}
END_BE_NAMESPACE

BE_INI_FUNCTION
{
  // Most of the Backend initialisation is only relevant for the first parameter point
  static bool scan_level = true;
  if (scan_level)
  {
    // Set the path to the likelihood data.
    std::string planck_path = runOptions->getValueOrDef<std::string>(PLC2_PATH,"plc_2.0_path");
    // Backward compability to previous implementation
    if (runOptions->hasKey("planck_path"))
      planck_path = runOptions->getValue<std::string>("planck_path");
    set_planck_path(planck_path);

    // Check whcih lilelihood is used and initialise it when needed.
    CHECK_AND_ACTIVATE(highl_TTTEEE_2015)
    CHECK_AND_ACTIVATE(highl_TT_2015)
    CHECK_AND_ACTIVATE(highl_TTTEEE_lite_2015)
    CHECK_AND_ACTIVATE(highl_TT_lite_2015)
    CHECK_AND_ACTIVATE(lowl_TEB_2015)
    CHECK_AND_ACTIVATE(lowl_TT_2015)
    CHECK_AND_ACTIVATE(lensing_2015)

    for (const auto& it: clik_map)
    {
      std::array<int,7> tmp{};

      clik_get_lmax(it.second.get(), tmp.data()+1, NULL);

      for (size_t i = 1; i < 7; ++i)
      {
        lmax_array[i] = std::max(lmax_array[i],tmp[i]);
      }
    }

    for (const auto& it: clik_lensing_map)
    {
      std::array<int,7> tmp{};

      clik_lensing_get_lmaxs(it.second.get(), tmp.data(), NULL);

      for (size_t i = 0; i < 7; ++i)
      {
        lmax_array[i] = std::max(lmax_array[i],tmp[i]);
      }
    }
  }
  scan_level = false;
}
END_BE_INI_FUNCTION

// undefine all macros
#undef CHECK_AND_ACTIVATE
#undef PLC_CLIK_LOGLIKE
#undef PLC_INITIALIZE

//#undef PLC2_PATH
