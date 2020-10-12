///  Functions of module NuclearBit
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Anders Kvellestad
///          (anders.kvellestad@fys.uio.no)
///  \date 2020 Sep
///
///  *********************************************

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/NuclearBit/NuclearBit_rollcall.hpp"
#include "gambit/Models/claw_singleton.hpp"

namespace Gambit
{

  namespace NuclearBit
  {
    using namespace LogTags;

    //************************************************************

    /// \name Module functions
    /// @{

    // Run the gledeli backend and retrieve all the results
    void getGledeliResults(std::map<std::string,double> &result)
    {
      using namespace Pipes::getGledeliResults;

      // The first time this function is run we let
      // gledeli know which GAMBIT models are in use.
      static bool first = true;
      if (first)
      {
        // Get a set<str> with the names of active GAMBIT models and cast it to a pybind11::list
        pybind11::list active_models_list = pybind11::cast(Models::ModelDB().get_activemodels());
        // Pass the list to the gledeli backend.
        BEreq::gledeliBE_set_model_names(active_models_list);
        first = false;
      }

      // Grab parameter values from the GAMBIT Param map (type map<str,double*>)
      // and populate a pybind11::dict that we pass to the gledeli backend
      pybind11::dict pars;
      for (auto& kv : Param)
      {
        pars[kv.first.c_str()] = *kv.second;
      }

      // Pass the parameters to the gledeli backend
      BEreq::gledeliBE_set_model_pars(pars);

      // Some dummy settings that can vary from point to point
      static int setting_A = 0;
      setting_A++;
      static int setting_B = setting_A * 10;
      setting_B++;

      pybind11::dict settings;
      settings["setting_A"] = setting_A;
      settings["setting_B"] = setting_B;

      // Now run the gledeli backend for this point with the above settings
      pybind11::dict gledeli_run_info = BEreq::gledeliBE_run(settings);

      // Exctract info from the gledeli_run_info dict
      if( !gledeli_run_info.contains("success") ) { NuclearBit_error().raise(LOCAL_INFO, "The expected 'success' key is not found in the pybind11:dict 'gledeli_run_info'"); }
      bool success = gledeli_run_info["success"].cast<bool>();

      if( !gledeli_run_info.contains("log_msg") ) { NuclearBit_error().raise(LOCAL_INFO, "The expected 'log_msg' key is not found in the pybind11:dict 'gledeli_run_info'"); }
      std::string log_msg = gledeli_run_info["log_msg"].cast<std::string>();

      // Pass log message to the GAMBIT logger
      logger() << log_msg << EOM;

      // If the calculation failed, invalidate this parameter point
      if(!success)
      {
          invalid_point().raise("The gledeli backend reported an unsuccessful calculation. Will invalidate this parameter point.");
      }

      // Get the result (map<str,double>) from the gledeli backend
      pybind11::dict gledeli_results = BEreq::gledeliBE_get_results();

      // std::map<std::string,double> gledeli_results& = BEreq::gledeliBE_get_results();
      // gledeli_results = BEreq::gledeliBE_get_results();
      // result = gledeli_results;

      // Fill the result map (map<str,double>) of this function
      for (auto& kv : gledeli_results)
      {
        result[kv.first.cast<std::string>()] = kv.second.cast<double>();
      }
    }


    // Extract the log-likelihood from the gledeli results
    // so it can be added to the GAMBIT scan likelihood
    void getGledeliLogLike(double &result)
    {
      using namespace Pipes::getGledeliLogLike;

      // Get the gledeli results
      const std::map<std::string,double>& gledeliBE_results = *Dep::gledeliResults;

      // Extract the loglike entry and return that as the result
      if(gledeliBE_results.count("loglike") > 0)
      {
        result = gledeliBE_results.at("loglike");
      }
      else
      {
        NuclearBit_error().raise(LOCAL_INFO, "The required 'loglike' key was not found in the provided map!");
      }

    }


    /// @}

  }

}

