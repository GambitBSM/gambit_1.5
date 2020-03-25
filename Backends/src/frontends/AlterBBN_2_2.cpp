//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for AlterBBN backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
/// \author Janina Renk
///         (janina.renk@fysik.su.se)
/// \date 2018 Jun
/// \date 2020 Jan
///
/// \author Patrick Stöcker
///         (stoecker@physik.rwth-achen.de)
/// \date 2019 Sep
///
/// \author Will Handley
///         (wh260@cam.ac.uk)
/// \date 2020 Mar
///
///  *********************************************

#include <sstream>
#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/AlterBBN_2_2.hpp"

#define NNUC 26 // number of element abundances computed in AlterBBN 2.2

// Initialisation
BE_INI_FUNCTION{


}
END_BE_INI_FUNCTION


// Convenience functions (definitions)
BE_NAMESPACE
{
  static int                 nucl_err_res=0;
  static map_str_dbl         prev_AlterBBN_input{};
  static std::vector<double> prev_ratioH(0., NNUC+1);
  static std::vector<double> prev_cov_ratioH(0., (NNUC+1)*(NNUC+1));

  /// string set containing the name of all members of the AlterBBN relicparam structures that can currently 
  /// be set with GAMBIT. If you add a new model and need to pass a different option to AlterBBN add e.g. "life_neutron" 
  /// here and in the function fill_cosmomodel below to modify the lifetime of the neutron
  std::set<std::string> known_relicparam_options = {"eta0", "Nnu", "dNnu", "err","failsafe"};

  /// pass all values of AlterBBN's relicparam structure that  
  void fill_cosmomodel(AlterBBN::AlterBBN_2_2::relicparam * input_relicparam, map_str_dbl & AlterBBN_input)
  {
    // check that only options that are known to the interface, i.e. that they are one of the options 
    // defined in the string set 'known_relicparam_options' above, are passed to AlterBBN
    static bool first_run = true;
    if(first_run)
    {
      std::map<std::string, double>::const_iterator it;
      for ( it = AlterBBN_input.begin(); it != AlterBBN_input.end(); it++ )
      {
        // throw error if an option is not currently being passed on to AlterBBN & tell user how to change that
        if(known_relicparam_options.count(it->first)==0)
        {
          std::ostringstream errormsg;
          errormsg << "Option '" << it->first << "' tried to pass to AlterBBN_" << VERSION<< " is not know to the relicparam strucutre "<<std::endl;
          errormsg << "Known options are:";
          for(auto item : known_relicparam_options)
          {
                  errormsg <<" '"<<item << "'";
          }
          errormsg << ". Check for typos or add new one to the string set 'known_relicparam_options' in Backends/src/frontends/AlterBBN_"<< VERSION<< ".cpp"<<std::endl;
          backend_error().raise(LOCAL_INFO,errormsg.str());
        }
      }
      first_run = false;
    }

    // if a value for a member of the relicparam structure has been added to the AlterBBN_input map 
    // set the value accordingly
    if (AlterBBN_input.count("eta0")){input_relicparam->eta0 = AlterBBN_input["eta0"];}
    if (AlterBBN_input.count("Nnu")){input_relicparam->Nnu = AlterBBN_input["Nnu"];}
    if (AlterBBN_input.count("dNnu")){input_relicparam->dNnu = AlterBBN_input["dNnu"];}
    // to overwrite the default value of another relicparam structure member, e.g. life_neutron,
    // to AlterBBN uncomment the line below
    // if (AlterBBN_input.count("life_neutron")){input_relicparam->life_neutron = AlterBBN_input["life_neutron"];}

    // set error handling related parameters
    if (AlterBBN_input.count("failsafe")){input_relicparam->failsafe = (int)AlterBBN_input["failsafe"];}
    if (AlterBBN_input.count("err")){input_relicparam->err = (int)AlterBBN_input["err"];}
  }
  
  /// calls the AlterBBN routine nucl_err with the filled relicparam structure. This will fill the array ratioH with 
  /// all computed element abundances, and cov_ratioH with their errors & covariances
  int call_nucl_err(map_str_dbl &AlterBBN_input, double* ratioH, double* cov_ratioH )
  {
    if (AlterBBN_input != prev_AlterBBN_input)
    {
      AlterBBN::AlterBBN_2_2::relicparam input_relicparam;
      Init_cosmomodel(&input_relicparam); // initialise valuse of relicparam structure to their defaults
      fill_cosmomodel(&input_relicparam, AlterBBN_input); // fill strucutre with values contained in AlerBBN_input map which is filled in CosmoBit.ccp for different models

      nucl_err_res = nucl_err(&input_relicparam, ratioH, cov_ratioH);
      //int bbn_exc = bbn_excluded_chi2(&input_relicparam); just for testing purposes to compare internal AlterBBN output to GAMBIT result
      prev_ratioH = std::vector<double>(ratioH, ratioH+NNUC+1);
      prev_cov_ratioH = std::vector<double>(cov_ratioH, cov_ratioH+(NNUC+1)*(NNUC+1));
      prev_AlterBBN_input = AlterBBN_input;
    }
    else
    {
      for (int i=0;i<NNUC+1;i++) ratioH[i] = prev_ratioH[i];
      for (int i=0;i<(NNUC+1)*(NNUC+1);i++) cov_ratioH[i] = prev_cov_ratioH[i];
    }
    return nucl_err_res;
  }
 
  /// return the NNUC -- global parameter of AlterBBN specifying the number of 
  /// elements for which abundances are calculated -> length of array ratioH is NNUC+1 (AlterBBN
  /// starts filling @ position 1)
  int get_NNUC(){return NNUC;}


  /// create a map that translates element name to position of element in ratioH vector 
  /// (holding the computed element abundances)
  /// BE convinience function just in case it changes with a different version of AlterBBN
  map_str_int get_abund_map_AlterBBN()
  {

    map_str_int abund_map;

    // maps elements to their position in 'ratioH' array in AlterBBN holding 
    // primordial element abundances relative to H abundance
    abund_map["H2"] = 3;
    abund_map["D"] = 3;
    abund_map["H3"] = 4;
    abund_map["He3"] = 5;
    abund_map["He4"] = 6;   
    abund_map["Yp"] = 6;
    abund_map["Li6"] = 7;
    abund_map["Li7"] = 8;
    abund_map["Be7"] = 9;
    abund_map["Li8"] = 10;

    return abund_map;

  }


}
END_BE_NAMESPACE
