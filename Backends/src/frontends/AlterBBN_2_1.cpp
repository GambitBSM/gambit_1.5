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
///  \author Pat Scott
///          (pat.scott@uq.edu.au)
///  \date 2020 Apr
///
///  *********************************************

#include <sstream>
#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/AlterBBN_2_1.hpp"

#define NNUC 26 // number of element abundances computed in AlterBBN 2.1

// Initialisation
BE_INI_FUNCTION{


}
END_BE_INI_FUNCTION


// Convenience functions (definitions)
BE_NAMESPACE
{

  /// string set containing the name of all members of the AlterBBN relicparam structures that can currently
  /// be set with GAMBIT. If you add a new model and need to pass a different option to AlterBBN add e.g. "life_neutron"
  /// here and in the function fill_cosmomodel below to modify the lifetime of the neutron
  std::set<std::string> known_relicparam_options = {"eta0", "Nnu", "dNnu", "err","failsafe"};

  /// pass all values of AlterBBN's relicparam structure that
  void fill_cosmomodel(AlterBBN_2_1::relicparam * input_relicparam, map_str_dbl & AlterBBN_input)
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
    AlterBBN_2_1::relicparam input_relicparam;
    Init_cosmomodel(&input_relicparam); // initialise valuse of relicparam structure to their defaults
    fill_cosmomodel(&input_relicparam, AlterBBN_input); // fill strucutre with values contained in AlerBBN_input map which is filled in CosmoBit.ccp for different models

    int nucl_err_res = nucl_err(&input_relicparam, ratioH, cov_ratioH);
    //int bbn_exc = bbn_excluded_chi2(&input_relicparam); just for testing purposes to compare internal AlterBBN output to GAMBIT result

    return nucl_err_res;
  }

  /// return the NNUC -- global parameter of AlterBBN specifying the number of
  /// elements for which abundances are calculated -> length of array ratioH is NNUC+1 (AlterBBN
  /// starts filling @ position 1)
  int get_NNUC() { return NNUC; }

  /// create a map that translates element name to position of element in ratioH vector
  /// (holding the computed element abundances)
  /// BE convinience function just in case it changes with a different version of AlterBBN
  map_str_int get_abund_map_AlterBBN()
  { return {{"H2",3},{"D",3},{"H3",4},{"He3",5},{"He4",6},{"Yp",6},{"Li6",7},{"Li7",8},{"Be7",9},{"Li8",10}}; }

}
END_BE_NAMESPACE
