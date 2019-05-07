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
///
///  *********************************************

#include <sstream>
#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/AlterBBN_2_1.hpp"
//#include "gambit/Backends/backend_types/AlterBBN_2_1/AlterBBN_2_1.hpp"

#define NNUC 26 // number of element abundances computed in AlterBBN 2.1

// Initialisation
BE_INI_FUNCTION{


}
END_BE_INI_FUNCTION


// Convenience functions (definitions)
BE_NAMESPACE
{
  std::set<std::string> known_relicparam_options = {"eta0", "Nnu", "dNnu", "life_neutron", "life_neutron_error",
                    "entropy_model", "energy_model", "relicmass","err","failsafe"};
        
  void fill_cosmomodel(AlterBBN::AlterBBN_2_1::relicparam * input_relicparam, map_str_dbl & AlterBBN_input)
  {
    static bool first_run = true;
    if(first_run)
    {
      std::map<std::string, double>::const_iterator it;
      for ( it = AlterBBN_input.begin(); it != AlterBBN_input.end(); it++ )
      {
        if(known_relicparam_options.count(it->first)==0)
        {
          std::ostringstream errormsg;
          errormsg << "Option '" << it->first << "' tried to pass to AlterBBN_" << VERSION<< " is not know to the relicparam strucutre "<<std::endl;
          errormsg << "Known options are:";
          for(auto item : known_relicparam_options)
          {
                  errormsg <<" '"<<item << "'";
          }
          errormsg << ". Check for typos or add new one."<<std::endl;
          backend_error().raise(LOCAL_INFO,errormsg.str());
        }
      }
      first_run = false;
    }

    if (AlterBBN_input.count("eta0")){input_relicparam->eta0 = AlterBBN_input["eta0"];}
  	if (AlterBBN_input.count("Nnu")){input_relicparam->Nnu = AlterBBN_input["Nnu"];}
  	if (AlterBBN_input.count("dNnu")){input_relicparam->dNnu = AlterBBN_input["dNnu"];}
  	
    if (AlterBBN_input.count("failsafe")){input_relicparam->dNnu = (int)AlterBBN_input["failsafe"];}
  	if (AlterBBN_input.count("err")){input_relicparam->err = (int)AlterBBN_input["err"];}
  }

  int call_nucl_err(map_str_dbl &AlterBBN_input, double* ratioH, double* cov_ratioH )
  {
    struct AlterBBN::AlterBBN_2_1::relicparam input_relicparam;
    Init_cosmomodel(&input_relicparam); // initialise valuse of relicparam structure to their defaults 
    fill_cosmomodel(&input_relicparam, AlterBBN_input); // fill strucutre with values contained in AlerBBN_input map which is filled in CosmoBit.ccp for different models

    int nucl_err_res = nucl_err(&input_relicparam, ratioH, cov_ratioH);

    return nucl_err_res;
  }

  int get_NNUC(){return NNUC;}

}
END_BE_NAMESPACE
