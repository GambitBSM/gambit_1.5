//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Source code for utilities needed in module CosmoBit.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2019 Mar
///  *********************************************

//#include <string>
//#include <iostream>
//#include <stdlib.h>     /* malloc, free, rand */
//#include <valarray>

#include "gambit/CosmoBit/CosmoBit_utils.hpp"
#include "gambit/CosmoBit/CosmoBit_types.hpp"
#include "gambit/Utils/numerical_constants.hpp"

namespace Gambit
{
  namespace CosmoBit
  {

    // Return entropy density of SM as function of Temperature T. By default: T interpreted to be in K, 
    // set T_in_eV = True if T is in eV
    double entropy_density_SM(double T, bool T_in_eV)
    {
        if(T_in_eV == true) {T = T/_kB_eV_over_K_;}

        return (2.*pow(pi,2)/45.) * (43./11.) * pow((_kB_eV_over_K_*T),3);
    }

    // convert neutrino masses into string compatible with CLASS input
    std::vector<double> m_ncdm_classInput(std::map<std::string,double> NuMasses_SM)
    {
      std::vector<double> numasses;

      if(NuMasses_SM["mNu1"]>0.)
        numasses.push_back(NuMasses_SM["mNu1"]);
      if(NuMasses_SM["mNu2"]>0.)
        numasses.push_back(NuMasses_SM["mNu2"]);
      if(NuMasses_SM["mNu3"]>0.)
        numasses.push_back(NuMasses_SM["mNu3"]);

      // Do a quick check if the size of numasses fits the expectation.
      // If this particular error is thrown I messed up big time.
      if (numasses.size() != NuMasses_SM["N_ncdm"])
      {
        std::ostringstream err;
        err << "Sonmething went wrong in \'m_ncdm_classInput\'. The size of the vector \'numasses\' (=" << numasses.size() <<") is not what is expected (=" << NuMasses_SM["N_ncdm"] <<").";
        CosmoBit_error().raise(LOCAL_INFO, err.str());
      }

      return numasses;
    }
  }
}
