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

#include <string>
#include <iostream>
#include <stdlib.h>     /* malloc, free, rand */
#include <valarray>

#include "gambit/CosmoBit/CosmoBit_utils.hpp"
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


    // Utility function to set the number of massive neutrino species and prepare the input for class
    std::vector<double> set_nu_masses(double mNu1, double mNu2, double mNu3, int& N_ncdm)
    {
      // !! masses mNu1, mNu2, mNu3 in eV unlike the definition in StandardMOdel_SLHA2 !!

      // Reset N_ncdm
      N_ncdm = 0;
      std::vector<double> neutrino_masses;

      // check for every mass if it is positive. If so add it to the vector and increase N_ncdm
      if (mNu1 > 0.)
      {
        N_ncdm++;
        neutrino_masses.push_back(mNu1);
      }
      if (mNu2 > 0.)
      {
        N_ncdm++;
        neutrino_masses.push_back(mNu2);
      }
      if (mNu3 > 0.)
      {
        N_ncdm++;
        neutrino_masses.push_back(mNu3);
      }

      return neutrino_masses;
    }


  }
}
