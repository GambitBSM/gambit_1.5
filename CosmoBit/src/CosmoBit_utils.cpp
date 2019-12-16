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
///  \date 20190 Mar, June
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

    // function to merge python dictionary b into a. If both dictionaries have the same key
    // the values of the keys will be concatenated (without duplicating an entry)
    // -> very specific for merging the classy input dictionaries
    // typical example for this would be 
    // dict a: 'output' : 'tCl nCl' and 
    // dict b: 'output' : 'tCl mPk' => mPk has
    // to be added to 'output' : 'tCl nCl mPk'
    void merge_pybind_dicts(pybind11::dict& a, pybind11::dict b) 
    {
      // loop through 2nd dict (better if this is the shorter one)
      static bool first_run = true;
      for (auto item : b)
      {  
        pybind11::str key = pybind11::str(item.first);
        pybind11::str arg = pybind11::str(item.second);
        
        // if item not contained in b but not a it will be added to a
        if(!a.attr("has_key")(key).cast<bool>())
        {
          a[key] = arg;
          //std::cout << "Adding key = " << std::string(pybind11::str(item.first)) << ", "<< "value=" << std::string(pybind11::str(item.second)) << std::endl;
        }
        // if item contained in both: split b by spaces and iterate through single entries 
        // to see if they are included in a
        else
        { 
          
          // python string.find("x") returns -1 if "x" not contained
          if(key.attr("find")(pybind11::str("output")).cast<int>()!=-1 ||
             key.attr("find")(pybind11::str("modes")).cast<int>()!=-1)
          {
            // split string by spaces into list
            pybind11::list list = b[key].attr("split")();
            for(auto it : list)
            { 
              std::string list_entry = std::string(pybind11::str(it));
              
              // python string.find("x") returns -1 if "x" not contained
              if(a[key].attr("find")(pybind11::str(list_entry)).cast<int>()==-1)
              { 
                // add part of b[key] string that is not contained in a[key] string to a[key]
                std::string new_arg=std::string(pybind11::str(a[key]))+std::string(" ")+std::string(list_entry);
                a[key]= new_arg;
              }
            }
          }
          else
          {
            if(first_run)
            {
              // (JR) TODO see what other combinations could go wrong here -> different likelihoods
              // asking for different redshift bins?
              // But need this check only on the first run (the inputs we are worries about are 
              // parameters specifying which output CLASS should produce so they won't change 
              // in-between the calculation of different points in one scan)
              std::cout <<"___________________________________________________________________"<<std::endl;
              std::cout <<"___________________________________________________________________"<<std::endl;
              std::cout << "Both dictionaries to merge contain key" << std::string(key) << "with entries "
                  << std::string(pybind11::str(a[key]))<< " and " << std::string(pybind11::str(b[key])) << ". Don't know how to deal with that, yet. Will be taken care of soon." << std::endl;
              std::cout <<"___________________________________________________________________"<<std::endl;
              std::cout <<"___________________________________________________________________"<<std::endl;
            }
          }
        }
      }
      first_run = false;
    }

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
