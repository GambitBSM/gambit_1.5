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

    /// function to merge python dictionary extra_dict into input_dict. If both dictionaries have the same key
    /// it has to be decided on a case-by-case basis how to deal with that
    /// some examples are already implemented, might have to be extended in the future
    /// if more CLASS features are used. 
    /// -> very specific for merging the classy input dictionaries
    /// typical example for this would be the key 'output': 
    ///    dict input_dict: 'output' : 'tCl nCl' and 
    ///    dict extra_dict: 'output' : 'tCl mPk' => mPk has to be added
    ///    => results in    'output' : 'tCl nCl mPk'
    void merge_pybind_dicts(pybind11::dict& input_dict, pybind11::dict extra_dict, bool first_run) 
    {
      
      // loop through extra_dict (should typically be the shorter one)
      for (auto item : extra_dict)
      {  
        pybind11::str key = pybind11::str(item.first);
        pybind11::str arg = pybind11::str(item.second);
        
        // if item not contained in extra_dict but not in input_dict it will be added to input_dict
        if(!input_dict.attr("has_key")(key).cast<bool>())
        {
          input_dict[key] = arg;
          //std::cout << "Adding key = " << std::string(pybind11::str(item.first)) << ", "<< "value=" << std::string(pybind11::str(item.second)) << std::endl;
        }

        // if item contained in both: have to decide on a case-by-case basis what to do! 
        else
        { 
          // if 'output' is defined twice the entries have to be merged
          // e.g. input_dict['output'] = 'A B C'
          //      extra_dict['output'] = 'A B X Y'
          // should result in input_dict['output'] = 'A B C X Y'
          // (python string.find("x") returns -1 if "x" not contained)
          if(key.attr("find")(pybind11::str("output")).cast<int>()!=-1 ||
             key.attr("find")(pybind11::str("modes")).cast<int>()!=-1) 
          {
            // split string of extra_dict['output'] by spaces into list 
            // (-> in the example above this would give list = ['A', 'B', 'X', 'Y'])
            pybind11::list list = extra_dict[key].attr("split")();
            // iterate through list and check if current item is also contained in input_dict['output']
            // if it is not it will be added
            for(auto it : list)
            { 
              std::string list_entry = std::string(pybind11::str(it));

              // add entry if it is not already in input_dict['output']
              if(input_dict[key].attr("find")(pybind11::str(list_entry)).cast<int>()==-1)
              { 
                // add part of extra_dict[key] string that is not contained in input_dict[key] string to input_dict[key]
                std::string new_arg=std::string(pybind11::str(input_dict[key]))+std::string(" ")+std::string(list_entry);
                input_dict[key]= new_arg;
              }
            }
          }

          // if 'l_max...' (scalars or tensors) is defined twice use the maximum requested value
          // e.g. input_dict['l_max_scalars'] = '500'
          //      extra_dict['l_max_scalars'] = '2500'
          // should result in input_dict['l_max_scalars'] = '2500'
          else if(key.attr("find")(pybind11::str("l_max")).cast<int>()!=-1)
          {
            // cast pybind11::detail::item_accessor to pybind11::str, to c++ string and then to int
            //  (I know... the problem is that the yaml file entries are all parsed as strings so 
            //  we don't have to distinguish between the different types of different CLASS input keys. 
            //  The python wrapper parses everything as string so it is fine. Only here, where we actually
            //  have to compare the passed numbers this causes troubles and we have to cast from a python 
            //  string to a c++ int to be able to tell which entry is the larger one. If you have an idea
            //  to make this nicer --> go for it!!)
            int lmax_input_dict = std::stoi(std::string(pybind11::str(input_dict[key])));
            int lmax_extra_dict = std::stoi(std::string(pybind11::str(extra_dict[key])));

            // if lmax_extra_dict is higher than the entry in the input_dict replace it
            if (lmax_input_dict < lmax_extra_dict){input_dict[key] = lmax_extra_dict;}
            // if not 'input_dict'already contains the higher value and there is nothing to do here. 
          }
          else
          {
            if(first_run)
            {
              // (JR) TODO see what other combinations could go wrong here -> different likelihoods
              // asking for different redshift bins?
              // But need this check only on the first run (the inputs we are worried about are 
              // parameters specifying which output CLASS should produce so they won't change 
              // in-between the calculation of different points in one scan)
              std::cout <<"___________________________________________________________________"<<std::endl;
              std::cout <<"___________________________________________________________________"<<std::endl;
              std::cout << "Both dictionaries to merge contain key" << std::string(key) << "with entries "
                  << std::string(pybind11::str(input_dict[key]))<< " and " << std::string(pybind11::str(extra_dict[key])) << ". Don't know how to deal with that, yet. Will be taken care of soon." << std::endl;
              std::cout <<"___________________________________________________________________"<<std::endl;
              std::cout <<"___________________________________________________________________"<<std::endl;
            }
          }
        }
      }
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
