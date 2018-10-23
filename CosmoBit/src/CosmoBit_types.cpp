//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Source code for types for module CosmoBit.
///  For instructions on adding new types, see
///  the corresponding header.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2017 Nov
///
///  \author Selim Hotinli
///	 \date 2018 Jan
///  *********************************************

#include <string>
#include <iostream>
#include <stdlib.h>     /* malloc, free, rand */

#include "gambit/CosmoBit/CosmoBit_types.hpp"

namespace Gambit
{
  namespace CosmoBit
  {

    BBN_container::BBN_container()
    { 
      int NNUC = 26; // amount of abundances, global param in AlterBBN
      BBN_abund.resize(NNUC+1, 0.);
      BBN_covmat.resize(NNUC+1, std::vector<double>(NNUC+1,0.));
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

      std::map<std::string,std::vector<double>> BBN_obs_dict;
    }
    
    std::map<std::string,int> BBN_container::get_map()
    {
      return abund_map;
    }   

    void BBN_container::fill_obs_dict(std::map<std::string,std::vector<double>> dict)
    {
      BBN_obs_dict = dict;
     }

    std::map<std::string,std::vector<double>> BBN_container::get_obs_dict()
    {
      return BBN_obs_dict;
    }  

    Class_container::Class_container() : lmax(2508)
    {
      input.clear();
      Cl_TT.resize(lmax+1, 0.);
      Cl_TE.resize(lmax+1, 0.);
      Cl_EE.resize(lmax+1, 0.);
      Cl_BB.resize(lmax+1, 0.);
      Cl_PhiPhi.resize(lmax+1, 0.);
      //std::cout << "Hello it's me. I am a Class_container" << std::endl;

      Pk_S.resize(lmax+1, 0.); // Primordial Scalar Power Spectrum
      Pk_T.resize(lmax+1, 0.); // Primordial Tensor Power Spectrum
      k_ar.resize(lmax+1, 0.); // Corresponding wavenumbers.
    }
/*
    Class_container::~Class_container()
    {
      //std::cout << "And I am out" << std::endl;
    }
*/
    void ClassInput::addEntry(std::string key, std::string val)
    {
      input_list[key] = val;
    }

    void ClassInput::addEntry(std::string key, double val)
    {
      std::ostringstream stringified_val;
      stringified_val << val;
      addEntry(key,stringified_val.str());
    }

    void ClassInput::addEntry(std::string key, int val)
    {
      std::ostringstream stringified_val;
      stringified_val << val;
      addEntry(key,stringified_val.str());
    }

    void ClassInput::clear()
    {
      input_list.clear();
    }

    std::map<std::string,std::string> ClassInput::get_map()
    {
      return input_list;
    }


  }
}
