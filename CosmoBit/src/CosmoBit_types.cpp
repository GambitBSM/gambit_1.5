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
///  *********************************************

#include <string>
#include <iostream>
#include <stdlib.h>     /* malloc, free, rand */

#include "gambit/CosmoBit/CosmoBit_types.hpp"

namespace Gambit
{
  namespace CosmoBit
  {
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
