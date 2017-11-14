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
    Class_container::Class_container() : already_calculated(false), is_initialized(false)
    {
      std::cout << "Hello it's me. I am a Class_container" << std::endl;
    }

    Class_container::~Class_container()
    {
      std::cout << "And I am out" << std::endl;
    }
  }
}
