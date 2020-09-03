//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Type definition header for utilities for module CosmoBit.
///
///  Compile-time registration of type definitions
///  required for the rest of the code to
///  communicate with CosmoBit.
///
///  Add to this if you want to define a new type
///  for the functions in CosmoBit to return, but
///  you don't expect that type to be needed by
///  any other modules.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2019 Mar
///  \date 2019 June
///
///  *********************************************


#ifndef __CosmoBit_utils_hpp__
#define __CosmoBit_utils_hpp__

namespace Gambit
{

  namespace CosmoBit
  {

    namespace CosmoBit_utils
    {

      // set value of Neff that is assumed by default. 
      // Note: the reason why it's not fixed here, is that 
      // the SM value is also needed in the CosmoBit type 'SM_time_evo'
      // where we don't have access to the result of a capability
      double set_Neff_SM_value();

    }
  }
}

#endif // defined __CosmoBit_types_hpp__
