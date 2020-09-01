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

#include "gambit/CosmoBit/CosmoBit_utils.hpp"
#include "gambit/CosmoBit/CosmoBit_types.hpp"
#include "gambit/Utils/numerical_constants.hpp"

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
      double set_Neff_SM_value()
      {
        return 3.045;
      }

    }
  }
}
