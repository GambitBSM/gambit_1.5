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
#include "gambit/Backends/frontends/AlterBBN_2_0.hpp"
#include "gambit/Backends/backend_types/AlterBBN.hpp"

#define NNUC 26 // number of element abundances computed in AlterBBN 2.0

// Initialisation
BE_INI_FUNCTION{


}
END_BE_INI_FUNCTION


// Convenience functions (definitions)
BE_NAMESPACE
{
  int get_NNUC()
  {
  	return NNUC;
  }

}
END_BE_NAMESPACE
