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

/// Number of observables the SuperIso returns for B0 -> K(*) mu mu
//#define Nobs_BKsll 30
//#define Nobs_BKll 2

#define LOWERR -0.2    // Minimum value of the lower relative error of the nuclear rates
#define UPPERR 1.      // Maximum value of the upper relative error of the nuclear rates

#define NNUCREAC 100  // Expanded from the original 88 reactions
#define NNUC 26

// Initialisation
BE_INI_FUNCTION{}
END_BE_INI_FUNCTION

// Convenience functions (definitions)
BE_NAMESPACE
{

}
END_BE_NAMESPACE
