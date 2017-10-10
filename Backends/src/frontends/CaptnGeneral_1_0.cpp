//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for Capt'n General 1.0 backend
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  Aaron Vincent
///  25/09/2017
///  *********************************************

#include "gambit/Backends/frontend_macros.hpp"
#include "gambit/Backends/frontends/CaptnGeneral_1_0.hpp"
//#include "gambit/Utils/file_lock.hpp"
//#include "gambit/Utils/mpiwrapper.hpp"

//#define square(x) ((x) * (x))  // square a number

// Initialisation function (definition)
BE_INI_FUNCTION
{
    static bool scan_level = true;
    if (scan_level)
    {
    cout << "Initializing CapGen, hold on to your hamsters. \n";
    const int clen = 300;
    char solarmodel[clen];
    Utils::strcpy2f(solarmodel, clen, runOptions->getValueOrDef<str>(backendDir +
                                                                    "/solarmodels/model_gs98_nohead.dat", "solarmodel"));

   cout << solarmodel << "\n";
  captn_init(solarmodel[0]);
    cout << "Hamsters held. \n";
    }
    scan_level = false;
}
END_BE_INI_FUNCTION

// Convenience functions (definitions)

