//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for AlterBBN v 2.0beta1
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Janina Renk
///				(janina.renk@fysik.su.se)
///  \date   2018 Jun
///
///  *********************************************


#define BACKENDNAME AlterBBN
#define BACKENDLANG CC
#define VERSION 2.0
#define SAFE_VERSION 2_0

LOAD_LIBRARY

BE_ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)
// TODO: add SMASH
BE_FUNCTION(Init_cosmomodel, void, (relicparam*), "Init_cosmomodel", "Init_cosmomodel")
BE_FUNCTION(nucl_err, int, (const relicparam*, double* , double* ), "nucl_err", "nucl_err")

//BE_FUNCTION(heff, double, (double, relicparam*), "heff", "heff")
//BE_FUNCTION(geff, double, (double, relicparam*), "geff", "geff")
//BE_FUNCTION(sgStar, double, (double, relicparam*), "sgStar", "sgStar")

BE_CONV_FUNCTION(get_NNUC, int, (), "get_NNUC")

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
