//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend for AlterBBN v 2.0
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Janina Renk
///             (janina.renk@fysik.su.se)
///  \date   2018 Jun
///
///  *********************************************


#define BACKENDNAME AlterBBN
#define BACKENDLANG CC
#define VERSION 2.0
#define SAFE_VERSION 2_0

//#include "gambit/Backends/backend_types/AlterBBN_2_0/identification_AlterBBN_2_0.hpp"

LOAD_LIBRARY


BE_ALLOW_MODELS(LCDM_dNeffCMB_dNeffBBN_etaBBN)

BE_FUNCTION(Init_cosmomodel, void, (AlterBBN::AlterBBN_2_0::relicparam*), "Init_cosmomodel", "Init_cosmomodel")
BE_FUNCTION(nucl_err, int, (AlterBBN::AlterBBN_2_0::relicparam*, double* , double* ), "nucl_err", "nucl_err")
//BE_FUNCTION(bbn_excluded_chi2, int, (AlterBBN::AlterBBN_2_0::relicparam*), "bbn_excluded_chi2", "bbn_excluded_chi2")

//BE_FUNCTION(heff, double, (double, relicparam*), "heff", "heff")
//BE_FUNCTION(geff, double, (double, relicparam*), "geff", "geff")
//BE_FUNCTION(sgStar, double, (double, relicparam*), "sgStar", "sgStar")

BE_CONV_FUNCTION(get_NNUC, int, (), "get_NNUC")
BE_CONV_FUNCTION(fill_cosmomodel, void, (AlterBBN::AlterBBN_2_0::relicparam*, map_str_dbl &), "Init_AlterBBN")
BE_CONV_FUNCTION(call_nucl_err, int, (map_str_dbl &, double* , double* ), "call_nucl_err")

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
