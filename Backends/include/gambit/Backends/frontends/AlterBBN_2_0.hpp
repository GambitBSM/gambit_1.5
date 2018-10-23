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

#define NNUC 26

LOAD_LIBRARY

BE_ALLOW_MODELS(LCDM, LCDM_dNeff_Smu, LCDM_dNeff_Smu_etaBBN)

BE_FUNCTION(Init_cosmomodel, void, (relicparam*), "Init_cosmomodel", "Init_cosmomodel")
BE_FUNCTION(nucl_err, int, (const relicparam*, double* , double* ), "nucl_err", "nucl_err")

//BE_FUNCTION(bbn_excluded_chi2, int, (const relicparam*), "bbn_excluded_chi2", "bbn_excluded_chi2")

//BE_FUNCTION(invert_matrix, int, (int , double*, double* ), "invert_matrix", "bbn_excluded_chi2")
// BE_FUNCTION(Init_modeleff, void, (int, relicparam*), "Init_modeleff", "Init_modeleff")
//BE_FUNCTION(Init_wimp, void, (double, int, int, int, int, int, double, relicparam*), "Init_wimp", "Init_wimp")
//BE_FUNCTION(Init_nonthermal, void, (double, double, double, relicparam*), "Init_nonthermal", "Init_nonthermal")
//BE_FUNCTION(Init_gravitino, void, (double, relicparam*), "Init_gravitino", "Init_gravitino")

// // mod QCD eg
// BE_FUNCTION(heff, double, (double, relicparam*), "heff", "heff")
// BE_FUNCTION(geff, double, (double, relicparam*), "geff", "geff")
// BE_FUNCTION(sgStar, double, (double, relicparam*), "sgStar", "sgStar")

// // calculate chi2 for Yp and 2H/H --> write own routine in gambit to be able to adopt measured values and which abundances are used 
BE_FUNCTION(bbn_excluded_chi2, int, (const relicparam*), "bbn_excluded_chi2", "bbn_excluded_chi2")



// Convenience functions:
//BE_CONV_FUNCTION(BKstarmumu_CONV, Flav_KstarMuMu_obs, (const parameters*, double, double), "BKstarmumu_CONV", (MSSM63atQ, MSSM63atMGUT, WC))

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
