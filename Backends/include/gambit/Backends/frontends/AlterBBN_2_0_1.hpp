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
///  \date   2018 Jun
///
///  *********************************************



#define BACKENDNAME AlterBBN
#define BACKENDLANG CC
#define VERSION 2.0.1
#define SAFE_VERSION 2_0_1

LOAD_LIBRARY

// BE_ALLOW_MODELS(MSSM63atQ, MSSM63atMGUT, WC)

BE_FUNCTION(Init_cosmomodel, void, (relicparam*), "Init_cosmomodel", "Init_cosmomodel")

BE_FUNCTION(Init_cosmomodel_param, void, (double, double, double, double, double, double, double, double, relicparam*), "Init_cosmomodel_param", "Init_cosmomodel_param")
BE_FUNCTION(Init_dark_density, void, (double, double, double, relicparam*), "Init_dark_density", "Init_dark_density")
BE_FUNCTION(Init_dark_density2, void, (double, double, double, relicparam*), "Init_dark_density2", "Init_dark_density2")
BE_FUNCTION(Init_dark_entropy, void, (double, double, double, relicparam*), "Init_dark_entropy", "Init_dark_entropy")
BE_FUNCTION(Init_dark_entropySigmaD, void, (double, double, double, relicparam*), "Init_dark_entropySigmaD", "Init_dark_entropySigmaD")
BE_FUNCTION(Init_entropySigmarad, void, (double, double, double, relicparam*), "Init_entropySigmarad", "Init_entropySigmarad")
BE_FUNCTION(Init_nonthermal, void, (double, double, double, relicparam*), "Init_nonthermal", "Init_nonthermal")
BE_FUNCTION(Init_gravitino, void, (double, relicparam*), "Init_gravitino", "Init_gravitino")
BE_FUNCTION(Init_scalarfield, void, (double, double, double, relicparam*), "Init_scalarfield", "Init_scalarfield")
BE_FUNCTION(Init_wimp, void, (double, int, int, int, int, int, double, relicparam*), "Init_wimp", "Init_wimp")
BE_FUNCTION(Init_modeleff, void, (int, relicparam*), "Init_modeleff", "Init_modeleff")


// compute nu density including effects on nu degeneracy (and derivative)
BE_FUNCTION(neutdens, void, (double, relicparam*), "neutdens", "neutdens")
BE_FUNCTION(neutdens_deriv, void, (double, double, double, relicparam*), "neutdens_deriv", "neutdens_deriv")

// mod QCD eg
BE_FUNCTION(heff, void, (double, relicparam*), "heff", "heff")
BE_FUNCTION(geff, void, (double, relicparam*), "geff", "geff")
BE_FUNCTION(sgStar, void, (double, relicparam*), "sgStar", "sgStar")



// Convenience functions:
//BE_CONV_FUNCTION(BKstarmumu_CONV, Flav_KstarMuMu_obs, (const parameters*, double, double), "BKstarmumu_CONV", (MSSM63atQ, MSSM63atMGUT, WC))

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
