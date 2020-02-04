//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  General model 'cosmo_nuisance' and 
///  
///  
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2019 June
///  *********************************************

#ifndef __CosmoNuisanceModels_hpp__
#define __CosmoNuisanceModels_hpp__


// Supernovae -- JLA 
#define MODEL cosmo_nuisance_JLA
  START_MODEL
  DEFINEPARS(alpha,beta,M,Delta_M)
#undef MODEL

/// Pantheon -> child of JLA (light curve params fitted with SaltMu2 therefore only 1 nuisance param)
/// use same model for JLA_simple likelihood
#define MODEL cosmo_nuisance_Pantheon
  #define PARENT cosmo_nuisance_JLA
    START_MODEL
    DEFINEPARS(M)
    INTERPRET_AS_PARENT_FUNCTION(cosmo_nuisance_Pantheon_to_cosmo_nuisance_JLA)
  #undef PARENT
#undef MODEL

/// nuisance params for bicep/keck array likelihood implemented in MontePython
#define MODEL cosmo_nuisance_BK14
    START_MODEL
    DEFINEPARS(BBdust, BBsync,BBalphadust,BBbetadust,BBTdust,BBalphasync,BBbetasync,BBdustsynccorr,EEtoBB_dust,EEtoBBsync)
#undef MODEL

/// nuisance params for CFHTLenS tomographic weak lensing likelihood implemented in MontePython
#define MODEL cosmo_nuisance_CFHTLens_correlation
    START_MODEL
    // parameter name in MP likelihood: epsilon renamed to epsilon_CFHT, have to translate back before passing to MP
    // (we can not have several parameters with the same name in GAMBIT)
    DEFINEPARS(epsilon_CFHT)
#undef MODEL

/// nuisance params for euclid lensing likelihood implemented in MontePython
#define MODEL cosmo_nuisance_euclid_lensing
    START_MODEL
    // parameter name in MP likelihood: epsilon renamed to epsilon_euclid, have to translate back before passing to MP
    // (we can not have several parameters with the same name in GAMBIT)
    DEFINEPARS(epsilon_euclid)
#undef MODEL

/// nuisance params for euclid galaxy clustering likelihood implemented in MontePython
#define MODEL cosmo_nuisance_euclid_pk
    START_MODEL
    // actual parameter names: beta_0^Euclid,beta_1^Euclid but macros don't like ^ -> have to take care of renaming when 
    // filling data.mcmc_parameters dictionary for MontePython (TODO)
    // parameter name in MP likelihood: sigma_euclid renamed to sigma_NL_euclid, have to translate back before passing to MP
    // (we can not have several parameters with the same name in GAMBIT)
    DEFINEPARS(sigma_NL_euclid,beta_0Euclid,beta_1Euclid,P_shot)
#undef MODEL

/// nuisance params for tomographic ISW likelihood implemented in MontePython
#define MODEL cosmo_nuisance_ISW
    START_MODEL
    DEFINEPARS(A_ISW,b0_sdss,b1_sdss,b2_sdss,b3_sdss,b4_sdss,b0_qso,b1_qso,b2_qso,b0_mpz,b1_mpz,b2_mpz,b0_wisc,b1_wisc,b2_wisc,b0_nvss)
#undef MODEL

/// nuisance params for KiDS weak lensing likelihood implemented in MontePython
#define MODEL cosmo_nuisance_kids450_qe_likelihood_public
    START_MODEL
    DEFINEPARS(m_corr,A_IA,exp_IA,A_bary,A_noise_z1,A_noise_z2,A_noise_z3)
    // D_z1,D_z2,D_z3 params were not used in public analysis, don't think they are taken into account in the public 
    // likelihood calulation either
#undef MODEL

// contains nuisance params from all ska likelihoods implemented in MontePython-- make sure to check which ones are needed by the specific 
// you use (and we have another epsilon here..)
#define MODEL cosmo_nuisance_ska
    START_MODEL
    // actual parameter names: beta_0^IM,beta_1^IM,beta_0^SKA1,beta_1^SKA1,beta_0^SKA2,beta_1^SKA2 but macros don't like ^ -> have to take care of renaming when 
    // filling data.mcmc_parameters dictionary for MontePython (TODO)
    // parameter name in MP likelihood: epsilon,sigma_NL renamed to epsilon_ska & sigma_NL_ska, have to translate back before passing to MP
    // (we can not have several parameters with the same name in GAMBIT)
    DEFINEPARS(sigma_NL_ska,beta_0IM,beta_1IM,Omega_HI0,alpha_HI,beta_0SKA1,beta_1SKA1,beta_0SKA2,beta_1SKA2,epsilon_ska)
#undef MODEL

/*  // todo consistent tag for comments with new features for others
#define MODEL cosmo_nuisance_FOR_YOUR_NEW_LIKE 
    START_MODEL
    DEFINEPARS(your_param1,your_param2,your_param3) 
#undef MODEL
*/

#endif
