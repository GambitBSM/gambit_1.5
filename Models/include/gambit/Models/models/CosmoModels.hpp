//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  LCDM model and data dependent nuisance parameter declarations.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Selim C. Hotinli
///          (selim.hotinli14@imperial.ac.uk)
///  \date 2016 Oct
///  \date 2017 Jan
///
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2018 Oct
///  *********************************************

#ifndef __CosmoModels_hpp__
#define __CosmoModels_hpp__

#define MODEL LCDM_dNeffCMB_dNeffBBN_etaBBN
  START_MODEL
  DEFINEPARS(omega_b,omega_cdm,H0,ln10A_s,n_s,tau_reio,dNeff,dNeff_BBN,eta_BBN)
#undef MODEL

#define MODEL LCDM_dNeffCMB_dNeffBBN
 #define PARENT LCDM_dNeffCMB_dNeffBBN_etaBBN
  START_MODEL
  DEFINEPARS(omega_b,omega_cdm,H0,ln10A_s,n_s,tau_reio,dNeff,dNeff_BBN)
  INTERPRET_AS_PARENT_FUNCTION(LCDM_dNeffCMB_dNeffBBN_to_LCDM_dNeffCMB_dNeffBBN_etaBBN)
  INTERPRET_AS_PARENT_DEPENDENCY(etaBBN, double)
 #undef PARENT
#undef MODEL

#define MODEL LCDM_dNeffCMB
 #define PARENT LCDM_dNeffCMB_dNeffBBN
  START_MODEL
  DEFINEPARS(omega_b,omega_cdm,H0,ln10A_s,n_s,tau_reio,dNeff)
  INTERPRET_AS_PARENT_FUNCTION(LCDM_dNeffCMB_to_LCDM_dNeffCMB_dNeffBBN)
 #undef PARENT
#undef MODEL

// get dNeff from external calculation (e.g. for ALPs)
#define MODEL LCDM_ExtdNeffCMB_ExtetaBBN
 #define PARENT LCDM_dNeffCMB_dNeffBBN_etaBBN
  START_MODEL
  DEFINEPARS(omega_b,omega_cdm,H0,ln10A_s,n_s,tau_reio)
  INTERPRET_AS_PARENT_FUNCTION(LCDM_ExtdNeffCMB_ExtetaBBN_to_LCDM_dNeffCMB_dNeffBBN_etaBBN)
  INTERPRET_AS_PARENT_DEPENDENCY(etaBBN, double)
  INTERPRET_AS_PARENT_DEPENDENCY(ExtdNeffCMB, double)
 #undef PARENT
#undef MODEL

#define MODEL LCDM
 #define PARENT LCDM_dNeffCMB
  START_MODEL
  DEFINEPARS(omega_b,omega_cdm,H0,ln10A_s,n_s,tau_reio)
  INTERPRET_AS_PARENT_FUNCTION(LCDM_to_LCDM_dNeffCMB)
 #undef PARENT
#undef MODEL


/*
#define MODEL rLCDM
START_MODEL
DEFINEPARS(omega_b,omega_cdm,H0,tau_reio)
#undef MODEL

#define MODEL rLCDMtensor
  START_MODEL
  DEFINEPARS(omega_b,omega_cdm,H0,tau_reio,r_tensor)
#undef MODEL
*/

#define MODEL LCDMtensor
START_MODEL
DEFINEPARS(omega_b,omega_cdm,H0,ln10A_s,n_s,tau_reio,r_tensor)
#undef MODEL

/*
#define MODEL inflation // a minimally defined general inflationary model with 3 sets of parameters.
  START_MODEL
  DEFINEPARS(phi0,dphi0,vparams1,vparams2,vparams3)
#undef MODEL
*/
// simplest 6 parameter cosmology+inflation model: 0.5 m^2 phi^2 --- quadratic inflation
// A_s, n_s and r are given by inflationary model
// parameters: N_piv, m^2.
#define MODEL inf_SR1quad_LCDMt
  START_MODEL
  DEFINEPARS(m2_inflaton,N_pivot,omega_b,omega_cdm,H0,tau_reio)
#undef MODEL

//  6 parameter cosmology+inflation model: 0.25 \lambda phi^4  --- quartic inflation
// A_s, n_s and r are given by inflationary model
// parameters: N_piv, lambda.
#define MODEL inf_1quarInf_LCDMt
START_MODEL
DEFINEPARS(lambda,N_pivot,omega_b,omega_cdm,H0,tau_reio)
#undef MODEL

//  6 parameter cosmology+inflation model: 2/3 \lambda phi^2/3 --- inflation
// A_s, n_s and r are given by inflationary model
// parameters: N_piv, lambda.
#define MODEL inf_1mono32Inf_LCDMt
START_MODEL
DEFINEPARS(lambda,N_pivot,omega_b,omega_cdm,H0,tau_reio)
#undef MODEL

//  6 parameter cosmology+inflation model: m phi --- linear inflation
// A_s, n_s and r are given by inflationary model
// parameters: N_piv, m^2.
#define MODEL inf_1linearInf_LCDMt
START_MODEL
DEFINEPARS(lambda,N_pivot,omega_b,omega_cdm,H0,tau_reio)
#undef MODEL

// simplest 8 parameter cosmology smash inflation model --- smash inflation
// A_s, n_s and r are given by inflationary model
// parameters: xi, m^2.
#define MODEL inf_smashInf_LCDMt
START_MODEL
DEFINEPARS(log10_xi,log10_beta,log10_lambda,N_pivot,omega_b,omega_cdm,H0,tau_reio)
#undef MODEL

#define MODEL inf_1naturalInf_LCDMt // N-flation (axions)
START_MODEL
DEFINEPARS(lambda,faxion,N_pivot,omega_b,omega_cdm,H0,tau_reio)
#undef MODEL

#define MODEL inf_1hilltopInf_LCDMt // Hilltop
START_MODEL
DEFINEPARS(lambda,mu,N_pivot,omega_b,omega_cdm,H0,tau_reio)
#undef MODEL


/*

#define MODEL inf_diff1 // Lambda^4 - mu*phi^4
  START_MODEL
  DEFINEPARS(phi0,dphi0,vparams1,vparams2,vparams3)
#undef MODEL

#define MODEL inf_exp // Product of exponentials
  START_MODEL
  DEFINEPARS(phi0,dphi0,vparams1,vparams2,vparams3)
#undef MODEL

#define MODEL inf_hybrid // Canonical two-field hybrid
  START_MODEL
  DEFINEPARS(phi0,dphi0,vparams1,vparams2,vparams3)
#undef MODEL

#define MODEL inf_offset // V0 + m_i^2 phi_i^2
  START_MODEL
  DEFINEPARS(phi0,dphi0,vparams1,vparams2,vparams3)
#undef MODEL

// N-quadratic w/one quartic interaction
// term phi_i^2 + phi_{lightest}^2*phi_i^2
#define MODEL inf_intrx
  START_MODEL
  DEFINEPARS(phi0,dphi0,vparams1,vparams2,vparams3)
#undef MODEL

// Mass matrix with diagonal terms = m_i^2
// Off-diagonal terms = \eps
#define MODEL inf_offdiag
  START_MODEL
  DEFINEPARS(phi0,dphi0,vparams1,vparams2,vparams3)
#undef MODEL

#define MODEL inf_step // Multifield step potential
  START_MODEL
  DEFINEPARS(phi0,dphi0,vparams1,vparams2,vparams3)
#undef MODEL

#define MODEL inf_monomial // (1/p) lambda_i |phi_i|^p --- N-monomial
  START_MODEL
  DEFINEPARS(phi0,dphi0,vparams1,vparams2,vparams3)
#undef MODEL

#define MODEL inf_gaxion // Generalized axions
  START_MODEL
  DEFINEPARS(phi0,dphi0,vparams1,vparams2,vparams3)
#undef MODEL

#define MODEL inf_smash // SMASH potential
  START_MODEL
  DEFINEPARS(phi0,dphi0,vparams1,vparams2,vparams3)
#undef MODEL
 */

#define MODEL plik_dx11dr2_HM_v18_TT
  START_MODEL
  DEFINEPARS(A_cib_217,cib_index,xi_sz_cib,A_sz,ps_A_100_100,ps_A_143_143,ps_A_143_217,ps_A_217_217,ksz_norm,gal545_A_100,gal545_A_143,gal545_A_143_217,gal545_A_217,calib_100T,calib_217T,A_planck)
#undef MODEL

#define MODEL lowl_SMW_70_dx11d_2014_10_03_v5c_Ap
  #define PARENT plik_dx11dr2_HM_v18_TT
    START_MODEL
    DEFINEPARS(A_planck)
  #undef PARENT
#undef MODEL

#define MODEL Planck_TTTEEE
  START_MODEL
    DEFINEPARS(A_cib_217,cib_index,xi_sz_cib,A_sz,ps_A_100_100,ps_A_143_143,ps_A_143_217,ps_A_217_217,ksz_norm,gal545_A_100,gal545_A_143,gal545_A_143_217,gal545_A_217,galf_EE_A_100,galf_EE_A_100_143,galf_EE_A_100_217,galf_EE_A_143,galf_EE_A_143_217,galf_EE_A_217,galf_EE_index,galf_TE_A_100,galf_TE_A_100_143,galf_TE_A_100_217,galf_TE_A_143,galf_TE_A_143_217,galf_TE_A_217,galf_TE_index,calib_100T,calib_217T,calib_100P,calib_143P,calib_217P,A_pol,A_planck)
#undef MODEL

#define MODEL Planck_TT
  #define PARENT Planck_TTTEEE
    START_MODEL
    DEFINEPARS(A_cib_217,cib_index,xi_sz_cib,A_sz,ps_A_100_100,ps_A_143_143,ps_A_143_217,ps_A_217_217,ksz_norm,gal545_A_100,gal545_A_143,gal545_A_143_217,gal545_A_217,calib_100T,calib_217T,A_planck)
  #undef PARENT
#undef MODEL

#define MODEL Planck_lite
  #define PARENT Planck_TT
    START_MODEL
    DEFINEPARS(A_planck)
  #undef PARENT
#undef MODEL

// model used for SNe likelihood implemented in CosmoBit (not for use with MontePython)
#define MODEL cosmo_nuisance_params 
  START_MODEL
  DEFINEPARS(M_AbsMag_SNe)
#undef MODEL

// Following: Nuisance parameter models for each MPlike likelihood 

// Supernovae -- JLA 
#define MODEL cosmo_nuisance_params_JLA
  START_MODEL
  DEFINEPARS(alpha,beta,M,Delta_M)
#undef MODEL

// Pantheon -> child of JLA (light curve params fitted with SaltMu2 therefore only 1 nuisance param)
// use same model for JLA_simple likelihood
#define MODEL cosmo_nuisance_params_Pantheon
  #define PARENT cosmo_nuisance_params_JLA
    START_MODEL
    DEFINEPARS(M)
    INTERPRET_AS_PARENT_FUNCTION(cosmo_nuisance_params_Pantheon_to_cosmo_nuisance_params_JLA)
  #undef PARENT
#undef MODEL

#define MODEL cosmo_nuisance_params_BK14
  START_MODEL
  DEFINEPARS(BBdust, BBsync,BBalphadust,BBbetadust,BBTdust,BBalphasync,BBbetasync,BBdustsynccorr,EEtoBB_dust,EEtoBBsync)
#undef MODEL

#define MODEL cosmo_nuisance_params_CFHTLens_correlation
  START_MODEL
  DEFINEPARS(epsilon)
#undef MODEL

/*  /\
    || (JR), what happens when two models have same param name?? Can't really chance this though..
    \/   otherwise we will have problems passing it to MP such that is understands which parameter to use with likelihood */

#define MODEL cosmo_nuisance_params_euclid_lensing
  START_MODEL
  DEFINEPARS(epsilon) 
#undef MODEL

#define MODEL cosmo_nuisance_params_euclid_pk
  START_MODEL
  // actual parameter names: beta_0^Euclid,beta_1^Euclid but macros don't like ^ -> have to take care of renaming when 
  // filling data.mcmc_parameters dictionary for MontePython (TODO)
  DEFINEPARS(sigma_NL,beta_0Euclid,beta_1Euclid,P_shot) 
#undef MODEL

#define MODEL cosmo_nuisance_params_ISW
  START_MODEL
  DEFINEPARS(A_ISW,b0_sdss,b1_sdss,b2_sdss,b3_sdss,b4_sdss,b0_qso,b1_qso,b2_qso,b0_mpz,b1_mpz,b2_mpz,b0_wisc,b1_wisc,b2_wisc,b0_nvss)
#undef MODEL

#define MODEL cosmo_nuisance_params_kids450_qe_likelihood_public
  START_MODEL
  DEFINEPARS(m_corr,A_IA,exp_IA,A_bary,A_noise_z1,A_noise_z2,A_noise_z3)
  // D_z1,D_z2,D_z3 params were not used in public analysis, don't think they are taken into account in the public 
  // likelihood calulation either
#undef MODEL

// contains nuisance params from all ska likelihoods -- make sure to check which ones are needed by the specific 
// you use (and we have another epsilon here..)
#define MODEL cosmo_nuisance_params_ska
  START_MODEL
  // actual parameter names: beta_0^IM,beta_1^IM,beta_0^SKA1,beta_1^SKA1,beta_0^SKA2,beta_1^SKA2 but macros don't like ^ -> have to take care of renaming when 
  // filling data.mcmc_parameters dictionary for MontePython (TODO)
  DEFINEPARS(sigma_NL,beta_0IM,beta_1IM,Omega_HI0,alpha_HI,beta_0SKA1,beta_1SKA1,beta_0SKA2,beta_1SKA2,epsilon)
#undef MODEL

#define MODEL cosmo_nuisance_params_wmap
  START_MODEL
  DEFINEPARS(A_SZ) // heads-up: this is not actually used in the wmap likelihood included in MontePyhton (3.1.0 at least)
#undef MODEL

/*
#define MODEL cosmo_nuisance_params_YOUR_NEW_LIKELIHOOD_NUISANCE_PRAMS_HERE
  START_MODEL
  DEFINEPARS(your_param1,your_param2,your_param3) 
#undef MODEL
*/


//#define MODEL inflation
//START_MODEL
//DEFINEPARS(num_inflaton, potential_choice, slowroll_infl_end, instreheat, vparam_rows, use_deltaN_SR, evaluate_modes, use_horiz_cross_approx, get_runningofrunning, ic_sampling, energy_scale, numb_samples, save_iso_N, N_iso_ref, param_sampling, vp_prior_min, vp_prior_max, varying_N_pivot, use_first_priorval, phi_init0, dphi_init0, vparams, N_pivot, k_pivot, dlnk, turning_choice  calc_full_pk,  steps,  kmin,  kmax,  phi0_priors_min,  phi0_priors_max,  dphi0_priors_min,  dphi0_priors_max,  N_pivot_prior_min,  N_pivot_prior_max)
//#undef MODEL

#endif
