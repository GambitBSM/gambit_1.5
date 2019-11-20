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
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2019 Feb, July
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2018 Oct
///  \date 2019 June
///
///  \author Sanjay Bloor
///          (sanjay.bloor12@imperial.ac.uk)
///   \date 2019 Nov
///
///  *********************************************

#ifndef __CosmoModels_hpp__
#define __CosmoModels_hpp__

// Vanilla ΛCDM.
#define MODEL LCDM
  START_MODEL
  DEFINEPARS(omega_b,omega_cdm,H0,ln10A_s,n_s,tau_reio)
#undef MODEL

// ΛCDM parameters without those relating to the primordial power spectrum (A_s, n_s)
// This model should be scanned alongside an inflationary model able to provide
// a primordial power spectrum. 
#define MODEL LCDM_no_primordial
  #define PARENT LCDM
  START_MODEL
  DEFINEPARS(omega_b,omega_cdm,H0,tau_reio)
  INTERPRET_AS_PARENT_FUNCTION(LCDM_to_LCDM_no_primordial)
  #undef PARENT
#undef MODEL

/* CMB + BBN */

// η (baryon-to-photon ratio) defined at BBN.
// r_CMB(BBN): Ratio of temperatures in non-cold DM compared to the SM with 3 massive neutrinos:
// r_CMB = T_v(BSM)/T_v(SM) at CMB (BBN)
// dNeff_CMB(_BBN): ΔN_effective defined at CMB(BBN) from additional radiation.
#define MODEL etaBBN_rBBN_rCMB_dNeffBBN_dNeffCMB
  START_MODEL
  DEFINEPARS(eta_BBN)
  MAP_TO_CAPABILITY(eta_BBN, etaBBN)
  DEFINEPARS(r_BBN,r_CMB)
  DEFINEPARS(dNeff_BBN,dNeff_CMB)
#undef MODEL

// No additional radiation or changes to the neutrino temperature.
// Just the baryon-to-photon ratio η at BBN.
#define MODEL etaBBN
 #define PARENT etaBBN_rBBN_rCMB_dNeffBBN_dNeffCMB
  START_MODEL
  DEFINEPARS(eta_BBN)
  INTERPRET_AS_PARENT_FUNCTION(etaBBN_to_etaBBN_rBBN_rCMB_dNeffBBN_dNeffCMB)
 #undef PARENT
#undef MODEL

// As etaBBN_rBBN_rCMB_dNeffBBN_dNeffCMB, but with the
// baryon-to-photon ratio η at BBN set equal to η today.
#define MODEL rBBN_rCMB_dNeffBBN_dNeffCMB
 #define PARENT etaBBN_rBBN_rCMB_dNeffBBN_dNeffCMB
  START_MODEL
  DEFINEPARS(r_BBN,r_CMB)
  DEFINEPARS(dNeff_BBN,dNeff_CMB)
  INTERPRET_AS_PARENT_FUNCTION(rBBN_rCMB_dNeffBBN_dNeffCMB_to_etaBBN_rBBN_rCMB_dNeffBBN_dNeffCMB)
  INTERPRET_AS_PARENT_DEPENDENCY(eta0, double)
 #undef PARENT
#undef MODEL

// As above, but with no additional radiation.
#define MODEL rBBN_rCMB
 #define PARENT rBBN_rCMB_dNeffBBN_dNeffCMB
  START_MODEL
  DEFINEPARS(r_BBN,r_CMB)
  INTERPRET_AS_PARENT_FUNCTION(rBBN_rCMB_to_rBBN_rCMB_dNeffBBN_dNeffCMB)
 #undef PARENT
#undef MODEL

// As above, but with neutrino temperature at BBN the same as at recombination.
#define MODEL rCMB
 #define PARENT rBBN_rCMB
  START_MODEL
  DEFINEPARS(r_CMB)
  INTERPRET_AS_PARENT_FUNCTION(rCMB_to_rBBN_rCMB)
 #undef PARENT
#undef MODEL

// As rBBN_rCMB_dNeffBBN_dNeffCMB, but with no changes to the neutrino temperature.
#define MODEL dNeffBBN_dNeffCMB
 #define PARENT rBBN_rCMB_dNeffBBN_dNeffCMB
  START_MODEL
  DEFINEPARS(dNeff_BBN,dNeff_CMB)
  INTERPRET_AS_PARENT_FUNCTION(dNeffBBN_dNeffCMB_to_rBBN_rCMB_dNeffBBN_dNeffCMB)
 #undef PARENT
#undef MODEL

// As above, but with additional radiation the same at BBN as at recombination.
#define MODEL dNeffCMB
 #define PARENT dNeffBBN_dNeffCMB
  START_MODEL
  DEFINEPARS(dNeff_CMB)
  INTERPRET_AS_PARENT_FUNCTION(dNeffCMB_to_dNeffBBN_dNeffCMB)
 #undef PARENT
#undef MODEL
  

/* INFLATION */

// inflationary models -- if one of them is in use you have to use the model LCDM_no_primordial 
// to scan over the four standard cosmological parameters (H0, omega_b, omega_cdm, tau_reio) and
// the shape of the primordial power spectrum will be determined by the inflation model in use

#define MODEL Inflation_tensor
START_MODEL
DEFINEPARS(ln10A_s,n_s,r_tensor)
#undef MODEL

/*
#define MODEL inflation // a minimally defined general inflationary model with 3 sets of parameters.
  START_MODEL
  DEFINEPARS(phi0,dphi0,vparams1,vparams2,vparams3)
#undef MODEL
*/

// Inflation model: 0.5 m^2 phi^2 --- quadratic inflation
// A_s, n_s and r are given by inflationary model
// parameters: N_piv, m^2.
#define MODEL Inflation_SR1quad
  START_MODEL
  DEFINEPARS(m2_inflaton,N_pivot)
#undef MODEL

// Inflation model: 0.25 \lambda phi^4  --- quartic inflation
// A_s, n_s and r are given by inflationary model
// parameters: N_piv, lambda.
#define MODEL Inflation_1quar
START_MODEL
DEFINEPARS(lambda,N_pivot)
#undef MODEL

// Inflation model: 2/3 \lambda phi^2/3 --- inflation
// A_s, n_s and r are given by inflationary model
// parameters: N_piv, lambda.
#define MODEL Inflation_1mono32Inf
START_MODEL
DEFINEPARS(lambda,N_pivot)
#undef MODEL

// Inflation model: m phi --- linear inflation
// A_s, n_s and r are given by inflationary model
// parameters: N_piv, m^2.
#define MODEL Inflation_1linearInf
START_MODEL
DEFINEPARS(lambda,N_pivot)
#undef MODEL

// simplest 8 parameter cosmology smash inflation model --- smash inflation
// A_s, n_s and r are given by inflationary model
// parameters: xi, m^2.
#define MODEL Inflation_smash
START_MODEL
DEFINEPARS(log10_xi,log10_beta,log10_lambda,N_pivot)
#undef MODEL

#define MODEL Inflation_1natural // N-flation (axions)
START_MODEL
DEFINEPARS(lambda,faxion,N_pivot)
#undef MODEL

#define MODEL Inflation_1hilltopInf // Hilltop
START_MODEL
DEFINEPARS(lambda,mu,N_pivot)
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

/* PLANCK NUISANCE PARAMETERS */

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
  //#define PARENT Planck_TTTEEE
    START_MODEL
    DEFINEPARS(A_cib_217,cib_index,xi_sz_cib,A_sz,ps_A_100_100,ps_A_143_143,ps_A_143_217,ps_A_217_217,ksz_norm,gal545_A_100,gal545_A_143,gal545_A_143_217,gal545_A_217,calib_100T,calib_217T,A_planck)
  //#undef PARENT
#undef MODEL

#define MODEL Planck_lite
  //#define PARENT Planck_TT
    START_MODEL
    DEFINEPARS(A_planck)
  //#undef PARENT
#undef MODEL

//#define MODEL inflation
//START_MODEL
//DEFINEPARS(num_inflaton, potential_choice, slowroll_infl_end, instreheat, vparam_rows, use_deltaN_SR, evaluate_modes, use_horiz_cross_approx, get_runningofrunning, ic_sampling, energy_scale, numb_samples, save_iso_N, N_iso_ref, param_sampling, vp_prior_min, vp_prior_max, varying_N_pivot, use_first_priorval, phi_init0, dphi_init0, vparams, N_pivot, k_pivot, dlnk, turning_choice  calc_full_pk,  steps,  kmin,  kmax,  phi0_priors_min,  phi0_priors_max,  dphi0_priors_min,  dphi0_priors_max,  N_pivot_prior_min,  N_pivot_prior_max)
//#undef MODEL

#endif
