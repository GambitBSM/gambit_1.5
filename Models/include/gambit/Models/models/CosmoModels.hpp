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
///  *********************************************


#ifndef __CosmoModels_hpp__
#define __CosmoModels_hpp__


#define MODEL LCDM
  START_MODEL
  DEFINEPARS(omega_b,omega_cdm,H0,ln10A_s,n_s,tau_reio)
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

/*
#define MODEL inf_SR1quad_LCDMt // m_i^2 phi_i^2 --- N-quadratic
 START_MODEL
 DEFINEPARS(phi0,dphi0,vparams1,vparams2,vparams3)
#undef MODEL

 #define MODEL inf_axions // N-flation (axions)
  START_MODEL
  DEFINEPARS(phi0,dphi0,vparams1,vparams2,vparams3)
#undef MODEL

#define MODEL inf_nquar // lambda_i phi_i^4 --- N-quartic
  START_MODEL
  DEFINEPARS(phi0,dphi0,vparams1,vparams2,vparams3)
#undef MODEL

#define MODEL inf_nlinear // lambda_i phi_i --- N-linear
  START_MODEL
  DEFINEPARS(phi0,dphi0,vparams1,vparams2,vparams3)
#undef MODEL

#define MODEL inf_power23 // lambda_i phi_i^(2/3)
  START_MODEL
  DEFINEPARS(phi0,dphi0,vparams1,vparams2,vparams3)
#undef MODEL

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


#endif
