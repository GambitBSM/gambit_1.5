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
///  \author Sebastian Hoof
///          (hoof@uni-goettingen.de)
///   \date 2020 Mar
///
///  *********************************************

#ifndef __CosmoModels_hpp__
#define __CosmoModels_hpp__

// Vanilla ΛCDM.
#define MODEL LCDM
  START_MODEL
  DEFINEPARS(T_cmb,omega_b,omega_cdm,H0,ln10A_s,n_s,tau_reio)
  MAP_TO_CAPABILITY(T_cmb,T_cmb)
  MAP_TO_CAPABILITY(H0, H0)
#undef MODEL

// ΛCDM parameters without those relating to the primordial power spectrum (A_s, n_s)
// This model should be scanned alongside an inflationary model able to provide
// a primordial power spectrum.
#define MODEL LCDM_no_primordial
  #define PARENT LCDM
  START_MODEL
  DEFINEPARS(T_cmb,omega_b,omega_cdm,H0,tau_reio)
  INTERPRET_AS_PARENT_FUNCTION(LCDM_to_LCDM_no_primordial)
  #undef PARENT
#undef MODEL

// Vanilla ΛCDM.
#define MODEL LCDM_theta
  START_MODEL
  DEFINEPARS(T_cmb,omega_b,omega_cdm,100theta_s,ln10A_s,n_s,tau_reio)
  MAP_TO_CAPABILITY(T_cmb,T_cmb)
#undef MODEL

// ΛCDM parameters without those relating to the primordial power spectrum (A_s, n_s)
// This model should be scanned alongside an inflationary model able to provide
// a primordial power spectrum.
#define MODEL LCDM_theta_no_primordial
  #define PARENT LCDM_theta
  START_MODEL
  DEFINEPARS(T_cmb,omega_b,omega_cdm,100theta_s,tau_reio)
  INTERPRET_AS_PARENT_FUNCTION(LCDM_theta_to_LCDM_theta_no_primordial)
  #undef PARENT
#undef MODEL



/* CMB + BBN */

// η (baryon-to-photon ratio) defined at BBN.
// r_CMB(BBN): Ratio of temperatures in non-cold DM compared to the SM with 3 massive neutrinos:
// r_CMB = T_v(BSM)/T_v(SM) at CMB (BBN)
// dNurCMB(_BBN): ΔN_effective defined at CMB(BBN) from additional radiation.
#define MODEL etaBBN_rBBN_rCMB_dNurBBN_dNurCMB
  START_MODEL
  DEFINEPARS(eta_BBN)
  MAP_TO_CAPABILITY(eta_BBN, etaBBN)
  DEFINEPARS(r_BBN,r_CMB)
  DEFINEPARS(dNur_BBN,dNur_CMB)
#undef MODEL

// No additional radiation or changes to the neutrino temperature.
// Just the baryon-to-photon ratio η at BBN.
#define MODEL etaBBN
 #define PARENT etaBBN_rBBN_rCMB_dNurBBN_dNurCMB
  START_MODEL
  DEFINEPARS(eta_BBN)
  INTERPRET_AS_PARENT_FUNCTION(etaBBN_to_etaBBN_rBBN_rCMB_dNurBBN_dNurCMB)
 #undef PARENT
#undef MODEL

// As etaBBN_rBBN_rCMB_dNurBBN_dNurCMB, but with the
// baryon-to-photon ratio η at BBN set equal to η today.
#define MODEL rBBN_rCMB_dNurBBN_dNurCMB
 #define PARENT etaBBN_rBBN_rCMB_dNurBBN_dNurCMB
  START_MODEL
  DEFINEPARS(r_BBN,r_CMB)
  DEFINEPARS(dNur_BBN,dNur_CMB)
  INTERPRET_AS_PARENT_FUNCTION(rBBN_rCMB_dNurBBN_dNurCMB_to_etaBBN_rBBN_rCMB_dNurBBN_dNurCMB)
  INTERPRET_AS_PARENT_DEPENDENCY(eta0, double)
 #undef PARENT
#undef MODEL

// As above, but with no additional radiation.
#define MODEL rBBN_rCMB
 #define PARENT rBBN_rCMB_dNurBBN_dNurCMB
  START_MODEL
  DEFINEPARS(r_BBN,r_CMB)
  INTERPRET_AS_PARENT_FUNCTION(rBBN_rCMB_to_rBBN_rCMB_dNurBBN_dNurCMB)
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

// As rBBN_rCMB_dNurBBN_dNurCMB, but with no changes to the neutrino temperature.
#define MODEL dNurBBN_dNurCMB
 #define PARENT rBBN_rCMB_dNurBBN_dNurCMB
  START_MODEL
  DEFINEPARS(dNur_BBN,dNur_CMB)
  INTERPRET_AS_PARENT_FUNCTION(dNurBBN_dNurCMB_to_rBBN_rCMB_dNurBBN_dNurCMB)
 #undef PARENT
#undef MODEL

// As above, but with additional radiation the same at BBN as at recombination.
#define MODEL dNurCMB
 #define PARENT dNurBBN_dNurCMB
  START_MODEL
  DEFINEPARS(dNur_CMB)
  INTERPRET_AS_PARENT_FUNCTION(dNurCMB_to_dNurBBN_dNurCMB)
 #undef PARENT
#undef MODEL


/* INFLATION */

// Inflationary models -- if one of them is in use you have to use the model LCDM_theta_no_primordial
// to scan over the four standard cosmological parameters (H0 or theta, omega_b, omega_cdm, tau_reio).
// The shape of the (either parameterised or full) primordial power spectrum will be determined by
// the inflation model in use.

// Single field, monomic inflation with exponent 2/3 (assuming instant reheating)
// Potential: V(phi) = 1.5 lambda phi^(2/3)
#define MODEL Inflation_InstReh_1mono23
  START_MODEL
  DEFINEPARS(lambda)
#undef MODEL

// Single field, quadratic inflation (assuming instant reheating)
// Potential: V(phi) = lambda phi
#define MODEL Inflation_InstReh_1linear
  START_MODEL
  DEFINEPARS(lambda)
#undef MODEL

// Single field, quadratic inflation (assuming instant reheating)
// Potential: V(phi) = 0.5 m^2 phi^2
#define MODEL Inflation_InstReh_1quadratic
  START_MODEL
  DEFINEPARS(m_phi)
#undef MODEL

// Single field, quartic inflation (assuming instant reheating)
// Potential: V(phi) = 0.25 lambda phi^4
#define MODEL Inflation_InstReh_1quartic
  START_MODEL
  DEFINEPARS(lambda)
#undef MODEL

// Single field, natural inflation (assuming instant reheating)
// Potential: V(phi) = Lambda^4 [1 + cos(phi/f)]
#define MODEL Inflation_InstReh_1natural
  START_MODEL
  DEFINEPARS(Lambda, f_phi)
#undef MODEL

// Single field, Starobinsky - aka R^2 - inflation (assuming instant reheating)
// Potential: V(phi) = ...
#define MODEL Inflation_InstReh_1Starobinsky
  START_MODEL
  DEFINEPARS(Lambda)
#undef MODEL

#endif
