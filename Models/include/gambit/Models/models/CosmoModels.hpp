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

#define MODEL LCDMtensor //question: +tensor is a parent? probably not qualified.
  START_MODEL
  DEFINEPARS(omega_b,omega_cdm,H0,ln10A_s,n_s,tau_reio,r_tensor)
#undef MODEL

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

//#define MODEL inflation
//START_MODEL
//DEFINEPARS(num_inflaton, potential_choice, slowroll_infl_end, instreheat, vparam_rows, use_deltaN_SR, evaluate_modes, use_horiz_cross_approx, get_runningofrunning, ic_sampling, energy_scale, numb_samples, save_iso_N, N_iso_ref, param_sampling, vp_prior_min, vp_prior_max, varying_N_pivot, use_first_priorval, phi_init0, dphi_init0, vparams, N_pivot, k_pivot, dlnk, turning_choice  calc_full_pk,  steps,  kmin,  kmax,  phi0_priors_min,  phi0_priors_max,  dphi0_priors_min,  dphi0_priors_max,  N_pivot_prior_min,  N_pivot_prior_max)
//#undef MODEL

#endif
