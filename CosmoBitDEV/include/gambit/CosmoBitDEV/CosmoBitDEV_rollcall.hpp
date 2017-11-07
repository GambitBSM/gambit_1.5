//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Rollcall header for the developement
///  version of the CosmoBit module.
///
///  Compile-time registration of available
///  functions and variables from this backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Selim C. Hotinli
///          (selim.hotinli14@imperial.ac.uk)
///  \date 2017 Jul
///
///  *********************************************

#ifndef __CosmoBitDEV_rollcall_hpp__
#define __CosmoBitDEV_rollcall_hpp__

#define MODULE CosmoBitDEV
START_MODULE

  #define CAPABILITY compute_vanilla_lowp_TT_loglike
  START_CAPABILITY
    #define FUNCTION function_vanilla_lowp_TT_loglike
	START_FUNCTION(double)
	  ALLOW_MODELS(LCDM,plik_dx11dr2_HM_v18_TT)
	  BACKEND_REQ(clik_compute_loglike, (clik_tag), double ,(clik_object*,double*,clik_error**))
	  BACKEND_REQ(clik_initialize, (clik_tag), clik_object* ,(char*,clik_error**))
	  BACKEND_REQ(clik_initialize_error, (clik_tag), clik_error* ,())
	  BACKEND_REQ(clik_get_list_cls, (clik_tag), void , (clik_object*,int,clik_error**))
	  BACKEND_REQ(clik_get_extra_parameter_names, (clik_tag), int ,(clik_object*,parname**,clik_error**))
	  BACKEND_REQ(clik_get_lmax, (clik_tag), void ,(clik_object*,int,clik_error**))
	  BACKEND_REQ(clik_cleanup, (clik_tag), void ,(clik_object**))
	  BACKEND_REQ(clik_get_version, (clik_tag), char* , (clik_object*,clik_error**))
	  BACKEND_REQ(data_initialize, (clik_tag), void, ())
	  BACKEND_REQ(return_high_TT,(clik_tag),clik_object*,())
	  BACKEND_REQ(return_lowp_TT,(clik_tag),clik_object*,())
	  BACKEND_REQ(class_input_initialize,(class_tag),int,(struct file_content*, struct precision*, struct background*, struct thermo*, struct perturbs*, struct transfers*, struct primordial*, struct spectra*, struct nonlinear*, struct lensing*, struct output*, char*))
	  BACKEND_REQ(class_background_initialize,(class_tag),int,(struct precision*, struct background*))
      BACKEND_REQ(class_thermodynamics_initialize,(class_tag),int,(struct precision*, struct background*, struct thermo*))
      BACKEND_REQ(class_perturb_initialize,(class_tag),int,(struct precision*, struct background*, struct thermo*, struct perturbs*))
	  BACKEND_REQ(class_primordial_initialize,(class_tag),int,(struct precision*, struct perturbs*, struct primordial*))
	  BACKEND_REQ(class_nonlinear_initialize,(class_tag),int,(struct precision*, struct background*, struct thermo*, struct perturbs*, struct primordial*, struct nonlinear*))
  	  BACKEND_REQ(class_transfer_initialize,(class_tag),int,(struct precision*, struct background*, struct thermo*, struct perturbs*, struct nonlinear*, struct transfers*))
	  BACKEND_REQ(class_spectra_initialize,(class_tag),int,(struct precision*, struct background*, struct perturbs*, struct primordial*, struct nonlinear*, struct transfers*, struct spectra*))
  	  BACKEND_REQ(class_lensing_initialize,(class_tag),int,(struct precision*, struct perturbs*,struct spectra*, struct nonlinear*,struct lensing*))
      BACKEND_REQ(class_output_initialize,(class_tag),int,(struct background*,struct thermo*, struct perturbs*, struct primordial*, struct nonlinear*, struct transfers*, struct spectra*, struct nonlinear*, struct lensing*, struct output*))
	  BACKEND_REQ(class_parser_initialize,(class_tag),int,(struct file_content*,int ,char*, char*))
	  BACKEND_REQ(class_output_total_cl_at_l,(class_tag),int, (struct spectra*, struct lensing* , struct output*, int, double* ))
      BACKEND_REQ(class_lensing_free,(class_tag),int, (struct lensing*))
      BACKEND_REQ(class_spectra_free,(class_tag),int, (struct spectra*))
      BACKEND_REQ(class_transfer_free,(class_tag),int, (struct transfers*))
      BACKEND_REQ(class_nonlinear_free,(class_tag),int, (struct nonlinear*))
	  BACKEND_REQ(class_primordial_free,(class_tag),int, (struct primordial*))
      BACKEND_REQ(class_perturb_free,(class_tag),int, (struct perturbs*))
      BACKEND_REQ(class_thermodynamics_free,(class_tag),int, (struct thermo*))
	  BACKEND_REQ(class_background_free,(class_tag),int, (struct background*))
	  #undef FUNCTION
	#undef CAPABILITY


  #define CAPABILITY compute_LCDMtensor_lowp_TT_loglike
  START_CAPABILITY
    #define FUNCTION function_LCDMtensor_lowp_TT_loglike
    START_FUNCTION(double)
      ALLOW_MODELS(LCDMtensor,plik_dx11dr2_HM_v18_TT)
      BACKEND_REQ(clik_compute_loglike, (clik_tag), double ,(clik_object*,double*,clik_error**))
      BACKEND_REQ(clik_initialize, (clik_tag), clik_object* ,(char*,clik_error**))
      BACKEND_REQ(clik_initialize_error, (clik_tag), clik_error* ,())
      BACKEND_REQ(clik_get_list_cls, (clik_tag), void , (clik_object*,int,clik_error**))
      BACKEND_REQ(clik_get_extra_parameter_names, (clik_tag), int ,(clik_object*,parname**,clik_error**))
      BACKEND_REQ(clik_get_lmax, (clik_tag), void ,(clik_object*,int,clik_error**))
      BACKEND_REQ(clik_cleanup, (clik_tag), void ,(clik_object**))
      BACKEND_REQ(clik_get_version, (clik_tag), char* , (clik_object*,clik_error**))
      BACKEND_REQ(data_initialize, (clik_tag), void, ())
      BACKEND_REQ(return_high_TT,(clik_tag),clik_object*,())
      BACKEND_REQ(return_lowp_TT,(clik_tag),clik_object*,())
      BACKEND_REQ(class_input_initialize,(class_tag),int,(struct file_content*, struct precision*, struct background*, struct thermo*, struct perturbs*, struct transfers*, struct primordial*, struct spectra*, struct nonlinear*, struct lensing*, struct output*, char*))
      BACKEND_REQ(class_background_initialize,(class_tag),int,(struct precision*, struct background*))
      BACKEND_REQ(class_thermodynamics_initialize,(class_tag),int,(struct precision*, struct background*, struct thermo*))
      BACKEND_REQ(class_perturb_initialize,(class_tag),int,(struct precision*, struct background*, struct thermo*, struct perturbs*))
      BACKEND_REQ(class_primordial_initialize,(class_tag),int,(struct precision*, struct perturbs*, struct primordial*))
      BACKEND_REQ(class_nonlinear_initialize,(class_tag),int,(struct precision*, struct background*, struct thermo*, struct perturbs*, struct primordial*, struct nonlinear*))
      BACKEND_REQ(class_transfer_initialize,(class_tag),int,(struct precision*, struct background*, struct thermo*, struct perturbs*, struct nonlinear*, struct transfers*))
      BACKEND_REQ(class_spectra_initialize,(class_tag),int,(struct precision*, struct background*, struct perturbs*, struct primordial*, struct nonlinear*, struct transfers*, struct spectra*))
      BACKEND_REQ(class_lensing_initialize,(class_tag),int,(struct precision*, struct perturbs*,struct spectra*, struct nonlinear*,struct lensing*))
      BACKEND_REQ(class_output_initialize,(class_tag),int,(struct background*,struct thermo*, struct perturbs*, struct primordial*, struct nonlinear*, struct transfers*, struct spectra*, struct nonlinear*, struct lensing*, struct output*))
      BACKEND_REQ(class_parser_initialize,(class_tag),int,(struct file_content*,int ,char*, char*))
      BACKEND_REQ(class_output_total_cl_at_l,(class_tag),int, (struct spectra*, struct lensing* , struct output*, int, double* ))
      BACKEND_REQ(class_lensing_free,(class_tag),int, (struct lensing*))
      BACKEND_REQ(class_spectra_free,(class_tag),int, (struct spectra*))
      BACKEND_REQ(class_transfer_free,(class_tag),int, (struct transfers*))
      BACKEND_REQ(class_nonlinear_free,(class_tag),int, (struct nonlinear*))
      BACKEND_REQ(class_primordial_free,(class_tag),int, (struct primordial*))
      BACKEND_REQ(class_perturb_free,(class_tag),int, (struct perturbs*))
      BACKEND_REQ(class_thermodynamics_free,(class_tag),int, (struct thermo*))
      BACKEND_REQ(class_background_free,(class_tag),int, (struct background*))
      #undef FUNCTION
    #undef CAPABILITY


  #define CAPABILITY compute_LCDMtensor_inflation_lowp_TT_loglike
  START_CAPABILITY
    #define FUNCTION function_LCDMtensor_inflation_lowp_TT_loglike
    START_FUNCTION(double)
      ALLOW_MODELS(LCDMtensor,plik_dx11dr2_HM_v18_TT)
      BACKEND_REQ(multimodecode_gambit_driver,(modecode_tag), void, (gambit_inflation_observables*,int& ,int& ,  int& ,  int& ,int& ,  int& ,  int& ,  int& ,  int& ,int& ,double& ,int& ,  int& ,double& ,int& ,double*,double*,  int& ,  int& ,double*,double*,double*,double& ,double& ,double& ,  int& ,int& ,double& ,double*,double*,double*,double*,double& ,double&))
      BACKEND_REQ(multimodecode_gambit_driver_test,(modecode_tag), void, (gambit_inflation_observables*,int& ,int& ,  int& ,  int& ,int& ,  int& ,  int& ,  int& ,  int& ,int& ,double& ,int& ,  int& ,double& ,int& ,  int& ,  int& ,double& ,double& ,double& ,  int& ,int& ,double& ,double& ,double&))
      BACKEND_REQ(clik_compute_loglike, (clik_tag), double ,(clik_object*,double*,clik_error**))
      BACKEND_REQ(clik_initialize, (clik_tag), clik_object* ,(char*,clik_error**))
      BACKEND_REQ(clik_initialize_error, (clik_tag), clik_error* ,())
      BACKEND_REQ(clik_get_list_cls, (clik_tag), void , (clik_object*,int,clik_error**))
      BACKEND_REQ(clik_get_extra_parameter_names, (clik_tag), int ,(clik_object*,parname**,clik_error**))
      BACKEND_REQ(clik_get_lmax, (clik_tag), void ,(clik_object*,int,clik_error**))
      BACKEND_REQ(clik_cleanup, (clik_tag), void ,(clik_object**))
      BACKEND_REQ(clik_get_version, (clik_tag), char* , (clik_object*,clik_error**))
      BACKEND_REQ(data_initialize, (clik_tag), void, ())
      BACKEND_REQ(return_high_TT,(clik_tag),clik_object*,())
      BACKEND_REQ(return_lowp_TT,(clik_tag),clik_object*,())
      BACKEND_REQ(class_input_initialize,(class_tag),int,(struct file_content*, struct precision*, struct background*, struct thermo*, struct perturbs*, struct transfers*, struct primordial*, struct spectra*, struct nonlinear*, struct lensing*, struct output*, char*))
      BACKEND_REQ(class_background_initialize,(class_tag),int,(struct precision*, struct background*))
      BACKEND_REQ(class_thermodynamics_initialize,(class_tag),int,(struct precision*, struct background*, struct thermo*))
      BACKEND_REQ(class_perturb_initialize,(class_tag),int,(struct precision*, struct background*, struct thermo*, struct perturbs*))
      BACKEND_REQ(class_primordial_initialize,(class_tag),int,(struct precision*, struct perturbs*, struct primordial*))
      BACKEND_REQ(class_nonlinear_initialize,(class_tag),int,(struct precision*, struct background*, struct thermo*, struct perturbs*, struct primordial*, struct nonlinear*))
      BACKEND_REQ(class_transfer_initialize,(class_tag),int,(struct precision*, struct background*, struct thermo*, struct perturbs*, struct nonlinear*, struct transfers*))
      BACKEND_REQ(class_spectra_initialize,(class_tag),int,(struct precision*, struct background*, struct perturbs*, struct primordial*, struct nonlinear*, struct transfers*, struct spectra*))
      BACKEND_REQ(class_lensing_initialize,(class_tag),int,(struct precision*, struct perturbs*,struct spectra*, struct nonlinear*,struct lensing*))
      BACKEND_REQ(class_output_initialize,(class_tag),int,(struct background*,struct thermo*, struct perturbs*, struct primordial*, struct nonlinear*, struct transfers*, struct spectra*, struct nonlinear*, struct lensing*, struct output*))
      BACKEND_REQ(class_parser_initialize,(class_tag),int,(struct file_content*,int ,char*, char*))
      BACKEND_REQ(class_output_total_cl_at_l,(class_tag),int, (struct spectra*, struct lensing* , struct output*, int, double* ))
      BACKEND_REQ(class_lensing_free,(class_tag),int, (struct lensing*))
      BACKEND_REQ(class_spectra_free,(class_tag),int, (struct spectra*))
      BACKEND_REQ(class_transfer_free,(class_tag),int, (struct transfers*))
      BACKEND_REQ(class_nonlinear_free,(class_tag),int, (struct nonlinear*))
      BACKEND_REQ(class_primordial_free,(class_tag),int, (struct primordial*))
      BACKEND_REQ(class_perturb_free,(class_tag),int, (struct perturbs*))
      BACKEND_REQ(class_thermodynamics_free,(class_tag),int, (struct thermo*))
      BACKEND_REQ(class_background_free,(class_tag),int, (struct background*))
      #undef FUNCTION
    #undef CAPABILITY


#undef MODULE

#endif /* defined __CosmoBitDEV_rollcall_hpp__ */
