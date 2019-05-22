//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend header for the class 2.6.3 backend.
///
///  Compile-time registration of available
///  functions and variables from this backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Selim C. Hotinli
///          (selim.hotinli14@pimperial.ac.uk)
///  \date 2016 Oct
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2017 Nov
///  \date 2018 Apr, May
///  
///	 \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2018 Oct
///
///  *********************************************

#define BACKENDNAME class
#define BACKENDLANG CC
#define VERSION 2.6.3
#define SAFE_VERSION 2_6_3

// Load it
LOAD_LIBRARY

BE_FUNCTION(class_input_initialize, int,(Class::file_content*, Class::precision*, Class::background*, Class::thermo*, Class::perturbs*, Class::transfers*, Class::primordial*, Class::spectra*, Class::nonlinear*, Class::lensing*, Class::output*, char*),"input_init","class_input_initialize")
BE_FUNCTION(class_background_initialize, int,(Class::precision*, Class::background*), "background_init","class_background_initialize")
BE_FUNCTION(class_thermodynamics_initialize, int,(Class::precision*, Class::background*, Class::thermo*), "thermodynamics_init", "class_thermodynamics_initialize")
BE_FUNCTION(class_perturb_initialize,int,(Class::precision*, Class::background*, Class::thermo*, Class::perturbs*), "perturb_init","class_perturb_initialize")
BE_FUNCTION(class_primordial_initialize,int,(Class::precision*, Class::perturbs*, Class::primordial*), "primordial_init","class_primordial_initialize")
BE_FUNCTION(class_nonlinear_initialize,int,(Class::precision*, Class::background*, Class::thermo*, Class::perturbs*, Class::primordial*, Class::nonlinear*), "nonlinear_init","class_nonlinear_initialize")
BE_FUNCTION(class_transfer_initialize,int,(Class::precision*, Class::background*, Class::thermo*, Class::perturbs*, Class::nonlinear*, Class::transfers*), "transfer_init","class_transfer_initialize")
BE_FUNCTION(class_spectra_initialize,int,(Class::precision*, Class::background*, Class::perturbs*, Class::primordial*, Class::nonlinear*, Class::transfers*, Class::spectra*), "spectra_init","class_spectra_initialize")
BE_FUNCTION(class_lensing_initialize,int,(Class::precision*, Class::perturbs*,Class::spectra*, Class::nonlinear*,Class::lensing*), "lensing_init", "class_lensing_initialize")
BE_FUNCTION(class_output_initialize,int,(Class::background*,Class::thermo*, Class::perturbs*, Class::primordial*, Class::nonlinear*, Class::transfers*, Class::spectra*, Class::nonlinear*, Class::lensing*, Class::output*), "output_init", "class_output_initialize")
BE_FUNCTION(class_parser_initialize, int, (Class::file_content*,int ,char*, char*), "parser_init", "class_parser_initialize")
BE_FUNCTION(class_output_total_cl_at_l,int,(Class::spectra*, Class::lensing* , Class::output*, int, double* ), "output_total_cl_at_l" , "class_output_total_cl_at_l")
BE_FUNCTION(class_lensing_free, int, (Class::lensing*), "lensing_free", "class_lensing_free")
BE_FUNCTION(class_spectra_free, int, (Class::spectra*), "spectra_free", "class_spectra_free")
BE_FUNCTION(class_transfer_free, int, (Class::transfers*), "transfer_free", "class_transfer_free")
BE_FUNCTION(class_nonlinear_free, int, (Class::nonlinear*), "nonlinear_free", "class_nonlinear_free")
BE_FUNCTION(class_primordial_free, int, (Class::primordial*), "primordial_free", "class_primordial_free")
BE_FUNCTION(class_perturb_free, int, (Class::perturbs*), "perturb_free", "class_perturb_free")
BE_FUNCTION(class_thermodynamics_free, int, (Class::thermo*), "thermodynamics_free", "class_thermodynamics_free")
BE_FUNCTION(class_background_free, int, (Class::background*), "background_free", "class_background_free")


BE_FUNCTION(background_tau_of_z, int, (Class::background*,double,double*), "background_tau_of_z", "class_background_tau_of_z")
BE_FUNCTION(background_at_tau, int, (Class::background*,double,short,short,int*,double*), "background_at_tau", "class_background_at_tau")
BE_FUNCTION(spectra_sigma, int, (Class::background*,Class::primordial*,Class::spectra*,double,double,double*), "spectra_sigma", "class_spectra_sigma")

BE_CONV_FUNCTION(class_get_cl, std::vector<double>, (str), "class_get_cl")
BE_CONV_FUNCTION(class_get_Da, double, (double), "class_get_Da")
BE_CONV_FUNCTION(class_get_Dl, double, (double), "class_get_Dl")
BE_CONV_FUNCTION(class_get_scale_independent_growth_factor, double, (double), "class_get_scale_independent_growth_factor")
BE_CONV_FUNCTION(class_get_scale_independent_growth_factor_f, double, (double), "class_get_scale_independent_growth_factor_f")
BE_CONV_FUNCTION(class_get_z_of_r, std::vector<std::vector<double>>, (std::vector<double>), "class_get_z_of_r")
BE_CONV_FUNCTION(class_get_Hz, double, (double), "class_get_Hz")
BE_CONV_FUNCTION(class_get_rs, double, (), "class_get_rs")
BE_CONV_FUNCTION(class_get_Omega_m, double, (), "class_get_Omega_m")
BE_CONV_FUNCTION(class_get_sigma8, double, (double), "class_get_sigma8")

//BE_CONV_FUNCTION(class_2_6_3_run, void, (), "class_run_internal")
//BE_CONV_FUNCTION(class_2_6_3_free, void, (), "class_free_internal")
BE_CONV_FUNCTION(get_class_2_6_3, CosmoBit::Class_container, (), "get_ptr_to_class")
//BE_CONV_FUNCTION(class_2_6_3_set_parameter, void, (CosmoBit::ClassInput), "class_set_parameter_internal")

BE_INI_DEPENDENCY(class_set_parameter,CosmoBit::Class_container)

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
