//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend header for the class backend.
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
///
///  *********************************************

#define BACKENDNAME class
#define VERSION 2.6.1
#define SAFE_VERSION 2_6_1

// Load it
LOAD_LIBRARY

BE_FUNCTION(class_input_initialize, int,(struct class_2_6::file_content*, struct class_2_6::precision*, struct class_2_6::background*, struct class_2_6::thermo*, struct class_2_6::perturbs*, struct class_2_6::transfers*, struct class_2_6::primordial*, struct class_2_6::spectra*, struct class_2_6::nonlinear*, struct class_2_6::lensing*, struct class_2_6::output*, char*),"input_init","class_input_initialize")

BE_FUNCTION(class_background_init, int,(struct class_2_6::precision*, struct class_2_6::background*), "background_init","class_background_initialize")

BE_FUNCTION(class_thermodynamics_init, int,(struct class_2_6::precision*, struct class_2_6::background*, struct class_2_6::thermo*), "thermodynamics_init", "class_thermodynamics_initialize")

BE_FUNCTION(class_perturb_init,int,(struct class_2_6::precision*, struct class_2_6::background*, struct class_2_6::thermo*, struct class_2_6::perturbs*), "perturb_init","class_perturb_initialize")

BE_FUNCTION(class_primordial_init,int,(struct class_2_6::precision*, struct class_2_6::perturbs*, struct class_2_6::primordial*), "primordial_init","class_primordial_initialize")

BE_FUNCTION(class_nonlinear_init,int,(struct class_2_6::precision*, struct class_2_6::background*, struct class_2_6::thermo*, struct class_2_6::perturbs*, struct class_2_6::primordial*, struct class_2_6::nonlinear*), "nonlinear_init","class_nonlinear_initialize")

BE_FUNCTION(class_transfer_init,int,(struct class_2_6::precision*, struct class_2_6::background*, struct class_2_6::thermo*, struct class_2_6::perturbs*, struct class_2_6::nonlinear*, struct class_2_6::transfers*), "transfer_init","class_transfer_initialize")

BE_FUNCTION(class_spectra_init,int,(struct class_2_6::precision*, struct class_2_6::background*, struct class_2_6::perturbs*, struct class_2_6::primordial*, struct class_2_6::nonlinear*, struct class_2_6::transfers*, struct class_2_6::spectra*), "spectra_init","class_spectra_initialize")

BE_FUNCTION(class_lensing_init,int,(struct class_2_6::precision*, struct class_2_6::perturbs*,struct class_2_6::spectra*, struct class_2_6::nonlinear*,struct class_2_6::lensing*), "lensing_init", "class_lensing_initialize")

BE_FUNCTION(class_output_init,int,(struct class_2_6::background*,struct class_2_6::thermo*, struct class_2_6::perturbs*, struct class_2_6::primordial*, struct class_2_6::nonlinear*, struct class_2_6::transfers*, struct class_2_6::spectra*, struct class_2_6::nonlinear*, struct class_2_6::lensing*, struct class_2_6::output*), "output_init", "class_output_initialize")

BE_FUNCTION(class_parser_init, int, (struct class_2_6::file_content*,int ,char*, char*), "parser_init", "class_parser_initialize")

BE_FUNCTION(class_output_total_cl_at_l,int,(struct class_2_6::spectra*, struct class_2_6::lensing* , struct class_2_6::output*, int, double* ), "output_total_cl_at_l" , "class_output_total_cl_at_l")

BE_FUNCTION(class_lensing_free, int, (struct class_2_6::lensing*), "lensing_free", "class_lensing_free")

BE_FUNCTION(class_spectra_free, int, (struct class_2_6::spectra*), "spectra_free", "class_spectra_free")

BE_FUNCTION(class_transfer_free, int, (struct class_2_6::transfers*), "transfer_free", "class_transfer_free")

BE_FUNCTION(class_nonlinear_free, int, (struct class_2_6::nonlinear*), "nonlinear_free", "class_nonlinear_free")

BE_FUNCTION(class_primordial_free, int, (struct class_2_6::primordial*), "primordial_free", "class_primordial_free")

BE_FUNCTION(class_perturb_free, int, (struct class_2_6::perturbs*), "perturb_free", "class_perturb_free")

BE_FUNCTION(class_thermodynamics_free, int, (struct class_2_6::thermo*), "thermodynamics_free", "class_thermodynamics_free")

BE_FUNCTION(class_background_free, int, (struct class_2_6::background*), "background_free", "class_background_free")

namespace Gambit
{
	namespace Backends
	{
		namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
		{
			
			/* Convenience functions go here */
			
		} /* end namespace BACKENDNAME_SAFE_VERSION */
	} /* end namespace Backends */
} /* end namespace Gambit */


//BE_CONV_FUNCTION(awesomenessByAnders, double, "awesomeness")

BE_INI_FUNCTION{}
END_BE_INI_FUNCTION

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
