//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend header for the classy backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2019 June
///
///  \author Sanjay Bloor
///          (sanjay.bloor12@imperial.ac.uk)
///  \date 2019 June
///
///  *********************************************

#define BACKENDNAME classy
#define BACKENDLANG Python2
#define VERSION 2.6.3
#define SAFE_VERSION 2_6_3

LOAD_LIBRARY

//BE_FUNCTION(init, void, (), "init", "MontePythonLike_init")

BE_CONV_FUNCTION(path_to_classy, std::string, (), "path_to_classy")
//BE_CONV_FUNCTION(classy_create_class_instance, void, (pybind11::object&), 				 "classy_create_class_instance")
//BE_CONV_FUNCTION(classy_compute, 	 void, (CosmoBit::Classy_cosmo_container&), "classy_compute")
//BE_CONV_FUNCTION(classy_2_6_3_set_parameter, void, (pybind11::dict, pybind11::dict), "classy_2_6_3_set_parameter")

BE_CONV_FUNCTION(get_classy_cosmo_object, pybind11::object, (), "get_classy_cosmo_object")

//BE_CONV_FUNCTION(class_get_cl, std::vector<double>, (str), "class_get_cl")
BE_CONV_FUNCTION(class_get_Da, double, (double), "class_get_Da")
BE_CONV_FUNCTION(class_get_Dl, double, (double), "class_get_Dl")
BE_CONV_FUNCTION(class_get_scale_independent_growth_factor, double, (double), "class_get_scale_independent_growth_factor")
BE_CONV_FUNCTION(class_get_scale_independent_growth_factor_f, double, (double), "class_get_scale_independent_growth_factor_f")
BE_CONV_FUNCTION(class_get_Hz, double, (double), "class_get_Hz")
BE_CONV_FUNCTION(class_get_rs, double, (), "class_get_rs")
BE_CONV_FUNCTION(class_get_Omega_m, double, (), "class_get_Omega_m")
BE_CONV_FUNCTION(class_get_Omega0_nu, double, (), "class_get_Omega0_nu")
BE_CONV_FUNCTION(class_get_Omega0_Lambda, double, (), "class_get_Omega0_Lambda")
BE_CONV_FUNCTION(class_get_sigma8, double, (), "class_get_sigma8")
BE_CONV_FUNCTION(class_get_Neff, double, (), "class_get_Neff")

BE_INI_DEPENDENCY(get_Classy_cosmo_container,pybind11::dict)
// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
