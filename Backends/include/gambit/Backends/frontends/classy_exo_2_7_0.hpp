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
#define VERSION exo_2.7.0
#define SAFE_VERSION exo_2_7_0

LOAD_LIBRARY

//BE_FUNCTION(init, void, (), "init", "MontePythonLike_init")

BE_CONV_FUNCTION(path_to_classy, std::string, (), "path_to_classy")
BE_CONV_FUNCTION(classy_create_class_instance, void, (pybind11::object&), 				 "classy_create_class_instance")
BE_CONV_FUNCTION(classy_compute, 	 void, (CosmoBit::Classy_cosmo_container&), "classy_compute")

BE_CONV_FUNCTION(classy_set_energy_injection_efficiency_input, map_str_dblptr, (), "classy_set_energy_injection_efficiency_input")

BE_INI_CONDITIONAL_DEPENDENCY(energy_injection_efficiency, DarkAges::fz_table, TestDecayingDM)
// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
