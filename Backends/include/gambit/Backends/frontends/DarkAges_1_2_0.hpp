//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend header for the DarkAges backend
///
///  Compile-time registration of available
///  functions and variables from this backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Patrick Stoecker
///          (stoecker@physik.rwth-aachen.de)
///  \date 2019 Oct
///
///  *********************************************


#define BACKENDNAME DarkAges
#define BACKENDLANG Python
#define VERSION 1.2.0
#define SAFE_VERSION 1_2_0

LOAD_LIBRARY

/* Include an extra header of pybind11.
 * Needed for casting of python list into std::vector and vice versa
 */
#include <pybind11/numpy.h>

/* Syntax for BE_FUNCTION (same as for any other backend):
 * BE_FUNCTION([choose function name], [type], [arguement types], "[exact symbol name]", "[choose capability name]")
 */

BE_FUNCTION(initialize, void, (), "initialize", "DA_initialize")
BE_FUNCTION(calc_f_decay, void, (pybind11::array_t<double>, pybind11::array_t<double>, pybind11::array_t<double>, double, double), "calculate_f_for_decay", "DA_calc_f")
BE_FUNCTION(get_result, pybind11::array_t<double>, (str), "get_result","DA_get_result")

/* Syntax for BE_VARIABLE:
 * BE_VARIABLE([name], [type], "[exact symbol name]", "[choose capability name]")
 * */

BE_VARIABLE(already_calculated, bool, "alreadyCalculated", "DA_alreadyCalculated")

/* We use BE_INI_DEPENDENCY, since DarkAges needs the spectra of injected electrons, positrons and photons
 * to calculate f(z)
 * */

BE_INI_DEPENDENCY(injection_spectrum,DarkAges::injectionSpectrum)
BE_INI_DEPENDENCY(DM_mass,double)
BE_INI_DEPENDENCY(lifetime,double)

/* Now register any convenience functions and wrap them in functors.
 *
 * Syntax for BE_CONV_FUNCTION:
 * BE_CONV_FUNCTION([function name], type, (arguments), "[choose capability name]") */

BE_CONV_FUNCTION(gather_results, DarkAges::fz_table, (), "DA_efficiency_function")

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
