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
///  \date 2018 Mar
///
///  *********************************************


#define BACKENDNAME DarkAges
#define BACKENDLANG Python
#define VERSION 1.0
#define SAFE_VERSION 1_0

LOAD_LIBRARY

/* Next we use macros BE_VARIABLE and BE_FUNCTION to extract pointers
 * to the variables and functions within the Python module.
 *
 * The macros create functors that wrap the library variables and functions.
 * These are used by the Core for dependency resolution and to set up a suitable
 * interface to the library functions/variables at module level. */

/* Syntax for BE_FUNCTION (same as for any other backend):
 * BE_FUNCTION([choose function name], [type], [arguement types], "[exact symbol name]", "[choose capability name]")
 */

BE_FUNCTION(initialize, void, (int), "initialize", "DA_initialize")
BE_FUNCTION(multiplyToArray, void, (double), "multiplyToArray", "DA_multiply")
BE_FUNCTION(returnArray, pybind11::object, (), "returnArray","DA_returnArray")
BE_FUNCTION(returnSumOfArray, double, (), "returnSumOfArray","DA_returnSumOfArray")

/* Syntax for BE_VARIABLE:
 * BE_VARIABLE([name], [type], "[exact symbol name]", "[choose capability name]")
 * */

BE_VARIABLE(arrayLen, int, "arrayLen", "arrayLen")
BE_VARIABLE(someFactor, double, "someFactor", "SomeFactor")

/* At this point we have a minimal interface to the loaded library.
 * Any additional convenience functions could be constructed below
 * using the available pointers. All convenience functions must be
 * registred/wrapped via the macro BE_CONV_FUNCTION (see below). */

BE_NAMESPACE
{
  /* Convenience functions go here */
  void awesomenessNeitherByAndersNorByPat(std::vector<double>& result)
  {
    logger().send("Message from 'awesomenessNeitherByAndersNorByPat' backend convenience function in DarkAges v1.0 wrapper",LogTags::info);
    multiplyToArray(*someFactor);
    pybind11::object tmp_result;
    tmp_result = returnArray();
    /* Mapping from numpy array onto std::vector goes here */
    result.resize(10,2.5); // this is a dummy such that somethign happens.
  }

  void greatnessNeitherByAndersNorByPat(double& result)
  {
    logger().send("Message from 'greatnessNeitherByAndersNorByPat' backend convenience function in DarkAges v1.0 wrapper",LogTags::info);
    multiplyToArray(*someFactor);
    result = returnSumOfArray();
  }
}
END_BE_NAMESPACE

/* Note that BE_NAMESPACE is just
 * namespace Gambit
 * {
 *   namespace Backends
 *   {
 *     namespace CAT_3(BACKENDNAME,_,SAFE_VERSION)
 * and END_BE_NAMESPACE is just
 *   }
 * }
 */

/* Now register any convenience functions and wrap them in functors.
 *
 * Syntax for BE_CONV_FUNCTION:
 * BE_CONV_FUNCTION([function name], type, (arguments), "[choose capability name]") */

BE_CONV_FUNCTION(awesomenessNeitherByAndersNorByPat, void, (std::vector<double>&), "DA_awesomeness")
BE_CONV_FUNCTION(greatnessNeitherByAndersNorByPat, void, (double&), "DA_greatness")

BE_INI_FUNCTION
{
  static bool scan_level = true;
  if (scan_level)
  {
    *arrayLen = runOptions->getValueOrDef<int>(10, "arrayLen");
    initialize(*arrayLen);
  }
  scan_level = false;
  *someFactor = runOptions->getValueOrDef<double>(2.5, "someFactor");
}
END_BE_INI_FUNCTION

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
