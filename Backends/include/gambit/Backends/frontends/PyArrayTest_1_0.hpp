//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Frontend header for the PyArrayTest backend
///  as a playground example for the Binding:
///
///   std::vector<T>  --> pybind11::list
///
///  and vice versa. 
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


#define BACKENDNAME PyArrayTest
#define BACKENDLANG Python
#define VERSION 1.0
#define SAFE_VERSION 1_0

LOAD_LIBRARY

/* Include an extra header of pybind11.
 * Needed for casting of python list into std::vector and vice versa
 */
#include <pybind11/stl.h> 

/* Next we use macros BE_VARIABLE and BE_FUNCTION to extract pointers
 * to the variables and functions within the Python module.
 *
 * The macros create functors that wrap the library variables and functions.
 * These are used by the Core for dependency resolution and to set up a suitable
 * interface to the library functions/variables at module level. */

/* Syntax for BE_FUNCTION (same as for any other backend):
 * BE_FUNCTION([choose function name], [type], [arguement types], "[exact symbol name]", "[choose capability name]")
 */

BE_FUNCTION(initialize, void, (int), "initialize", "PyArrayTest_initialize")
BE_FUNCTION(multiplyToArray, void, (double), "multiplyToArray", "PyArrayTest_multiply")
BE_FUNCTION(returnArray, pybind11::list, (), "returnArray","returnPythonArray")
BE_FUNCTION(readArray, double, (pybind11::list), "readArray","readPythonArray")

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
  std::vector<double> awesomeMapping_Python_to_Gambit()
  {
    logger().send("Message from 'AwesomeMapping_Python_to_Gambit' backend convenience function in PyArrayTest v1.0 wrapper",LogTags::info);
    multiplyToArray(*someFactor);
    pybind11::list tmp_result = returnArray();
    /* Mapping from numpy-array (or list) onto std::vector goes here */
    std::vector<double> result = pybind11::cast< std::vector<double> >(tmp_result); 
    return result;
  }

  double awesomeMapping_Gambit_to_Python(int len)
  {
    logger().send("Message from 'AwesomeMapping_Gambit_to_Python' backend convenience function in PyArrayTest v1.0 wrapper",LogTags::info);
    std::vector<double> tmp_vec;
    /* Populate the vector with some entries */
    for (int i=0; i<len; i++)
    {
      tmp_vec.push_back(1./pow(2.,i));
    }
    /* Map std::vector onto a python-list*/
    pybind11::list someExternalArray = pybind11::cast(tmp_vec);
    double result = readArray(someExternalArray);
    return result;
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

BE_CONV_FUNCTION(awesomeMapping_Python_to_Gambit, std::vector<double>, (), "PyArrayTest_Py_to_cpp")
BE_CONV_FUNCTION(awesomeMapping_Gambit_to_Python, double, (int), "PyArrayTest_cpp_to_Py")

BE_INI_FUNCTION
{
  static bool scan_level = true;
  if (scan_level)
  {
    *arrayLen = runOptions->getValueOrDef<int>(10, "arrayLen");
    *someFactor = runOptions->getValueOrDef<double>(2.5, "someFactor");
  }
  scan_level = false;
  initialize(*arrayLen);
}
END_BE_INI_FUNCTION

// Undefine macros to avoid conflict with other backends
#include "gambit/Backends/backend_undefs.hpp"
