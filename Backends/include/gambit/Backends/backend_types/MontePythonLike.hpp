//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Definitions of container classes
///  for the DarkAges backend.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Janina Renk
///          (janina.renk@fysik.su.se)
///  \date 2019 Jun
///
///  *********************************************

#ifndef __MontePythonLike_types_hpp__
#define __MontePythonLike_types_hpp__

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

namespace Gambit
{

  /// Shorthand for a string to pybing object map map to avoid commas in macros
  typedef std::map<std::string,pybind11::object> map_str_pyobj;
  
}

#endif // defined __DarkAges_types_hpp__
