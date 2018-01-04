//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Class python_variable, needed to overload
///  constructor and assignment operators to send
///  messages through pybind11.
///
///  ***********************************************
///
///  Authors
///  =======
///
///  (add name and date if you modify)
///
///  \author Pat Scott
///          (p.scott@imperial.ac.uk)
///  \date 2017 Dec
///
///  ***********************************************

#ifndef __python_variable_hpp__
#define __python_variable_hpp__

#include "gambit/Elements/ini_catch.hpp"
#include "gambit/Backends/backend_singleton.hpp"
#include "gambit/cmake/cmake_variables.hpp"

#include <pybind11/pybind11.h>

namespace Gambit
{

  namespace Backends
  {
    /// Holds the info about a python backend variable, and defines conversion functions.
    template <typename TYPE>
    class python_variable
    {

      private:

        /// Wrapper to the variable
        pybind11::object var;

        /// Indication of whether or not the function has been successfully loaded
        bool handle_works;

      public:

        /// Constructor
        python_variable(const str& be, const str& ver, const str& symbol) : handle_works(false)
        {
          try
          {
            // Extract the backend module pointer from the backendInfo object
            pybind11::object* mod;
            if (backendInfo().works.at(be+ver))
            {
              mod = backendInfo().loaded_python_backends.at(be+ver);
            }
            else
            {
              mod = (pybind11::module*)0;
              return;
            }

            // Extract the pointer to the variable from the module
            try
            {
              var = mod->attr(symbol.c_str());
              handle_works = true;
            }
            catch (std::exception& e)
            {
              std::ostringstream err;
              err << "Failed to retrieve handle to variable " << symbol << " from Python module for " << be+ver << endl
                  << "The backend variable associated with this symbol will be disabled (i.e. get status = -2)." << endl
                  << "Python error was: " << e.what() << endl;
              backend_warning().raise(LOCAL_INFO, err.str());
              backendInfo().dlerrors[be+ver] = symbol;
            }
          }
          catch (std::exception& e) { ini_catch(e); }
        }

        /// Assignment operator for python_variable from equivalent C++ type
        python_variable& operator=(const TYPE& val)
        {
          if (not handle_works) backend_error().raise(LOCAL_INFO, "Attempted to use a Python backend variable that was not successfully loaded.");
          var = pybind11::cast(val);
          return *this;
        }

        /// Cast operator from python_variable to equivalent C++ type
        operator TYPE const()
        {
          if (not handle_works) backend_error().raise(LOCAL_INFO, "Attempted to use a Python backend variable that was not successfully loaded.");
          return var.cast<TYPE>();
        }

    };

  }

}

#endif /* __python_variable_hpp__ */
