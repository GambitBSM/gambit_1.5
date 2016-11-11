//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Declarations of the mathematica wrapper functions
///
///  ***********************************************
///
///  Authors
///  =======
///
///  (add name and date if you modify)
///
///  \author Tomas Gonzalo
///          (t.e.gonzalo@fys.uio.no)
///  \date 2016 Nov
///
///  ***********************************************

#include "gambit/Utils/util_types.hpp"
#include "gambit/cmake/cmake_variables.hpp"

#ifdef HAVE_MATHEMATICA
#include MATHEMATICA_WSTP_H

#ifndef __mathematica_variable_hpp__
#define __mathematica_variable_hpp__

namespace Gambit
{
  namespace Backends
  {
  
    template <typename TYPE>
    class mathematica_variable : public TYPE
    {
      public:

       // Constructor
       mathematica_variable();

       // Assignment operator with class
       mathematica_variable& operator=(const mathematica_variable&);

       // Assignment operator with TYPE
       mathematica_variable& operator=(const TYPE&);

    };


  }
}
#endif /* __mathematica_variable_hpp__ */

#endif /* HAVE_MATHEMATICA */
