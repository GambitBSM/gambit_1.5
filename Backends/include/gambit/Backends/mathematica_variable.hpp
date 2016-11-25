//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Declarations of the mathematica wrapper variables 
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
    class mathematica_variable
    {

      private:
        TYPE _value;
        WSLINK _WSlink;
        str _symbol;

      public:

       // Constructor
       mathematica_variable(WSLINK, str);

       // Assignment operator with TYPE
       mathematica_variable& operator=(const TYPE&);

       // Cast operator for type TYPE
       operator TYPE const();

      // Overloaded functions to get data through WSTP
      int WSGetVariable(WSLINK, int*);
      int WSGetVariable(WSLINK, float*);
      int WSGetVariable(WSLINK, double*);
      int WSGetVariable(WSLINK, bool*);
      int WSGetVariable(WSLINK, char*);
      int WSGetVariable(WSLINK, str*);

      // Overloaded functions to put data through WSTP
      int WSPutVariable(WSLINK, int);
      int WSPutVariable(WSLINK, float);
      int WSPutVariable(WSLINK, double);
      int WSPutVariable(WSLINK, bool);
      int WSPutVariable(WSLINK, char);
      int WSPutVariable(WSLINK, str);


    };


  }
}
#endif /* __mathematica_variable_hpp__ */

#endif /* HAVE_MATHEMATICA */
