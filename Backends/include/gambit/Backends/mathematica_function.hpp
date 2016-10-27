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
///  \date 2016 Oct
///
///  ***********************************************

#include "gambit/Utils/util_types.hpp"
#include "gambit/cmake/cmake_variables.hpp"

#ifdef HAVE_MATHEMATICA
#include MATHEMATICA_WSTP_H

#ifndef __mathematica_function_hpp__
#define __mathematica_function_hpp__

namespace Gambit
{
  namespace Backends
  {
  
    template<typename TYPE, typename... ARGS>
    TYPE some_function(ARGS...);

    template <typename TYPE, typename... ARGS>
    class mathematica_function : std::function<TYPE(ARGS...)>
    {
      public:

       // Constructor
       mathematica_function(void *, str);

       // Getter for the function name
       str get_function_name();

       // Wrapper to execute the function and provide the return value
       TYPE operator()(ARGS&&...);

      protected:

        // Pointer to the WSTP link stablished during loading 
        WSLINK WSlink;

        // Name of the function to the called through WSTP
        str function_name;

    };


  }
}

///
/// Definitions of the mathematica wrapper functions
///
/// **********************************************************
///

namespace Gambit
{
  namespace Backends
  {

    // Test
    template <typename TYPE, typename... ARGS>
    TYPE some_function(ARGS... args)
    {
      cout << "function" << endl;
    }

    // Constructor
    template <typename TYPE, typename... ARGS>
    mathematica_function<TYPE, ARGS...>::mathematica_function(void *pHandle, str symbol_name) :
       WSlink((WSLINK) pHandle),
       function_name(symbol_name)
    {
    }

    // Getter for the function name
    template <typename TYPE, typename... ARGS>
    str mathematica_function<TYPE, ARGS...>::get_function_name()
    {
      return function_name;
    }

    // Actual function that takes care of the wrapping
    template <typename TYPE, typename... ARGS>
    TYPE mathematica_function<TYPE, ARGS...>::operator()(ARGS&&... args)
    {

      try
      {

        TYPE return_value;

        if(!WSPutFunction(WSlink, "CalculateSquare", 1)
            or !WSPutInteger(WSlink, 7)
            or !WSEndPacket(WSlink))
          cout << "Error sending packet" << endl;

        int pkt;
        while( (pkt = WSNextPacket(WSlink), pkt) && pkt != RETURNPKT)
        {
          WSNewPacket(WSlink);
          if (WSError(WSlink)) cout << "Error reading package" << endl;
        }

        double square;
        if (WSGetReal(WSlink, &square))
        {
          cout << "Calculate square of " << 7 << " is " << square << endl;
        }
        else
        {
          cout << "Error" << endl;;
        }

        return return_value;

      }
      catch(std::exception &e)
      {
        //TODO: Some backend error handling
      }
    }
  }
}

#endif /* __mathematica_function_hpp__ */

#endif /* HAVE_MATHEMATICA */
