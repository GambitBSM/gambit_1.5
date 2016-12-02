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
        TYPE _var;
        WSLINK _WSlink;
        str _symbol;

      public:

        // Constructor
        mathematica_variable(WSLINK WSlink, str symbol) :  _WSlink(WSlink), _symbol(symbol)
        {
        //using namespace CAT_3(BACKENDNAME,_,SAFE_VERSION);

          /* If TYPE is a numeric type, send N first */
          if(boost::is_same<TYPE, int>::value or
             boost::is_same<TYPE, float>::value or
             boost::is_same<TYPE, double>::value)
            if(!WSPutFunction(WSlink, "N", 1))
              backend_warning().raise(LOCAL_INFO, "Error sending packet through WSTP");

          /* Send the variable symbol */
          if(!WSPutSymbol(WSlink, symbol.c_str()))
            backend_warning().raise(LOCAL_INFO, "Error sending packet through WSTP");

          /* Mark the end of the message */
          if(!WSEndPacket(WSlink))
            backend_warning().raise(LOCAL_INFO, "Error sending packet through WSTP");

          /* Wait to receive a packet from the kernel */
          int pkt;
          while( (pkt = WSNextPacket(WSlink), pkt) && pkt != RETURNPKT)
          {
            WSNewPacket(WSlink);
            if (WSError(WSlink))
              backend_warning().raise(LOCAL_INFO, "Error reading packet from WSTP");
          }

          /* Read the received packet into the return value, unless it's void */
          if(!boost::is_same<TYPE, void>::value)
          {
            if(!boost::is_same<TYPE, int>::value and
               !boost::is_same<TYPE, float>::value and
               !boost::is_same<TYPE, double>::value and
               !boost::is_same<TYPE, bool>::value and
               !boost::is_same<TYPE, char>::value and
               !boost::is_same<TYPE, str>::value)
              backend_warning().raise(LOCAL_INFO, "Error, WSTP type nor recognised");
            else if(!WSGetVariable(WSlink, &_var))
              backend_warning().raise(LOCAL_INFO, "Error reading packet from WSTP");

          }
        }


        // Assignment operator with TYPE
        mathematica_variable& operator=(const TYPE& val)
        {
          /* Send the expression SYMBOL = val */
          if(!WSPutSymbol(_WSlink, _symbol.c_str()) or
             !WSPutSymbol(_WSlink, "="))
            backend_warning().raise(LOCAL_INFO, "Error sending packet through WSTP");

          if(!boost::is_same<TYPE, void>::value)
          {
            if(!boost::is_same<TYPE, int>::value and
               !boost::is_same<TYPE, float>::value and
               !boost::is_same<TYPE, double>::value and
               !boost::is_same<TYPE, bool>::value and
               !boost::is_same<TYPE, char>::value and
               !boost::is_same<TYPE, str>::value)
              backend_warning().raise(LOCAL_INFO, "Error, WSTP type nor recognised");
            else if(!WSPutVariable(_WSlink, val))
              backend_warning().raise(LOCAL_INFO, "Error reading packet from WSTP");

          }

          /* Mark the end of the message */
          if(!WSEndPacket(_WSlink))
            backend_warning().raise(LOCAL_INFO, "Error sending packet through WSTP");

          return *this;

        }


       // Cast operator for type TYPE
       operator TYPE const() { return _var; }

       // Overloaded functions to get data through WSTP
       int WSGetVariable(WSLINK WSlink, int* val) { return WSGetInteger(WSlink, val); }
       int WSGetVariable(WSLINK WSlink, float* val) { return WSGetReal32(WSlink, val); }
       int WSGetVariable(WSLINK WSlink, double* val) { return WSGetReal64(WSlink, val); } 
       int WSGetVariable(WSLINK WSlink, bool* val) { return WSGetInteger8(WSlink, (char *)val); }
       int WSGetVariable(WSLINK WSlink, char* val) { return WSGetInteger8(WSlink, val); }
       int WSGetVariable(WSLINK WSlink, str* val) {  return WSGetString(WSlink, (char **)val); } 

       // Overloaded functions to put data through WSTP
       int WSPutVariable(WSLINK WSlink, int val) { return WSPutInteger32(WSlink, val); }
       int WSPutVariable(WSLINK WSlink, float val) { return WSPutReal32(WSlink, val); }
       int WSPutVariable(WSLINK WSlink, double val) { return WSPutReal64(WSlink, val); }
       int WSPutVariable(WSLINK WSlink, bool val) { return WSPutInteger8(WSlink, (char)val); }
       int WSPutVariable(WSLINK WSlink, char val) {  return WSPutInteger8(WSlink, val); }
       int WSPutVariable(WSLINK WSlink, str val) { return WSPutString(WSlink, (char *)val); }

    };

  }
}
#endif /* __mathematica_variable_hpp__ */

#endif /* HAVE_MATHEMATICA */
