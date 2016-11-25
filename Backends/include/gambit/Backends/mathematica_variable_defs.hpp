//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Definitions of the mathematica wrapper variables
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
#include "gambit/Backends/mathematica_variable.hpp"

namespace Gambit
{
  namespace Backends
  {

    // Constructor  
    template <typename T>                                                                       
    mathematica_variable<T>::mathematica_variable(WSLINK WSlink, str symbol) : 
    _WSlink(WSlink), _symbol(symbol)
    {
      //using namespace CAT_3(BACKENDNAME,_,SAFE_VERSION);

      /* If TYPE is a numeric type, send N first */
      if(boost::is_same<T, int>::value or 
         boost::is_same<T, float>::value or 
         boost::is_same<T, double>::value)
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
      if(!boost::is_same<T, void>::value)
      {
        if(!boost::is_same<T, int>::value and
           !boost::is_same<T, float>::value and
           !boost::is_same<T, double>::value and
           !boost::is_same<T, bool>::value and
           !boost::is_same<T, char>::value and
           !boost::is_same<T, str>::value)
          backend_warning().raise(LOCAL_INFO, "Error, WSTP type nor recognised");
        else if(!WSGetVariable(WSlink, &_value))
          backend_warning().raise(LOCAL_INFO, "Error reading packet from WSTP");

      } 
    }

    // Assignment operator with TYPE
    template <typename T>
    mathematica_variable<T>& mathematica_variable<T>::operator=(const T& val)
    {
      /* Send the expression SYMBOL = val */
      if(!WSPutSymbol(_WSlink, _symbol.c_str()) or
         !WSPutSymbol(_WSlink, "="))
        backend_warning().raise(LOCAL_INFO, "Error sending packet through WSTP");

      if(!boost::is_same<T, void>::value)
      {
        if(!boost::is_same<T, int>::value and
           !boost::is_same<T, float>::value and
           !boost::is_same<T, double>::value and
           !boost::is_same<T, bool>::value and
           !boost::is_same<T, char>::value and
           !boost::is_same<T, str>::value)
          backend_warning().raise(LOCAL_INFO, "Error, WSTP type nor recognised");
        else if(!WSPutVariable(_WSlink, val))
          backend_warning().raise(LOCAL_INFO, "Error reading packet from WSTP");

      }

      /* Mark the end of the message */
      if(!WSEndPacket(_WSlink))
        backend_warning().raise(LOCAL_INFO, "Error sending packet through WSTP");

    }

    // Cast operator for type TYPE
    template <typename T> 
    mathematica_variable<T>::operator T const()
    {
      return _value;
    }

    // Overloaded functions to get data through WSTP
    template <typename T>
    int mathematica_variable<T>::WSGetVariable(WSLINK WSlink, int* val)
    {
      return WSGetInteger32(WSlink, val);
    }

    template <typename T>
    int mathematica_variable<T>::WSGetVariable(WSLINK WSlink, float* val)
    {
      return WSGetReal32(WSlink, val);
    }

    template <typename T>
    int mathematica_variable<T>::WSGetVariable(WSLINK WSlink, double* val)
    {
      return WSGetReal64(WSlink, val);
    }

    template <typename T>
    int mathematica_variable<T>::WSGetVariable(WSLINK WSlink, bool* val)
    {
      return WSGetInteger8(WSlink, (char *)val);
    }

    template <typename T>
    int mathematica_variable<T>::WSGetVariable(WSLINK WSlink, char* val)
    {
      return WSGetInteger8(WSlink, val);
    }

    template <typename T>
    int mathematica_variable<T>::WSGetVariable(WSLINK WSlink, str* val)
    {
      return WSGetString(WSlink, (char **)val);
    }

    // Overloaded functions to put data through WSTP
    template <typename T>
    int mathematica_variable<T>::WSPutVariable(WSLINK WSlink, int val)
    {
      return WSPutInteger32(WSlink, val);
    }

    template <typename T>
    int mathematica_variable<T>::WSPutVariable(WSLINK WSlink, float val)
    {
      return WSPutReal32(WSlink, val);
    }

    template <typename T>
    int mathematica_variable<T>::WSPutVariable(WSLINK WSlink, double val)
    {
      return WSPutReal64(WSlink, val);
    }

    template <typename T>
    int mathematica_variable<T>::WSPutVariable(WSLINK WSlink, bool val)
    {
      return WSPutInteger8(WSlink, (char)val);
    }

    template <typename T>
    int mathematica_variable<T>::WSPutVariable(WSLINK WSlink, char val)
    {
      return WSPutInteger8(WSlink, val);
    }

    template <typename T>
    int mathematica_variable<T>::WSPutVariable(WSLINK WSlink, str val)
    {
      return WSPutString(WSlink, (char *)val);
    }


  }
}
