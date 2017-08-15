//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Rollcall header for NeutrinoBit.
///
///  Compile-time registration of available
///  observables and likelihoods for neutrino observables.
///
///  Don't put typedefs or other type definitions
///  in this file; see
///  Core/include/types_rollcall.hpp for further
///  instructions on how to add new types.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Tomas Gonzalo
///          (t.e.gonzalo@fys.uio.no)
///  \date 2017 July
///
///  *********************************************


#ifndef __NeutrinoBit_rollcall_hpp__
#define __NeutrinoBit_rollcall_hpp__

#define MODULE NeutrinoBit
START_MODULE

  // Neutrino mass matrix
  #define CAPABILITY m_nu
  START_CAPABILITY
    // Neutrino masss matrix
    #define FUNCTION M_nu
    START_FUNCTION(Eigen::Matrix3cd)
    ALLOW_MODELS(StandardModel_Neutrinos)
    #undef FUNCTION

  #undef CAPABILITY

  // Construct the neutrino PMNS matrix
  #define CAPABILITY UPMNS
  START_CAPABILITY
    // Neutrino PMNS matrix in a parametrization where the charged Yukawas are diagonal
    #define FUNCTION UPMNS
    START_FUNCTION(Eigen::Matrix3cd)
    ALLOW_MODELS(StandardModel_Neutrinos)
    #undef FUNCTION
  #undef CAPABILITY


  // Construct the neutrino Seesaw I type Theta matrix
  #define CAPABILITY SeesawI_Theta
  START_CAPABILITY
    // Theta in the Casas-Ibarra parametrization
    #define FUNCTION CI_Theta
    START_FUNCTION(Eigen::Matrix3cd)
    DEPENDENCY(m_nu, Eigen::Matrix3cd)
    DEPENDENCY(UPMNS, Eigen::Matrix3cd)
    DEPENDENCY(SMINPUTS, SMInputs)
    ALLOW_MODELS(StandardModel_SLHA2, StandardModel_Higgs, SN_dev)
    #undef FUNCTION
  #undef CAPABILITY

  // Corrected UMPNS for Seesaw type I
  #define CAPABILITY SeesawI_UPMNS
  START_CAPABILITY
    #define FUNCTION V_nu
    START_FUNCTION(Eigen::Matrix3cd)
    ALLOW_MODELS(SN_dev)
    #undef FUNCTION
  #undef CAPABILITY

#undef MODULE


#endif /* defined(__NeutrinoBit_rollcall_hpp__) */


