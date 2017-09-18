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
    ALLOW_MODELS(StandardModel_SLHA2, StandardModel_Neutrinos)
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

  #define CAPABILITY SeesawI_Vnu
  START_CAPABILITY
    #define FUNCTION Vnu
    START_FUNCTION(Eigen::Matrix3cd)
    DEPENDENCY(UPMNS, Eigen::Matrix3cd)
    ALLOW_MODELS(SN_dev)
    #undef FUNCTION
  #undef CAPABILITY

  // Unitarity checks for the mixing matrix
  #define CAPABILITY Unitarity
  START_CAPABILITY
    // Unitarity of the PMNS matrix
    #define FUNCTION Unitarity_UPMNS
    START_FUNCTION(bool)
    DEPENDENCY(UPMNS, Eigen::Matrix3cd)
    DEPENDENCY(m_nu, Eigen::Matrix3cd)
    #undef FUNCTION

    // Unitarity of the mixing matrix in seesaw I
    #define FUNCTION Unitarity_SeesawI
    START_FUNCTION(bool)
    DEPENDENCY(SeesawI_Theta, Eigen::Matrix3cd)
    DEPENDENCY(SeesawI_Vnu, Eigen::Matrix3cd)
    DEPENDENCY(m_nu, Eigen::Matrix3cd)
    ALLOW_MODEL(SN_dev)
    #undef FUNCTION
    
  #undef CAPABILITY

#undef MODULE


#endif /* defined(__NeutrinoBit_rollcall_hpp__) */


