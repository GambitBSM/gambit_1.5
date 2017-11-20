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
///  \author Suraj Krishnamurthy
///          (S.Krishnamurthy@uva.nl)
///  \date 2017 Feb
///
///  \author Tomas Gonzalo
///          (t.e.gonzalo@fys.uio.no)
///  \date 2017 July
///
///  *********************************************


#ifndef __NeutrinoBit_rollcall_hpp__
#define __NeutrinoBit_rollcall_hpp__

#include <Eigen/Sparse>
#include <Eigen/Dense>

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


  #define CAPABILITY Theta_sq
  START_CAPABILITY

    #define FUNCTION CI_param
    START_FUNCTION(Eigen::Matrix3d)
    ALLOW_MODEL(SN_dev)
    DEPENDENCY(SeesawI_Theta, Eigen::Matrix3cd)
    #undef FUNCTION
 
  #undef CAPABILITY

  #define CAPABILITY Ue1
  START_CAPABILITY

    #define FUNCTION printable_Ue1
    START_FUNCTION(double)
    DEPENDENCY(Theta_sq, Eigen::Matrix3d)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY Um1
  START_CAPABILITY

    #define FUNCTION printable_Um1
    START_FUNCTION(double)
    DEPENDENCY(Theta_sq, Eigen::Matrix3d)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY Ut1
  START_CAPABILITY

    #define FUNCTION printable_Ut1
    START_FUNCTION(double)
    DEPENDENCY(Theta_sq, Eigen::Matrix3d)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY Ue2
  START_CAPABILITY

    #define FUNCTION printable_Ue2
    START_FUNCTION(double)
    DEPENDENCY(Theta_sq, Eigen::Matrix3d)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY Um2
  START_CAPABILITY

    #define FUNCTION printable_Um2
    START_FUNCTION(double)
    DEPENDENCY(Theta_sq, Eigen::Matrix3d)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY Ut2
  START_CAPABILITY

    #define FUNCTION printable_Ut2
    START_FUNCTION(double)
    DEPENDENCY(Theta_sq, Eigen::Matrix3d)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY Ue3
  START_CAPABILITY

    #define FUNCTION printable_Ue3
    START_FUNCTION(double)
    DEPENDENCY(Theta_sq, Eigen::Matrix3d)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY Um3
  START_CAPABILITY

    #define FUNCTION printable_Um3
    START_FUNCTION(double)
    DEPENDENCY(Theta_sq, Eigen::Matrix3d)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY Ut3
  START_CAPABILITY

    #define FUNCTION printable_Ut3
    START_FUNCTION(double)
    DEPENDENCY(Theta_sq, Eigen::Matrix3d)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY bbn_lifetime
  START_CAPABILITY
    #define FUNCTION SN_bbn_lifetime
    START_FUNCTION(std::vector<double>)
    DEPENDENCY(SMINPUTS, SMInputs)
    DEPENDENCY(Theta_sq, Eigen::Matrix3d)
    ALLOW_MODEL(SN_dev)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY lnLbbn
  START_CAPABILITY
    #define FUNCTION lnL_bbn
    START_FUNCTION(double)
    DEPENDENCY(bbn_lifetime, std::vector<double>)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY R_pi
  START_CAPABILITY
    #define FUNCTION SN_R_pi
    START_FUNCTION(double)
    DEPENDENCY(SMINPUTS, SMInputs)
    DEPENDENCY(Theta_sq, Eigen::Matrix3d)
    ALLOW_MODELS(SN_dev)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY R_K
  START_CAPABILITY
    #define FUNCTION SN_R_K
    START_FUNCTION(double)
    DEPENDENCY(SMINPUTS, SMInputs)
    DEPENDENCY(Theta_sq, Eigen::Matrix3d)
    ALLOW_MODEL(SN_dev)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY R_tau
  START_CAPABILITY
    #define FUNCTION SN_R_tau
    START_FUNCTION(double)
    DEPENDENCY(SMINPUTS, SMInputs)
    DEPENDENCY(Theta_sq, Eigen::Matrix3d)
    ALLOW_MODEL(SN_dev)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY lnLlepuniv
  START_CAPABILITY
    #define FUNCTION lnL_lepuniv
    START_FUNCTION(double)
    DEPENDENCY(R_pi, double)
    DEPENDENCY(R_K, double)
    DEPENDENCY(R_tau, double)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY m_GERDA
  START_CAPABILITY
    #define FUNCTION SN_m_GERDA
    START_FUNCTION(double)
    DEPENDENCY(m_nu, Eigen::Matrix3cd)
    DEPENDENCY(UPMNS, Eigen::Matrix3cd)
    DEPENDENCY(SeesawI_Theta, Eigen::Matrix3cd)
    ALLOW_MODEL(SN_dev)
    #undef FUNCTION
  #undef CAPABILITY 

  #define CAPABILITY m_Kam
  START_CAPABILITY
    #define FUNCTION SN_m_Kam
    START_FUNCTION(double)
    DEPENDENCY(m_nu, Eigen::Matrix3cd)
    DEPENDENCY(UPMNS, Eigen::Matrix3cd)
    DEPENDENCY(SeesawI_Theta, Eigen::Matrix3cd)
    ALLOW_MODEL(SN_dev)
    #undef FUNCTION
  #undef CAPABILITY 

  #define CAPABILITY lnL0nubb
  START_CAPABILITY
    #define FUNCTION lnL_0nubb
    START_FUNCTION(double)
    DEPENDENCY(m_GERDA, double)
    DEPENDENCY(m_Kam, double)
   #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY ckm_V_ud 
  START_CAPABILITY
    #define FUNCTION SN_ckm_V_ud
    START_FUNCTION(double)
    ALLOW_MODEL(SN_dev)
    DEPENDENCY(SMINPUTS, SMInputs)
    DEPENDENCY(Theta_sq, Eigen::Matrix3d)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY lnLckm
  START_CAPABILITY
    #define FUNCTION lnL_ckm
    START_FUNCTION(double)
    ALLOW_MODEL(SN_dev)
    DEPENDENCY(SMINPUTS, SMInputs)
    DEPENDENCY(Theta_sq, Eigen::Matrix3d)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY lnLpienu
  START_CAPABILITY
    #define FUNCTION lnL_pienu
    START_FUNCTION(double)
    ALLOW_MODEL(SN_dev)
    DEPENDENCY(Ue1, double)
    DEPENDENCY(Ue2, double)
    DEPENDENCY(Ue3, double)
    #undef FUNCTION
  #undef CAPABILITY

  #define CAPABILITY lnLps191e
  START_CAPABILITY

    #define FUNCTION lnL_ps191_e
    START_FUNCTION(double)
    ALLOW_MODEL(SN_dev)
    DEPENDENCY(Ue1, double)
    DEPENDENCY(Ue2, double)
    DEPENDENCY(Ue3, double)
    DEPENDENCY(Um1, double)
    DEPENDENCY(Um2, double)
    DEPENDENCY(Um3, double)
    DEPENDENCY(Ut1, double)
    DEPENDENCY(Ut2, double)
    DEPENDENCY(Ut3, double)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY lnLps191mu
  START_CAPABILITY

    #define FUNCTION lnL_ps191_mu
    START_FUNCTION(double)
    ALLOW_MODEL(SN_dev)
    DEPENDENCY(Ue1, double)
    DEPENDENCY(Ue2, double)
    DEPENDENCY(Ue3, double)
    DEPENDENCY(Um1, double)
    DEPENDENCY(Um2, double)
    DEPENDENCY(Um3, double)
    DEPENDENCY(Ut1, double)
    DEPENDENCY(Ut2, double)
    DEPENDENCY(Ut3, double)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY lnLcharme
  START_CAPABILITY

    #define FUNCTION lnL_charm_e
    START_FUNCTION(double)
    ALLOW_MODEL(SN_dev)
    DEPENDENCY(Ue1, double)
    DEPENDENCY(Ue2, double)
    DEPENDENCY(Ue3, double)
    DEPENDENCY(Um1, double)
    DEPENDENCY(Um2, double)
    DEPENDENCY(Um3, double)
    DEPENDENCY(Ut1, double)
    DEPENDENCY(Ut2, double)
    DEPENDENCY(Ut3, double)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY lnLcharmmu
  START_CAPABILITY

    #define FUNCTION lnL_charm_mu
    START_FUNCTION(double)
    ALLOW_MODEL(SN_dev)
    DEPENDENCY(Ue1, double)
    DEPENDENCY(Ue2, double)
    DEPENDENCY(Ue3, double)
    DEPENDENCY(Um1, double)
    DEPENDENCY(Um2, double)
    DEPENDENCY(Um3, double)
    DEPENDENCY(Ut1, double)
    DEPENDENCY(Ut2, double)
    DEPENDENCY(Ut3, double)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY lnLdelphi
  START_CAPABILITY

    #define FUNCTION lnL_delphi
    START_FUNCTION(double)
    ALLOW_MODEL(SN_dev)
    DEPENDENCY(Ue1, double)
    DEPENDENCY(Ue2, double)
    DEPENDENCY(Ue3, double)
    DEPENDENCY(Um1, double)
    DEPENDENCY(Um2, double)
    DEPENDENCY(Um3, double)
    DEPENDENCY(Ut1, double)
    DEPENDENCY(Ut2, double)
    DEPENDENCY(Ut3, double)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY lnLatlase
  START_CAPABILITY

    #define FUNCTION lnL_atlas_e
    START_FUNCTION(double)
    ALLOW_MODEL(SN_dev)
    DEPENDENCY(Ue1, double)
    DEPENDENCY(Ue2, double)
    DEPENDENCY(Ue3, double)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY lnLatlasmu
  START_CAPABILITY

    #define FUNCTION lnL_atlas_mu
    START_FUNCTION(double)
    ALLOW_MODEL(SN_dev)
    DEPENDENCY(Um1, double)
    DEPENDENCY(Um2, double)
    DEPENDENCY(Um3, double)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY lnLe949
  START_CAPABILITY

    #define FUNCTION lnL_e949
    START_FUNCTION(double)
    ALLOW_MODEL(SN_dev)
    DEPENDENCY(Um1, double)
    DEPENDENCY(Um2, double)
    DEPENDENCY(Um3, double)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY lnLnutev
  START_CAPABILITY

    #define FUNCTION lnL_nutev
    START_FUNCTION(double)
    ALLOW_MODEL(SN_dev)
    DEPENDENCY(Um1, double)
    DEPENDENCY(Um2, double)
    DEPENDENCY(Um3, double)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY lnLtau
  START_CAPABILITY

    #define FUNCTION lnL_tau
    START_FUNCTION(double)
    ALLOW_MODEL(SN_dev)
    DEPENDENCY(Ut1, double)
    DEPENDENCY(Ut2, double)
    DEPENDENCY(Ut3, double)
    #undef FUNCTION

  #undef CAPABILITY

  #define CAPABILITY U_ps191e
  START_CAPABILITY

    #define FUNCTION printable_ps191e
    START_FUNCTION(double)
    ALLOW_MODEL(SN_dev)
    DEPENDENCY(Ue1, double)
    DEPENDENCY(Ue2, double)
//    DEPENDENCY(Ue3, double)
    DEPENDENCY(Um1, double)
    DEPENDENCY(Um2, double)
//    DEPENDENCY(Um3, double)
    DEPENDENCY(Ut1, double)
    DEPENDENCY(Ut2, double)
//    DEPENDENCY(Ut3, double)
    #undef FUNCTION

  #undef CAPABILITY



#undef MODULE


#endif /* defined(__NeutrinoBit_rollcall_hpp__) */


