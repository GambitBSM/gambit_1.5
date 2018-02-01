//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Sterile neutrino scan; using Casas-Ibarra parameterization
///
///  *********************************************
/// 
///  Authors
/// 
///  \author Suraj Krishnamurthy
///          (S.Krishnamurthy@uva.nl) 
///  \date 2017 February
///
///  \author Tomas Gonzalo
///          (t.e.gonzalo@fys.uio.no)
///  \date 2017 Oct
///
///  *********************************************

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/NeutrinoBit/NeutrinoBit_rollcall.hpp"
#include "gambit/Utils/numerical_constants.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include "gambit/NeutrinoBit/spline.h"
#include "gambit/Utils/statistics.hpp"

using namespace Eigen;

namespace Gambit
{
  namespace NeutrinoBit
  {

    // BBN constraint: lifetime must be less than 0.1s [arXiv:1202.2841] 
    void SN_bbn_lifetime(std::vector<double>& result_lifetime)
    {
      using namespace Pipes::SN_bbn_lifetime;
      SMInputs sminputs = *Dep::SMINPUTS;
      static double conv_fact = 6.58e-16;  // conversion factor from ev^-1 to s
      static double G_F_sq = pow(sminputs.GF, 2.0);  // GeV^-4
      // TODO (CW): Should come from SM input file
      static double g_L_twid_sq = 0.0771;  // g_L_twid = -0.5 + s_W_sq
      static double g_R_sq = 0.0494;  // g_R = s_W^2
      static double g_L_sq = 0.5217;  // g_L = 0.5 + s_W^2
      //double temp_bbn;
      std::vector<double> lifetime(3), M(3);
      M[0] = *Param["M_1"];
      M[1] = *Param["M_2"];
      M[2] = *Param["M_3"];
      Matrix3d Usq = Dep::SeesawI_Theta->cwiseAbs2(); // |\Theta_{ij}|^2

      for (int i=0; i<3; i++)
      {
        lifetime[i] = (96*pow(pi,3.0)*1e-9*conv_fact) / (G_F_sq*pow(M[i],5.0))*( ((1 + g_L_twid_sq + g_R_sq)*(Usq(1,i) + Usq(2,i))) + ((1 + g_L_sq + g_R_sq)*Usq(0,i)) );
      }
      result_lifetime = lifetime;
    }
 
    // BBN constraint likelihood : lifetime must be less than 0.1s
    // [arXiv:1202.2841] Since the limit is approximate, we simply implement it
    // as a hard ~5 sigma cut.
    void lnL_bbn(double& result_bbn)
    {
      using namespace Pipes::lnL_bbn;
      std::vector<double> lifetime = *Dep::bbn_lifetime;
      result_bbn = 0.0;
      for(int i=0; i<3; i++)
      {
        if(lifetime[i]>0.1)
        {
          result_bbn = -12.5;
          break;
        }
      }
    }
   
    // Lepton universality constraint: R_(e,mu)_pi/R_(e,mu)_K should be within experimental limits [R_pi_SM, R_K_SM: Phys. Rev. Lett 99, 231801; R_tau_SM: Int. J. Mod. Phys. A 24, 715, 2009; R_pi experimental limits: Phys. Rev. Lett. 70, 17; iR_K experimental limits (NA62): Phys. Lett. B 719 (2013), 326; R_tau experimental limits: Phys. Rev. D 86, 010001]
    void SN_R_pi(double& R_pi)
    {
      using namespace Pipes::SN_R_pi;
      SMInputs sminputs = *Dep::SMINPUTS;
      static double m_pi = meson_masses.pi_plus;
      static double R_pi_SM = 1.2354e-4;
      static double r_e_pi = pow(sminputs.mE,2)/pow(m_pi,2);
      static double r_mu_pi = pow(sminputs.mMu,2)/pow(m_pi,2);
      double e_f_pi, mu_f_pi, d_r_pi;
      std::vector<double> M(3), r_I_pi(3), G_e_pi(3), G_mu_pi(3), e_fac_pi(3), mu_fac_pi(3);
      Matrix3d Usq = Dep::SeesawI_Theta->cwiseAbs2();

      e_f_pi = 0.0;
      mu_f_pi = 0.0;
      M[0] = *Param["M_1"];
      M[1] = *Param["M_2"];
      M[2] = *Param["M_3"];

      for (int i=0; i<3; i++)
      {
        e_fac_pi[i] = 0;
        mu_fac_pi[i] = 0;
        r_I_pi[i] = pow(M[i], 2.0)/pow(m_pi, 2.0);
 
        if(M[i] + sminputs.mMu < m_pi)
        {
          G_mu_pi[i] = (r_mu_pi + r_I_pi[i] - pow((r_mu_pi - r_I_pi[i]), 2.0) * sqrt(1.0 - 2.0*(r_mu_pi + r_I_pi[i]) + pow(r_mu_pi - r_I_pi[i], 2.0))) / (r_mu_pi * pow((1.0 - r_mu_pi), 2.0));
          mu_fac_pi[i] = Usq(1,i) * (G_mu_pi[i] - 1.0);
        } 
        if(M[i] + sminputs.mE < m_pi)
        {
          G_e_pi[i] = (r_e_pi + r_I_pi[i] - pow((r_e_pi - r_I_pi[i]), 2.0) * sqrt(1.0 - 2.0*(r_e_pi + r_I_pi[i]) + pow((r_e_pi - r_I_pi[i]), 2.0))) / (r_e_pi * pow((1.0 - r_e_pi), 2.0));
          e_fac_pi[i] = Usq(0,i) * (G_e_pi[i] - 1.0); 
        }
        e_f_pi += e_fac_pi[i];
        mu_f_pi += mu_fac_pi[i];
      }
      d_r_pi = ((1.0 + e_f_pi)/(1.0 + mu_f_pi)) - 1.0;
      R_pi = R_pi_SM * (1.0 + d_r_pi);
 
    }

    void SN_R_K(double& R_K)
    {
      using namespace Pipes::SN_R_K;
      SMInputs sminputs = *Dep::SMINPUTS;
      static double m_pi = meson_masses.pi_plus;
      static double m_K = meson_masses.kaon_plus; 
      static double R_K_SM = 2.477e-5;
      static double r_e_K = pow(sminputs.mE,2)/pow(m_K,2);
      static double r_mu_K = pow(sminputs.mMu,2)/pow(m_K,2);
      double e_f_K, mu_f_K, d_r_K;
      std::vector<double> M(3), r_I_K(3), G_e_K(3), G_mu_K(3), e_fac_K(3), mu_fac_K(3);
      Matrix3d Usq = Dep::SeesawI_Theta->cwiseAbs2();

      e_f_K = 0.0;
      mu_f_K = 0.0;
      M[0] = *Param["M_1"];
      M[1] = *Param["M_2"];
      M[2] = *Param["M_3"];

      for (int i=0; i<3; i++)
      {
        e_fac_K[i] = 0;
        mu_fac_K[i] = 0;
        r_I_K[i] = pow(M[i], 2.0)/pow(m_K,2.0);

        if(M[i] + sminputs.mMu < m_K and M[i] + sminputs.mMu > m_pi)
        {
          G_mu_K[i] = (r_mu_K + r_I_K[i] - pow((r_mu_K - r_I_K[i]), 2.0) * sqrt(1.0 - 2.0*(r_mu_K + r_I_K[i]) + pow(r_mu_K - r_I_K[i], 2.0))) / (r_mu_K * pow((1.0 - r_mu_K), 2.0));
          mu_fac_K[i] = Usq(1,i) * (G_mu_K[i] - 1.0);
        } 
        if(M[i] + sminputs.mE < m_K and M[i] + sminputs.mE > m_pi)
        {
          G_e_K[i] = (r_e_K + r_I_K[i] - pow((r_e_K - r_I_K[i]), 2.0) * sqrt(1.0 - 2.0*(r_e_K + r_I_K[i]) + pow((r_e_K - r_I_K[i]), 2.0))) / (r_e_K * pow((1.0 - r_e_K), 2.0));
          e_fac_K[i] = Usq(0,i) * (G_e_K[i] - 1.0);
        }
          
        e_f_K += e_fac_K[i];
        mu_f_K += mu_fac_K[i];
      }
      d_r_K = ((1.0 + e_f_K)/(1.0 + mu_f_K)) - 1.0;
      R_K = R_K_SM * (1.0 + d_r_K);
    }

    void SN_R_tau(double& R_tau)
    {
      using namespace Pipes::SN_R_tau;
      SMInputs sminputs = *Dep::SMINPUTS;
      static double m_tau = sminputs.mTau;  // GeV
      static double R_tau_SM = 0.973;
      double e_f_tau, mu_f_tau, d_r_tau;
      std::vector<double> M(3), e_fac_tau(3), mu_fac_tau(3);
      Matrix3d Usq = Dep::SeesawI_Theta->cwiseAbs2();

      e_f_tau = 0.0;
      mu_f_tau = 0.0;
      M[0] = *Param["M_1"];
      M[1] = *Param["M_2"];
      M[2] = *Param["M_3"];

      for (int i=0; i<3; i++)
      {
        if(M[i] + sminputs.mE > m_tau)
          e_fac_tau[i] = Usq(0,i);
        if(M[i] + sminputs.mMu > m_tau)
          mu_fac_tau[i] = Usq(1,i);

        e_f_tau += e_fac_tau[i];
        mu_f_tau += mu_fac_tau[i];
      }
      d_r_tau = ((1.0 - e_f_tau)/(1.0 - mu_f_tau)) - 1.0;
      R_tau = R_tau_SM * (1.0 + d_r_tau);
    }

    void lnL_lepuniv(double& result_lepuniv)
    {
      using namespace Pipes::lnL_lepuniv;
      double R_pi = *Dep::R_pi;
      double R_K = *Dep::R_K;
      double R_tau = *Dep::R_tau;

      // TODO: change to 1sigma
      double R_pi_exp = 1.23e-4;
      double R_pi_err = 0.012e-4;
      double R_K_exp = 2.488e-5;
      double R_K_err = 0.03e-5;
      double R_tau_exp = 0.9764;
      double R_tau_err = 0.009;

      result_lepuniv = 0.0;
      result_lepuniv += Stats::gaussian_loglikelihood(R_pi, R_pi_exp, 0.0, R_pi_err, false);
      result_lepuniv += Stats::gaussian_loglikelihood(R_K, R_K_exp, 0.0, R_K_err, false);
      result_lepuniv += Stats::gaussian_loglikelihood(R_tau, R_tau_exp, 0.0, R_tau_err, false);
    }

    // Neutrinoless double-beta decay constraint: m_bb should be less than the experimentally determined limits in GERDA and KamLAND-Zen [GERDA: Phys. Rev. Lett. 111 (2013) 122503; KamLAND-Zen: Phys. Rev. Lett 117 (2016) 082503]
    void SN_m_GERDA(double &m_GERDA)
    {
      using namespace Pipes::SN_m_GERDA;
      std::vector<double> M(3);
      std::complex<double> m_temp_GERDA = {0.0,0.0};

      Matrix3cd m_light = *Dep::m_nu;
      Matrix3cd U_light = *Dep::UPMNS;
      Matrix3cd theta = *Dep::SeesawI_Theta;

      m_GERDA = 0.0;
      M[0] = *Param["M_1"];
      M[1] = *Param["M_2"];
      M[2] = *Param["M_3"];

      for (int i=0; i<3; i++)
        m_temp_GERDA += pow(U_light(0,i),2)*m_light(i,i) + pow(theta(0,i),2)*M[i]*(pow(*Param["L_Ge"], 2.0)/(pow(*Param["L_Ge"], 2.0)+pow(M[i], 2.0)));

      m_GERDA = abs(m_temp_GERDA);
    }

    // Calculate 0nubb decay rate [1/s] for 136Xe 0nubb detector, for sterile
    // neutrino model
    void Gamma_0nubb_Xe_SN(double& result)
    {
      using namespace Pipes::Gamma_0nubb_Xe_SN;
      double mp, Gamma_0nu, A_0nubb_Xe, p2_0nubb_Xe, prefactor;
      std::vector<double> M(3);
      std::complex<double> sum = {0.0,0.0};

      // Relevant model parameters
      Matrix3cd m_light = *Dep::m_nu;
      Matrix3cd U_light = *Dep::UPMNS;
      Matrix3cd theta = *Dep::SeesawI_Theta;
      M[0] = *Param["M_1"];
      M[1] = *Param["M_2"];
      M[2] = *Param["M_3"];

      // NOTE: For the time being, we retreive nuisance parameters as yaml file options for the
      // A_0nubb_Xe = *Param["A_0nubb_Xe"];  // Range: 4.41 - 19.7 [1e-10 1/yr]
      // p2_0nubb_Xe = pow(*Param["p_0nubb_Xe"], 2.0);  // Range: 178.0 - 211.0 [MeV]

      // Nuisance parameters following the definitions in Faessler et al. 2014 (1408.6077)
      A_0nubb_Xe = runOptions->getValueOrDef<double>(8.74, "A");
      p2_0nubb_Xe = pow(runOptions->getValueOrDef<double>(183.0, "p"), 2.0);
      p2_0nubb_Xe *= 1e-3;  // MeV --> GeV
      mp = 0.938;  // [GeV] (PDG 2014)

      // Lifetime equation is adopted from Faessler+14, Eq. (13)
      prefactor = A_0nubb_Xe*mp*mp/p2_0nubb_Xe/p2_0nubb_Xe;
        for (int i=0; i<3; i++)
        {
          sum += pow(U_light(0,i),2)*m_light(i,i) + pow(theta(0,i),2)*M[i]*p2_0nubb_Xe/(p2_0nubb_Xe+pow(M[i], 2.0));
        }
      result = prefactor * abs(sum) * abs(sum);
    }

    void lnL_0nubb_KamLAND_Zen(double& result)
    {
      using namespace Pipes::lnL_0nubb_KamLAND_Zen;
      double tau_limit = 1.07e26*3.156e7;  // [s] 90% CL

      double Gamma = *Dep::Gamma_0nubb_Xe;

      // Factor 1.28155 corresponds to one-sided UL at 90% CL
      result = Stats::gaussian_loglikelihood(Gamma, 0., 0., 1./tau_limit/1.28155, false);
    }

    void SN_m_Kam(double& m_Kam)
    {
      using namespace Pipes::SN_m_Kam;
      std::vector<double> M(3);
      std::complex<double> m_temp_Kam = {0.0,0.0};

      Matrix3cd m_light = *Dep::m_nu;
      Matrix3cd U_light = *Dep::UPMNS;
      Matrix3cd theta = *Dep::SeesawI_Theta;

      m_Kam = 0.0;
      M[0] = *Param["M_1"];
      M[1] = *Param["M_2"];
      M[2] = *Param["M_3"];

       for (int i=0; i<3; i++)
        m_temp_Kam += pow(U_light(0,i),2)*m_light(i,i) + pow(theta(0,i),2)*M[i]*(pow(*Param["L_Xe"], 2.0)/(pow(*Param["L_Xe"], 2.0)+pow(M[i], 2.0)));

      m_Kam = abs(m_temp_Kam);
    }

    void lnL_0nubb(double& result_0nubb)
    {
      using namespace Pipes::lnL_0nubb;
      static double m_bb_GERDA = 4e-10;  // GeV
      static double m_bb_Kam = 1.65e-10;  // GeV
      double m_GERDA = *Dep::m_GERDA;
      double m_Kam = *Dep::m_Kam;

      result_0nubb += 
        Stats::gaussian_loglikelihood(m_GERDA, 0., 0., m_bb_GERDA, false)+
        Stats::gaussian_loglikelihood(m_Kam, 0., 0., m_bb_Kam, false);

      /*
      if ((m_GERDA < m_bb_GERDA) && (m_Kam < m_bb_Kam))
      {
        result_0nubb = 0.0;

      }
      else
      {
        result_0nubb = -1.0E+05;
      }
      */
    }

    // CKM unitarity constraint: V_ud should lie within 3sigma of the world average [PDG 2016]
    void SN_ckm_V_ud(double &V_ud)
    {
      /*
      using namespace Pipes::SN_ckm_V_ud;
      SMInputs sminputs = *Dep::SMINPUTS;
      Matrix3cd Theta = *Dep::SeesawI_Theta;
      static double G_F_sq = pow(sminputs.GF, 2.0);  // GeV^-4
      static double V_us_exp[7] = {0.2235,0.2227,0.2244,0.2238,0.2242,0.2262,0.2214};
      static double err[7] = {0.0006,0.0013,0.0008,0.0006,0.0011,0.0013,0.0022};
      static double G_mu_sq = 1.3605e-10;  // GeV^-4
      static double V_ud_0 = 0.97427;
      static double err_0 = 0.00015;
      static double err_factor = 5.3779e7;
      double V_ud_exp[7], f[7];
      double V_ud_sq;

      Matrix3d ThetaNorm = (Theta * Theta.adjoint()).real();

      for (int i=0;i<7;i++)
      {
        V_ud_exp[i] = sqrt(1 - pow(V_us_exp[i], 2.0));
      }
      f[0] = (G_F_sq/G_mu_sq)*(1 - ThetaNorm(0,0));
      f[3] = (G_F_sq/G_mu_sq)*(1 - ThetaNorm(1,1));
      f[5] = 1 + ThetaNorm(1,1);
      f[6] = 1 + ThetaNorm(0,0) + ThetaNorm(1,1) + ThetaNorm(2,2);
      f[1] = f[0];
      f[2] = f[0];
      f[4] = f[3];
      V_ud_sq = 0.0;
      for (int j=0; j<7; j++)
      {
        V_ud_sq += pow(V_ud_exp[j], 2.0)/(pow(err[j], 2.0)*f[j]);
      }
      V_ud_sq += pow(V_ud_0, 2.0)/(pow(err_0, 2.0)*(1 + ThetaNorm(0,0)));
      V_ud_sq /= err_factor;

      V_ud = sqrt(V_ud_sq);*/
    }

    void lnL_ckm(double& result_ckm)
    {
      using namespace Pipes::lnL_ckm;
      SMInputs sminputs = *Dep::SMINPUTS;
      Matrix3cd Theta = *Dep::SeesawI_Theta;
      double G_mu = *Dep::Gmu;
      double V_us = *Param["CKM_lambda"];
 
      // Experimental values determined for K and tau decays. From table 1 in 1502.00477
      double V_us_exp[] = {0.2163, 0.2166, 0.2155, 0.2160, 0.2158, 0.2262, 0.2214, 0.2173};
      double err_V_us_exp[] = {0.0006, 0.0006, 0.0013, 0.0011, 0.0014, 0.0013, 0.0022, 0.0022};
      double f_plus = 0.959;
      double err_f_plus = 0.005;
      for(int i=0; i<5; i++)
      {
        V_us_exp[i] /= f_plus;
        err_V_us_exp[i] = sqrt(pow(err_V_us_exp[i] / f_plus,2) + pow(V_us_exp[i] * err_f_plus / f_plus, 2));
      }  

      // Superallowed beta decays and more, from 1509.0474
      // TODO: Wait for Marcin to check these out, only include pion beta decays for now
      //static double V_ud_exp[] = {0.97417, 0.9754, 0.9734, 0.9718, 0.9749};
      //static double err_V_ud_exp[] = {0.00021, 0.0014, 0.0027, 0.0017, 0.0026};
      static double V_ud_exp = 0.9749;
      static double err_V_ud_exp = 0.0026;

      double f[8];
      Matrix3d ThetaNorm = (Theta * Theta.adjoint()).real();

      f[0] = pow(sminputs.GF/G_mu,2)*(1 - ThetaNorm(0,0));
      f[1] = f[0];
      f[2] = f[0];
      f[3] = pow(sminputs.GF/G_mu,2)*(1 - ThetaNorm(1,1));
      f[4] = f[3];
      f[5] = 1 + ThetaNorm(1,1);
      f[6] = 1 + ThetaNorm(0,0) + ThetaNorm(1,1) - ThetaNorm(2,2);
      f[7] = 1 + 0.2*ThetaNorm(0,0) - 0.9*ThetaNorm(1,1) - 0.2*ThetaNorm(2,2);

      double chi2 = 0.0;
      for (int i=0; i<7; i++)
        chi2 += pow( (sqrt(pow(V_us,2)*f[i]) - V_us_exp[i]) / err_V_us_exp[i], 2.0);

      // According to 1407.6607 the correction for Vud is the same as K->pi e nu (f[0])
      chi2 += pow( (sqrt((1 - pow(V_us,2))*f[0]) - V_ud_exp)/ err_V_ud_exp, 2.0);

      result_ckm = -0.5*chi2;
    }

    // Likelihood contribution from PIENU; searched for extra peaks in the spectrum of pi -> mu + nu. Constrains |U_ei|^2 in the mass range 60-129 MeV. [Phys. Rev. D, 84(5), 2011]
    void lnL_pienu(double& result_pienu)
    {
      using namespace Pipes::lnL_pienu;
      static bool read_table_pienu = true;
      static tk::spline s_pienu;
      static std::vector<double> M_temp_pienu(140), U_temp_pienu(140);
      double M_1, M_2, M_3;
      std::vector<double> U_pienu(3), mixing_sq(3);

      mixing_sq[0] = *Dep::Ue1;
      mixing_sq[1] = *Dep::Ue2;
      mixing_sq[2] = *Dep::Ue3;
      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];

      if (read_table_pienu)
      {
        double array_pienu[140][2];
        std::ifstream f_pienu("NeutrinoBit/data/pienu.csv");
        for (int row=0; row<140; ++row)
        {
          std::string line_pienu;
          getline(f_pienu, line_pienu);
          if (!f_pienu.good())
            break;
          std::stringstream iss_pienu(line_pienu);
          for (int col=0; col<2; ++col)
          {
            std::string val_pienu;
            getline(iss_pienu, val_pienu, ',');
            if (!iss_pienu)
              break;
            std::stringstream conv_pienu(val_pienu);
            conv_pienu >> array_pienu[row][col];
          }
        }

        for (int i=0; i<140; i++)
        {
          M_temp_pienu[i] = array_pienu[i][0];
          U_temp_pienu[i] = array_pienu[i][1];
        }
        s_pienu.set_points(M_temp_pienu, U_temp_pienu);
        read_table_pienu = false;
      }

      U_pienu[0] = s_pienu(M_1);
      U_pienu[1] = s_pienu(M_2);
      U_pienu[2] = s_pienu(M_3);

      //TODO: Replace with Gaussian likelihood
      result_pienu = -2.44*((mixing_sq[0]/U_pienu[0]) + (mixing_sq[1]/U_pienu[1]) + (mixing_sq[2]/U_pienu[2]));
    }

    // Likelihood contribution from PS191, electron sector; looked for charged tracks originating from RHN decays: nu_r -> l(-) + l(+) + nu / l + pi / e + pi(+) + pi(0). Constrains |U_ei|^2 in the mass range 20-450 MeV. Function also incorporates a later re-interpretation of the data to account for neutral current interaction (ignored in original) as well as the RHNs' Majorana nature. [Original: Phys. Lett. B, 203(3):332-334, 1988][Re-interp.: JHEP, 2012(6):1-27]
    void lnL_ps191_e(double& result_ps191e)
    {
      using namespace Pipes::lnL_ps191_e;
      static bool read_table_ps191e = true;
      static tk::spline s_ps191e;
      double M_1, M_2, M_3;
      static std::vector<double> M_temp_ps191e(116), U_temp_ps191e(116);
      std::vector<double> U_ps191e(3), mixing_sq_ps191e(3);
//      double mixing_sq;
      double c_e = 0.5711;
      double c_mu = 0.1265;
      double c_tau = 0.1265;

      mixing_sq_ps191e[0] = *Dep::Ue1 * ((c_e * *Dep::Ue1) + (c_mu * *Dep::Um1) + (c_tau * *Dep::Ut1));
      mixing_sq_ps191e[1] = *Dep::Ue2 * ((c_e * *Dep::Ue2) + (c_mu * *Dep::Um2) + (c_tau * *Dep::Ut2));
      mixing_sq_ps191e[2] = *Dep::Ue3 * ((c_e * *Dep::Ue3) + (c_mu * *Dep::Um3) + (c_tau * *Dep::Ut3));
//      mixing_sq = mixing_sq_ps191e[0] + mixing_sq_ps191e[1];

     if (read_table_ps191e)
      {
        double array_ps191e[116][2];
        std::ifstream f_ps191e("NeutrinoBit/data/ps191_e.csv");
        for (int row=0; row<116; ++row)
        {
          std::string line_ps191e;
          getline(f_ps191e, line_ps191e);
          if (!f_ps191e.good())
            break;
          std::stringstream iss_ps191e(line_ps191e);
          for (int col=0; col<2; ++col)
          {
            std::string val_ps191e;
            getline(iss_ps191e, val_ps191e, ',');
            if (!iss_ps191e)
              break;
            std::stringstream conv_ps191e(val_ps191e);
            conv_ps191e >> array_ps191e[row][col];
          }
        }

        for (int i=0; i<116; i++)
        {
          M_temp_ps191e[i] = array_ps191e[i][0];
          U_temp_ps191e[i] = array_ps191e[i][1];
        }
        s_ps191e.set_points(M_temp_ps191e, U_temp_ps191e);
        read_table_ps191e = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U_ps191e[0] = s_ps191e(M_1);
      U_ps191e[1] = s_ps191e(M_2);
      U_ps191e[2] = s_ps191e(M_3);
//      result_ps191e = -2.44*((mixing_sq_ps191e[0]/pow(U_ps191e[0], 2.0)) + (mixing_sq_ps191e[1]/pow(U_ps191e[1], 2.0)) + (mixing_sq_ps191e[2]/pow(U_ps191e[2], 2.0)));
//      result_ps191e = -2.44*(mixing_sq/pow(U_ps191e[0], 2.0));
      result_ps191e = -2.44*((sqrt(mixing_sq_ps191e[0])/U_ps191e[0]) + (sqrt(mixing_sq_ps191e[1])/U_ps191e[1]) + (sqrt(mixing_sq_ps191e[2])/U_ps191e[2]));
    }

    // Likelihood contribution from PS191, muon sector. Constrains |U_(mu,i)|^2 in the mass range 20-450 MeV. Description & references above.
    void lnL_ps191_mu(double& result_ps191mu)
    {
      using namespace Pipes::lnL_ps191_mu;
      static bool read_table_ps191mu = true;
      static tk::spline s_ps191mu;
      double M_1, M_2, M_3;
      static std::vector<double> M_temp_ps191mu(102), U_temp_ps191mu(102);
      std::vector<double> U_ps191mu(3), mixing_sq_ps191mu(3);
      double c_e = 0.5711;
      double c_mu = 0.1265;
      double c_tau = 0.1265;

      mixing_sq_ps191mu[0] = *Dep::Um1 * ((c_e * *Dep::Ue1) + (c_mu * *Dep::Um1) + (c_tau * *Dep::Ut1));
      mixing_sq_ps191mu[1] = *Dep::Um2 * ((c_e * *Dep::Ue2) + (c_mu * *Dep::Um2) + (c_tau * *Dep::Ut2));
      mixing_sq_ps191mu[2] = *Dep::Um3 * ((c_e * *Dep::Ue3) + (c_mu * *Dep::Um3) + (c_tau * *Dep::Ut3));

      if (read_table_ps191mu)
      {
        double array_ps191mu[102][2];
        std::ifstream f_ps191mu("NeutrinoBit/data/ps191_mu.csv");
        for (int row=0; row<102; ++row)
        {
          std::string line_ps191mu;
          getline(f_ps191mu, line_ps191mu);
          if (!f_ps191mu.good())
            break;
          std::stringstream iss_ps191mu(line_ps191mu);
          for (int col=0; col<2; ++col)
          {
            std::string val_ps191mu;
            getline(iss_ps191mu, val_ps191mu, ',');
            if (!iss_ps191mu)
              break;
            std::stringstream conv_ps191mu(val_ps191mu);
            conv_ps191mu >> array_ps191mu[row][col];
          }
        }  

        for (int i=0; i<102; i++)
        {
          M_temp_ps191mu[i] = array_ps191mu[i][0];
          U_temp_ps191mu[i] = array_ps191mu[i][1];
        }
        s_ps191mu.set_points(M_temp_ps191mu, U_temp_ps191mu);
        read_table_ps191mu = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U_ps191mu[0] = s_ps191mu(M_1);
      U_ps191mu[1] = s_ps191mu(M_2);
      U_ps191mu[2] = s_ps191mu(M_3);
      result_ps191mu = -2.44*((mixing_sq_ps191mu[0]/pow(U_ps191mu[0], 2.0)) + (mixing_sq_ps191mu[1]/pow(U_ps191mu[1], 2.0)) + (mixing_sq_ps191mu[2]/pow(U_ps191mu[2], 2.0)));
    }

    // Likelihood contribution from CHARM, electron sector; searched for charged and neutral current decays of RHNs. Constrains |U_ei|^2 in the mass range 0.5-2.8 GeV. [Phys. Lett. B, 166(4):473-478, 1986]
    void lnL_charm_e(double& result_charme)
    {
      using namespace Pipes::lnL_charm_e;
      static bool read_table_charme = true;
      static tk::spline s_charme;
      double M_1, M_2, M_3;
      static std::vector<double> M_temp_charme(56), U_temp_charme(56);
      std::vector<double> U_charme(3), mixing_sq(3);
      double c_e = 0.5711;
      double c_mu = 0.1265;
      double c_tau = 0.1265;

      mixing_sq[0] = *Dep::Ue1 * ((c_e * *Dep::Ue1) + (c_mu * *Dep::Um1) + (c_tau * *Dep::Ut1));
      mixing_sq[1] = *Dep::Ue2 * ((c_e * *Dep::Ue2) + (c_mu * *Dep::Um2) + (c_tau * *Dep::Ut2));
      mixing_sq[2] = *Dep::Ue3 * ((c_e * *Dep::Ue3) + (c_mu * *Dep::Um3) + (c_tau * *Dep::Ut3));

      if (read_table_charme)
      {
        double array_charme[56][2];
        std::ifstream f_charme("NeutrinoBit/data/charm_e.csv");
        for (int row=0; row<56; ++row)
        {
          std::string line_charme;
          getline(f_charme, line_charme);
          if (!f_charme.good())
            break;
          std::stringstream iss_charme(line_charme);
          for (int col=0; col<2; ++col)
          {
            std::string val_charme;
            getline(iss_charme, val_charme, ',');
            if (!iss_charme)
              break;
            std::stringstream conv_charme(val_charme);
            conv_charme >> array_charme[row][col];
          }
        }

        for (int i=0; i<56; i++)
        {
          M_temp_charme[i] = array_charme[i][0];
          U_temp_charme[i] = array_charme[i][1];
        }
        s_charme.set_points(M_temp_charme, U_temp_charme);
        read_table_charme = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U_charme[0] = s_charme(M_1);
      U_charme[1] = s_charme(M_2);
      U_charme[2] = s_charme(M_3);
      result_charme = -2.44*((mixing_sq[0]/pow(U_charme[0], 2.0)) + (mixing_sq[1]/pow(U_charme[1], 2.0)) + (mixing_sq[2]/pow(U_charme[2], 2.0)));
    }

    // Likelihood contribution from CHARM, muon sector. Constrains |U_(mu,i)|^2 in the mass range 0.5-2.8 GeV. Description & references above.
    void lnL_charm_mu(double& result_charmmu)
    {
      using namespace Pipes::lnL_charm_mu;
      static bool read_table_charmmu = true;
      static tk::spline s_charmmu;
      double M_1, M_2, M_3;
      static std::vector<double> M_temp_charmmu(34), U_temp_charmmu(34);
      std::vector<double> U_charmmu(3), mixing_sq(3);
      double c_e = 0.5711;
      double c_mu = 0.1265;
      double c_tau = 0.1265;

      mixing_sq[0] = *Dep::Um1 * ((c_e * *Dep::Ue1) + (c_mu * *Dep::Um1) + (c_tau * *Dep::Ut1));
      mixing_sq[1] = *Dep::Um2 * ((c_e * *Dep::Ue2) + (c_mu * *Dep::Um2) + (c_tau * *Dep::Ut2));
      mixing_sq[2] = *Dep::Um3 * ((c_e * *Dep::Ue3) + (c_mu * *Dep::Um3) + (c_tau * *Dep::Ut3));

      if (read_table_charmmu)
      {
        double array_charmmu[34][2];
        std::ifstream f_charmmu("NeutrinoBit/data/charm_mu.csv");
        for (int row=0; row<34; ++row)
        {
          std::string line_charmmu;
          getline(f_charmmu, line_charmmu);
          if (!f_charmmu.good())
            break;
          std::stringstream iss_charmmu(line_charmmu);
          for (int col=0; col<2; ++col)
          {
            std::string val_charmmu;
            getline(iss_charmmu, val_charmmu, ',');
            if (!iss_charmmu)
              break;
            std::stringstream conv_charmmu(val_charmmu);
            conv_charmmu >> array_charmmu[row][col];
          }
        }

        for (int i=0; i<34; i++)
        {
          M_temp_charmmu[i] = array_charmmu[i][0];
          U_temp_charmmu[i] = array_charmmu[i][1];
        }
        s_charmmu.set_points(M_temp_charmmu, U_temp_charmmu);
        read_table_charmmu = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U_charmmu[0] = s_charmmu(M_1);
      U_charmmu[1] = s_charmmu(M_2);
      U_charmmu[2] = s_charmmu(M_3);
      result_charmmu = -2.44*((mixing_sq[0]/pow(U_charmmu[0], 2.0)) + (mixing_sq[1]/pow(U_charmmu[1], 2.0)) + (mixing_sq[2]/pow(U_charmmu[2], 2.0)));
    }

    // Likelihood contribution from DELPHI; searched for charged and neutral current decays of RHNs. Constrains |U_ei|^2, |U_(mu,i)|^2 as well as |U_(tau,i)|^2 in the mass range 3.5-50 GeV. [Z. Phys. C, 74(1):57-71, 1997]
    void lnL_delphi(double& result_delphi)
    {
      using namespace Pipes::lnL_delphi;
      static bool read_table_delphi = true;
      static tk::spline s_delphi;
      double M_1, M_2, M_3;
      static std::vector<double> M_temp_delphi(180), U_temp_delphi(180);
      std::vector<double> U_delphi(3), mixing_sq(9);

      mixing_sq[0] = *Dep::Ue1;  // This is |U_{e1}|^2 etc
      mixing_sq[1] = *Dep::Ue2;
      mixing_sq[2] = *Dep::Ue3;
      mixing_sq[3] = *Dep::Um1;
      mixing_sq[4] = *Dep::Um2;
      mixing_sq[5] = *Dep::Um3;
      mixing_sq[6] = *Dep::Ut1;
      mixing_sq[7] = *Dep::Ut2;
      mixing_sq[8] = *Dep::Ut3;

      if (read_table_delphi)
      {
        double array_delphi[180][2];
        std::ifstream f_delphi("NeutrinoBit/data/delphi.csv");
        for (int row=0; row<180; ++row)
        {
          std::string line_delphi;
          getline(f_delphi, line_delphi);
          if (!f_delphi.good())
            break;
          std::stringstream iss_delphi(line_delphi);
          for (int col=0; col<2; ++col)
          {
            std::string val_delphi;
            getline(iss_delphi, val_delphi, ',');
            if (!iss_delphi)
              break;
            std::stringstream conv_delphi(val_delphi);
            conv_delphi >> array_delphi[row][col];
          }
        }

        for (int i=0; i<180; i++)
        {
          M_temp_delphi[i] = array_delphi[i][0];
          U_temp_delphi[i] = array_delphi[i][1];
        }
        s_delphi.set_points(M_temp_delphi, U_temp_delphi);
        read_table_delphi = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U_delphi[0] = s_delphi(M_1);
      U_delphi[1] = s_delphi(M_2);
      U_delphi[2] = s_delphi(M_3);

      // Assume scaling with |U|^4, zero bkg, number of events at 95% CL is
      // reverse engineered.  We assume that lnL = mu_sig is a faithful
      // approximation to the true Poisson likelihood.
      result_delphi = -3.09*
         (pow(mixing_sq[0]/U_delphi[0], 2.0) +
          pow(mixing_sq[1]/U_delphi[1], 2.0) + pow(mixing_sq[2]/U_delphi[2], 2.0) +
          pow(mixing_sq[3]/U_delphi[0], 2.0) + pow(mixing_sq[4]/U_delphi[1], 2.0) +
          pow(mixing_sq[5]/U_delphi[2], 2.0) + pow(mixing_sq[6]/U_delphi[0], 2.0) +
          pow(mixing_sq[7]/U_delphi[1], 2.0) + pow(mixing_sq[8]/U_delphi[2], 2.0));

    }

    // Likelihood contribution from ATLAS, electron sector; looked at the production and decay chain: pp -> W*(+-) -> l(+-) + nu_r. nu_r then decays into an on-shell W and a lepton; the W decays primarily into a qq pair. Constrains |U_ei|^2 in the mass range 50-500 GeV. [JHEP, 07:162, 2015]
    void lnL_atlas_e(double& result_atlase)
    {
      using namespace Pipes::lnL_atlas_e;

      static bool read_table_atlase = true;
      static tk::spline s_atlase;
      double M_1, M_2, M_3;
      static std::vector<double> M_temp_atlase(87), U_temp_atlase(87);
      std::vector<double> U_atlase(3), mixing_sq_atlase(3);

      mixing_sq_atlase[0] = *Dep::Ue1;
      mixing_sq_atlase[1] = *Dep::Ue2;
      mixing_sq_atlase[2] = *Dep::Ue3;

      if (read_table_atlase)
      {
        double array_atlase[87][2];
        std::ifstream f_atlase("NeutrinoBit/data/atlas_e.csv");
        for (int row=0; row<87; ++row)
        {
          std::string line_atlase;
          getline(f_atlase, line_atlase);
          if (!f_atlase.good())
            break;
          std::stringstream iss_atlase(line_atlase);
          for (int col=0; col<2; ++col)
          {
            std::string val_atlase;
            getline(iss_atlase, val_atlase, ',');
            if (!iss_atlase)
              break;
            std::stringstream conv_atlase(val_atlase);
            conv_atlase >> array_atlase[row][col];
          }
        }

        for (int i=0; i<87; i++)
        {
          M_temp_atlase[i] = array_atlase[i][0];
          U_temp_atlase[i] = array_atlase[i][1];
        }
        s_atlase.set_points(M_temp_atlase, U_temp_atlase);
        read_table_atlase = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U_atlase[0] = s_atlase(M_1);
      U_atlase[1] = s_atlase(M_2);
      U_atlase[2] = s_atlase(M_3);
      result_atlase = -3.09*(pow((mixing_sq_atlase[0]/U_atlase[0]), 2.0) + pow((mixing_sq_atlase[1]/U_atlase[1]), 2.0) + pow((mixing_sq_atlase[2]/U_atlase[2]), 2.0));
    }

    // Likelihood contribution from ATLAS, muon sector. Constrains |U_(mu,i)|^2 in the mass range 50-500 GeV. Description & references above.
    void lnL_atlas_mu(double& result_atlasmu)
    {
      using namespace Pipes::lnL_atlas_mu;
      static bool read_table_atlasmu = true;
      static tk::spline s_atlasmu;
      double M_1, M_2, M_3;
      static std::vector<double> M_temp_atlasmu(87), U_temp_atlasmu(87);
      std::vector<double> U_atlasmu(3), mixing_sq_atlasmu(3);

      mixing_sq_atlasmu[0] = *Dep::Um1;
      mixing_sq_atlasmu[1] = *Dep::Um2;
      mixing_sq_atlasmu[2] = *Dep::Um3;

      if (read_table_atlasmu)
      {
        double array_atlasmu[87][2];
        std::ifstream f_atlasmu("NeutrinoBit/data/atlas_mu.csv");
        for (int row=0; row<87; ++row)
        {
          std::string line_atlasmu;
          getline(f_atlasmu, line_atlasmu);
          if (!f_atlasmu.good())
            break;
          std::stringstream iss_atlasmu(line_atlasmu);
          for (int col=0; col<2; ++col)
          {
            std::string val_atlasmu;
            getline(iss_atlasmu, val_atlasmu, ',');
            if (!iss_atlasmu)
              break;
            std::stringstream conv_atlasmu(val_atlasmu);
            conv_atlasmu >> array_atlasmu[row][col];
          }
        }

        for (int i=0; i<87; i++)
        {
          M_temp_atlasmu[i] = array_atlasmu[i][0];
          U_temp_atlasmu[i] = array_atlasmu[i][1];
        }
        s_atlasmu.set_points(M_temp_atlasmu, U_temp_atlasmu);
        read_table_atlasmu = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U_atlasmu[0] = s_atlasmu(M_1);
      U_atlasmu[1] = s_atlasmu(M_2);
      U_atlasmu[2] = s_atlasmu(M_3);
      result_atlasmu = -3.09*(pow((mixing_sq_atlasmu[0]/U_atlasmu[0]), 2.0) + pow((mixing_sq_atlasmu[1]/U_atlasmu[1]), 2.0) + pow((mixing_sq_atlasmu[2]/U_atlasmu[2]), 2.0));
    }

    // Likelihood contribution from E949; used the kaon decay: K(+) -> mu(+) + nu_r. Constrains |U_(mu,i)|^2 in the mass range 175-300 MeV. [Phys. Rev. D, 91, 2015]
    void lnL_e949(double& result_e949)
    {
      using namespace Pipes::lnL_e949;
      static bool read_table_e949 = true;
      static tk::spline s_e949;
      double M_1, M_2, M_3;
      static std::vector<double> M_temp_e949(112), U_temp_e949(112);
      std::vector<double> U_e949(3), mixing_sq_e949(3);

      mixing_sq_e949[0] = *Dep::Um1;
      mixing_sq_e949[1] = *Dep::Um2;
      mixing_sq_e949[2] = *Dep::Um3;

      if (read_table_e949)
      {
        double array_e949[112][2];
        std::ifstream f_e949("NeutrinoBit/data/e949.csv");
        for (int row=0; row<112; ++row)
        {
          std::string line_e949;
          getline(f_e949, line_e949);
          if (!f_e949.good())
            break;
          std::stringstream iss_e949(line_e949);
          for (int col=0; col<2; ++col)
          {
            std::string val_e949;
            getline(iss_e949, val_e949, ',');
            if (!iss_e949)
              break;
            std::stringstream conv_e949(val_e949);
            conv_e949 >> array_e949[row][col];
          }
        }

        for (int i=0; i<112; i++)
        {
          M_temp_e949[i] = array_e949[i][0];
          U_temp_e949[i] = array_e949[i][1];
        }
        s_e949.set_points(M_temp_e949, U_temp_e949);
        read_table_e949 = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U_e949[0] = s_e949(M_1);
      U_e949[1] = s_e949(M_2);
      U_e949[2] = s_e949(M_3);
      result_e949 = -2.44*((mixing_sq_e949[0]/U_e949[0]) + (mixing_sq_e949[1]/U_e949[1]) + (mixing_sq_e949[2]/U_e949[2]));
    }

    // Likelihood contribution from NuTeV; used RHN decays into muonic final states (mu + mu + nu / mu + e + nu / mu + pi / mu + rho). Constrains |U_(mu,i)|^2 in the mass range 0.25-2 GeV. [Phys. Rev. Lett., 83:4943-4946, 1999]
    void lnL_nutev(double& result_nutev)
    {
      using namespace Pipes::lnL_nutev;
      static bool read_table_nutev = true;
      static tk::spline s_nutev;
      double M_1, M_2, M_3;
      static std::vector<double> M_temp_nutev(249), U_temp_nutev(249);
      std::vector<double> U_nutev(3), mixing_sq(3);

      mixing_sq[0] = *Dep::Um1;
      mixing_sq[1] = *Dep::Um2;
      mixing_sq[2] = *Dep::Um3;

      if (read_table_nutev)
      {
        double array_nutev[249][2];
        std::ifstream f_nutev("NeutrinoBit/data/nutev.csv");
        for (int row=0; row<249; ++row)
        {
          std::string line_nutev;
          getline(f_nutev, line_nutev);
          if (!f_nutev.good())
            break;
          std::stringstream iss_nutev(line_nutev);
          for (int col=0; col<2; ++col)
          {
            std::string val_nutev;
            getline(iss_nutev, val_nutev, ',');
            if (!iss_nutev)
              break;
            std::stringstream conv_nutev(val_nutev);
            conv_nutev >> array_nutev[row][col];
          }
        }

        for (int i=0; i<249; i++)
        {
          M_temp_nutev[i] = array_nutev[i][0];
          U_temp_nutev[i] = array_nutev[i][1];
        }
        s_nutev.set_points(M_temp_nutev, U_temp_nutev);
        read_table_nutev = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U_nutev[0] = s_nutev(M_1);
      U_nutev[1] = s_nutev(M_2);
      U_nutev[2] = s_nutev(M_3);
      result_nutev = -2.44*(pow(mixing_sq[0]/U_nutev[0], 2.0) + pow(mixing_sq[1]/U_nutev[1], 2.0) + pow(mixing_sq[2]/U_nutev[2], 2.0));
    }

    // Likelihood contribution from a re-interpretation of CHARM data; assumes tau mixing is dominant. Constrains |U_(tau,i)|^2 in the mass range 10-290 MeV. [Phys. Lett. B, 550(1-2):8-15, 2002]
    void lnL_charm_tau(double& result_tau)
    {
      using namespace Pipes::lnL_charm_tau;
      static bool read_table_tau = true;
      static tk::spline s_tau;
      double M_1, M_2, M_3;
      static std::vector<double> M_temp_tau(172), U_temp_tau(172);
      std::vector<double> U_tau(3), mixing_sq(3);

      mixing_sq[0] = *Dep::Ut1;
      mixing_sq[1] = *Dep::Ut2;
      mixing_sq[2] = *Dep::Ut3;

      if (read_table_tau)
      {
        double array_tau[172][2];
        std::ifstream f_tau("NeutrinoBit/data/tau.csv");
        for (int row=0; row<172; ++row)
        {
          std::string line_tau;
          getline(f_tau, line_tau);
          if (!f_tau.good())
            break;
          std::stringstream iss_tau(line_tau);
          for (int col=0; col<2; ++col)
          {
            std::string val_tau;
            getline(iss_tau, val_tau, ',');
            if (!iss_tau)
              break;
            std::stringstream conv_tau(val_tau);
            conv_tau >> array_tau[row][col];
          }
        } 

        for (int i=0; i<172; i++)
        {
          M_temp_tau[i] = array_tau[i][0];
          U_temp_tau[i] = array_tau[i][1];
        }
        s_tau.set_points(M_temp_tau, U_temp_tau);
        read_table_tau = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U_tau[0] = s_tau(M_1);
      U_tau[1] = s_tau(M_2);
      U_tau[2] = s_tau(M_3);
      result_tau = -2.44*((mixing_sq[0]/U_tau[0]) + (mixing_sq[1]/U_tau[1]) + (mixing_sq[2]/U_tau[2]));
    }

    void Ue1(double& Ue1_sq)
    {
      using namespace Pipes::Ue1;
      Ue1_sq = ((*Dep::SeesawI_Theta).cwiseAbs2())(0,0);
    }

    void Um1(double& Um1_sq)
    {
      using namespace Pipes::Um1;
      Um1_sq = (Dep::SeesawI_Theta->cwiseAbs2())(1,0);
    }

    void Ut1(double& Ut1_sq)
    {
      using namespace Pipes::Ut1;
      Ut1_sq = (Dep::SeesawI_Theta->cwiseAbs2())(2,0);
    }

    void Ue2(double& Ue2_sq)
    {
      using namespace Pipes::Ue2;
      Ue2_sq = (Dep::SeesawI_Theta->cwiseAbs2())(0,1);
    }

    void Um2(double& Um2_sq)
    {
      using namespace Pipes::Um2;
      Um2_sq = (Dep::SeesawI_Theta->cwiseAbs2())(1,1);
    }

    void Ut2(double& Ut2_sq)
    {
      using namespace Pipes::Ut2;
      Ut2_sq = (Dep::SeesawI_Theta->cwiseAbs2())(2,1);
    }

    void Ue3(double& Ue3_sq)
    {
      using namespace Pipes::Ue3;
      Ue3_sq = (Dep::SeesawI_Theta->cwiseAbs2())(0,2);
    }

    void Um3(double& Um3_sq)
    {
      using namespace Pipes::Um3;
      Um3_sq = (Dep::SeesawI_Theta->cwiseAbs2())(1,2);
    }

    void Ut3(double& Ut3_sq)
    {
      using namespace Pipes::Ut3;
      Ut3_sq = (Dep::SeesawI_Theta->cwiseAbs2())(2,2);
    }

    void printable_ps191e(double& U_ps191e)
    {
      using namespace Pipes::printable_ps191e;
      double c_e = 0.5711;
      double c_mu = 0.1265;
      double c_tau = 0.1265;
      U_ps191e = (*Dep::Ue1 * ((c_e * *Dep::Ue1) + (c_mu * *Dep::Um1) + (c_tau * *Dep::Ut1))) + (*Dep::Ue2 * ((c_e * *Dep::Ue2) + (c_mu * *Dep::Um2) + (c_tau * *Dep::Ut2)));
    }

    void perturbativity_likelihood(double &lnL)
    {
      using namespace Pipes::perturbativity_likelihood;
      Matrix3d Usq = Dep::SeesawI_Theta->cwiseAbs2();
      
      if( *Param["M_1"] and 
          Usq(0,0) < *Param["M_2"] / *Param["M_1"] * Usq(0,1) + *Param["M_3"] / *Param["M_1"] * Usq(0,2) and
          Usq(1,0) < *Param["M_2"] / *Param["M_1"] * Usq(1,1) + *Param["M_3"] / *Param["M_1"] * Usq(1,2) and
          Usq(2,0) < *Param["M_2"] / *Param["M_1"] * Usq(2,1) + *Param["M_3"] / *Param["M_1"] * Usq(2,2) )
        lnL = 0;
      else
        lnL = -1E10;
    }

  }

}
