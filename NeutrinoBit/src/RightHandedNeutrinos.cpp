//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Right handed neutrino scan; using Casas-Ibarra parameterization
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
    void RHN_bbn_lifetime(std::vector<double>& result_lifetime)
    {
      using namespace Pipes::RHN_bbn_lifetime;
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
   
    // Lepton universality constraint: R_(e,mu)_pi/R_(e,mu)_K should be within experimental limits [R_pi_SM, R_K_SM: Phys. Rev. Lett 99, 231801; R_tau_SM: Int. J. Mod. Phys. A 24, 715, 2009; R_pi experimental limits: Phys. Rev. Lett. 70, 17; R_K experimental limits (NA62): Phys. Lett. B 719 (2013), 326; R_tau experimental limits: Phys. Rev. D 86, 010001]
    void RHN_R_pi(double& R_pi)
    {
      using namespace Pipes::RHN_R_pi;
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

    void RHN_R_K(double& R_K)
    {
      using namespace Pipes::RHN_R_K;
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

    void RHN_R_tau(double& R_tau)
    {
      using namespace Pipes::RHN_R_tau;
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
      d_r_tau = ((1.0 - mu_f_tau)/(1.0 - e_f_tau)) - 1.0;
      R_tau = R_tau_SM * (1.0 + d_r_tau);
    }

    // Lepton universality from W decays
    // 0: R(W->mu nu/W->e nu) from LHCb 1608.01484
    // 1: R(W->tau nu/W->e nu) from LEP 1302.3415
    // 2: R(W->tau nu/W->mu nu) from LEP 1302.3415
    void RHN_R_W(std::vector<double> &R_W)
    {
      using namespace Pipes::RHN_R_W;
      Matrix3d ThetaNorm = (*Dep::SeesawI_Theta * Dep::SeesawI_Theta->adjoint()).real();

      if(*Param["M_1"] < Dep::mw->central)
      {
        R_W.push_back(sqrt((1.0 - ThetaNorm(1,1))/(1.0 - ThetaNorm(0,0))));
        R_W.push_back(sqrt((1.0 - ThetaNorm(2,2))/(1.0 - ThetaNorm(0,0))));
        R_W.push_back(sqrt((1.0 - ThetaNorm(2,2))/(1.0 - ThetaNorm(1,1))));
      }
      else
      {
        R_W = {1.0, 1.0, 1.0};
      }
    }

    void lnL_lepuniv(double& result_lepuniv)
    {
      using namespace Pipes::lnL_lepuniv;
      double R_pi = *Dep::R_pi;
      double R_K = *Dep::R_K;
      double R_tau = *Dep::R_tau;
      std::vector<double> R_W = *Dep::R_W;

      double R_pi_exp = 1.23e-4; // Phys.Rev.Lett. 70 (1993) 17-20  
      double R_pi_err = 0.005e-4;
      double R_K_exp = 2.488e-5; // 1212.4012
      double R_K_err = 0.010e-5;
      double R_tau_exp = 0.9762; // 1612.07233 
      double R_tau_err = 0.0028;
      std::vector<double> R_W_exp = {0.980, 1.063, 1.070};
      std::vector<double> R_W_err = {0.018, 0.027, 0.026};

      result_lepuniv = 0.0;
      result_lepuniv += Stats::gaussian_loglikelihood(R_pi, R_pi_exp, 0.0, R_pi_err, false);
      result_lepuniv += Stats::gaussian_loglikelihood(R_K, R_K_exp, 0.0, R_K_err, false);
      result_lepuniv += Stats::gaussian_loglikelihood(R_tau, R_tau_exp, 0.0, R_tau_err, false);
      result_lepuniv += Stats::gaussian_loglikelihood(R_W[0], R_W_exp[0], 0.0, R_W_err[0], false);
      result_lepuniv += Stats::gaussian_loglikelihood(R_W[1], R_W_exp[1], 0.0, R_W_err[1], false);
      result_lepuniv += Stats::gaussian_loglikelihood(R_W[2], R_W_exp[2], 0.0, R_W_err[2], false);
    }

    // Neutrinoless double-beta decay constraint: m_bb should be less than the experimentally determined limits in GERDA and KamLAND-Zen [GERDA: Phys. Rev. Lett. 111 (2013) 122503; KamLAND-Zen: Phys. Rev. Lett 117 (2016) 082503]
    void RHN_m_GERDA(double &m_GERDA)
    {
      using namespace Pipes::RHN_m_GERDA;
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

    // Calculate 0nubb decay rate [1/s] for 136Xe 0nubb detector, for right-handed
    // neutrino model
    void RHN_Gamma_0nubb_Xe(double& result)
    {
      using namespace Pipes::RHN_Gamma_0nubb_Xe;
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

    void RHN_m_Kam(double& m_Kam)
    {
      using namespace Pipes::RHN_m_Kam;
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

    }

    // CKM unitarity constraint: V_ud should lie within 3sigma of the world average [PDG 2016]
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

    // Function to fill a spline object from a file
    tk::spline fill_spline(std::string file)
    {
      tk::spline s;
      std::vector<double> M_temp, U_temp;

      std::vector<std::pair<double,double> > array;
      std::ifstream f(file);
      while(f.good())
      {
        std::string line;
        getline(f, line);
        if (!f.good())
          break;
        std::stringstream iss(line);
        std::pair<double,double> point;
        iss >> point.first;
        iss.ignore();
        iss >> point.second;
        array.push_back(point);
      }

      for (unsigned int i=0; i<array.size(); i++)
      {
        M_temp.push_back(array[i].first);
        U_temp.push_back(array[i].second);
      }
      s.set_points(M_temp, U_temp);

      return s;
    }

    // Likelihood contribution from PIENU; searched for extra peaks in the spectrum of pi -> mu + nu. Constrains |U_ei|^2 at 90% in the mass range 60-129 MeV. [Phys. Rev. D, 84(5), 2011]
    void lnL_pienu(double& result)
    {
      using namespace Pipes::lnL_pienu;
      static bool read_table = true;
      static tk::spline s;
      double M_1, M_2, M_3;
      std::vector<double> U(3), mixing_sq(3);

      mixing_sq[0] = *Dep::Ue1;
      mixing_sq[1] = *Dep::Ue2;
      mixing_sq[2] = *Dep::Ue3;
      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];

      if (read_table)
      {
        s = fill_spline("NeutrinoBit/data/pienu.csv");
        read_table = false;
      }

      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);

      // Assume Gaussian errors with zero mean and that limits scale as |U|^2.
      result = 0;
      for(int i=0; i<3; i++)
        result += Stats::gaussian_upper_limit(mixing_sq[i]/U[i], 0, 0, 1/1.28, false);  // exp_error = abs(exp_value - 90CL_value)/, exp_value = 0. 1.28: 90% CL limit for half-Gaussian.

    }

    // Likelihood contribution from PS191, electron sector; looked for charged tracks originating from RHN decays: nu_r -> l(-) + l(+) + nu / l + pi / e + pi(+) + pi(0). Constrains |U_ei|^2 at 90% in the mass range 20-450 MeV. Function also incorporates a later re-interpretation of the data to account for neutral current interaction (ignored in original) as well as the RHNs' Majorana nature. [Original: Phys. Lett. B, 203(3):332-334, 1988][Re-interp.: JHEP, 2012(6):1-27]
    void lnL_ps191_e(double& result)
    {
      using namespace Pipes::lnL_ps191_e;
      static bool read_table = true;
      static tk::spline s;
      double M_1, M_2, M_3;
      std::vector<double> U(3), mixing_sq(3);
      double c_e = 0.5711;
      double c_mu = 0.1265;
      double c_tau = 0.1265;

      mixing_sq[0] = *Dep::Ue1 * ((c_e * *Dep::Ue1) + (c_mu * *Dep::Um1) + (c_tau * *Dep::Ut1));
      mixing_sq[1] = *Dep::Ue2 * ((c_e * *Dep::Ue2) + (c_mu * *Dep::Um2) + (c_tau * *Dep::Ut2));
      mixing_sq[2] = *Dep::Ue3 * ((c_e * *Dep::Ue3) + (c_mu * *Dep::Um3) + (c_tau * *Dep::Ut3));

      if (read_table)
      {
        s = fill_spline("NeutrinoBit/data/ps191_e.csv");
        read_table = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      // Division by sqrt(2) to correct for Majorana nature of RHNs.
      U[0] = s(M_1)/sqrt(2);
      U[1] = s(M_2)/sqrt(2);
      U[2] = s(M_3)/sqrt(2);

      // Assume scaling with |U|^4, zero bkg, number of events at 90% CL is
      // reverse engineered.  We assume that lnL = mu_sig is a faithful
      // approximation to the true Poisson likelihood.
      result = -2.44*((mixing_sq[0]/pow(U[0], 2.0)) + (mixing_sq[1]/pow(U[1], 2.0)) + (mixing_sq[2]/pow(U[2], 2.0)));
    }

    // Likelihood contribution from PS191, muon sector. Constrains |U_(mu,i)|^2 at 90% in the mass range 20-450 MeV. Description & references above.
    void lnL_ps191_mu(double& result)
    {
      using namespace Pipes::lnL_ps191_mu;
      static bool read_table = true;
      static tk::spline s;
      double M_1, M_2, M_3;
      std::vector<double> U(3), mixing_sq(3);
      double c_e = 0.5711;
      double c_mu = 0.1265;
      double c_tau = 0.1265;

      mixing_sq[0] = *Dep::Um1 * ((c_e * *Dep::Ue1) + (c_mu * *Dep::Um1) + (c_tau * *Dep::Ut1));
      mixing_sq[1] = *Dep::Um2 * ((c_e * *Dep::Ue2) + (c_mu * *Dep::Um2) + (c_tau * *Dep::Ut2));
      mixing_sq[2] = *Dep::Um3 * ((c_e * *Dep::Ue3) + (c_mu * *Dep::Um3) + (c_tau * *Dep::Ut3));

      if (read_table)
      {
        s = fill_spline("NeutrinoBit/data/ps191_mu.csv");
        read_table = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      // Division by sqrt(2) to correct for Majorana nature of RHNs.
      U[0] = s(M_1)/sqrt(2);
      U[1] = s(M_2)/sqrt(2);
      U[2] = s(M_3)/sqrt(2);

      // Assume scaling with |U|^4, zero bkg, number of events at 90% CL is
      // reverse engineered.  We assume that lnL = mu_sig is a faithful
      // approximation to the true Poisson likelihood.
      result = -2.44*((mixing_sq[0]/pow(U[0], 2.0)) + (mixing_sq[1]/pow(U[1], 2.0)) + (mixing_sq[2]/pow(U[2], 2.0)));
    }

    // Likelihood contribution from CHARM, electron sector; searched for charged and neutral current decays of RHNs. Constrains |U_ei|^2 at 90% in the mass range 0.5-2.8 GeV. [Phys. Lett. B, 166(4):473-478, 1986]
    void lnL_charm_e(double& result)
    {
      using namespace Pipes::lnL_charm_e;
      static bool read_table = true;
      static tk::spline s;
      double M_1, M_2, M_3;
      std::vector<double> U(3), mixing_sq(3);
      double c_e = 0.5711;
      double c_mu = 0.1265;
      double c_tau = 0.1265;

      mixing_sq[0] = *Dep::Ue1 * ((c_e * *Dep::Ue1) + (c_mu * *Dep::Um1) + (c_tau * *Dep::Ut1));
      mixing_sq[1] = *Dep::Ue2 * ((c_e * *Dep::Ue2) + (c_mu * *Dep::Um2) + (c_tau * *Dep::Ut2));
      mixing_sq[2] = *Dep::Ue3 * ((c_e * *Dep::Ue3) + (c_mu * *Dep::Um3) + (c_tau * *Dep::Ut3));

      if (read_table)
      {
        s = fill_spline("NeutrinoBit/data/charm_e.csv");
        read_table = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);

      // Assume scaling with |U|^4, zero bkg, number of events at 90% CL is
      // reverse engineered.  We assume that lnL = mu_sig is a faithful
      // approximation to the true Poisson likelihood.
      result = -2.44*((mixing_sq[0]/pow(U[0], 2.0)) + (mixing_sq[1]/pow(U[1], 2.0)) + (mixing_sq[2]/pow(U[2], 2.0)));
    }

    // Likelihood contribution from CHARM, muon sector. Constrains |U_(mu,i)|^2 at 90% in the mass range 0.5-2.8 GeV. Description & references above.
    void lnL_charm_mu(double& result)
    {
      using namespace Pipes::lnL_charm_mu;
      static bool read_table = true;
      static tk::spline s;
      double M_1, M_2, M_3;
      std::vector<double> U(3), mixing_sq(3);
      double c_e = 0.5711;
      double c_mu = 0.1265;
      double c_tau = 0.1265;

      mixing_sq[0] = *Dep::Um1 * ((c_e * *Dep::Ue1) + (c_mu * *Dep::Um1) + (c_tau * *Dep::Ut1));
      mixing_sq[1] = *Dep::Um2 * ((c_e * *Dep::Ue2) + (c_mu * *Dep::Um2) + (c_tau * *Dep::Ut2));
      mixing_sq[2] = *Dep::Um3 * ((c_e * *Dep::Ue3) + (c_mu * *Dep::Um3) + (c_tau * *Dep::Ut3));

      if (read_table)
      {
        s = fill_spline("NeutrinoBit/data/charm_mu.csv");
        read_table = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);

      // Assume scaling with |U|^4, zero bkg, number of events at 90% CL is
      // reverse engineered.  We assume that lnL = mu_sig is a faithful
      // approximation to the true Poisson likelihood.
      result = -2.44*((mixing_sq[0]/pow(U[0], 2.0)) + (mixing_sq[1]/pow(U[1], 2.0)) + (mixing_sq[2]/pow(U[2], 2.0)));
    }

    // Likelihood contribution from DELPHI; searched for charged and neutral current decays of RHNs. Constrains |U_ei|^2, |U_(mu,i)|^2 as well as |U_(tau,i)|^2 at 95% in the mass range 3.5-50 GeV. [Z. Phys. C, 74(1):57-71, 1997]
    void lnL_delphi(double& result)
    {
      using namespace Pipes::lnL_delphi;
      static bool read_table = true;
      static tk::spline s;
      double M_1, M_2, M_3;
      std::vector<double> U(3), mixing_sq(9);

      mixing_sq[0] = *Dep::Ue1;  // This is |U_{e1}|^2 etc
      mixing_sq[1] = *Dep::Ue2;
      mixing_sq[2] = *Dep::Ue3;
      mixing_sq[3] = *Dep::Um1;
      mixing_sq[4] = *Dep::Um2;
      mixing_sq[5] = *Dep::Um3;
      mixing_sq[6] = *Dep::Ut1;
      mixing_sq[7] = *Dep::Ut2;
      mixing_sq[8] = *Dep::Ut3;

      if (read_table)
      {
        s = fill_spline("NeutrinoBit/data/delphi.csv");
        read_table = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);

      // Assume scaling with |U|^4, zero bkg, number of events at 95% CL is
      // reverse engineered.  We assume that lnL = mu_sig is a faithful
      // approximation to the true Poisson likelihood.
      result = -3.09*
         (pow(mixing_sq[0]/U[0], 2.0) + pow(mixing_sq[1]/U[1], 2.0) + pow(mixing_sq[2]/U[2], 2.0) +
          pow(mixing_sq[3]/U[0], 2.0) + pow(mixing_sq[4]/U[1], 2.0) + pow(mixing_sq[5]/U[2], 2.0) + 
          pow(mixing_sq[6]/U[0], 2.0) + pow(mixing_sq[7]/U[1], 2.0) + pow(mixing_sq[8]/U[2], 2.0));
    }

    // Likelihood contribution from ATLAS, electron sector; looked at the production and decay chain: pp -> W*(+-) -> l(+-) + nu_r. nu_r then decays into an on-shell W and a lepton; the W decays primarily into a qq pair. Constrains |U_ei|^2 at 95% in the mass range 50-500 GeV. [JHEP, 07:162, 2015]
    void lnL_atlas_e(double& result)
    {
      using namespace Pipes::lnL_atlas_e;
      static bool read_table = true;
      static tk::spline s;
      double M_1, M_2, M_3;
      std::vector<double> U(3), mixing_sq(3);

      mixing_sq[0] = *Dep::Ue1;
      mixing_sq[1] = *Dep::Ue2;
      mixing_sq[2] = *Dep::Ue3;

      if (read_table)
      {
        s = fill_spline("NeutrinoBit/data/atlas_e.csv");
        read_table = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);

      // Assume Gaussian errors with zero mean and that limits scale as |U|^4.
      result = 0;
      for(int i=0; i<3; i++)
        result += Stats::gaussian_upper_limit(pow((mixing_sq[i]/U[i]), 2.0), 0, 0, 1/1.64, false);  // exp_error = abs(exp_value - 95CL_value)/1.64, exp_value = 0. 1.64: 95% CL limit for half-Gaussian.
    }

    // Likelihood contribution from ATLAS, muon sector. Constrains |U_(mu,i)|^2 at 95% in the mass range 50-500 GeV. Description & references above.
    void lnL_atlas_mu(double& result)
    {
      using namespace Pipes::lnL_atlas_mu;
      static bool read_table = true;
      static tk::spline s;
      double M_1, M_2, M_3;
      std::vector<double> U(3), mixing_sq(3);

      mixing_sq[0] = *Dep::Um1;
      mixing_sq[1] = *Dep::Um2;
      mixing_sq[2] = *Dep::Um3;

      if (read_table)
      {
        s = fill_spline("NeutrinoBit/data/atlas_mu.csv");
        read_table = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);

      // Assume Gaussian errors with zero mean and that limits scale as |U|^4.
      result = 0;
      for(int i=0; i<3; i++)
        result += Stats::gaussian_upper_limit(pow((mixing_sq[i]/U[i]), 2.0), 0, 0, 1/1.64, false);  // exp_error = abs(exp_value - 95CL_digitized)/1.64, exp_value = 0. 1.64: 95% CL limit for half-Gaussian.
    }

    // Likelihood contribution from E949; used the kaon decay: K(+) -> mu(+) + nu_r. Constrains |U_(mu,i)|^2 at 90% in the mass range 175-300 MeV. [Phys. Rev. D, 91, 2015]
    void lnL_e949(double& result)
    {
      using namespace Pipes::lnL_e949;
      static bool read_table = true;
      static tk::spline s;
      double M_1, M_2, M_3;
      std::vector<double> U(3), mixing_sq(3);

      mixing_sq[0] = *Dep::Um1;
      mixing_sq[1] = *Dep::Um2;
      mixing_sq[2] = *Dep::Um3;

      if (read_table)
      {
        s = fill_spline("NeutrinoBit/data/e949.csv");
        read_table = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      // Division by sqrt(2) to correct for Majorana nature of RHNs.
      U[0] = s(M_1)/sqrt(2);
      U[1] = s(M_2)/sqrt(2);
      U[2] = s(M_3)/sqrt(2);

      // Assume Gaussian errors with zero mean and that limits scale as |U|^2.
      result = 0;
      for(int i=0; i<3; i++)
        result += Stats::gaussian_upper_limit(mixing_sq[i]/U[i], 0, 0, 1/1.28, false);  // exp_error = abs(exp_value - 90CL_value)/1.28, exp_value = 0. 1.28: 90% CL limit for half-Gaussian.
    }

    // Likelihood contribution from NuTeV; used RHN decays into muonic final states (mu + mu + nu / mu + e + nu / mu + pi / mu + rho). Constrains |U_(mu,i)|^2 at 90% CL in the mass range 0.25-2 GeV. [Phys. Rev. Lett., 83:4943-4946, 1999]
    void lnL_nutev(double& result)
    {
      using namespace Pipes::lnL_nutev;
      static bool read_table = true;
      static tk::spline s;
      double M_1, M_2, M_3;
      std::vector<double> U(3), mixing_sq(3);

      mixing_sq[0] = *Dep::Um1;
      mixing_sq[1] = *Dep::Um2;
      mixing_sq[2] = *Dep::Um3;

      if (read_table)
      {
        s = fill_spline("NeutrinoBit/data/nutev.csv");
        read_table = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);

      // Assume scaling with |U|^4, zero bkg, number of events at 90% CL is
      // reverse engineered.  We assume that lnL = mu_sig is a faithful
      // approximation to the true Poisson likelihood.
      result = -2.44*(pow(mixing_sq[0]/U[0], 2.0) + pow(mixing_sq[1]/U[1], 2.0) + pow(mixing_sq[2]/U[2], 2.0));
    }

    // Likelihood contribution from a re-interpretation of CHARM data; assumes tau mixing is dominant. Constrains |U_(tau,i)|^2 at 90% CL in the mass range 10-290 MeV. [Phys. Lett. B, 550(1-2):8-15, 2002]
    void lnL_charm_tau(double& result)
    {
      using namespace Pipes::lnL_charm_tau;
      static bool read_table = true;
      static tk::spline s;
      double M_1, M_2, M_3;
      std::vector<double> U(3), mixing_sq(3);

      mixing_sq[0] = *Dep::Ut1;
      mixing_sq[1] = *Dep::Ut2;
      mixing_sq[2] = *Dep::Ut3;

      if (read_table)
      {
        s = fill_spline("NeutrinoBit/data/tau.csv");
        read_table = false;
      }

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      // Division by sqrt(2) to correct for Majorana nature of RHNs.
      U[0] = s(M_1)/sqrt(2);
      U[1] = s(M_2)/sqrt(2);
      U[2] = s(M_3)/sqrt(2);

      // Assume Gaussian errors with zero mean and that limits scale as |U|^4.
      result = 0;
      for(int i=0; i<3; i++)
        result += Stats::gaussian_upper_limit(pow((mixing_sq[i]/U[i]), 2.0), 0, 0, 1/1.28, false);  // exp_error = abs(exp_value - 90CL_value)/1.28, but exp_value = 0. 1.28: 90% CL limit for half-Gaussian.
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
