// GAMBIT: Global and Modular BSM Inference Tool
//
// *********************************************
//
// Sterile neutrino scan; using Casas-Ibarra parameterization
//
// *********************************************
//
// Authors
//
// Suraj Krishnamurthy
// 2017 February
//
// *********************************************

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/DarkBit/DarkBit_rollcall.hpp"
#include "gambit/Utils/numerical_constants.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "gambit/DarkBit/spline.h"
#include <vector>
using namespace Eigen;

namespace Gambit
{
  namespace DarkBit
  {
    // Gets the matrix Theta in the C-I parametrization from NeutrinoBit and returns its squared absolute |Theta|^2
    void CI_param(Matrix3d& result) 
    {
      using namespace Pipes::CI_param;
      //Matrix3d t_sq;
      Matrix3cd t;

      t = *Dep::SeesawI_Theta;
      result = t.cwiseAbs2();
    }

    // BBN constraint: lifetime must be less than 0.1s [arXiv:1202.2841]
    void lnL_bbn(double& result_bbn)
    {
      using namespace Pipes::lnL_bbn;
      SMInputs sminputs = *Dep::SMINPUTS;
      static double conv_fact = 6.58e-16;  // conversion factor from ev^-1 to s
      static double G_F_sq = pow(sminputs.GF, 2.0);  // GeV^-4
      // TODO (CW): Should come from SM input file
      static double g_L_twid_sq = 0.0771;  // g_L_twid = -0.5 + s_W_sq
      static double g_R_sq = 0.0494;  // g_R = s_W^2
      static double g_L_sq = 0.5217;  // g_L = 0.5 + s_W^2
      double temp_bbn;
      std::vector<double> lifetime(3), M(3);
      M[0] = *Param["M_1"];
      M[1] = *Param["M_2"];
      M[2] = *Param["M_3"];
      Matrix3d t_sq = *Dep::Theta_sq;

      for (int i=0; i<3; i++)
      {
        lifetime[i] = (96*pow(pi,3.0)*1e-9*conv_fact) / (G_F_sq*pow(M[i],5.0))*( ((1 + g_L_twid_sq + g_R_sq)*(t_sq(1,i) + t_sq(2,i))) + ((1 + g_L_sq + g_R_sq)*t_sq(0,i)) );
        if(lifetime[i]<0.1)
        {
          temp_bbn = 0.0;
        }
        else
        {
          temp_bbn = -100.0;
        }
      }
      result_bbn = temp_bbn;
//      result_bbn = lifetime[0];
    }

    // Lepton universality constraint: R_(e,mu)_pi/R_(e,mu)_K should be within experimental limits [R_pi_SM, R_K_SM: Phys. Rev. Lett 99, 231801; R_tau_SM: Int. J. Mod. Phys. A 24, 715, 2009; R_pi experimental limits: Phys. Rev. Lett. 70, 17; iR_K experimental limits (NA62): Phys. Lett. B 719 (2013), 326; R_tau experimental limits: Phys. Rev. D 86, 010001]
    void lnL_lepuniv(double& result_lepuniv)
    {
      using namespace Pipes::lnL_lepuniv;
      SMInputs sminputs = *Dep::SMINPUTS;
      static double m_pi = meson_masses.pi_plus;
      static double m_K = meson_masses.kaon_plus; 
      static double m_tau = sminputs.mTau;  // GeV
      static double R_pi_SM = 1.2354e-4;
      static double R_K_SM = 2.477e-5;
      static double R_tau_SM = 0.973;
      static double r_e_pi = pow(sminputs.mE,2)/pow(m_pi,2);
      static double r_mu_pi = pow(sminputs.mMu,2)/pow(m_pi,2);
      static double r_e_K = pow(sminputs.mE,2)/pow(m_K,2);
      static double r_mu_K = pow(sminputs.mMu,2)/pow(m_K,2);
      double e_f_pi, mu_f_pi, e_f_K, mu_f_K, e_f_tau, mu_f_tau, d_r_pi, d_r_K, d_r_tau, R_pi, R_K, R_tau, temp_lepuniv;
      std::vector<double> M(3), r_I_pi(3), G_e_pi(3), G_mu_pi(3), e_fac_pi(3), mu_fac_pi(3), r_I_K(3), G_e_K(3), G_mu_K(3), e_fac_K(3), mu_fac_K(3), e_fac_tau(3), mu_fac_tau(3);
      Matrix3d t_sq = *Dep::Theta_sq;

      e_f_pi = 0.0;
      mu_f_pi = 0.0;
      e_f_K = 0.0;
      mu_f_K = 0.0;
      e_f_tau = 0.0;
      mu_f_tau = 0.0;
      M[0] = *Param["M_1"];
      M[1] = *Param["M_2"];
      M[2] = *Param["M_3"];

      for (int i=0; i<3; i++)
      {
        if(M[i] < m_K)
        {
          if(M[i] < m_pi)
          {
            r_I_pi[i] = pow(M[i], 2.0)/pow(m_pi, 2.0);
            G_e_pi[i] = (r_e_pi + r_I_pi[i] - pow((r_e_pi - r_I_pi[i]), 2.0) * sqrt(1.0 - (2.0*pow((r_e_pi + r_I_pi[i]), 2.0)) + pow((r_e_pi - r_I_pi[i]), 2.0))) / (r_e_pi * pow((1.0 - r_e_pi), 2.0));
            G_mu_pi[i] = (r_mu_pi + r_I_pi[i] - pow((r_mu_pi - r_I_pi[i]), 2.0) * sqrt(1.0 - (2.0*pow((r_mu_pi + r_I_pi[i]), 2.0)) + pow((r_mu_pi - r_I_pi[i]), 2.0))) / (r_mu_pi * pow((1.0 - r_mu_pi), 2.0));
            e_fac_pi[i] = t_sq(0,i) * (G_e_pi[i] - 1.0);
            mu_fac_pi[i] = t_sq(1,i) * (G_mu_pi[i] - 1.0);
          }
          else
          {
            r_I_K[i] = pow(M[i], 2.0)/pow(m_K, 2.0);
            G_e_K[i] = (r_e_K + r_I_K[i] - pow((r_e_K - r_I_K[i]), 2.0) * sqrt(1.0 - (2.0*pow((r_e_K + r_I_K[i]), 2.0)) + pow((r_e_K - r_I_K[i]), 2.0))) / (r_e_pi * pow((1.0 - r_e_pi), 2.0));
            G_mu_K[i] = (r_mu_K + r_I_K[i] - pow((r_mu_K - r_I_K[i]), 2.0) * sqrt(1.0 - (2.0*pow((r_mu_K + r_I_K[i]), 2.0)) + pow((r_mu_K - r_I_K[i]), 2.0))) / (r_mu_K * pow((1.0 - r_mu_K), 2.0));
            e_fac_K[i] = t_sq(0,i) * (G_e_pi[i] - 1.0);
            mu_fac_K[i] = t_sq(1,i) * (G_mu_pi[i] - 1.0);
          }
        } 
        else if(M[i] > m_tau)
        {
          e_fac_tau[i] = t_sq(0,i);
          mu_fac_tau[i] = t_sq(1,i);
        }
        else
        {
          temp_lepuniv = 0.0;
        }
        e_f_pi += e_fac_pi[i];
        mu_f_pi += mu_fac_pi[i];
        e_f_K += e_fac_K[i];
        mu_f_K += mu_fac_K[i];
        e_f_tau += e_fac_tau[i];
        mu_f_tau += mu_fac_tau[i];
      }
      d_r_pi = ((1.0 + e_f_pi)/(1.0 + mu_f_pi)) - 1.0;
      d_r_K = ((1.0 + e_f_K)/(1.0 + mu_f_K)) - 1.0;
      d_r_tau = ((1.0 - e_f_tau)/(1.0 - mu_f_tau)) - 1.0;
      R_pi = R_pi_SM * (1.0 + d_r_pi);
      R_K = R_K_SM * (1.0 + d_r_K);
      R_tau = R_tau_SM * (1.0 + d_r_tau);
      if (((1.218e-4<R_pi) && (R_pi<1.242e-4)) || ((2.458e-5<R_K) && (R_K<2.518e-5)) || ((0.9674<R_tau) && (R_tau<0.9854)))
      {
        temp_lepuniv = 0.0;
      }
      else
      {
        temp_lepuniv = -100.0;
      }
      result_lepuniv = temp_lepuniv;
//      result_lepuniv = R_pi;
    }

    // Neutrinoless double-beta decay constraint: m_bb should be less than the experimentally determined limits in GERDA and KamLAND-Zen [GERDA: Phys. Rev. Lett. 111 (2013) 122503; KamLAND-Zen: Phys. Rev. Lett 117 (2016) 082503]
    void lnL_0nubb(double& result_0nubb)
    {
      using namespace Pipes::lnL_0nubb;
      static double m_bb_GERDA = 4e-10;  // GeV
      static double m_bb_Kam = 1.65e-10;  // GeV
      Matrix3cd m_light, U_light;
      Matrix3d U_light_sq, t_sq;
      std::vector<double> M(3), m_temp_GERDA(3), m_temp_Kam(3);
      double m_GERDA, m_Kam;
      m_light = *Dep::m_nu;
      U_light = *Dep::UPMNS;
      U_light_sq = U_light.cwiseAbs2();
      t_sq = *Dep::Theta_sq;
      m_GERDA = 0.0;
      m_Kam = 0.0;
      M[0] = *Param["M_1"];
      M[1] = *Param["M_2"];
      M[2] = *Param["M_3"];

      for (int i=0; i<3; i++)
      {
        m_temp_GERDA[i] = U_light_sq(0,i)*abs(m_light(i,i)) + t_sq(0,i)*M[i]*(pow(*Param["L_Ge"], 2.0)/(pow(*Param["L_Ge"], 2.0)+pow(M[i], 2.0)));
        m_GERDA += m_temp_GERDA[i];
        m_temp_Kam[i] = U_light_sq(0,i)*abs(m_light(i,i)) + t_sq(0,i)*M[i]*(pow(*Param["L_Xe"], 2.0)/(pow(*Param["L_Xe"], 2.0)+pow(M[i], 2.0)));
        m_Kam +=m_temp_Kam[i];
      }
//      if ((m_GERDA < m_bb_GERDA) && (m_Kam < m_bb_Kam))
//      if (m_Kam < m_bb_Kam)
//      {
//        result_0nubb = 0.0;
//      }
//      else
//      {
//        result_0nubb = -100.0;
//      }
      result_0nubb = m_GERDA;
    }

    // CKM unitarity constraint: V_ud should lie within 3sigma of the world average [PDG 2016]
    void lnL_ckm(double& result_ckm)
    {
      using namespace Pipes::lnL_ckm;
      SMInputs sminputs = *Dep::SMINPUTS;
      static double G_F_sq = pow(sminputs.GF, 2.0);  // GeV^-4
      static double V_us_exp[7] = {0.2235,0.2227,0.2244,0.2238,0.2242,0.2262,0.2214};
      static double err[7] = {0.0006,0.0013,0.0008,0.0006,0.0011,0.0013,0.0022};
      static double G_mu_sq = 1.3605e-10;  // GeV^-4
      static double V_ud_0 = 0.97427;
      static double err_0 = 0.00015;
      static double err_factor = 5.3779e7;
      double V_ud_exp[7], f[7];
      double V_ud_sq, chi2;
      Matrix3d t_sq = *Dep::Theta_sq;

      for (int i=0;i<7;i++)
      {
        V_ud_exp[i] = sqrt(1 - pow(V_us_exp[i], 2.0));
      }
      f[0] = (G_F_sq/G_mu_sq)*(1 - t_sq(0,0));
      f[3] = (G_F_sq/G_mu_sq)*(1 - t_sq(1,1));
      f[5] = 1 + t_sq(1,1);
      f[6] = 1 + t_sq(0,0) + t_sq(1,1) + t_sq(2,2);
      f[1] = f[0];
      f[2] = f[0];
      f[4] = f[3];
      V_ud_sq = 0.0;
      for (int j=0; j<7; j++)
      {
        V_ud_sq += pow(V_ud_exp[j], 2.0)/(pow(err[j], 2.0)*f[j]);
      }
      V_ud_sq += pow(V_ud_0, 2.0)/(pow(err_0, 2.0)*(1 + t_sq(0,0)));
      V_ud_sq /= err_factor;
      chi2 = 0.0;
      for (int k=0; k<7; k++)
      {
        chi2 += pow(((sqrt((1 - V_ud_sq)*f[k]) - V_us_exp[k])/err[k]), 2.0);
      }
      chi2 += pow((((sqrt(V_ud_sq)*(1 + t_sq(0,0))) - V_ud_0)/err_0), 2.0);
      if (chi2 < 23.5744)
      { 
        result_ckm = 0.0;
      }
      else
      {
        result_ckm = -100.0;
      }
//      result_ckm = V_ud_sq;
    }

    // Likelihood contribution from PIENU; searched for extra peaks in the spectrum of pi -> mu + nu. Constrains |U_ei|^2 in the mass range 60-129 MeV. [Phys. Rev. D, 84(5), 2011]
    void lnL_pienu(double& result_pienu)
    {
      using namespace Pipes::lnL_pienu;
      static bool read_table_pienu = true;
      static tk::spline s_pienu;
      static std::vector<double> M_temp_pienu(140), U_temp_pienu(140);
      double M_1, M_2, M_3;
      std::vector<double> U_pienu(3), mixing_sq_pienu(3);

      mixing_sq_pienu[0] = *Dep::Ue1;
      mixing_sq_pienu[1] = *Dep::Ue2;
      mixing_sq_pienu[2] = *Dep::Ue3;
      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];

      if (read_table_pienu)
      {
        double array_pienu[140][2];
        std::ifstream f_pienu("DarkBit/data/pienu.csv");
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
      result_pienu = -2.44*((mixing_sq_pienu[0]/U_pienu[0]) + (mixing_sq_pienu[1]/U_pienu[1]) + (mixing_sq_pienu[2]/U_pienu[2]));
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
        std::ifstream f_ps191e("DarkBit/data/ps191_e.csv");
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
        std::ifstream f_ps191mu("DarkBit/data/ps191_mu.csv");
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
      std::vector<double> U_charme(3), mixing_sq_charme(3);
      double c_e = 0.5711;
      double c_mu = 0.1265;
      double c_tau = 0.1265;

      mixing_sq_charme[0] = *Dep::Ue1 * ((c_e * *Dep::Ue1) + (c_mu * *Dep::Um1) + (c_tau * *Dep::Ut1));
      mixing_sq_charme[1] = *Dep::Ue2 * ((c_e * *Dep::Ue2) + (c_mu * *Dep::Um2) + (c_tau * *Dep::Ut2));
      mixing_sq_charme[2] = *Dep::Ue3 * ((c_e * *Dep::Ue3) + (c_mu * *Dep::Um3) + (c_tau * *Dep::Ut3));

      if (read_table_charme)
      {
        double array_charme[56][2];
        std::ifstream f_charme("DarkBit/data/charm_e.csv");
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
      result_charme = -2.44*((mixing_sq_charme[0]/pow(U_charme[0], 2.0)) + (mixing_sq_charme[1]/pow(U_charme[1], 2.0)) + (mixing_sq_charme[2]/pow(U_charme[2], 2.0)));
    }

    // Likelihood contribution from CHARM, muon sector. Constrains |U_(mu,i)|^2 in the mass range 0.5-2.8 GeV. Description & references above.
    void lnL_charm_mu(double& result_charmmu)
    {
      using namespace Pipes::lnL_charm_mu;
      static bool read_table_charmmu = true;
      static tk::spline s_charmmu;
      double M_1, M_2, M_3;
      static std::vector<double> M_temp_charmmu(34), U_temp_charmmu(34);
      std::vector<double> U_charmmu(3), mixing_sq_charmmu(3);
      double c_e = 0.5711;
      double c_mu = 0.1265;
      double c_tau = 0.1265;

      mixing_sq_charmmu[0] = *Dep::Um1 * ((c_e * *Dep::Ue1) + (c_mu * *Dep::Um1) + (c_tau * *Dep::Ut1));
      mixing_sq_charmmu[1] = *Dep::Um2 * ((c_e * *Dep::Ue2) + (c_mu * *Dep::Um2) + (c_tau * *Dep::Ut2));
      mixing_sq_charmmu[2] = *Dep::Um3 * ((c_e * *Dep::Ue3) + (c_mu * *Dep::Um3) + (c_tau * *Dep::Ut3));

      if (read_table_charmmu)
      {
        double array_charmmu[34][2];
        std::ifstream f_charmmu("DarkBit/data/charm_mu.csv");
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
      result_charmmu = -2.44*((mixing_sq_charmmu[0]/pow(U_charmmu[0], 2.0)) + (mixing_sq_charmmu[1]/pow(U_charmmu[1], 2.0)) + (mixing_sq_charmmu[2]/pow(U_charmmu[2], 2.0)));
    }

    // Likelihood contribution from DELPHI; searched for charged and neutral current decays of RHNs. Constrains |U_ei|^2, |U_(mu,i)|^2 as well as |U_(tau,i)|^2 in the mass range 3.5-50 GeV. [Z. Phys. C, 74(1):57-71, 1997]
    void lnL_delphi(double& result_delphi)
    {
      using namespace Pipes::lnL_delphi;
      static bool read_table_delphi = true;
      static tk::spline s_delphi;
      double M_1, M_2, M_3;
      static std::vector<double> M_temp_delphi(180), U_temp_delphi(180);
      std::vector<double> U_delphi(3), mixing_sq_delphi(9);

      mixing_sq_delphi[0] = *Dep::Ue1;
      mixing_sq_delphi[1] = *Dep::Ue2;
      mixing_sq_delphi[2] = *Dep::Ue3;
      mixing_sq_delphi[3] = *Dep::Um1;
      mixing_sq_delphi[4] = *Dep::Um2;
      mixing_sq_delphi[5] = *Dep::Um3;
      mixing_sq_delphi[6] = *Dep::Ut1;
      mixing_sq_delphi[7] = *Dep::Ut2;
      mixing_sq_delphi[8] = *Dep::Ut3;

      if (read_table_delphi)
      {
        double array_delphi[180][2];
        std::ifstream f_delphi("DarkBit/data/delphi.csv");
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
      result_delphi = -0.1*((mixing_sq_delphi[0]/U_delphi[0]) + (mixing_sq_delphi[1]/U_delphi[1]) + (mixing_sq_delphi[2]/U_delphi[2]) + (mixing_sq_delphi[3]/U_delphi[0]) + (mixing_sq_delphi[4]/U_delphi[1]) + (mixing_sq_delphi[5]/U_delphi[2]) + (mixing_sq_delphi[6]/U_delphi[0]) + (mixing_sq_delphi[7]/U_delphi[1]) + (mixing_sq_delphi[8]/U_delphi[2]));
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
        std::ifstream f_atlase("DarkBit/data/atlas_e.csv");
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
        std::ifstream f_atlasmu("DarkBit/data/atlas_mu.csv");
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
        std::ifstream f_e949("DarkBit/data/e949.csv");
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
      std::vector<double> U_nutev(3), mixing_sq_nutev(3);

      mixing_sq_nutev[0] = *Dep::Um1;
      mixing_sq_nutev[1] = *Dep::Um2;
      mixing_sq_nutev[2] = *Dep::Um3;

      if (read_table_nutev)
      {
        double array_nutev[249][2];
        std::ifstream f_nutev("DarkBit/data/nutev.csv");
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
      result_nutev = -2.44*((mixing_sq_nutev[0]/U_nutev[0]) + (mixing_sq_nutev[1]/U_nutev[1]) + (mixing_sq_nutev[2]/U_nutev[2]));
    }

    // Likelihood contribution from a re-interpretation of CHARM data; assumes tau mixing is dominant. Constrains |U_(tau,i)|^2 in the mass range 10-290 MeV. [Phys. Lett. B, 550(1-2):8-15, 2002]
    void lnL_tau(double& result_tau)
    {
      using namespace Pipes::lnL_tau;
      static bool read_table_tau = true;
      static tk::spline s_tau;
      double M_1, M_2, M_3;
      static std::vector<double> M_temp_tau(172), U_temp_tau(172);
      std::vector<double> U_tau(3), mixing_sq_tau(3);

      mixing_sq_tau[0] = *Dep::Ut1;
      mixing_sq_tau[1] = *Dep::Ut2;
      mixing_sq_tau[2] = *Dep::Ut3;

      if (read_table_tau)
      {
        double array_tau[172][2];
        std::ifstream f_tau("DarkBit/data/tau.csv");
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
      result_tau = -2.44*((mixing_sq_tau[0]/U_tau[0]) + (mixing_sq_tau[1]/U_tau[1]) + (mixing_sq_tau[2]/U_tau[2]));
    }

    void printable_Ue1(double& Ue1_sq)
    {
      namespace myPipe2 = Pipes::printable_Ue1;
      Matrix3d t_1(*myPipe2::Dep::Theta_sq);
      Ue1_sq = t_1(0,0);
    }

    void printable_Um1(double& Um1_sq)
    {
      namespace myPipe3 = Pipes::printable_Um1;
      Matrix3d t_2(*myPipe3::Dep::Theta_sq);
      Um1_sq = t_2(1,0);
    }

    void printable_Ut1(double& Ut1_sq)
    {
      namespace myPipe4 = Pipes::printable_Ut1;
      Matrix3d t_3(*myPipe4::Dep::Theta_sq);
      Ut1_sq = t_3(2,0);
    }

    void printable_Ue2(double& Ue2_sq)
    {
      namespace myPipe5 = Pipes::printable_Ue2;
      Matrix3d t_4(*myPipe5::Dep::Theta_sq);
      Ue2_sq = t_4(0,1);
    }

    void printable_Um2(double& Um2_sq)
    {
      namespace myPipe6 = Pipes::printable_Um2;
      Matrix3d t_5(*myPipe6::Dep::Theta_sq);
      Um2_sq = t_5(1,1);
    }

    void printable_Ut2(double& Ut2_sq)
    {
      namespace myPipe7 = Pipes::printable_Ut2;
      Matrix3d t_6(*myPipe7::Dep::Theta_sq);
      Ut2_sq = t_6(2,1);
    }

    void printable_Ue3(double& Ue3_sq)
    {
      namespace myPipe8 = Pipes::printable_Ue3;
      Matrix3d t_7(*myPipe8::Dep::Theta_sq);
      Ue3_sq = t_7(0,2);
    }

    void printable_Um3(double& Um3_sq)
    {
      namespace myPipe9 = Pipes::printable_Um3;
      Matrix3d t_8(*myPipe9::Dep::Theta_sq);
      Um3_sq = t_8(1,2);
    }

    void printable_Ut3(double& Ut3_sq)
    {
      namespace myPipe10 = Pipes::printable_Ut3;
      Matrix3d t_9(*myPipe10::Dep::Theta_sq);
      Ut3_sq = t_9(2,2);
    }

    void printable_ps191e(double& U_ps191e)
    {
      using namespace Pipes::printable_ps191e;
      double c_e = 0.5711;
      double c_mu = 0.1265;
      double c_tau = 0.1265;
      U_ps191e = (*Dep::Ue1 * ((c_e * *Dep::Ue1) + (c_mu * *Dep::Um1) + (c_tau * *Dep::Ut1))) + (*Dep::Ue2 * ((c_e * *Dep::Ue2) + (c_mu * *Dep::Um2) + (c_tau * *Dep::Ut2)));
    }

  }

}
