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
//#include "gambit/Elements/numerical_constants.hpp"
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
      typedef std::complex<double> dcomp;
      SMInputs sminputs = *Dep::SMINPUTS;
      static double conv_fact = 6.58e-16;  // conversion factor from eV^-1 to s, for lifetime
      static double G_F_sq = pow(sminputs.GF, 2.0);  // GeV^-4
      static double g_L_twid_sq = 0.0771;  // g_L_twid = -0.5 + s_W_sq
      static double g_R_sq = 0.0494;  // g_R = s_W^2
      static double g_L_sq = 0.5217;  // g_L = 0.5 + s_W^2
      static double R_pi_SM = 1.2354e-4;
      static double R_K_SM = 2.477e-5;
      static double R_tau_SM = 0.973;
      static double r_e_pi = 1.3399e-5;  // r_e_pi = m_e^2/m_pi^2
      static double r_mu_pi = 0.5733;  // r_mu_pi = m_mu^2/m_pi^2
      static double r_e_K = 1.0713e-6;  // r_e_K = m_e^2/m_K^2
      static double r_mu_K = 0.0458;  // r_mu_K = m_mu^2/m_K^2
      static double m_pi = 0.1396; // GeV
      static double m_K = 0.4937;  // GeV
      static double m_tau = 1.7768;  // GeV
      static double m_bb_GERDA = 2e-10;  // GeV
      static double m_bb_Kam = 1.61e-10;  // GeV
      static double G_mu_sq = 1.3605e-10;  // GeV^-4
      dcomp I(0.0, 1.0);
      Matrix3d M_I, t_sq, result_temp_bbn, result_temp_lepuniv, result_temp_0nubb, result_temp_ckm, U_light_sq;  // M_I not complex; circumvents type mismatch in l(M)
      Matrix3cd m_light, U_light ,t;
      std::vector<double> lifetime(3);
      std::vector<double> r_I_pi(3), G_e_pi(3), G_mu_pi(3), e_fac_pi(3), mu_fac_pi(3);
      std::vector<double> r_I_K(3), G_e_K(3), G_mu_K(3), e_fac_K(3), mu_fac_K(3);
      std::vector<double> e_fac_tau(3), mu_fac_tau(3);
      std::vector<double> m_temp_GERDA(3), m_temp_Kam(3);
      double V_ud_exp[7], f[7];
      double e_f_pi, mu_f_pi, e_f_K, mu_f_K, e_f_tau, mu_f_tau;
      double d_r_pi, d_r_K, d_r_tau, R_pi, R_K, R_tau;
      double m_GERDA, m_Kam;
      double total_err, V_ud_sq, chi2;

      M_I << *Param["M_1"], 0.0, 0.0,
             0.0, *Param["M_2"], 0.0,
             0.0, 0.0, *Param["M_3"];

      t = *Dep::SeesawI_Theta;
      t_sq = t.cwiseAbs2();

      // BBN constraint: lifetime must be less than 0.1s
      result_temp_bbn << 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0;
      for (int i=0; i<3; i++)
      {
        lifetime[i] = (96*pow(pi,3.0)*1e-9*conv_fact) / (G_F_sq*pow(M_I(i,i),5.0))*( ((1 + g_L_twid_sq + g_R_sq)*(t_sq(1,i) + t_sq(2,i))) + ((1 + g_L_sq + g_R_sq)*t_sq(0,i)) );
        if(lifetime[i]<0.1)
        {
          result_temp_bbn(0,i) = t_sq(0,i);
          result_temp_bbn(1,i) = t_sq(1,i);
          result_temp_bbn(2,i) = t_sq(2,i);
        }
      }

      // Lepton universality constraint: R_eu_pi/R_eu_K should be within experimental limits (TODO: implement 3-sigma limits as likelihood?)
      result_temp_lepuniv << 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0;
      e_f_pi = 0.0;
      mu_f_pi = 0.0;
      e_f_K = 0.0;
      mu_f_K = 0.0;
      e_f_tau = 0.0;
      mu_f_tau = 0.0;
      for (int j=0; j<3; j++)
      {
        if(M_I(j,j)<m_K)
        {
          if(M_I(j,j)<m_pi)
          {
            r_I_pi[j] = pow(M_I(j,j), 2.0)/pow(m_pi, 2.0);
            G_e_pi[j] = (r_e_pi + r_I_pi[j] - pow((r_e_pi - r_I_pi[j]), 2.0) * sqrt(1.0 - (2.0*pow((r_e_pi + r_I_pi[j]), 2.0)) + pow((r_e_pi - r_I_pi[j]), 2.0))) / (r_e_pi * pow((1.0 - r_e_pi), 2.0));
            G_mu_pi[j] = (r_mu_pi + r_I_pi[j] - pow((r_mu_pi - r_I_pi[j]), 2.0) * sqrt(1.0 - (2.0*pow((r_mu_pi + r_I_pi[j]), 2.0)) + pow((r_mu_pi - r_I_pi[j]), 2.0))) / (r_mu_pi * pow((1.0 - r_mu_pi), 2.0));
            e_fac_pi[j] = result_temp_bbn(0,j) * (G_e_pi[j] - 1.0);
            mu_fac_pi[j] = result_temp_bbn(1,j) * (G_mu_pi[j] - 1.0);
          }
          else
          {
            r_I_K[j] = pow(M_I(j,j), 2.0)/pow(m_K, 2.0);
            G_e_K[j] = (r_e_K + r_I_K[j] - pow((r_e_K - r_I_K[j]), 2.0) * sqrt(1.0 - (2.0*pow((r_e_K + r_I_K[j]), 2.0)) + pow((r_e_K - r_I_K[j]), 2.0))) / (r_e_pi * pow((1.0 - r_e_pi), 2.0));
            G_mu_K[j] = (r_mu_K + r_I_K[j] - pow((r_mu_K - r_I_K[j]), 2.0) * sqrt(1.0 - (2.0*pow((r_mu_K + r_I_K[j]), 2.0)) + pow((r_mu_K - r_I_K[j]), 2.0))) / (r_mu_K * pow((1.0 - r_mu_K), 2.0));
            e_fac_K[j] = result_temp_bbn(0,j) * (G_e_pi[j] - 1.0);
            mu_fac_K[j] = result_temp_bbn(1,j) * (G_mu_pi[j] - 1.0);
          }
        } 
        else if(M_I(j,j)>m_tau)
        {
          e_fac_tau[j] = result_temp_bbn(0,j);
          mu_fac_tau[j] = result_temp_bbn(1,j);
        }
        else
        {
          for (int k=0; k<3; k++)
          {
            result_temp_lepuniv(0,k) = result_temp_bbn(0,k);
            result_temp_lepuniv(1,k) = result_temp_bbn(1,k);
            result_temp_lepuniv(2,k) = result_temp_bbn(2,k);
          }
        }
        e_f_pi += e_fac_pi[j];
        mu_f_pi += mu_fac_pi[j];
        e_f_K += e_fac_K[j];
        mu_f_K += mu_fac_K[j];
        e_f_tau += e_fac_tau[j];
        mu_f_tau += mu_fac_tau[j];
      }
      d_r_pi = ((1.0 + e_f_pi)/(1.0 + mu_f_pi)) - 1.0;
      d_r_K = ((1.0 + e_f_K)/(1.0 + mu_f_K)) - 1.0;
      d_r_tau = ((1.0 - e_f_tau)/(1.0 - mu_f_tau)) - 1.0;
      R_pi = R_pi_SM * (1.0 + d_r_pi);
      R_K = R_K_SM * (1.0 + d_r_K);
      R_tau = R_tau_SM * (1.0 + d_r_tau);
      if (((1.226e-4<R_pi) && (R_pi<1.234e-4)) || ((2.478e-5<R_K) && (R_K<2.498e-5)) || ((0.9734<R_tau) && (R_tau<0.9794)))
      {
        for (int k=0; k<3; k++)
        {
          result_temp_lepuniv(0,k) = result_temp_bbn(0,k);
          result_temp_lepuniv(1,k) = result_temp_bbn(1,k);
          result_temp_lepuniv(2,k) = result_temp_bbn(2,k);
        }
      }

      // Neutrinoless double-beta decay constraint: m_bb should be less than the experimentally determined limits in GERDA and KamLAND-Zen
      result_temp_0nubb << 0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0,
                           0.0, 0.0, 0.0;
      m_light = *Dep::m_nu;
      U_light = *Dep::UPMNS;
      U_light_sq = U_light.cwiseAbs2();
      m_GERDA = 0.0;
      m_Kam = 0.0;
      for (int p=0; p<3; p++)
      {
        m_temp_GERDA[p] = U_light_sq(0,p)*abs(m_light(0,p)) + result_temp_lepuniv(0,p)*M_I(p,p)*(pow(*Param["L_Ge"], 2.0)/(pow(*Param["L_Ge"], 2.0)+pow(M_I(p,p), 2.0)));
        m_GERDA +=m_temp_GERDA[p];
        m_temp_Kam[p] = U_light_sq(0,p)*abs(m_light(0,p)) + result_temp_lepuniv(0,p)*M_I(p,p)*(pow(*Param["L_Xe"], 2.0)/(pow(*Param["L_Xe"], 2.0)+pow(M_I(p,p), 2.0)));
        m_Kam +=m_temp_Kam[p];
      }
      if ((m_GERDA < m_bb_GERDA) && (m_Kam < m_bb_Kam))
      {
        for (int q=0; q<3; q++)
        {
          result_temp_0nubb(0,q) = result_temp_lepuniv(0,q);
          result_temp_0nubb(1,q) = result_temp_lepuniv(1,q);
          result_temp_0nubb(2,q) = result_temp_lepuniv(2,q);

        }
      }

      // CKM unitarity constraint: V_ud should lie within 3sigma of the world average (TODO: implement as likelihood)
      result_temp_ckm << 0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0,
                         0.0, 0.0, 0.0;
      static double V_us_exp[7] = {0.2235,0.2227,0.2244,0.2238,0.2242,0.2262,0.2214};
      static double err[7] = {0.0006,0.0013,0.0008,0.0006,0.0011,0.0013,0.0022};
      total_err = 0.0;
      for (int a=0; a<7; a++)
      {
        total_err += pow(err[a], 2.0);
      }
      for (int i=0;i<7;i++)
      {
        V_ud_exp[i] = sqrt(1 - pow(V_us_exp[i], 2.0));
      }
      f[0] = (G_F_sq/G_mu_sq)*(1 - result_temp_0nubb(0,0) - result_temp_0nubb(0,1) - result_temp_0nubb(0,2));
      f[3] = (G_F_sq/G_mu_sq)*(1 - result_temp_0nubb(1,0) - result_temp_0nubb(1,1) - result_temp_0nubb(1,2));
      f[5] = 1 + result_temp_0nubb(1,0) + result_temp_0nubb(1,1) + result_temp_0nubb(1,2);
      f[6] = 1 + result_temp_0nubb(0,0) + result_temp_0nubb(0,1) + result_temp_0nubb(0,2) + result_temp_0nubb(1,0) + result_temp_0nubb(1,1) + result_temp_0nubb(1,2) + result_temp_0nubb(2,0) + result_temp_0nubb(2,1) + result_temp_0nubb(2,2);
      f[1] = f[0];
      f[2] = f[0];
      f[4] = f[3];
      V_ud_sq = 0.0;
      for (int j=0; j<7; j++)
      {
        V_ud_sq += (pow(V_ud_exp[j], 2.0)/(pow(err[j], 2.0)*f[j])) / total_err;
      }
      chi2 = 0.0;
      for (int k=0; k<7; k++)
      {
        chi2 += pow((V_ud_exp[k] - V_ud_sq), 2.0) / pow(err[k], 2.0);
      }

      result = t_sq;
    }

    void lnL(double& lnLike)
    {
      using namespace Pipes::lnL;
      lnLike = *Dep::lnLpienu + *Dep::lnLps191e + *Dep::lnLps191mu + *Dep::lnLcharme + *Dep::lnLcharmmu + *Dep::lnLdelphi + *Dep::lnLatlase + *Dep::lnLatlasmu + *Dep::lnLe949 + *Dep::lnLnutev + *Dep::lnLtau;
    }

    // Likelihood contribution from PIENU; searched for extra peaks in the spectrum of pi -> mu + nu. Constrains |U_ei|^2 in the mass range 60-129 MeV. [Phys. Rev. D, 84(5), 2011]
    void lnL_pienu(double& result)
    {
      using namespace Pipes::lnL_pienu;
      double array[140][2], M_1, M_2, M_3;
      std::vector<double> M_temp(140), U_temp(140), U(3), mixing_sq(3);

      mixing_sq[0] = *Dep::Ue1;
      mixing_sq[1] = *Dep::Ue2;
      mixing_sq[2] = *Dep::Ue3;

      std::ifstream f("DarkBit/data/pienu.csv");
      for (int row=0; row<140; ++row)
      {
        std::string line;
        getline(f, line);
        if (!f.good())
          break;
        std::stringstream iss(line);
        for (int col=0; col<2; ++col)
        {
          std::string val;
          getline(iss, val, ',');
          if (!iss)
            break;
          std::stringstream conv(val);
          conv >> array[row][col];
        }
      }

      for (int i=0; i<140; i++)
      {
        M_temp[i] = array[i][0];
        U_temp[i] = array[i][1];
      }
      tk::spline s;
      s.set_points(M_temp, U_temp);

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);
      result = -2.44*((mixing_sq[0]/U[0]) + (mixing_sq[1]/U[1]) + (mixing_sq[2]/U[2]));
    }

    // Likelihood contribution from PS191, electron sector; looked for charged tracks originating from RHN decays: nu_r -> l(-) + l(+) + nu / l + pi / e + pi(+) + pi(0). Constrains |U_ei|^2 in the mass range 20-450 MeV. Function also incorporates a later re-interpretation of the data to account for neutral current interaction (ignored in original) as well as the RHNs' Majorana nature. [Original: Phys. Lett. B, 203(3):332-334, 1988][Re-interp.: JHEP, 2012(6):1-27]
    void lnL_ps191_e(double& result)
    {
      using namespace Pipes::lnL_ps191_e;
      double array[116][2], M_1, M_2, M_3;
      std::vector<double> M_temp(116), U_temp(116), U(3), mixing_sq(3);
      double c_e = 0.5711;
      double c_mu = 0.1265;
      double c_tau = 0.1265;

      mixing_sq[0] = *Dep::Ue1 * ((c_e * *Dep::Ue1) + (c_mu * *Dep::Um1) + (c_tau * *Dep::Ut1));
      mixing_sq[1] = *Dep::Ue2 * ((c_e * *Dep::Ue2) + (c_mu * *Dep::Um2) + (c_tau * *Dep::Ut2));
      mixing_sq[2] = *Dep::Ue3 * ((c_e * *Dep::Ue3) + (c_mu * *Dep::Um3) + (c_tau * *Dep::Ut3));

      std::ifstream f("DarkBit/data/ps191_e.csv");
      for (int row=0; row<116; ++row)
      {
        std::string line;
        getline(f, line);
        if (!f.good())
          break;
        std::stringstream iss(line);
        for (int col=0; col<2; ++col)
        {
          std::string val;
          getline(iss, val, ',');
          if (!iss)
            break;
          std::stringstream conv(val);
          conv >> array[row][col];
        }
      }

      for (int i=0; i<116; i++)
      {
        M_temp[i] = array[i][0];
        U_temp[i] = array[i][1];
      }
      tk::spline s;
      s.set_points(M_temp, U_temp);

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);
      result = -2.44*((mixing_sq[0]/pow(U[0], 2.0)) + (mixing_sq[1]/pow(U[1], 2.0)) + (mixing_sq[2]/pow(U[2], 2.0)));
    }

    // Likelihood contribution from PS191, muon sector. Constrains |U_(mu,i)|^2 in the mass range 20-450 MeV. Description & references above.
    void lnL_ps191_mu(double& result)
    {
      using namespace Pipes::lnL_ps191_mu;
      double array[102][2], M_1, M_2, M_3;
      std::vector<double> M_temp(102), U_temp(102), U(3), mixing_sq(3);
      double c_e = 0.5711;
      double c_mu = 0.1265;
      double c_tau = 0.1265;

      mixing_sq[0] = *Dep::Um1 * ((c_e * *Dep::Ue1) + (c_mu * *Dep::Um1) + (c_tau * *Dep::Ut1));
      mixing_sq[1] = *Dep::Um2 * ((c_e * *Dep::Ue2) + (c_mu * *Dep::Um2) + (c_tau * *Dep::Ut2));
      mixing_sq[2] = *Dep::Um3 * ((c_e * *Dep::Ue3) + (c_mu * *Dep::Um3) + (c_tau * *Dep::Ut3));

      std::ifstream f("DarkBit/data/ps191_mu.csv");
      for (int row=0; row<102; ++row)
      {
        std::string line;
        getline(f, line);
        if (!f.good())
          break;
        std::stringstream iss(line);
        for (int col=0; col<2; ++col)
        {
          std::string val;
          getline(iss, val, ',');
          if (!iss)
            break;
          std::stringstream conv(val);
          conv >> array[row][col];
        }
      }

      for (int i=0; i<102; i++)
      {
        M_temp[i] = array[i][0];
        U_temp[i] = array[i][1];
      }
      tk::spline s;
      s.set_points(M_temp, U_temp);

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);
      result = -2.44*((mixing_sq[0]/pow(U[0], 2.0)) + (mixing_sq[1]/pow(U[1], 2.0)) + (mixing_sq[2]/pow(U[2], 2.0)));
    }

    // Likelihood contribution from CHARM, electron sector; searched for charged and neutral current decays of RHNs. Constrains |U_ei|^2 in the mass range 0.5-2.8 GeV. [Phys. Lett. B, 166(4):473-478, 1986]
    void lnL_charm_e(double& result)
    {
      using namespace Pipes::lnL_charm_e;
      double array[56][2], M_1, M_2, M_3;
      std::vector<double> M_temp(56), U_temp(56), U(3), mixing_sq(3);
      double c_e = 0.5711;
      double c_mu = 0.1265;
      double c_tau = 0.1265;

      mixing_sq[0] = *Dep::Ue1 * ((c_e * *Dep::Ue1) + (c_mu * *Dep::Um1) + (c_tau * *Dep::Ut1));
      mixing_sq[1] = *Dep::Ue2 * ((c_e * *Dep::Ue2) + (c_mu * *Dep::Um2) + (c_tau * *Dep::Ut2));
      mixing_sq[2] = *Dep::Ue3 * ((c_e * *Dep::Ue3) + (c_mu * *Dep::Um3) + (c_tau * *Dep::Ut3));

      std::ifstream f("DarkBit/data/charm_e.csv");
      for (int row=0; row<56; ++row)
      {
        std::string line;
        getline(f, line);
        if (!f.good())
          break;
        std::stringstream iss(line);
        for (int col=0; col<2; ++col)
        {
          std::string val;
          getline(iss, val, ',');
          if (!iss)
            break;
          std::stringstream conv(val);
          conv >> array[row][col];
        }
      }

      for (int i=0; i<56; i++)
      {
        M_temp[i] = array[i][0];
        U_temp[i] = array[i][1];
      }
      tk::spline s;
      s.set_points(M_temp, U_temp);

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);
      result = -2.44*((mixing_sq[0]/pow(U[0], 2.0)) + (mixing_sq[1]/pow(U[1], 2.0)) + (mixing_sq[2]/pow(U[2], 2.0)));
    }

    // Likelihood contribution from CHARM, muon sector. Constrains |U_(mu,i)|^2 in the mass range 0.5-2.8 GeV. Description & references above.
    void lnL_charm_mu(double& result)
    {
      using namespace Pipes::lnL_charm_mu;
      double array[34][2], M_1, M_2, M_3;
      std::vector<double> M_temp(34), U_temp(34), U(3), mixing_sq(3);
      double c_e = 0.5711;
      double c_mu = 0.1265;
      double c_tau = 0.1265;

      mixing_sq[0] = *Dep::Um1 * ((c_e * *Dep::Ue1) + (c_mu * *Dep::Um1) + (c_tau * *Dep::Ut1));
      mixing_sq[1] = *Dep::Um2 * ((c_e * *Dep::Ue2) + (c_mu * *Dep::Um2) + (c_tau * *Dep::Ut2));
      mixing_sq[2] = *Dep::Um3 * ((c_e * *Dep::Ue3) + (c_mu * *Dep::Um3) + (c_tau * *Dep::Ut3));

      std::ifstream f("DarkBit/data/charm_mu.csv");
      for (int row=0; row<34; ++row)
      {
        std::string line;
        getline(f, line);
        if (!f.good())
          break;
        std::stringstream iss(line);
        for (int col=0; col<2; ++col)
        {
          std::string val;
          getline(iss, val, ',');
          if (!iss)
            break;
          std::stringstream conv(val);
          conv >> array[row][col];
        }
      }

      for (int i=0; i<34; i++)
      {
        M_temp[i] = array[i][0];
        U_temp[i] = array[i][1];
      }
      tk::spline s;
      s.set_points(M_temp, U_temp);

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);
      result = -2.44*((mixing_sq[0]/pow(U[0], 2.0)) + (mixing_sq[1]/pow(U[1], 2.0)) + (mixing_sq[2]/pow(U[2], 2.0)));
    }

    // Likelihood contribution from DELPHI; searched for charged and neutral current decays of RHNs. Constrains |U_ei|^2, |U_(mu,i)|^2 as well as |U_(tau,i)|^2 in the mass range 3.5-50 GeV. [Z. Phys. C, 74(1):57-71, 1997]
    void lnL_delphi(double& result)
    {
      using namespace Pipes::lnL_delphi;
      double array[180][2], M_1, M_2, M_3;
      std::vector<double> M_temp(180), U_temp(180), U(3), mixing_sq(9);

      mixing_sq[0] = *Dep::Ue1;
      mixing_sq[1] = *Dep::Ue2;
      mixing_sq[2] = *Dep::Ue3;
      mixing_sq[3] = *Dep::Um1;
      mixing_sq[4] = *Dep::Um2;
      mixing_sq[5] = *Dep::Um3;
      mixing_sq[6] = *Dep::Ut1;
      mixing_sq[7] = *Dep::Ut2;
      mixing_sq[8] = *Dep::Ut3;

      std::ifstream f("DarkBit/data/delphi.csv");
      for (int row=0; row<180; ++row)
      {
        std::string line;
        getline(f, line);
        if (!f.good())
          break;
        std::stringstream iss(line);
        for (int col=0; col<2; ++col)
        {
          std::string val;
          getline(iss, val, ',');
          if (!iss)
            break;
          std::stringstream conv(val);
          conv >> array[row][col];
        }
      }

      for (int i=0; i<180; i++)
      {
        M_temp[i] = array[i][0];
        U_temp[i] = array[i][1];
      }
      tk::spline s;
      s.set_points(M_temp, U_temp);

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);
      result = -3.09*((mixing_sq[0]/U[0]) + (mixing_sq[1]/U[1]) + (mixing_sq[2]/U[2]) + (mixing_sq[3]/U[0]) + (mixing_sq[4]/U[1]) + (mixing_sq[5]/U[2]) + (mixing_sq[6]/U[0]) + (mixing_sq[7]/U[1]) + (mixing_sq[8]/U[2]));
    }

    // Likelihood contribution from ATLAS, electron sector; looked at the production and decay chain: pp -> W*(+-) -> l(+-) + nu_r. nu_r then decays into an on-shell W and a lepton; the W decays primarily into a qq pair. Constrains |U_ei|^2 in the mass range 50-500 GeV. [JHEP, 07:162, 2015]
    void lnL_atlas_e(double& result)
    {
      using namespace Pipes::lnL_atlas_e;
      double array[87][2], M_1, M_2, M_3;
      std::vector<double> M_temp(87), U_temp(87), U(3), mixing_sq(3);

      mixing_sq[0] = *Dep::Ue1;
      mixing_sq[1] = *Dep::Ue2;
      mixing_sq[2] = *Dep::Ue3;

      std::ifstream f("DarkBit/data/atlas_e.csv");
      for (int row=0; row<87; ++row)
      {
        std::string line;
        getline(f, line);
        if (!f.good())
          break;
        std::stringstream iss(line);
        for (int col=0; col<2; ++col)
        {
          std::string val;
          getline(iss, val, ',');
          if (!iss)
            break;
          std::stringstream conv(val);
          conv >> array[row][col];
        }
      }

      for (int i=0; i<87; i++)
      {
        M_temp[i] = array[i][0];
        U_temp[i] = array[i][1];
      }
      tk::spline s;
      s.set_points(M_temp, U_temp);

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);
      result = -3.09*(pow((mixing_sq[0]/U[0]), 2.0) + pow((mixing_sq[1]/U[1]), 2.0) + pow((mixing_sq[2]/U[2]), 2.0));
    }

    // Likelihood contribution from ATLAS, muon sector. Constrains |U_(mu,i)|^2 in the mass range 50-500 GeV. Description & references above.
    void lnL_atlas_mu(double& result)
    {
      using namespace Pipes::lnL_atlas_mu;
      double array[87][2], M_1, M_2, M_3;
      std::vector<double> M_temp(87), U_temp(87), U(3), mixing_sq(3);

      mixing_sq[0] = *Dep::Um1;
      mixing_sq[1] = *Dep::Um2;
      mixing_sq[2] = *Dep::Um3;

      std::ifstream f("DarkBit/data/atlas_mu.csv");
      for (int row=0; row<87; ++row)
      {
        std::string line;
        getline(f, line);
        if (!f.good())
          break;
        std::stringstream iss(line);
        for (int col=0; col<2; ++col)
        {
          std::string val;
          getline(iss, val, ',');
          if (!iss)
            break;
          std::stringstream conv(val);
          conv >> array[row][col];
        }
      }

      for (int i=0; i<87; i++)
      {
        M_temp[i] = array[i][0];
        U_temp[i] = array[i][1];
      }
      tk::spline s;
      s.set_points(M_temp, U_temp);

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);
      result = -3.09*(pow((mixing_sq[0]/U[0]), 2.0) + pow((mixing_sq[1]/U[1]), 2.0) + pow((mixing_sq[2]/U[2]), 2.0));
    }

    // Likelihood contribution from E949; used the kaon decay: K(+) -> mu(+) + nu_r. Constrains |U_(mu,i)|^2 in the mass range 175-300 MeV. [Phys. Rev. D, 91, 2015]
    void lnL_e949(double& result)
    {
      using namespace Pipes::lnL_e949;
      double array[112][2], M_1, M_2, M_3;
      std::vector<double> M_temp(112), U_temp(112), U(3), mixing_sq(3);

      mixing_sq[0] = *Dep::Um1;
      mixing_sq[1] = *Dep::Um2;
      mixing_sq[2] = *Dep::Um3;

      std::ifstream f("DarkBit/data/e949.csv");
      for (int row=0; row<112; ++row)
      {
        std::string line;
        getline(f, line);
        if (!f.good())
          break;
        std::stringstream iss(line);
        for (int col=0; col<2; ++col)
        {
          std::string val;
          getline(iss, val, ',');
          if (!iss)
            break;
          std::stringstream conv(val);
          conv >> array[row][col];
        }
      }

      for (int i=0; i<112; i++)
      {
        M_temp[i] = array[i][0];
        U_temp[i] = array[i][1];
      }
      tk::spline s;
      s.set_points(M_temp, U_temp);

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);
      result = -2.44*((mixing_sq[0]/U[0]) + (mixing_sq[1]/U[1]) + (mixing_sq[2]/U[2]));
    }

    // Likelihood contribution from NuTeV; used RHN decays into muonic final states (mu + mu + nu / mu + e + nu / mu + pi / mu + rho). Constrains |U_(mu,i)|^2 in the mass range 0.25-2 GeV. [Phys. Rev. Lett., 83:4943-4946, 1999]
    void lnL_nutev(double& result)
    {
      using namespace Pipes::lnL_nutev;
      double array[249][2], M_1, M_2, M_3;
      std::vector<double> M_temp(249), U_temp(249), U(3), mixing_sq(3);

      mixing_sq[0] = *Dep::Um1;
      mixing_sq[1] = *Dep::Um2;
      mixing_sq[2] = *Dep::Um3;

      std::ifstream f("DarkBit/data/nutev.csv");
      for (int row=0; row<249; ++row)
      {
        std::string line;
        getline(f, line);
        if (!f.good())
          break;
        std::stringstream iss(line);
        for (int col=0; col<2; ++col)
        {
          std::string val;
          getline(iss, val, ',');
          if (!iss)
            break;
          std::stringstream conv(val);
          conv >> array[row][col];
        }
      }

      for (int i=0; i<249; i++)
      {
        M_temp[i] = array[i][0];
        U_temp[i] = array[i][1];
      }
      tk::spline s;
      s.set_points(M_temp, U_temp);

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);
      result = -2.44*((mixing_sq[0]/U[0]) + (mixing_sq[1]/U[1]) + (mixing_sq[2]/U[2]));
    }

    // Likelihood contribution from a re-interpretation of CHARM data; assumes tau mixing is dominant. Constrains |U_(tau,i)|^2 in the mass range 10-290 MeV. [Phys. Lett. B, 550(1-2):8-15, 2002]
    void lnL_tau(double& result)
    {
      using namespace Pipes::lnL_tau;
      double array[172][2], M_1, M_2, M_3;
      std::vector<double> M_temp(172), U_temp(172), U(3), mixing_sq(3);

      mixing_sq[0] = *Dep::Ut1;
      mixing_sq[1] = *Dep::Ut2;
      mixing_sq[2] = *Dep::Ut3;

      std::ifstream f("DarkBit/data/tau.csv");
      for (int row=0; row<172; ++row)
      {
        std::string line;
        getline(f, line);
        if (!f.good())
          break;
        std::stringstream iss(line);
        for (int col=0; col<2; ++col)
        {
          std::string val;
          getline(iss, val, ',');
          if (!iss)
            break;
          std::stringstream conv(val);
          conv >> array[row][col];
        }
      }

      for (int i=0; i<172; i++)
      {
        M_temp[i] = array[i][0];
        U_temp[i] = array[i][1];
      }
      tk::spline s;
      s.set_points(M_temp, U_temp);

      M_1 = *Param["M_1"];
      M_2 = *Param["M_2"];
      M_3 = *Param["M_3"];
      U[0] = s(M_1);
      U[1] = s(M_2);
      U[2] = s(M_3);
      result = -2.44*((mixing_sq[0]/U[0]) + (mixing_sq[1]/U[1]) + (mixing_sq[2]/U[2]));
    }
                                                                                                                                                                                                                 void printable_Ue1(double& Ue1_sq)
    {
      namespace myPipe2 = Pipes::printable_Ue1;
      Matrix3d t_1(*myPipe2::Dep::SN_stuff);
      Ue1_sq = t_1(0,0);
    }

    void printable_Um1(double& Um1_sq)
    {
      namespace myPipe3 = Pipes::printable_Um1;
      Matrix3d t_2(*myPipe3::Dep::SN_stuff);
      Um1_sq = t_2(1,0);
    }

    void printable_Ut1(double& Ut1_sq)
    {
      namespace myPipe4 = Pipes::printable_Ut1;
      Matrix3d t_3(*myPipe4::Dep::SN_stuff);
      Ut1_sq = t_3(2,0);
    }

    void printable_Ue2(double& Ue2_sq)
    {
      namespace myPipe5 = Pipes::printable_Ue2;
      Matrix3d t_4(*myPipe5::Dep::SN_stuff);
      Ue2_sq = t_4(0,1);
    }

    void printable_Um2(double& Um2_sq)
    {
      namespace myPipe6 = Pipes::printable_Um2;
      Matrix3d t_5(*myPipe6::Dep::SN_stuff);
      Um2_sq = t_5(1,1);
    }

    void printable_Ut2(double& Ut2_sq)
    {
      namespace myPipe7 = Pipes::printable_Ut2;
      Matrix3d t_6(*myPipe7::Dep::SN_stuff);
      Ut2_sq = t_6(2,1);
    }

    void printable_Ue3(double& Ue3_sq)
    {
      namespace myPipe8 = Pipes::printable_Ue3;
      Matrix3d t_7(*myPipe8::Dep::SN_stuff);
      Ue3_sq = t_7(0,2);
    }

    void printable_Um3(double& Um3_sq)
    {
      namespace myPipe9 = Pipes::printable_Um3;
      Matrix3d t_8(*myPipe9::Dep::SN_stuff);
      Um3_sq = t_8(1,2);
    }

    void printable_Ut3(double& Ut3_sq)
    {
      namespace myPipe10 = Pipes::printable_Ut3;
      Matrix3d t_9(*myPipe10::Dep::SN_stuff);
      Ut3_sq = t_9(2,2);
    }

  }

}
