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
#include "gambit/Elements/numerical_constants.hpp"
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
    double l_M(double M)
    {
      const double m_Z = 91.1876;  // GeV
      const double m_H = 125.09;  // GeV
      return 1.0/pow(4.0*pi, 2.0) * ( (3.0*log(pow(M/m_Z, 2.0)))/((pow(M/m_Z, 2.0)) - 1.0) + (log(pow(M/m_H, 2.0)))/((pow(M/m_H, 2.0)) - 1.0));
    }

    void CI_param(Matrix3d& result)   
    {
      using namespace Pipes::CI_param;
      typedef std::complex<double> dcomp;
      static double v = 246.0;  // GeV
      static double md21 = 7.5e-23;  // GeV^2
      static double md31 = 2.457e-21;  // GeV^2
      static double md23 = 2.449e-21;  // GeV^2
      static double theta23 = 0.7382;  // rad
      static double theta13 = 0.1483;  // rad
      static double theta12 = 0.5843;  // rad
      static double conv_fact = 6.58e-16;  // conversion factor from eV^-1 to s, for lifetime
      static double G_F_sq = 1.3604e-10;  // GeV^-4
      static double g_L_twid_sq = 0.0771;  // g_L_twid = -0.5 + s_W_sq
      static double g_R_sq = 0.0494;  // g_R = s_W^2
      static double g_L_sq = 0.5217;  // g_L = 0.5 + s_W^2
      dcomp I(0.0, 1.0);
      Matrix3d M_I;  // not complex; circumvents type mismatch in l(M)
      Matrix3cd M_twid_temp;
      Matrix3cd M_twid;
      Matrix3cd R_23;
      Matrix3cd R_13;
      Matrix3cd R_12;
      Matrix3cd R;
      Matrix3cd m_nu;
      Matrix3cd V_23;
      Matrix3cd V_13;
      Matrix3cd V_12;
      Matrix3cd U_pd;
      Matrix3cd U_nd;
      Matrix3cd Maj_phase;
      Matrix3cd U_nu;
      Matrix3cd t;
      Matrix3d t_sq;
      Matrix3d result_temp;
      std::vector<double> lifetime(3);
      double x23 = *Param["ReOm23"];
      double y23 = *Param["ImOm23"];
      double x13 = *Param["ReOm13"];
      double y13 = *Param["ImOm13"];
      double x12 = *Param["ReOm12"];
      double y12 = *Param["ImOm12"];
      double a1 = *Param["alpha1"];
      double a2 = *Param["alpha2"];
      double d = *Param["delta"];
      int o = *Param["ordering"];
      int m_min = *Param["min_mass"];

      M_I << *Param["M_1"], 0.0, 0.0,
             0.0, *Param["M_2"], 0.0,
             0.0, 0.0, *Param["M_3"];
      M_twid_temp(0,0) = M_I(0,0)  * (1.0 - (pow(M_I(0,0),2.0)*l_M(M_I(0,0))/pow(v,2.0)));
      M_twid_temp(0,1) = 0.0;
      M_twid_temp(0,2) = 0.0;
      M_twid_temp(1,0) = 0.0;
      M_twid_temp(1,1) = M_I(1,1)  * (1.0 - (pow(M_I(1,1),2.0)*l_M(M_I(1,1))/pow(v,2.0)));
      M_twid_temp(1,2) = 0.0;
      M_twid_temp(2,0) = 0.0;
      M_twid_temp(2,1) = 0.0;
      M_twid_temp(2,2) = M_I(2,2)  * (1.0 - (pow(M_I(2,2),2.0)*l_M(M_I(2,2))/pow(v,2.0)));
      M_twid = M_twid_temp.sqrt();

      R_23(0,0) = 1.0;
      R_23(0,1) = 0.0;
      R_23(0,2) = 0.0;
      R_23(1,0) = 0.0;
      R_23(1,1) = cos(x23)*cosh(y23) - I*sin(x23)*sinh(y23);
      R_23(1,2) = sin(x23)*cosh(y23) + I*cos(x23)*sinh(y23);
      R_23(2,0) = 0.0;
      R_23(2,1) = -sin(x23)*cosh(y23) - I*cos(x23)*sinh(y23);
      R_23(2,2) = cos(x23)*cosh(y23) - I*sin(x23)*sinh(y23);
      R_13(0,0) = cos(x13)*cosh(y13) - I*sin(x13)*sinh(y13);;
      R_13(0,1) = 0.0;
      R_13(0,2) = sin(x13)*cosh(y13) + I*cos(x13)*sinh(y13);
      R_13(1,0) = 0.0;
      R_13(1,1) = 1.0;
      R_13(1,2) = 0.0;
      R_13(2,0) = -sin(x13)*cosh(y13) - I*cos(x13)*sinh(y13);
      R_13(2,1) = 0.0;
      R_13(2,2) = cos(x13)*cosh(y13) - I*sin(x13)*sinh(y13);
      R_12(0,0) = cos(x12)*cosh(y12) - I*sin(x12)*sinh(y12);
      R_12(0,1) = sin(x12)*cosh(y12) + I*cos(x12)*sinh(y12);
      R_12(0,2) = 0.0;
      R_12(1,0) = -sin(x12)*cosh(y12) - I*cos(x12)*sinh(y12);
      R_12(1,1) = cos(x12)*cosh(y12) - I*sin(x12)*sinh(y12);
      R_12(1,2) = 0.0;
      R_12(2,0) = 0.0;
      R_12(2,1) = 0.0;
      R_12(2,2) = 1.0;
      R = R_23 * R_13 * R_12;

      m_nu(0,1) = 0.0;
      m_nu(0,2) = 0.0;
      m_nu(1,0) = 0.0;
      m_nu(1,2) = 0.0;
      m_nu(2,0) = 0.0;
      m_nu(2,1) = 0.0;
      if(o == 1)
      {
        if(m_min == 0)
          {
            m_nu(0,0) = 0.0;
            m_nu(1,1) = sqrt(md21);
            m_nu(2,2) = sqrt(md31);
          }
        else if(m_min == 1)
          {
            m_nu(0,0) = 2.3e-10;
            m_nu(1,1) = sqrt(pow(m_nu(0,0), 2.0) + md21);
            m_nu(2,2) = sqrt(pow(m_nu(0,0), 2.0) + md31);
          }
      }
      else if(o == 0)
      {
        if(m_min == 0)
          {
            m_nu(2,2) = 0.0;
            m_nu(1,1) = sqrt(md23);
            m_nu(0,0) = sqrt(pow(m_nu(1,1), 2.0) - md21);
          }
        else if(m_min == 0.23)
          {
            m_nu(2,2) = 2.3e-10;
            m_nu(1,1) = sqrt(pow(m_nu(2,2), 2.0) + md23);
            m_nu(0,0) = sqrt(pow(m_nu(1,1), 2.0) - md21);
          }
      }

      V_23 << 1.0, 0.0, 0.0,
              0.0, cos(theta23), sin(theta23),
              0.0, -sin(theta23), cos(theta23);
      V_13 << cos(theta13), 0.0, sin(theta13),
              0.0, 1.0, 0.0,
              -sin(theta13), 0.0, cos(theta13);
      V_12 << cos(theta12), sin(theta12), 0.0,
              -sin(theta12), cos(theta12), 0.0,
              0.0, 0.0, 1.0;
      U_pd << exp(-I*d/2.0), 0.0, 0.0,
              0.0, 1.0, 0.0,
              0.0, 1.0, exp(I*d/2.0);
      U_nd << exp(I*d/2.0), 0.0, 0.0,
              0.0, 1.0, 0.0,
              0.0, 1.0, exp(-I*d/2.0);
      Maj_phase << exp(I*a1/2.0), 0.0, 0.0,
                   0.0, exp(I*a2/2.0), 0.0,
                   0.0, 0.0, 1.0;
      U_nu = V_23 * U_pd * V_13 * U_nd* V_12 * Maj_phase;

      t = I * U_nu * m_nu.sqrt() * R * M_twid.inverse();
      t_sq = t.cwiseAbs2();
      result_temp << 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0;
      for (int i=0; i<3; i++)
      {
        lifetime[i] = (96*pow(pi,3.0)*1e-9*conv_fact) / (G_F_sq*pow(M_I(i,i),5.0))*( ((1 + g_L_twid_sq + g_R_sq)*(t_sq(1,i) + t_sq(2,i))) + ((1 + g_L_sq + g_R_sq)*t_sq(0,i)) );
        if(lifetime[i]<0.1)
        {
          result_temp(0,i) = t_sq(0,i);
          result_temp(1,i) = t_sq(1,i);
          result_temp(2,i) = t_sq(2,i);
        }
      }
      result = result_temp;
    }

    void lnL(double& lnLike)
    {
      namespace myPipe1 = Pipes::lnL;
      double array_nutev[249][2] = {};
      double array_delphi[180][2] = {};
      double array_pienu[140][2] = {};
      double array_ps191_e[116][2] = {};
      double array_ps191_mu[102][2] = {};
      double array_charm_e[56][2] = {};
      double array_charm_mu[34][2] = {};
      double array_atlas_e[87][2] = {};
      double array_atlas_mu[87][2] = {};
      double array_e949[112][2] = {};
      double array_tau[172][2] = {};
      double M_1, M_2, M_3;
      double c_e = 0.5711;
      double c_mu = 0.1265;
      double c_tau = 0.1265;
      std::vector<double> M_temp_nutev(249), U_temp_nutev(249), M_temp_delphi(180), U_temp_delphi(180), M_temp_pienu(140), U_temp_pienu(140), M_temp_ps191_e(116), U_temp_ps191_e(116), M_temp_ps191_mu(102), U_temp_ps191_mu(102), M_temp_charm_e(56), U_temp_charm_e(56), M_temp_charm_mu(34), U_temp_charm_mu(34), M_temp_atlas_e(87), U_temp_atlas_e(87), M_temp_atlas_mu(87), U_temp_atlas_mu(87), M_temp_e949(112), U_temp_e949(112), M_temp_tau(172), U_temp_tau(172);
      std::vector<double> U_pienu(3), U_ps191_e(3), U_ps191_mu(3), U_charm_e(3), U_charm_mu(3), U_delphi(3), U_atlas_e(3), U_atlas_mu(3), U_e949(3), U_nutev(3), U_tau(3);
      std::vector<double> mixing_ps191_e(3), mixing_ps191_mu(3), mixing_charm_e(3), mixing_charm_mu(3);

      std::ifstream f_pienu("DarkBit/data/pienu.csv");
      for (int row1=0; row1<140; ++row1)
      {
        std::string line1;
        getline(f_pienu, line1);
        if (!f_pienu.good())
          break;
        std::stringstream iss1(line1);
        for (int col1=0; col1<2; ++col1)
        {
          std::string val1;
          getline(iss1, val1, ',');
          if (!iss1)
            break;
          std::stringstream conv1(val1);
          conv1 >> array_pienu[row1][col1];
        }
      }
      std::ifstream f_ps191_e("DarkBit/data/ps191_e.csv");
      for (int row2=0; row2<116; ++row2)
      {
        std::string line2;
        getline(f_ps191_e, line2);
        if (!f_ps191_e.good())
          break;
        std::stringstream iss2(line2);
        for (int col2=0; col2<2; ++col2)
        {
          std::string val2;
          getline(iss2, val2, ',');
          if (!iss2)
            break;
          std::stringstream conv2(val2);
          conv2 >> array_ps191_e[row2][col2];
        }
      }
      std::ifstream f_ps191_mu("DarkBit/data/ps191_mu.csv");
      for (int row3=0; row3<102; ++row3)
      {
        std::string line3;
        getline(f_ps191_mu, line3);
        if (!f_ps191_mu.good())
          break;
        std::stringstream iss3(line3);
        for (int col3=0; col3<2; ++col3)
        {
          std::string val3;
          getline(iss3, val3, ',');
          if (!iss3)
            break;
          std::stringstream conv3(val3);
          conv3 >> array_ps191_mu[row3][col3];
        }
      }
      std::ifstream f_charm_e("DarkBit/data/charm_e.csv");
      for (int row4=0; row4<56; ++row4)
      {
        std::string line4;
        getline(f_charm_e, line4);
        if (!f_charm_e.good())
          break;
        std::stringstream iss4(line4);
        for (int col4=0; col4<2; ++col4)
        {
          std::string val4;
          getline(iss4, val4, ',');
          if (!iss4)
            break;
          std::stringstream conv4(val4);
          conv4 >> array_charm_e[row4][col4];
        }
      }
      std::ifstream f_charm_mu("DarkBit/data/charm_mu.csv");
      for (int row5=0; row5<34; ++row5)
      {
        std::string line5;
        getline(f_charm_mu, line5);
        if (!f_charm_mu.good())
          break;
        std::stringstream iss5(line5);
        for (int col5=0; col5<2; ++col5)
        {
          std::string val5;
          getline(iss5, val5, ',');
          if (!iss5)
            break;
          std::stringstream conv5(val5);
          conv5 >> array_charm_mu[row5][col5];
        }
      }
      std::ifstream f_delphi("DarkBit/data/delphi.csv");
      for (int row6=0; row6<180; ++row6)
      {
        std::string line6;
        getline(f_delphi, line6);
        if (!f_delphi.good())
          break;
        std::stringstream iss6(line6);
        for (int col6=0; col6<2; ++col6)
        {
          std::string val6;
          getline(iss6, val6, ',');
          if (!iss6)
            break;
          std::stringstream conv6(val6);
          conv6 >> array_delphi[row6][col6];
        }
      }
      std::ifstream f_atlas_e("DarkBit/data/atlas_e.csv");
      for (int row7=0; row7<87; ++row7)
      {
        std::string line7;
        getline(f_atlas_e, line7);
        if (!f_atlas_e.good())
          break;
        std::stringstream iss7(line7);
        for (int col7=0; col7<2; ++col7)
        {
          std::string val7;
          getline(iss7, val7, ',');
          if (!iss7)
            break;
          std::stringstream conv7(val7);
          conv7 >> array_atlas_e[row7][col7];
        }
      }
      std::ifstream f_atlas_mu("DarkBit/data/atlas_mu.csv");
      for (int row8=0; row8<87; ++row8)
      {
        std::string line8;
        getline(f_atlas_mu, line8);
        if (!f_atlas_mu.good())
          break;
        std::stringstream iss8(line8);
        for (int col8=0; col8<2; ++col8)
        {
          std::string val8;
          getline(iss8, val8, ',');
          if (!iss8)
            break;
          std::stringstream conv8(val8);
          conv8 >> array_atlas_mu[row8][col8];
        }
      }
      std::ifstream f_e949("DarkBit/data/e949.csv");
      for (int row9=0; row9<112; ++row9)
      {
        std::string line9;
        getline(f_e949, line9);
        if(!f_e949.good())
          break;
        std::stringstream iss9(line9);
        for (int col9=0; col9<2; ++col9)
        {
          std::string val9;
          getline(iss9, val9, ',');
          if (!iss9)
            break;
          std::stringstream conv9(val9);
          conv9 >> array_e949[row9][col9];
        }
      }
      std::ifstream f_nutev("DarkBit/data/nutev.csv");
      for (int row10=0; row10<249; ++row10)
      {
        std::string line10;
        getline(f_nutev, line10);
        if (!f_nutev.good())
          break;
        std::stringstream iss10(line10);
        for (int col10=0; col10<2; ++col10)
        {
          std::string val10;
          getline(iss10, val10, ',');
          if (!iss10)
            break;
          std::stringstream conv10(val10);
          conv10 >> array_nutev[row10][col10];
        }
      }
      std::ifstream f_tau("DarkBit/data/tau.csv");
      for (int row11=0; row11<172; ++row11)
      {
        std::string line11;
        getline(f_tau, line11);
        if(!f_tau.good())
          break;
        std::stringstream iss11(line11);
        for (int col11=0; col11<2; ++col11)
        {
          std::string val11;
          getline(iss11, val11, ',');
          if(!iss11)
            break;
          std::stringstream conv11(val11);
          conv11 >> array_tau[row11][col11];
        }
      }

      for (int i=0; i<140; i++)
      {
        M_temp_pienu[i] = array_pienu[i][0];
        U_temp_pienu[i] = array_pienu[i][1];
      }
      tk::spline s1;
      s1.set_points(M_temp_pienu, U_temp_pienu);
      for (int j=0; j<116; j++)
      {
        M_temp_ps191_e[j] = array_ps191_e[j][0];
        U_temp_ps191_e[j] = array_ps191_e[j][1];
      }
      tk::spline s2;
      s2.set_points(M_temp_ps191_e, U_temp_ps191_e);
      for (int k=0; k<102; k++)
      {
        M_temp_ps191_mu[k] = array_ps191_mu[k][0];
        U_temp_ps191_mu[k] = array_ps191_mu[k][1];
      }
      tk::spline s3;
      s3.set_points(M_temp_ps191_mu, U_temp_ps191_mu);
      for (int l=0; l<56; l++)
      {
        M_temp_charm_e[l] = array_charm_e[l][0];
        U_temp_charm_e[l] = array_charm_e[l][1];
      }
      tk::spline s4;
      s4.set_points(M_temp_charm_e, U_temp_charm_e);
      for (int n=0; n<34; n++)
      {
        M_temp_charm_mu[n] = array_charm_mu[n][0];
        U_temp_charm_mu[n] = array_charm_mu[n][1];
      }
      tk::spline s5;
      s5.set_points(M_temp_charm_mu, U_temp_charm_mu);
      for (int o=0; o<180; o++)
      {
        M_temp_delphi[o] = array_delphi[o][0];
        U_temp_delphi[o] = array_delphi[o][1];
      }
      tk::spline s6;
      s6.set_points(M_temp_delphi, U_temp_delphi);
      for (int p=0; p<87; p++)
      {
        M_temp_atlas_e[p] = array_atlas_e[p][0];
        U_temp_atlas_e[p] = array_atlas_e[p][1];
      }
      tk::spline s7;
      s7.set_points(M_temp_atlas_e, U_temp_atlas_e);
      for (int q=0; q<87; q++)
      {
        M_temp_atlas_mu[q] = array_atlas_mu[q][0];
        U_temp_atlas_mu[q] = array_atlas_mu[q][1];
      }
      tk::spline s8;
      s8.set_points(M_temp_atlas_mu, U_temp_atlas_mu);
      for (int r=0; r<112; r++)
      {
        M_temp_e949[r] = array_e949[r][0];
        U_temp_e949[r] = array_e949[r][1];
      }
      tk::spline s9;
      s9.set_points(M_temp_e949, U_temp_e949);
      for (int s=0; s<249; s++)
      {
        M_temp_nutev[s] = array_nutev[s][0];
        U_temp_nutev[s] = array_nutev[s][1];
      }
      tk::spline s10;
      s10.set_points(M_temp_nutev, U_temp_nutev);
      for (int t=0; t<172; t++)
      {
        M_temp_tau[t] = array_tau[t][0];
        U_temp_tau[t] = array_tau[t][1];
      }
      tk::spline s11;
      s11.set_points(M_temp_tau, U_temp_tau);

      std::vector<double> mixing_sq(9);
      Matrix3d m_sq(*myPipe1::Dep::SN_stuff);
      mixing_sq[0] = m_sq(0,0);
      mixing_sq[1] = m_sq(0,1);
      mixing_sq[2] = m_sq(0,2);
      mixing_sq[3] = m_sq(1,0);
      mixing_sq[4] = m_sq(1,1);
      mixing_sq[5] = m_sq(1,2);
      mixing_sq[6] = m_sq(2,0);
      mixing_sq[7] = m_sq(2,1);
      mixing_sq[8] = m_sq(2,2);
      M_1 = *myPipe1::Param["M_1"];
      M_2 = *myPipe1::Param["M_2"];
      M_3 = *myPipe1::Param["M_3"];
      U_pienu[0] = s1(M_1);
      U_pienu[1] = s1(M_2);
      U_pienu[2] = s1(M_3);
      U_ps191_e[0] = s2(M_1);
      U_ps191_e[1] = s2(M_2);
      U_ps191_e[2] = s2(M_3);
      U_ps191_mu[0] = s3(M_1);
      U_ps191_mu[1] = s3(M_2);
      U_ps191_mu[2] = s3(M_3);
      U_charm_e[0] = s4(M_1);
      U_charm_e[1] = s4(M_2);
      U_charm_e[2] = s4(M_3);
      U_charm_mu[0] = s5(M_1);
      U_charm_mu[1] = s5(M_2);
      U_charm_mu[2] = s5(M_3);
      U_delphi[0] = s6(M_1);
      U_delphi[1] = s6(M_2);
      U_delphi[2] = s6(M_3);
      U_atlas_e[0] = s7(M_1);
      U_atlas_e[1] = s7(M_2);
      U_atlas_e[2] = s7(M_3);
      U_atlas_mu[0] = s8(M_1);
      U_atlas_mu[1] = s8(M_2);
      U_atlas_mu[2] = s8(M_3);
      U_e949[0] = s9(M_1);
      U_e949[1] = s9(M_2);
      U_e949[2] = s9(M_3);
      U_nutev[0] = s10(M_1);
      U_nutev[1] = s10(M_2);
      U_nutev[2] = s10(M_3);
      U_tau[0] = s11(M_1);
      U_tau[1] = s11(M_2);
      U_tau[2] = s11(M_3);
//      mixing_pienu[i] = mixing_sq[i];
      mixing_ps191_e[0] = mixing_sq[0]*((c_e*mixing_sq[0])+(c_mu*mixing_sq[3])+(c_tau*mixing_sq[6]));
      mixing_ps191_e[1] = mixing_sq[1]*((c_e*mixing_sq[1])+(c_mu*mixing_sq[4])+(c_tau*mixing_sq[7]));
      mixing_ps191_e[2] = mixing_sq[2]*((c_e*mixing_sq[2])+(c_mu*mixing_sq[5])+(c_tau*mixing_sq[8]));
      mixing_ps191_mu[0] = mixing_sq[3]*((c_e*mixing_sq[0])+(c_mu*mixing_sq[3])+(c_tau*mixing_sq[6]));
      mixing_ps191_mu[1] = mixing_sq[4]*((c_e*mixing_sq[1])+(c_mu*mixing_sq[4])+(c_tau*mixing_sq[7]));
      mixing_ps191_mu[2] = mixing_sq[5]*((c_e*mixing_sq[2])+(c_mu*mixing_sq[5])+(c_tau*mixing_sq[8]));
      mixing_charm_e[0] = mixing_sq[0]*((c_e*mixing_sq[0])+(c_mu*mixing_sq[3])+(c_tau*mixing_sq[6]));
      mixing_charm_e[1] = mixing_sq[1]*((c_e*mixing_sq[1])+(c_mu*mixing_sq[4])+(c_tau*mixing_sq[7]));
      mixing_charm_e[2] = mixing_sq[2]*((c_e*mixing_sq[2])+(c_mu*mixing_sq[5])+(c_tau*mixing_sq[8]));
      mixing_charm_mu[0] = mixing_sq[3]*((c_e*mixing_sq[0])+(c_mu*mixing_sq[3])+(c_tau*mixing_sq[6]));
      mixing_charm_mu[1] = mixing_sq[4]*((c_e*mixing_sq[1])+(c_mu*mixing_sq[4])+(c_tau*mixing_sq[7]));
      mixing_charm_mu[2] = mixing_sq[5]*((c_e*mixing_sq[2])+(c_mu*mixing_sq[5])+(c_tau*mixing_sq[8]));
//      mixing_delphi[i] = mixing_sq[i]
//      mixing_nutev[i] = mixing_sq[i+3]
//      mixing_atlas_e[i] = mixing_sq[i]
//      mixing_atlas_mu[i] = mixing_sq[i+3]
//      mixing_e949[i] = mixing_sq[i+3]
//      mixing_tau[i] = mixing_sq[i+6]

      lnLike = -(2.44*mixing_sq[0])/U_pienu[0] -(2.44*mixing_sq[1])/U_pienu[1] -(2.44*mixing_sq[2])/U_pienu[2] -(2.44*mixing_ps191_e[0])/pow(U_ps191_e[0], 2.0) -(2.44*mixing_ps191_e[1])/pow(U_ps191_e[1], 2.0) -(2.44*mixing_ps191_e[2])/pow(U_ps191_e[2], 2.0) -(2.44*mixing_ps191_mu[0])/pow(U_ps191_mu[0], 2.0) -(2.44*mixing_ps191_mu[1])/pow(U_ps191_mu[1], 2.0) -(2.44*mixing_ps191_mu[2])/pow(U_ps191_mu[2], 2.0) -(2.44*mixing_charm_e[0])/pow(U_charm_e[0], 2.0) -(2.44*mixing_charm_e[1])/pow(U_charm_e[1], 2.0) -(2.44*mixing_charm_e[2])/pow(U_charm_e[2], 2.0) -(2.44*mixing_charm_mu[0])/pow(U_charm_mu[0], 2.0) -(2.44*mixing_charm_mu[1])/pow(U_charm_mu[1], 2.0) -(2.44*mixing_charm_mu[2])/pow(U_charm_mu[2], 2.0) -(3.09*mixing_sq[0])/U_delphi[0] -(3.09*mixing_sq[1])/U_delphi[1] -(3.09*mixing_sq[2])/U_delphi[2] - (3.09*mixing_sq[3])/U_delphi[0] -(3.09*mixing_sq[4])/U_delphi[1] -(3.09*mixing_sq[5])/U_delphi[2] -(3.09*mixing_sq[6])/U_delphi[0] -(3.09*mixing_sq[7])/U_delphi[1] -(3.09*mixing_sq[8])/U_delphi[2] -(3.09*mixing_sq[0]*mixing_sq[0])/pow(U_atlas_e[0], 2.0) -(3.09*mixing_sq[1]*mixing_sq[1])/pow(U_atlas_e[1], 2.0) -(3.09*mixing_sq[2]*mixing_sq[2])/pow(U_atlas_e[2], 2.0) -(3.09*mixing_sq[3]*mixing_sq[3])/pow(U_atlas_mu[0], 2.0) -(3.09*mixing_sq[4]*mixing_sq[4])/pow(U_atlas_e[1], 2.0) -(3.09*mixing_sq[5]*mixing_sq[5])/pow(U_atlas_e[2], 2.0) -(2.44*mixing_sq[3])/U_e949[0] -(2.44*mixing_sq[4])/U_e949[1] -(2.44*mixing_sq[5])/U_e949[2] -(2.44*mixing_sq[3])/U_nutev[0] -(2.44*mixing_sq[4])/U_nutev[1] -(2.44*mixing_sq[5])/U_nutev[2] - (2.44*mixing_sq[6])/U_tau[0] -(2.44*mixing_sq[7])/U_tau[1] -(2.44*mixing_sq[8])/U_tau[2];
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
