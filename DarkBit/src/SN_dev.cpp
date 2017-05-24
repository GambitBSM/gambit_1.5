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
      std::vector<double> temp1(3);
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
      result = t.cwiseAbs2();
    }

    void lnL(double& lnLike)
    {
      using namespace Pipes::lnL;
      static double m = 1e-7;
      static double s = 2e-18;
      double array_bebc[88][2] = {};
      double array_nutev[249][2] = {};
      double array_delphi[180][2] = {};
      double M, U_delphi;
    double U_bebc, U_nutev;
      std::vector<double> M_temp_bebc(88), U_temp_bebc(88), M_temp_nutev(249), U_temp_nutev(249;
      std::vector<double>  M_temp_delphi(180), U_temp_delphi(180);
      std::ifstream file("DarkBit/data/bebc.csv");      
      for (int row=0; row<88; ++row)
      {
        std::string line;
        getline(file, line);
        if (!file.good())
          break;
        std::stringstream iss(line);
        for (int col=0; col<2; ++col)
        {
          std::string val;
          getline(iss, val, ',');
          if (!iss)
            break;
          std::stringstream conv(val);
          conv >> array_bebc[row][col];
        }
      }
      std::ifstream file2("DarkBit/data/nutev.csv");
      for (int row2=0; row2<249; ++row2)
      {
        std::string line2;
        getline(file2, line2);
        if (!file.good())
          break;
        std::stringstream iss2(line2);
        for (int col2=0; col2<2; ++col2)
        {
          std::string val2;
          getline(iss2, val2, ',');
          if (!iss2)
            break;
          std::stringstream conv2(val2);
          conv2 >> array_nutev[row2][col2];
        }
      }
      std::ifstream file3("DarkBit/data/delphi.csv");
      for (int row3=0; row3<180; ++row3)
      {
        std::string line3;
        getline(file3, line3);
        if (!file3.good())
          break;
        std::stringstream iss3(line3);
        for (int col3=0; col3<2; ++col3)
        {
          std::string val3;
          getline(iss3, val3, ',');
          if (!iss3)
            break;
          std::stringstream conv3(val3);
          conv3 >> array_delphi[row3][col3];
        }
      }
      for (int i=0; i<88; i++)
      {
        M_temp_bebc[i] = array_bebc[i][0];
        U_temp_bebc[i] = array_bebc[i][1];
      }
      tk::spline s;
      s.set_points(M_temp_bebc, U_temp_bebc);
      for (int j=0; j<249; j++)
      {
        M_temp_nutev[j] = array_nutev[j][0];
        U_temp_nutev[j] = array_nutev[j][1];
      }
      tk::spline s2;
      s2.set_points(M_temp_nutev, U_temp_nutev);
      for (int k=0; k<180; k++)
      {
        M_temp_delphi[k] = array_delphi[k][0];
        U_temp_delphi[k] = array_delphi[k][1];
      }
      tk::spline s3;
      s3.set_points(M_temp_delphi, U_temp_delphi);

      double mixing_sq(*Dep::SN_stuff);
      M = *Param["M_2"];
      U_bebc = s(M);
      U_nutev = s2(M);
      U_delphi = s3(M);
      lnLike = -(2.44*mixing_sq)/U_bebc -(2.3*mixing_sq)/U_nutev -(5.14*mixing_sq)/U_delphi;
//      lnLike = -(5.14*mixing_sq*mixing_sq)/(U_delphi*U_delphi);
    }

    void printable_CI(double& Theta_sq)
    {
      namespace myPipe1 = Pipes::printable_CI;
      std::vector<double> temp1(9);
      Matrix3d t_sq(*myPipe1::Dep::SN_stuff);
      temp1[0] = t_sq(0,0);
      temp1[1] = t_sq(0,1);
      temp1[2] = t_sq(0,2);
      temp1[3] = t_sq(1,0);
      temp1[4] = t_sq(1,1);
      temp1[5] = t_sq(1,2);
      temp1[6] = t_sq(2,0);
      temp1[7] = t_sq(2,1);
      temp1[8] = t_sq(2,2);
      Theta_sq = temp1[0]+temp1[1]+temp1[2];
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

  }

}
