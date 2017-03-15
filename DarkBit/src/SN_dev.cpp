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
      static double m = 1e-7;
      static double s = 2e-18;
      namespace lnLPipe = Pipes::lnL;
      double mixing_sq(*lnLPipe::Dep::SN_stuff);
      lnLike = -(pow((mixing_sq - m), 2.0))/s;
    }

    void printable_CI(double& Theta_sq)
    {
      namespace myPipe = Pipes::printable_CI;
      std::vector<double> temp1(9);
      Matrix3d t_sq(*myPipe::Dep::SN_stuff);
      temp1[0] = t_sq(0,0);
      temp1[1] = t_sq(0,1);
      temp1[2] = t_sq(0,2);
      temp1[3] = t_sq(1,0);
      temp1[4] = t_sq(1,1);
      temp1[5] = t_sq(1,2);
      temp1[6] = t_sq(2,0);
      temp1[7] = t_sq(2,1);
      temp1[8] = t_sq(2,2);
      Theta_sq = temp1[0] + temp1[3] + temp1[6];
    }

  }

}
