//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Function definitions of NeutrinoBit.
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Tomas Gonzalo
///          (t.e.gonzalo@fys.uio.no)
///  \date 2017 Juky
///
///  *********************************************

#include "gambit/Elements/gambit_module_headers.hpp"
#include "gambit/NeutrinoBit/NeutrinoBit_rollcall.hpp"
#include <unsupported/Eigen/MatrixFunctions>

namespace Gambit
{

  namespace NeutrinoBit
  {

    using namespace LogTags;

    // Module functions

    // Neutrino mass matrix from true SM neutrino model
    void M_nu(Eigen::Matrix3cd& m_nu)
    {
      using namespace Pipes::M_nu;

//      int ordering = *Param["ordering"];
//      double m_min = *Param["min_mass"];
//      double md21 = *Param["md21"];
//      double md31 = *Param["md31"];
//      double md23 = *Param["md23"];    

      double mnu1 = *Param["mnu1"];
      double mnu2 = *Param["mnu2"];
      double mnu3 = *Param["mnu3"];
       
      m_nu(0,1) = 0.0;
      m_nu(0,2) = 0.0;
      m_nu(1,0) = 0.0;
      m_nu(1,2) = 0.0;
      m_nu(2,0) = 0.0;
      m_nu(2,1) = 0.0;

/*      if(ordering == 1) // Normal hierarchy
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
      else if(ordering == 0) // Inverted hierarchy
      {
        if(m_min == 0)
        {
          m_nu(2,2) = 0.0;
          m_nu(1,1) = sqrt(md23);
          m_nu(0,0) = sqrt(pow(m_nu(1,1), 2.0) - md21);
        }
        else if(m_min == 1)
        {
          m_nu(2,2) = 2.3e-10;
          m_nu(1,1) = sqrt(pow(m_nu(2,2), 2.0) + md23);
          m_nu(0,0) = sqrt(pow(m_nu(1,1), 2.0) - md21);
        }
      }
*/
      m_nu(0,0) = mnu1;
      m_nu(1,1) = mnu2;
      m_nu(2,2) = mnu3;
    }

 
    // PMNS matrix in the Casas-Ibarra paramtetrization
    void UPMNS(Eigen::Matrix3cd& U_nu)
    {
      using namespace Pipes::UPMNS;
     
      Eigen::Matrix3cd V_23, V_13, V_12, U_pd, U_nd, Maj_phase;
      // TODO: change parameter names when the conlfict of models is resolved
      double theta23 = *Param["theta23"];
      double theta12 = *Param["theta12"];
      double theta13 = *Param["theta13"];
      double delta = *Param["delta13"];
      double alpha1 = *Param["alpha1"];
      double alpha2 = *Param["alpha2"];
      std::complex<double> I(0.0, 1.0);

      V_23 << 1.0, 0.0, 0.0,
              0.0, cos(theta23), sin(theta23),
              0.0, -sin(theta23), cos(theta23);
      V_13 << cos(theta13), 0.0, sin(theta13),
              0.0, 1.0, 0.0,
              -sin(theta13), 0.0, cos(theta13);
      V_12 << cos(theta12), sin(theta12), 0.0,
              -sin(theta12), cos(theta12), 0.0,
              0.0, 0.0, 1.0;
      U_pd << exp(-I*delta/2.0), 0.0, 0.0,
              0.0, 1.0, 0.0,
              0.0, 0.0, exp(I*delta/2.0);
      U_nd << exp(I*delta/2.0), 0.0, 0.0,
              0.0, 1.0, 0.0,
              0.0, 0.0, exp(-I*delta/2.0);
      Maj_phase << exp(I*alpha1/2.0), 0.0, 0.0,
                   0.0, exp(I*alpha2/2.0), 0.0,
                   0.0, 0.0, 1.0;
      U_nu = V_23 * U_pd * V_13 * U_nd * V_12 * Maj_phase;

    }

    // Helper function for the heavy neutrino masses
    double l_M(double M, const double m_Z, const double m_H)
    {
      return 1.0/pow(4.0*pi, 2.0) * ( (3.0*log(pow(M/m_Z, 2.0)))/((pow(M/m_Z, 2.0)) - 1.0) + (log(pow(M/m_H, 2.0)))/((pow(M/m_H, 2.0)) - 1.0));
    }

    // Theta matrix in Seesaw I in the Casas-Ibarra parametrization
    void CI_Theta(Eigen::Matrix3cd& Theta)  // capability: SeesawI_Theta
    {
      using namespace Pipes::CI_Theta;
      SMInputs sminputs = *Dep::SMINPUTS;

      std::complex<double> I(0.0, 1.0);

      Eigen::Matrix3d M_I;  // M_I not complex; circumvents type mismatch in l(M)
      Eigen::Matrix3cd M_twid_temp, M_twid, R_23, R_13, R_12, R;

      double mZ = sminputs.mZ;
      double mH = *Param["mH"];
      double vev = 1. / sqrt(sqrt(2.)*sminputs.GF);

      double x23 = *Param["ReOm23"];
      double y23 = *Param["ImOm23"];
      double x13 = *Param["ReOm13"];
      double y13 = *Param["ImOm13"];
      double x12 = *Param["ReOm12"];
      double y12 = *Param["ImOm12"];

      M_I << *Param["M_1"], 0.0, 0.0,
             0.0, *Param["M_2"], 0.0,
             0.0, 0.0, *Param["M_3"];

      M_twid_temp(0,0) = M_I(0,0)  * (1.0 - (pow(M_I(0,0),2.0)*l_M(M_I(0,0),mZ,mH)/pow(vev,2.0)));
      M_twid_temp(0,1) = 0.0;
      M_twid_temp(0,2) = 0.0;
      M_twid_temp(1,0) = 0.0;
      M_twid_temp(1,1) = M_I(1,1)  * (1.0 - (pow(M_I(1,1),2.0)*l_M(M_I(1,1),mZ,mH)/pow(vev,2.0)));
      M_twid_temp(1,2) = 0.0;
      M_twid_temp(2,0) = 0.0;
      M_twid_temp(2,1) = 0.0;
      M_twid_temp(2,2) = M_I(2,2)  * (1.0 - (pow(M_I(2,2),2.0)*l_M(M_I(2,2),mZ,mH)/pow(vev,2.0)));
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

      Theta = I * *Dep::UPMNS * Dep::m_nu->sqrt() * R * M_twid.inverse();

      // This parametrisation is not valid when |Theta|^2_ij > 1, so invalidate those points
      Eigen::Matrix3d ThetaNorm = (Theta.adjoint() * Theta).real();
      Eigen::Matrix3d ThetaNorm2 = (Theta * Theta.adjoint()).real();
      for(int i=1; i<3; i++)
        for(int j=1; j<3; j++)
        {
          if(ThetaNorm(i,j) > 1 or ThetaNorm2(i,j) > 1 or abs(Theta(i,j)) > 1)
          {
            std::ostringstream msg;
            msg << "Casas-Ibarra parametrization breaks down for parameter point";
            logger() << msg.str() << EOM;
            invalid_point().raise(msg.str());
          }
        }
    }


    void Vnu(Eigen::Matrix3cd &V)
    {
      using namespace Pipes::Vnu;
      Eigen::Matrix3cd Theta = *Dep::SeesawI_Theta;
      Eigen::Matrix3cd U = *Dep::UPMNS;

      V = U - 0.5*Theta*Theta.adjoint()*U;

    }

    void Unitarity_UPMNS(bool &unitarity)
    {
      using namespace Pipes::Unitarity_UPMNS;

      Eigen::Matrix3cd Id;
      Id << 1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0;
      Eigen::Matrix3d Epsilon;
      Epsilon << 1e-08, 1e-08, 1e-08,
                 1e-08, 1e-08, 1e-08,
                 1e-08, 1e-08, 1e-08;

      Eigen::Matrix3cd Norm = Dep::UPMNS->adjoint() * *Dep::UPMNS;
      unitarity = true;
      for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
          if(std::abs(Norm(i,j) - Id(i,j)) > Epsilon(i,j))
            unitarity = false;

      if(!unitarity)
        return ;

      for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
          if(std::real((*Dep::m_nu)(i,j)*pow((*Dep::UPMNS)(i,j),2)) > Epsilon(i,j))
            unitarity = false;

   
    }

    void Unitarity_SeesawI(bool &unitarity)
    {
      using namespace Pipes::Unitarity_SeesawI;

      Eigen::Matrix3cd Id;
      Id << 1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0;
      Eigen::Matrix3d Epsilon;
      Epsilon << 1e-08, 1e-08, 1e-08,
                 1e-08, 1e-08, 1e-08,
                 1e-08, 1e-08, 1e-08;

      Eigen::Matrix3cd Vnu = *Dep::SeesawI_Vnu;
      Eigen::Matrix3cd Theta = *Dep::SeesawI_Theta;
      Eigen::Matrix3cd m_nu = *Dep::m_nu;

      Eigen::Matrix3cd Norm = Vnu.adjoint() * Vnu + Theta.adjoint() * Theta;
      unitarity = true;
      for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
          if(std::norm(Norm(i,j) - Id(i,j)) > Epsilon(i,j))
            unitarity = false;

      if(!unitarity)
        return ;

      unitarity = true;
      Eigen::Matrix3d MN;
      MN << *Param["M_1"], 0.0, 0.0,
            0.0, *Param["M_2"], 0.0,
            0.0, 0.0, *Param["M_3"];
      Eigen::Matrix3cd Unit;
      for(int i = 0; i < 3; i++)
      {
        Unit(i,i) = 0;
        for(int j = 0; j < 3; j++)
          Unit(i,i) += m_nu(j,j)*pow(Vnu(i,j),2) + MN(j,j) * pow(Theta(i,j),2);
        if(std::real(Unit(i,i)) > Epsilon(i,i))
          unitarity = false;
      }
    }

  }
}
