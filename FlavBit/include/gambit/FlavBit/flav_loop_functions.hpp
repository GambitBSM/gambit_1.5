//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Loop functions for flavour violating decays of charged leptons (from hep-ph/9403398)
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Tomas Gonzalo
///          (t.e.gonzalo@fys.uio.no)
///  \date 2017 Aug
///
///  *********************************************

#ifndef __flav_loop_functions_hpp
#define __flav_loop_functions_hpp

namespace Gambit
{

  namespace FlavBit
  {

    // Generic loop functions

    double Fgamma(const double x)
    {
      if(x == 0)
        return 0; 
      else if(x == 1)
        return -25.0/72.0;
      else
        return (-12.0*x - pow(x,2) + 7.0*pow(x,3))/(12.0*pow(1.0 - x,3)) - ((12.0*pow(x,2) - 10.0*pow(x,3) + pow(x,4))*log(x))/(6.0*pow(1.0 - x,4));
    }

    double Ggamma(const double x)
    {
      if(x == 0)
        return 0;
      else if(x == 1)
        return 1.0/8.0;
      else
        return -(-x + 5.0*pow(x,2) + 2.0*pow(x,3))/(4.0*pow(1.0 - x,3)) - (3.0*pow(x,3)*log(x))/(2.0*pow(1.0 - x,4));
    }

    double FZ(const double x)
    {
      if(x == 0)
        return 0;
      else if(x == 1)
        return -5.0/4.0;
      else
        return (-5.0*x)/(2.0*(1.0 - x)) - (5.0*pow(x,2)*log(x))/(2.0*pow(1.0 - x,2));
    }

    double GZ(const double x, const double y)
    {
      if(x == 0 and y == 0)
        return 0;
      else if(x == 1 or y == 1)
        return 0.5;
      else if(y == 0)
        return -(x*log(x))/(2.0*(1.0 - x));
      else if(x == 0)
        return GZ(y,x);
      else if(x == y)
        return -x/2.0 - (x*log(x))/(1.0 - x);
      else
        return -((pow(x,2)*(1.0 - y)*log(x))/(1.0 - x) - ((1.0 - x)*pow(y,2)*log(y))/(1.0 - y))/(2.0*(x - y));
    }

    double HZ(const double x, const double y)
    {
      if(x == 0 or y == 0)
        return 0;
      else if(x == 1 and y == 1)
        return 1.0/8.0;
      else if(y == 1)
        return sqrt(x)/4.0 * (3.0/(1.0 - x) - (pow(x,2) - 4*x)*log(x)/pow(1.0 - x,2));
      else if(x == 1)
        return HZ(y,x);
      else if(x == y)
        return 3.0/4.0 - x/4.0 - 3.0/(4.0*(1.0 - x)) - (pow(x,3) - 2.0*pow(x,2) + 4.0*x)*log(x)/(4.0*pow(1.0 - x,2));
      else
        return (sqrt(x*y)*(((-4.0*x + pow(x,2))*log(x))/(1.0 - x) - ((-4.0*y + pow(y,2))*log(y))/(1.0 - y)))/(4.0*(x - y));
    }

    double Fbox(const double x, const double y)
    {

      if(x == 0 and y == 0)
        return 1;
      else if(x == 1 and y == 1)
        return 3.0/4.0;
      else if( (x == 1 and y == 0) or (x == 0 and y == 0) )
        return 1.0/2.0;
      else if(y == 0)
        return 1.0/(1.0 - x) + x*log(x)/pow(1.0 - x,2);
      else if(x == 0)
        return Fbox(y,x);
      else if(y == 1)
        return - (5.0*pow(x,3) - 8.0*pow(x,2) + 7.0*x - 4.0)/(8.0*pow(1.0 - x,3)) - (pow(x,3) - 4.0*pow(x,2))/(4.0*pow(1.0 - x,3))*log(x);
      else if(x == 1)
        return Fbox(y,x);
      else if(x == y)
        return -(pow(x,4) - 16.0*pow(x,3) + 19.0*pow(x,2) - 4.0)/(4.0*pow(1.0 - x,3)) - (3.0*pow(x,3) + 4.0*pow(x,2) - 4.0*x)/(2.0*pow(1.0 - x,3))*log(x);
      else
        return 1.0/(x - y)*((1.0 + (x*y)/4.0)*(1.0/(1.0 - x) + (x*x*log(x))/pow(1.0 - x,2) - 1.0/(1.0 - y) - (y*y*log(y))/pow(1.0 - y,2)) - 2.0*x*y*(1.0/(1.0 - x) + (x*log(x))/pow(1.0 - x,2) - 1.0/(1.0 - y) - (y*log(y))/pow(1.0 - y,2)));

    }

    double Gbox(const double x, const double y)
    {
      if(x == 0 or y == 0)
        return 0;
      else if(x == 1 and y == 1)
        return 3.0/2.0;
      else if(y == 1)
        return - sqrt(x)*((pow(x,3) - 2.0*x*x + 7.0*x - 6.0)/(2.0*pow(1.0 - x,2)) + (x*x - 4.0*x)/pow(1.0 - x,3)*log(x));
      else if(x == 1)
        return Gbox(y,x);
      else if(x == y)
        return (2*pow(x,4) - 4.0*pow(x,3) + 8.0*x*x - 6.0*x)/pow(1.0 - x,3) - (pow(x,4) + pow(x,3) + 4.0*x)/pow(1 - x,3)*log(x);
      else
        return -((sqrt(x*y)*((4.0 + x*y)*(1.0/(1.0 - x) - 1.0/(1.0 - y) + (x*log(x))/pow(1.0 - x,2) - (y*log(y))/pow(1.0 - y,2)) - 2.0*(1.0/(1.0 - x) - 1.0/(1.0 - y) + (pow(x,2)*log(x))/pow(1.0 - x,2) - (pow(y,2)*log(y))/pow(1.0 - y,2))))/(x - y));

    }

    // Specific loop functions for Sterile Neutrinos

    std::complex<double> Fgamma(int l, int lp, Eigen::Matrix3cd Vnu, Eigen::Matrix3cd Theta, std::vector<double> x)
    {
      std::complex<double> fgamma = {0,0};

      for(int i = 0; i < 3; i++)
        fgamma += Vnu.adjoint()(i,l)*Vnu(lp,i)*Fgamma(x[i]);
      for(int i = 0; i < 3; i++)
        fgamma += Theta.adjoint()(i,l)*Theta(lp,i)*Fgamma(x[i+3]);

      return fgamma;
    }

    std::complex<double> Ggamma(int l, int lp, Eigen::Matrix3cd Vnu, Eigen::Matrix3cd Theta, std::vector<double> x)
    {
      std::complex<double> ggamma = {0,0};

      for(int i = 0; i < 3; i++)
        ggamma += Vnu.adjoint()(i,l)*Vnu(lp,i)*Ggamma(x[i]);
      for(int i = 0; i < 3; i++)
        ggamma += Theta.adjoint()(i,l)*Theta(lp,i)*Ggamma(x[i+3]);

      return ggamma;
    }

    std::complex<double> FZ(int l, int lp, Eigen::Matrix3cd Vnu, Eigen::Matrix3cd Theta, std::vector<double> x)
    {
      std::complex<double> fz = {0,0};
      Eigen::Matrix<std::complex<double>,3,6> U;
 
       for(int i = 0; i < 3; i++)
       {
        for(int j = 0; j < 3; j++)
        {
          U(i,j) = Vnu(i,j);
          U(i,j+3) = Theta(i,j);
        }
      }

      for(int i = 0; i < 6; i++)
      {
        fz += U.adjoint()(i,l)*U(lp,i)*FZ(x[i]);
        for(int j = 0; j < 6; j++)
          for(int k = 0; k < 3; k++)
            fz += U.adjoint()(i,l)*U(lp,j)*(U(k,i)*U.adjoint()(j,k)*GZ(x[i],x[j]) + U.adjoint()(i,k)*U(k,j)*HZ(x[i],x[j])); 
      }

      return fz;

    }
 
    std::complex<double> Fbox(int l, int lp, int l1, int l2, Eigen::Matrix3cd Vnu, Eigen::Matrix3cd Theta, std::vector<double> x)
    {
  
      std::complex<double> fbox = {0,0};
      Eigen::Matrix<std::complex<double>,3,6> U;

      for(int i = 0; i < 3; i++)
      {
        for(int j = 0; j < 3; j++)
        {
          U(i,j) = Vnu(i,j);
          U(i,j+3) = Theta(i,j);
        }
      }

      for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
          fbox += U.adjoint()(i,l)*U.adjoint()(j,l2)*(U(lp,i)*U(l1,j) + U(l1,i)*U(lp,j))*Fbox(x[i],x[j]) + U.adjoint()(i,l)*U.adjoint()(i,l2)*U(lp,j)*U(l1,j)*Gbox(x[i],x[j]);

      return fbox;
    }

  }
}

#endif //#defined __flav_loop_functions_hpp__
