//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  The \f$Z\f$-boson observables at two-loop from
///
///  Complete electroweak two-loop corrections to \f$Z\f$ boson production and
///  decay
///  Ievgen Dubovyk, Ayres Freitas, Janusz Gluza, Tord Riemann, Johann Usovitsch
///  arXiv:1804.10236
///
///  I refer to tables and equations in v1 (https://arxiv.org/pdf/1804.10236v1).
///
///  \example SM_Z.cpp
///
///  *********************************************
///
///  Authors (add name and date if you modify):
///
///  \author Andrew Fowlie
///          (andrew.fowlie@monash.edu)
///  \date 2018 May
///
///  *********************************************

#ifndef DECAYBIT_INCLUDE_GAMBIT_DECAYBIT_SM_Z_HPP_
#define DECAYBIT_INCLUDE_GAMBIT_DECAYBIT_SM_Z_HPP_

#include <cmath>

namespace SM_Z {

/** 
   @brief <ahref="
   http://pdglive.lbl.gov/BranchingRatio.action?desig=9&parCode=S044
   ">PDG</a> measurement of invisible width of \f$Z\f$ boson in MeV
*/
constexpr struct {
  const double mu = 499.0;
  const double sigma = 1.5;
} gamma_invisible;


/** @brief The central values of nuisances from eq. 13 */
constexpr struct {
  const double mh_OS = 125.7;  // GeV
  const double mt_OS = 173.2;
  const double MZ_OS = 91.1876;
  const double alpha_s_MSbar_MZ = 0.1184;
  const double delta_alpha_OS = 0.059;
} hat;

constexpr int kRows = 12;
constexpr int kCols = 9;
/** @brief Coefficient data in Table 5 */
constexpr double table_5[kRows][kCols] = {
  {83.983, -0.061, 0.810, -0.096, -0.01, 0.25, -1.1, 286, 0.001},
  {83.793, -0.060, 0.810, -0.095, -0.01, 0.25, -1.1, 285., 0.001},
  {167.176, -0.071, 1.26, -0.19 , -0.02, 0.36, -0.1, 504., 0.001},
  {299.993, -0.38, 4.08, 14.27, 1.6, 1.8 , -11.1, 1253, 0.002},
  {299.916, -0.38, 4.08, 14.27, 1.6, 1.8 , -11.1, 1253., 0.002},
  {382.828, -0.39 , 3.83, 10.20, -2.4 , 0.67, -10.1, 1470., 0.002},
  {375.889, -0.36, -2.14, 10.53, -2.4 , 1.2 , -10.1, 1459., 0.006},
  {2494.74, -2.3 , 19.9, 58.61, -4.0 , 8.0 , -56.0, 9273., 0.012},
  {20751.6, -7.8 , -37., 732.3 , -44 , 5.5 , -358, 11696., 0.1 },
  {172.22, -0.031 , 1.0 , 2.3 , 1.3 , 0.38, -1.2, 37., 0.01},
  {215.85, 0.029, -2.92, -1.32, -0.84, 0.032, 0.72 , -18., 0.01},
  {41489.6, 1.6 , 60.0, -579.6, 38., 7.3 , 85., 0.1},
};


/**
   @brief  Data in Table 6, though re-arranged to match columns in Table 5
   
   The final entry isn't in the table and instead comes from the text below
   eq. 13.
*/
constexpr double table_6[kRows] =
  {0.018, 0.018, 0.016, 0.11, 0.11, 0.08, 0.18, 0.4, 6.e-3, 5.e-5, 1.e-4, 6.};

class TwoLoop {
  /**
     @brief \f$Z\f$-boson observables at two-loop and the residual theory errors
  */
 public:
  // Partial widths in MeV
  double gamma_e() const {return observable(0);}
  double gamma_mu() const {return gamma_e();}
  double gamma_tau() const {return observable(1);}
  double gamma_nu_e() const {return observable(2);}
  double gamma_nu_mu() const {return gamma_nu_e();}
  double gamma_invisible() const {return 3. * gamma_nu_e();}
  double gamma_nu_tau() const {return observable(2);}
  double gamma_u() const {return observable(3);}
  double gamma_c() const {return observable(4);}
  double gamma_t() const {return 0.;}
  double gamma_d() const {return observable(5);}
  double gamma_s() const {return gamma_s();}
  double gamma_b() const {return observable(6);}
  double gamma_total() const {return observable(7);}

  // Residual theory error in partial widths in MeV
  double error_gamma_e() const {return error(0);}
  double error_gamma_mu() const {return error_gamma_e();}
  double error_gamma_tau() const {return error(1);}
  double error_gamma_nu_e() const {return error(2);}
  double error_gamma_nu_mu() const {return error_gamma_nu_e();}
  double error_gamma_invisible() const {return 3. * error_gamma_nu_e();}
  double error_gamma_nu_tau() const {return error(2);}
  double error_gamma_u() const {return error(3);}
  double error_gamma_c() const {return error(4);}
  double error_gamma_t() const {return 0.;}
  double error_gamma_d() const {return error(5);}
  double error_gamma_s() const {return error_gamma_d();}
  double error_gamma_b() const {return error(6);}
  double error_gamma_total() const {return error(7);}

  // Branching ratios
  double BR_e() const {return BR(0);}
  double BR_mu() const {return BR_e();}
  double BR_tau() const {return BR(1);}
  double BR_nu_e() const {return BR(2);}
  double BR_nu_mu() const {return BR_nu_e();}
  double BR_invisible() const {return 3. * BR_nu_e();}
  double BR_nu_tau() const {return BR(2);}
  double BR_u() const {return BR(3);}
  double BR_c() const {return BR(4);}
  double BR_t() const {return 0.;}
  double BR_d() const {return BR(5);}
  double BR_s() const {return BR_d();}
  double BR_b() const {return BR(6);}

  // Residual theory error in branching ratios
  double error_BR_e() const {return error_BR(0);}
  double error_BR_mu() const {return error_BR_e();}
  double error_BR_tau() const {return error_BR(1);}
  double error_BR_nu_e() const {return error_BR(2);}
  double error_BR_nu_mu() const {return error_BR_nu_e();}
  double error_BR_invisible() const {return 3. * error_BR_nu_e();}
  double error_BR_nu_tau() const {return error_BR(2);}
  double error_BR_u() const {return error_BR(3);}
  double error_BR_c() const {return error_BR(4);}
  double error_BR_t() const {return 0.;}
  double error_BR_d() const {return error_BR(5);}
  double error_BR_s() const {return error_BR_d();}
  double error_BR_b() const {return error_BR(6);}

  // Ratios of partial widths, defined in eq. 27
  double Rl() const {return 1.e-3 * observable(8);}
  double Rc() const {return 1.e-3 * observable(9);}
  double Rb() const {return 1.e-3 * observable(10);}

  // Residual theory error in ratios of partial widths
  double error_Rl() const {return 1.e-3 * error(8);}
  double error_Rc() const {return 1.e-3 * error(9);}
  double error_Rb() const {return 1.e-3 * error(10);}

  // Hadronic peak cross section in pb, defined in eq. 10
  double sigma_0_had() const {return observable(11);}

  // Residual theory error in hadronic peak cross section in pb
  double error_sigma_0_had() const {return error(11);}

  // Nuisance parameters
  double mh_OS;
  double mt_OS;
  double MZ_OS;
  double alpha_s_MSbar_MZ;
  double delta_alpha_OS;

  TwoLoop(double mh_OS = hat.mh_OS,
          double mt_OS = hat.mt_OS,
          double MZ_OS = hat.MZ_OS,
          double alpha_s_MSbar_MZ = hat.alpha_s_MSbar_MZ,
          double delta_alpha_OS = hat.delta_alpha_OS):
    mh_OS{mh_OS},
    mt_OS{mt_OS},
    MZ_OS{MZ_OS},
    alpha_s_MSbar_MZ{alpha_s_MSbar_MZ},
    delta_alpha_OS{delta_alpha_OS}
    {
    /**
       @param mh_OS - Higgs mass in OS scheme
       @param mt_OS - Top quark mass in OS scheme
       @param MZ_OS - \f$Z\f$-mass in OS scheme
       @param alpha_s_MSbar_MZ - Strong coupling in MS-bar scheme at \f$Q = M_Z\f$
       @param delta_alpha_OS - \f$\Delta\alpha\f$ parameter in OS scheme. Defined on p9
    */
    L_H = std::log(mh_OS / hat.mh_OS);
    delta_H = mh_OS / hat.mh_OS - 1.;
    delta_t = pow(mt_OS / hat.mt_OS, 2) - 1.;
    delta_z = MZ_OS / hat.MZ_OS - 1.;
    delta_alpha_s = alpha_s_MSbar_MZ / hat.alpha_s_MSbar_MZ - 1.;
    delta_delta_alpha = delta_alpha_OS / hat.delta_alpha_OS - 1.;
  }

 private:
  double L_H;
  double delta_H;
  double delta_t;
  double delta_z;
  double delta_alpha_s;
  double delta_delta_alpha;

  double error(int row) const {
    /**
       @brief Error in observable calculated from eq. 13
       
       We add the error from the parametric formula and theory error in
       quadrature.
       
       @param row Row number of Table 5 corresponding to quantity
       @returns Error in quantity
    */
    return std::sqrt(pow(table_5[row][9], 2) + pow(table_6[row], 2));
  }

  double observable(int row) const {
    /**
       @returns The observable calculated from eq. 13
       @param row Row number of Table 5 corresponding to quantity
    */    
    return table_5[row][0] +
           table_5[row][1] * L_H +
           table_5[row][2] * delta_t +
           table_5[row][3] * delta_H +
           table_5[row][4] * delta_alpha_s +
           table_5[row][5] * pow(delta_alpha_s, 2) +
           table_5[row][6] * delta_alpha_s * delta_t +
           table_5[row][7] * delta_delta_alpha +
           table_5[row][8] * delta_z;
  }

  double BR(int row) const {
    /*
      @param row Row number of Table 7 corresponding to quantity
      @returns Branching ratio
    */
    return observable(row) / gamma_total();
  }

  double error_BR(int row) const {
    /**
       @warning We propagate an error in \f$f = x / y\f$. In fact, though, we
       should propagate an error in \f$f = x / (x + y)\f$, since the partial
       width in the numerator contributes to the total width. Thus this formula
       is reliable only for small branching ratios.

       @param row Row number of Table 7 corresponding to quantity
       @returns Error in a branching ratio found by propagating errors
    */
    const double frac_error_sq = pow(error_gamma_total() / gamma_total(), 2)
      + pow(error(row) / observable(row), 2);
    return std::sqrt(frac_error_sq) * BR(row);
  }
};

}  // namespace SM_Z

#endif  //  DECAYBIT_INCLUDE_GAMBIT_DECAYBIT_SM_Z_HPP_
