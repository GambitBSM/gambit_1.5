//   GAMBIT: Global and Modular BSM Inference Tool
//   *********************************************
///  \file
///
///  Higgs boson decays to supersymmetric particles at tree-level from
///
///  The Anatomy of Electro-Weak Symmetry Breaking. II:
///  The Higgs bosons in the Minimal Supersymmetric Model
///  Abdelhak Djouadi
///  arXiv:hep-ph/0503173
///
///  I refer to equations in v2 (https://arxiv.org/abs/hep-ph/0503173v2).
///
///  \example MSSM_H.cpp
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

#ifndef DECAYBIT_INCLUDE_GAMBIT_DECAYBIT_MSSM_H_HPP_
#define DECAYBIT_INCLUDE_GAMBIT_DECAYBIT_MSSM_H_HPP_

#include <cmath>
#include <array>

namespace MSSM_H {

double lambda(double x, double y, double z) {
  /**
     @brief Phase-space function (similar to Kallen function)

     \f[
     \lambda(x, y; z) = (1 - x / z - y / z)^2 - 4 x y / z^2
     \f]
  */
  return pow(1. - x / z - y / z, 2) - 4. * x * y / pow(z, 2);
}

double gamma_h_chi_0(int i,
  int j,
  std::array<double, 4> m,
  std::array<std::array<double, 4>, 4> Z,
  double alpha,
  double mh = 125.,
  double mw = 80.385,
  double GF = 1.1663787e-5,
  double sw2 = 0.22) {
  /**
     @brief Lightest Higgs boson decay to neutralinos at tree-level in GeV
     
     @warning Arguments \f$i\f$ and \f$j\f$ are zero-based
     @warning Tree-level formula
     
     @returns \f$\Gamma(h \to \chi_i \chi_j)\f$ in GeV
     @param i Neutralino \f$\chi_i\f$ in final state
     @param j Neutralino \f$\chi_j\f$ in final state
     @param m Neutralino mass matrix with phases
     @param Z Real neutralino mixing matrix
     @param alpha \f$\alpha\f$, Higgs mixing angle
     @param mh Lightest Higgs mass, \f$m_h\f$
     @param GF Fermi constant, \f$G_F\f$
     @param sw2 Weinberg angle, \f$\sin^2\theta_W\f$
  */  
  const double sw = std::sqrt(sw2);
  const double cw = std::sqrt(1. - sw2);
  const double tw = sw / cw;

  // Eq. 1.113
  const double e2 = -std::sin(alpha);
  const double d2 = -std::cos(alpha);
  const double eps_2 = 1.;

  // Phase-space
  double l = lambda(pow(m[i], 2), pow(m[j], 2), pow(mh, 2));
  if (l <= 0.) {
    return 0.;
  }

  // Eq. 1.112 (without \f$\sin\theta_W\f$ which cancels)
  const double gL = 0.5 * (
    (Z[j][1] - tw * Z[j][0]) * (e2 * Z[i][2] + d2 * Z[i][3]) +
    (Z[i][1] - tw * Z[i][0]) * (e2 * Z[j][2] + d2 * Z[j][3]));
  const double gR = eps_2 * gL;

  // Eq. 2.56 without common factor (without f$\sin\theta_W\f$ which cancels)
  const double delta = (i == j) ? 1. : 0.;
  const double gamma_no_prefactor = std::sqrt(l) / (1. + delta) * ((
    pow(gL, 2) + pow(gR, 2)) * (1. - (pow(m[i], 2) + pow(m[j], 2)) / pow(mh, 2))
    -4. * gL * gR * m[i] * m[j] / pow(mh, 2));

  // Eq. 2.56 with common factor
  return GF * pow(mw, 2) / (2. * M_SQRT2 * M_PI) * mh * gamma_no_prefactor;
}
}  // namespace MSSM_H

#endif  //  DECAYBIT_INCLUDE_GAMBIT_DECAYBIT_MSSM_H_HPP_
