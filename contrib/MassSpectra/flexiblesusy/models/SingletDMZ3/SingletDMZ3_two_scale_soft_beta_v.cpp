// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Fri 9 Jun 2017 16:14:12

#include "SingletDMZ3_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

namespace {

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(double n, const Eigen::MatrixBase<Derived>& m)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m - Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(double n, const Eigen::MatrixBase<Derived>& m)
{
   return - m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

} // anonymous namespace

/**
 * Calculates the one-loop beta function of v.
 *
 * @return one-loop beta function
 */
double SingletDMZ3_soft_parameters::calc_beta_v_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_v;

   beta_v = Re(oneOver16PiSqr*(0.6*v*Sqr(g1) + v*(-3*traceYdAdjYd -
      traceYeAdjYe - 3*traceYuAdjYu + 3*Sqr(g2))));


   return beta_v;
}

/**
 * Calculates the two-loop beta function of v.
 *
 * @return two-loop beta function
 */
double SingletDMZ3_soft_parameters::calc_beta_v_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_v;

   beta_v = Re(-0.00125*twoLoop*v*(1257*Power(g1,4) + 10*Sqr(g1)*(86*
      traceYdAdjYd + 162*traceYeAdjYe + 206*traceYuAdjYu - 9*Sqr(g2)) - 25*(343
      *Power(g2,4) - 84*(3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu)*Sqr(g2
      ) - 4*(2*(-3*(9*traceYdAdjYdYdAdjYd - 2*traceYdAdjYuYuAdjYd + 3*
      traceYeAdjYeYeAdjYe + 9*traceYuAdjYuYuAdjYu) + 80*(traceYdAdjYd +
      traceYuAdjYu)*Sqr(g3)) + 12*Sqr(LamH) + Sqr(LamSH)))));


   return beta_v;
}

/**
 * Calculates the three-loop beta function of v.
 *
 * @return three-loop beta function
 */
double SingletDMZ3_soft_parameters::calc_beta_v_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v;

   beta_v = 0;


   return beta_v;
}

} // namespace flexiblesusy
