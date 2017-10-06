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

// File generated at Sun 24 Sep 2017 15:56:06

#include "NSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda3.
 *
 * @return 1-loop beta function
 */
double NSM_susy_parameters::calc_beta_Lambda3_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambda3;

   beta_Lambda3 = Re(0.5*oneOver16PiSqr*Lambda3*(4*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu + 6*Lambda1 + 12*Lambda2 + 4*Lambda3) - 3*
      Sqr(g1) - 9*Sqr(g2)));


   return beta_Lambda3;
}

/**
 * Calculates the 2-loop beta function of Lambda3.
 *
 * @return 2-loop beta function
 */
double NSM_susy_parameters::calc_beta_Lambda3_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda3;

   beta_Lambda3 = Re(-0.020833333333333332*twoLoop*Lambda3*(-557*Quad(g1)
      - 2*Sqr(g1)*(50*traceYdAdjYd + 150*traceYeAdjYe + 170*traceYuAdjYu + 576
      *Lambda1 + 48*Lambda3 + 45*Sqr(g2)) + 3*(145*Quad(g2) - 12*(15*
      traceYdAdjYd + 5*traceYeAdjYe + 15*traceYuAdjYu + 96*Lambda1 + 8*Lambda3)
      *Sqr(g2) + 8*(27*traceYdAdjYdYdAdjYd + 42*traceYdAdjYuYuAdjYd + 9*
      traceYeAdjYeYeAdjYe + 27*traceYuAdjYuYuAdjYu + 144*traceYdAdjYd*Lambda1 +
      48*traceYeAdjYe*Lambda1 + 144*traceYuAdjYu*Lambda1 + 48*traceYdAdjYd*
      Lambda3 + 16*traceYeAdjYe*Lambda3 + 48*traceYuAdjYu*Lambda3 + 288*Lambda1
      *Lambda3 + 576*Lambda2*Lambda3 - 80*(traceYdAdjYd + traceYuAdjYu)*Sqr(g3)
      + 120*Sqr(Lambda1) + 960*Sqr(Lambda2) + 84*Sqr(Lambda3)))));


   return beta_Lambda3;
}

/**
 * Calculates the 3-loop beta function of Lambda3.
 *
 * @return 3-loop beta function
 */
double NSM_susy_parameters::calc_beta_Lambda3_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda3;

   beta_Lambda3 = 0;


   return beta_Lambda3;
}

} // namespace flexiblesusy
