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

// File generated at Sun 24 Sep 2017 15:56:07

#include "NSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda1.
 *
 * @return 1-loop beta function
 */
double NSM_susy_parameters::calc_beta_Lambda1_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda1;

   beta_Lambda1 = Re(oneOver16PiSqr*(0.375*Quad(g1) + 1.125*Quad(g2) - 9*
      Lambda1*Sqr(g2) + 0.75*Sqr(g1)*(-4*Lambda1 + Sqr(g2)) + 2*(-3*
      traceYdAdjYdYdAdjYd - traceYeAdjYeYeAdjYe - 3*traceYuAdjYuYuAdjYu + 6*
      traceYdAdjYd*Lambda1 + 2*traceYeAdjYe*Lambda1 + 6*traceYuAdjYu*Lambda1 +
      12*Sqr(Lambda1) + Sqr(Lambda3))));


   return beta_Lambda1;
}

/**
 * Calculates the 2-loop beta function of Lambda1.
 *
 * @return 2-loop beta function
 */
double NSM_susy_parameters::calc_beta_Lambda1_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYdYdAdjYdYdAdjYd =
      TRACE_STRUCT.traceYdAdjYdYdAdjYdYdAdjYd;
   const double traceYdAdjYdYdAdjYuYuAdjYd =
      TRACE_STRUCT.traceYdAdjYdYdAdjYuYuAdjYd;
   const double traceYdAdjYuYuAdjYdYdAdjYd =
      TRACE_STRUCT.traceYdAdjYuYuAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYuYuAdjYd =
      TRACE_STRUCT.traceYdAdjYuYuAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYeYeAdjYe =
      TRACE_STRUCT.traceYeAdjYeYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYuYuAdjYu =
      TRACE_STRUCT.traceYuAdjYuYuAdjYuYuAdjYu;


   double beta_Lambda1;

   beta_Lambda1 = Re(twoLoop*(30*traceYdAdjYdYdAdjYdYdAdjYd - 12*
      traceYdAdjYdYdAdjYuYuAdjYd + 6*traceYdAdjYuYuAdjYdYdAdjYd - 6*
      traceYdAdjYuYuAdjYuYuAdjYd + 10*traceYeAdjYeYeAdjYeYeAdjYe + 30*
      traceYuAdjYuYuAdjYuYuAdjYu - 312*Cube(Lambda1) - 16*Cube(Lambda3) - 3*
      traceYdAdjYdYdAdjYd*Lambda1 - 42*traceYdAdjYuYuAdjYd*Lambda1 -
      traceYeAdjYeYeAdjYe*Lambda1 - 3*traceYuAdjYuYuAdjYu*Lambda1 -
      7.895833333333333*Power6(g1) + 19.0625*Power6(g2) - 0.125*(18*
      traceYdAdjYd + 6*traceYeAdjYe + 18*traceYuAdjYu + 73*Lambda1)*Quad(g2) +
      1.5*Lambda1*(15*traceYdAdjYd + 5*traceYeAdjYe + 15*traceYuAdjYu + 72*
      Lambda1)*Sqr(g2) - 0.020833333333333332*Quad(g1)*(-60*traceYdAdjYd + 300*
      traceYeAdjYe + 228*traceYuAdjYu - 1258*Lambda1 + 559*Sqr(g2)) - 32*
      traceYdAdjYdYdAdjYd*Sqr(g3) - 32*traceYuAdjYuYuAdjYu*Sqr(g3) + 80*
      traceYdAdjYd*Lambda1*Sqr(g3) + 80*traceYuAdjYu*Lambda1*Sqr(g3) - 144*
      traceYdAdjYd*Sqr(Lambda1) - 48*traceYeAdjYe*Sqr(Lambda1) - 144*
      traceYuAdjYu*Sqr(Lambda1) - 0.020833333333333332*Sqr(g1)*(289*Quad(g2) -
      12*(18*traceYdAdjYd + 22*traceYeAdjYe + 42*traceYuAdjYu + 39*Lambda1)*Sqr
      (g2) - 8*(8*traceYdAdjYdYdAdjYd - 24*traceYeAdjYeYeAdjYe - 16*
      traceYuAdjYuYuAdjYu + 25*traceYdAdjYd*Lambda1 + 75*traceYeAdjYe*Lambda1 +
      85*traceYuAdjYu*Lambda1 + 216*Sqr(Lambda1))) - 20*Lambda1*Sqr(Lambda3)))
      ;


   return beta_Lambda1;
}

/**
 * Calculates the 3-loop beta function of Lambda1.
 *
 * @return 3-loop beta function
 */
double NSM_susy_parameters::calc_beta_Lambda1_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda1;

   beta_Lambda1 = 0;


   return beta_Lambda1;
}

} // namespace flexiblesusy
