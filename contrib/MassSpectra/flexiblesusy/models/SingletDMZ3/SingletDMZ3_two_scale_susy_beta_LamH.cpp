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

// File generated at Fri 9 Jun 2017 16:14:09

#include "SingletDMZ3_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of LamH.
 *
 * @return one-loop beta function
 */
double SingletDMZ3_susy_parameters::calc_beta_LamH_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_LamH;

   beta_LamH = Re(oneOver16PiSqr*(0.27*Power(g1,4) + 2.25*Power(g2,4) +
      12*LamH*traceYdAdjYd - 12*traceYdAdjYdYdAdjYd + 4*LamH*traceYeAdjYe - 4*
      traceYeAdjYeYeAdjYe + 12*LamH*traceYuAdjYu - 12*traceYuAdjYuYuAdjYu - 9*
      LamH*Sqr(g2) + 0.9*Sqr(g1)*(-2*LamH + Sqr(g2)) + 12*Sqr(LamH) + 0.5*Sqr(
      LamSH)));


   return beta_LamH;
}

/**
 * Calculates the two-loop beta function of LamH.
 *
 * @return two-loop beta function
 */
double SingletDMZ3_susy_parameters::calc_beta_LamH_two_loop(const Susy_traces& susy_traces) const
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


   double beta_LamH;

   beta_LamH = Re(twoLoop*(-3.411*Power(g1,6) + 38.125*Power(g2,6) - 78*
      Power(LamH,3) - Power(LamSH,3) - 3*LamH*traceYdAdjYdYdAdjYd + 60*
      traceYdAdjYdYdAdjYdYdAdjYd + 12*traceYdAdjYdYdAdjYuYuAdjYd - 42*LamH*
      traceYdAdjYuYuAdjYd - 24*traceYdAdjYuYuAdjYdYdAdjYd - 12*
      traceYdAdjYuYuAdjYuYuAdjYd - LamH*traceYeAdjYeYeAdjYe + 20*
      traceYeAdjYeYeAdjYeYeAdjYe - 0.125*Power(g2,4)*(73*LamH + 12*(3*
      traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu)) - 3*LamH*
      traceYuAdjYuYuAdjYu + 60*traceYuAdjYuYuAdjYuYuAdjYu + 1.5*LamH*(36*LamH +
      5*(3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu))*Sqr(g2) - 0.015*
      Power(g1,4)*(-629*LamH - 60*traceYdAdjYd + 300*traceYeAdjYe + 228*
      traceYuAdjYu + 559*Sqr(g2)) + 80*LamH*traceYdAdjYd*Sqr(g3) - 64*
      traceYdAdjYdYdAdjYd*Sqr(g3) + 80*LamH*traceYuAdjYu*Sqr(g3) - 64*
      traceYuAdjYuYuAdjYu*Sqr(g3) - 72*traceYdAdjYd*Sqr(LamH) - 24*traceYeAdjYe
      *Sqr(LamH) - 72*traceYuAdjYu*Sqr(LamH) - 0.025*Sqr(g1)*(289*Power(g2,4) -
      6*(39*LamH + 36*traceYdAdjYd + 44*traceYeAdjYe + 84*traceYuAdjYu)*Sqr(g2
      ) - 4*(5*LamH*(5*traceYdAdjYd + 15*traceYeAdjYe + 17*traceYuAdjYu) + 16*(
      traceYdAdjYdYdAdjYd - 3*traceYeAdjYeYeAdjYe - 2*traceYuAdjYuYuAdjYu) +
      108*Sqr(LamH))) - 2.5*LamH*Sqr(LamSH)));


   return beta_LamH;
}

/**
 * Calculates the three-loop beta function of LamH.
 *
 * @return three-loop beta function
 */
double SingletDMZ3_susy_parameters::calc_beta_LamH_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamH;

   beta_LamH = 0;


   return beta_LamH;
}

} // namespace flexiblesusy
