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

// File generated at Sun 24 Sep 2017 15:56:09

#include "NSM_soft_parameters.hpp"
#include "config.h"
#include "global_thread_pool.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces
#define TRACE_STRUCT_TYPE Soft_traces
#define CALCULATE_TRACES(l) calc_soft_traces(l);

const int NSM_soft_parameters::numberOfParameters;

NSM_soft_parameters::NSM_soft_parameters(const NSM_input_parameters& input_)
   : NSM_susy_parameters(input_)
{
   set_number_of_parameters(numberOfParameters);
}

NSM_soft_parameters::NSM_soft_parameters(
   const NSM_susy_parameters& susy_model
   , double Lambda5_, double Lambda4_, double mS2_, double mH2_, double vH_,
   double vS_

)
   : NSM_susy_parameters(susy_model)
   , Lambda5(Lambda5_), Lambda4(Lambda4_), mS2(mS2_), mH2(mH2_), vH(vH_), vS(
   vS_)

{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd NSM_soft_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

NSM_soft_parameters NSM_soft_parameters::calc_beta(int loops) const
{
   double beta_Lambda5 = 0.;
   double beta_Lambda4 = 0.;
   double beta_mS2 = 0.;
   double beta_mH2 = 0.;
   double beta_vH = 0.;
   double beta_vS = 0.;

   if (loops > 0) {
      const auto TRACE_STRUCT = CALCULATE_TRACES(loops);

      beta_Lambda5 += calc_beta_Lambda5_1_loop(TRACE_STRUCT);
      beta_Lambda4 += calc_beta_Lambda4_1_loop(TRACE_STRUCT);
      beta_mS2 += calc_beta_mS2_1_loop(TRACE_STRUCT);
      beta_mH2 += calc_beta_mH2_1_loop(TRACE_STRUCT);
      beta_vH += calc_beta_vH_1_loop(TRACE_STRUCT);
      beta_vS += calc_beta_vS_1_loop(TRACE_STRUCT);

      if (loops > 1) {
         beta_Lambda5 += calc_beta_Lambda5_2_loop(TRACE_STRUCT);
         beta_Lambda4 += calc_beta_Lambda4_2_loop(TRACE_STRUCT);
         beta_mS2 += calc_beta_mS2_2_loop(TRACE_STRUCT);
         beta_mH2 += calc_beta_mH2_2_loop(TRACE_STRUCT);
         beta_vH += calc_beta_vH_2_loop(TRACE_STRUCT);
         beta_vS += calc_beta_vS_2_loop(TRACE_STRUCT);

         if (loops > 2) {
         #ifdef ENABLE_THREADS
            {

            }
         #else
         #endif

         }
      }
   }


   const NSM_susy_parameters susy_betas(NSM_susy_parameters::calc_beta(loops));

   return NSM_soft_parameters(susy_betas, beta_Lambda5, beta_Lambda4, beta_mS2, beta_mH2, beta_vH, beta_vS);
}

NSM_soft_parameters NSM_soft_parameters::calc_beta() const
{
   return calc_beta(get_loops());
}

void NSM_soft_parameters::clear()
{
   NSM_susy_parameters::clear();

   Lambda5 = 0.;
   Lambda4 = 0.;
   mS2 = 0.;
   mH2 = 0.;
   vH = 0.;
   vS = 0.;

}

Eigen::ArrayXd NSM_soft_parameters::get() const
{
   Eigen::ArrayXd pars(NSM_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(33) = Lambda5;
   pars(34) = Lambda4;
   pars(35) = mS2;
   pars(36) = mH2;
   pars(37) = vH;
   pars(38) = vS;


   return pars;
}

void NSM_soft_parameters::print(std::ostream& ostr) const
{
   NSM_susy_parameters::print(ostr);
   ostr << "soft parameters at Q = " << get_scale() << ":\n";
   ostr << "Lambda5 = " << Lambda5 << '\n';
   ostr << "Lambda4 = " << Lambda4 << '\n';
   ostr << "mS2 = " << mS2 << '\n';
   ostr << "mH2 = " << mH2 << '\n';
   ostr << "vH = " << vH << '\n';
   ostr << "vS = " << vS << '\n';

}

void NSM_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   NSM_susy_parameters::set(pars);

   Lambda5 = pars(33);
   Lambda4 = pars(34);
   mS2 = pars(35);
   mH2 = pars(36);
   vH = pars(37);
   vS = pars(38);

}

NSM_soft_parameters::Soft_traces NSM_soft_parameters::calc_soft_traces(int loops) const
{
   Soft_traces soft_traces;

   if (loops > 0) {
      TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());

   }

   if (loops > 1) {
      TRACE_STRUCT.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint())
         .trace());
      TRACE_STRUCT.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint())
         .trace());
      TRACE_STRUCT.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint())
         .trace());
      TRACE_STRUCT.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint())
         .trace());

   }

   if (loops > 2) {

   }

   return soft_traces;
}

std::ostream& operator<<(std::ostream& ostr, const NSM_soft_parameters& soft_pars)
{
   soft_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
