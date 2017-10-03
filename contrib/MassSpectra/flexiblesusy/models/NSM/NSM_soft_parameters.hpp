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

#ifndef NSM_soft_parameters_H
#define NSM_soft_parameters_H

#include "NSM_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class NSM_soft_parameters : public NSM_susy_parameters {
public:
   explicit NSM_soft_parameters(const NSM_input_parameters& input_ = NSM_input_parameters());
   NSM_soft_parameters(const NSM_susy_parameters& , double Lambda5_, double Lambda4_, double mS2_, double mH2_, double vH_,
   double vS_
);
   NSM_soft_parameters(const NSM_soft_parameters&) = default;
   NSM_soft_parameters(NSM_soft_parameters&&) = default;
   virtual ~NSM_soft_parameters() = default;
   NSM_soft_parameters& operator=(const NSM_soft_parameters&) = default;
   NSM_soft_parameters& operator=(NSM_soft_parameters&&) = default;

   virtual Eigen::ArrayXd beta() const override;
   virtual Eigen::ArrayXd get() const override;
   virtual void print(std::ostream&) const override;
   virtual void set(const Eigen::ArrayXd&) override;

   NSM_soft_parameters calc_beta() const;
   NSM_soft_parameters calc_beta(int) const;
   virtual void clear() override;

   void set_Lambda5(double Lambda5_) { Lambda5 = Lambda5_; }
   void set_Lambda4(double Lambda4_) { Lambda4 = Lambda4_; }
   void set_mS2(double mS2_) { mS2 = mS2_; }
   void set_mH2(double mH2_) { mH2 = mH2_; }
   void set_vH(double vH_) { vH = vH_; }
   void set_vS(double vS_) { vS = vS_; }

   double get_Lambda5() const { return Lambda5; }
   double get_Lambda4() const { return Lambda4; }
   double get_mS2() const { return mS2; }
   double get_mH2() const { return mH2; }
   double get_vH() const { return vH; }
   double get_vS() const { return vS; }


protected:
   double Lambda5{};
   double Lambda4{};
   double mS2{};
   double mH2{};
   double vH{};
   double vS{};


private:
   static const int numberOfParameters = 39;

   struct Soft_traces {
      double traceYdAdjYd{};
      double traceYeAdjYe{};
      double traceYuAdjYu{};
      double traceYdAdjYdYdAdjYd{};
      double traceYdAdjYuYuAdjYd{};
      double traceYeAdjYeYeAdjYe{};
      double traceYuAdjYuYuAdjYu{};

   };
   Soft_traces calc_soft_traces(int) const;

   double calc_beta_Lambda5_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda5_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda5_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda4_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda4_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Lambda4_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mS2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mS2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mS2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mH2_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mH2_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mH2_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vH_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vH_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vH_3_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vS_1_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vS_2_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_vS_3_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const NSM_soft_parameters&);

} // namespace flexiblesusy

#undef TRACE_STRUCT_TYPE

#endif
