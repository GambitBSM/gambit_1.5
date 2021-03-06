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

// File generated at Thu 10 May 2018 14:34:56

/**
 * @file MSSMatMGUTEFTHiggs_mAmu_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solving EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Thu 10 May 2018 14:34:56 with FlexibleSUSY
 * 2.0.1 (git commit: unknown) and SARAH 4.12.2 .
 */

#ifndef MSSMatMGUTEFTHiggs_mAmu_TWO_SCALE_H
#define MSSMatMGUTEFTHiggs_mAmu_TWO_SCALE_H

#include "MSSMatMGUTEFTHiggs_mAmu_model.hpp"
#include "MSSMatMGUTEFTHiggs_mAmu_mass_eigenstates.hpp"

#include "model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class MSSMatMGUTEFTHiggs_mAmu<Two_scale>
 * @brief model class with routines for determining masses and mixings and EWSB
 */
template<>
class MSSMatMGUTEFTHiggs_mAmu<Two_scale> : public Model, public MSSMatMGUTEFTHiggs_mAmu_mass_eigenstates {
public:
   explicit MSSMatMGUTEFTHiggs_mAmu(const MSSMatMGUTEFTHiggs_mAmu_input_parameters& input_ = MSSMatMGUTEFTHiggs_mAmu_input_parameters());
   MSSMatMGUTEFTHiggs_mAmu(const MSSMatMGUTEFTHiggs_mAmu&) = default;
   MSSMatMGUTEFTHiggs_mAmu(MSSMatMGUTEFTHiggs_mAmu&&) = default;
   virtual ~MSSMatMGUTEFTHiggs_mAmu() = default;
   MSSMatMGUTEFTHiggs_mAmu& operator=(const MSSMatMGUTEFTHiggs_mAmu&) = default;
   MSSMatMGUTEFTHiggs_mAmu& operator=(MSSMatMGUTEFTHiggs_mAmu&&) = default;

   // interface functions
   virtual void calculate_spectrum() override;
   virtual void clear_problems() override;
   virtual std::string name() const override;
   virtual void run_to(double scale, double eps = -1.0) override;
   virtual void print(std::ostream& out = std::cerr) const override;
   virtual void set_precision(double) override;
};

std::ostream& operator<<(std::ostream&, const MSSMatMGUTEFTHiggs_mAmu<Two_scale>&);

} // namespace flexiblesusy

#endif
