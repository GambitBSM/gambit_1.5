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

// File generated at Sun 24 Sep 2017 15:52:33

/**
 * @file MSSM_mAmu_ewsb_solver_interface.hpp
 *
 * @brief contains interface for EWSB solver
 *
 * This file was generated at Sun 24 Sep 2017 15:52:33 with FlexibleSUSY
 * 2.0.0-dev (git commit: 4d4c39a2702e9a6604f84813ccb0b85d40987f3b) and SARAH 4.11.0 .
 */

#ifndef MSSM_mAmu_EWSB_SOLVER_INTERFACE_H
#define MSSM_mAmu_EWSB_SOLVER_INTERFACE_H

namespace flexiblesusy {

class MSSM_mAmu_mass_eigenstates;

/**
 * @class MSSM_mAmu_ewsb_solver_interface
 * @brief interface for EWSB solvers to be used to solve the EWSB equations
 */
class MSSM_mAmu_ewsb_solver_interface {
public:
   virtual ~MSSM_mAmu_ewsb_solver_interface() {}

   virtual void set_loop_order(int) = 0;
   virtual void set_number_of_iterations(int) = 0;
   virtual void set_precision(double) = 0;

   virtual int get_loop_order() const = 0;
   virtual int get_number_of_iterations() const = 0;
   virtual double get_precision() const = 0;

   virtual int solve(MSSM_mAmu_mass_eigenstates&) = 0;
};

} // namespace flexiblesusy

#endif
