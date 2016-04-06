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

// File generated at Mon 22 Feb 2016 17:30:33

#ifndef SingletDMZ3_UTILITIES_H
#define SingletDMZ3_UTILITIES_H

#include "SingletDMZ3_two_scale_model.hpp"
#include "SingletDMZ3_info.hpp"
#include "wrappers.hpp"

#include <Eigen/Core>
#include <string>
#include <vector>
#include <valarray>
#include <utility>

namespace flexiblesusy {

class SingletDMZ3_parameter_getter {
public:
   Eigen::ArrayXd get_parameters(const SingletDMZ3_mass_eigenstates& model) {
      return model.get();
   }
   std::vector<std::string> get_parameter_names(const SingletDMZ3_mass_eigenstates&) const {
      using namespace SingletDMZ3_info;
      return std::vector<std::string>(parameter_names,
                                      parameter_names + NUMBER_OF_PARAMETERS);
   }
};

class SingletDMZ3_spectrum_plotter {
public:
   SingletDMZ3_spectrum_plotter();
   ~SingletDMZ3_spectrum_plotter() {}

   void extract_spectrum(const SingletDMZ3_mass_eigenstates&);
   void write_to_file(const std::string&) const;

private:
   struct TParticle {
      std::string name;
      std::string latex_name;
      std::valarray<double> masses;
      TParticle(const std::string& name_, const std::string& latex_name_,
                const std::valarray<double>& masses_)
         : name(name_)
         , latex_name(latex_name_)
         , masses(masses_)
         {}
   };
   typedef std::vector<TParticle> TSpectrum;
   TSpectrum spectrum;
   double scale;
   unsigned width;

   void write_spectrum(const TSpectrum&, std::ofstream&) const;
   static std::valarray<double> to_valarray(double);
   template <class Scalar, int M, int N>
   static std::valarray<double> to_valarray(const Eigen::Array<Scalar, M, N>&);
};

template <class Scalar, int M, int N>
std::valarray<double> SingletDMZ3_spectrum_plotter::to_valarray(const Eigen::Array<Scalar, M, N>& v)
{
   return std::valarray<double>(v.data(), v.size());
}

} // namespace flexiblesusy

#endif