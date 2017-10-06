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

#include "sm_threeloop_as.hpp"
#include <cmath>
#include <ostream>

namespace flexiblesusy {
namespace sm_threeloop_as {

namespace {
   const double Pi = 3.1415926535897932384626433832795;
   template <typename T> T power2(T x)  { return x*x; }
   template <typename T> T power3(T x)  { return x*x*x; }
} // anonymous namespace

/**
 * 1-loop O(alpha_s) contributions to Delta alpha_s, Eq (23) of
 * [hep-ph/0004189]
 *
 * @param pars parameters
 *
 * @return 1-loop Delta alpha_s
 */
double delta_alpha_s_1loop_as(const Parameters& pars)
{
   const double as = pars.as;
   const double mt2 = power2(pars.mt);
   const double Q2 = power2(pars.Q);
   const double L = std::log(Q2/mt2);

   return as / Pi * (1./6. * L);
}

/**
 * 2-loop O(alpha_s^2) contributions to Delta alpha_s, Eq (23) of
 * [hep-ph/0004189]
 *
 * @param pars parameters
 *
 * @return 2-loop Delta alpha_s
 */
double delta_alpha_s_2loop_as_as(const Parameters& pars)
{
   const double as = pars.as;
   const double mt2 = power2(pars.mt);
   const double Q2 = power2(pars.Q);
   const double L = std::log(Q2/mt2);

   return power2(as / Pi) * (-11./72. + 11./24*L + 1./36. * power2(L));
}

/**
 * 3-loop O(alpha_s^3) contributions to Delta alpha_s, Eq (23) of
 * [hep-ph/0004189]
 *
 * @param pars parameters
 *
 * @return 3-loop Delta alpha_s
 */
double delta_alpha_s_3loop_as_as_as(const Parameters& pars)
{
   const double as = pars.as;
   const double mt2 = power2(pars.mt);
   const double Q2 = power2(pars.Q);
   const double L = std::log(Q2/mt2);
   const double L2 = power2(L);
   const double L3 = power3(L);
   const double nl = 5.;
   const double zeta3 = 1.202056903159594;

   return power3(as / Pi) * (
      - 564731./124416.
      + 82043./27648. * zeta3
      + 2645./1728. * L
      + 167./576. * L2
      + 1./216. * L3
      + nl * (2633./31104. - 67./576. * L + 1./36. * L2)
   );
}

std::ostream& operator<<(std::ostream& out, const Parameters& pars)
{
   out <<
      "Delta alpha_s(SM) 2L parameters:\n"
      "alpha_s(SM(5)) = " <<  pars.as   << '\n' <<
      "mt             = " <<  pars.mt   << '\n' <<
      "Q              = " <<  pars.Q    << '\n';

   return out;
}

} // namespace sm_threeloop_as
} // namespace flexiblesusy
