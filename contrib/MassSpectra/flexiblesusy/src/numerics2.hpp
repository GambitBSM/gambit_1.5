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

#ifndef NUMERICS_HPP
#define NUMERICS_HPP

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <type_traits>

namespace flexiblesusy {

template <typename T>
typename std::enable_if<std::is_unsigned<T>::value, bool>::type
is_zero(T a, T prec = std::numeric_limits<T>::epsilon())
{
   return a <= prec;
}

template <typename T>
typename std::enable_if<!std::is_unsigned<T>::value, bool>::type
is_zero(T a, T prec = std::numeric_limits<T>::epsilon())
{
   return std::abs(a) <= prec;
}

template <typename T>
bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon())
{
   return is_zero(a - b, prec);
}

template <typename T>
bool is_equal(std::complex<T> a, std::complex<T> b,
              T prec = std::numeric_limits<T>::epsilon())
{
   return (is_equal(a.real(), b.real(), prec)
           && is_equal(a.imag(), b.imag(), prec));
}

template <typename T>
typename std::enable_if<!std::is_unsigned<T>::value, bool>::type
is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon())
{
   if (is_equal(a, b, std::numeric_limits<T>::epsilon()))
      return true;

   if (std::abs(a) < std::numeric_limits<T>::epsilon() ||
       std::abs(b) < std::numeric_limits<T>::epsilon())
      return false;

   return std::abs((a - b)/a) < prec;
}

template <typename T>
typename std::enable_if<std::is_unsigned<T>::value, bool>::type
is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon())
{
   using ST = typename std::make_signed<T>::type;
   const auto sa = static_cast<ST>(a);
   const auto sb = static_cast<ST>(b);
   const auto sprec = static_cast<ST>(prec);

   return is_equal_rel(sa, sb, sprec);
}

bool is_finite(const double*, long length);

template <std::size_t N>
bool is_finite(const double v[N])
{
   bool is_finite = true;

   for (std::size_t i = 0; i < N; i++)
      is_finite = is_finite && std::isfinite(v[i]);

   return is_finite;
}

template <typename T, std::size_t N>
bool is_finite(const std::array<T, N>& v)
{
   return is_finite<N>(&v[0]);
}

template <class T>
std::complex<T> fast_log(const std::complex<T>& z)
{
   const T r = std::abs(z);
   return std::complex<T>(std::log(r), std::arg(z));
}

} // namespace flexiblesusy

#endif
