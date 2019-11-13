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

#include "precise.hpp"

#ifndef NUMERICS_HPP
#define NUMERICS_HPP

//#include <array>
#include <boost/array.hpp>
#include <cmath>
//#include <complex>
#include <limits>
#include <cstddef>
#include <cstdlib>
#include <type_traits>

namespace flexiblesusy {

template <typename T>
typename boost::enable_if_c<boost::is_unsigned<T>::value, bool>::type
is_zero(T a, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   return a <= prec;
}

template <typename T>
typename boost::enable_if_c<!std::is_unsigned<T>::value, bool>::type
is_zero(T a, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   return abs(a) <= prec;
}


bool is_equal(precise_real_type a, precise_real_type b,
  precise_real_type prec = std::numeric_limits<precise_real_type>::epsilon()) noexcept;

//template <typename T>
bool is_equal(precise_complex_type a, precise_complex_type b,
              precise_real_type prec = std::numeric_limits<precise_real_type>::epsilon()) noexcept;

template <typename T>
typename boost::enable_if_c<!std::is_unsigned<T>::value, bool>::type
is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   if (is_equal(a, b, std::numeric_limits<T>::epsilon()))
      return true;

   if (abs(a) < std::numeric_limits<T>::epsilon() ||
       abs(b) < std::numeric_limits<T>::epsilon())
      return false;

   return abs((a - b)/a) < prec;
}

template <typename T>
typename boost::enable_if_c<std::is_unsigned<T>::value, bool>::type
is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   using ST = typename std::make_signed<T>::type;
   const auto sa = static_cast<ST>(a);
   const auto sb = static_cast<ST>(b);
   const auto sprec = static_cast<ST>(prec);

   return is_equal_rel(sa, sb, sprec);
}

bool is_finite(const precise_real_type*, long length);

template <std::size_t N>
bool is_finite(const precise_real_type v[N]) noexcept
{
   bool is_finite = true;

   for (std::size_t i = 0; i < N; i++)
      is_finite = is_finite && isfinite(v[i]);

   return is_finite;
}

template <typename T, std::size_t N>
bool is_finite(const boost::array<T, N>& v) noexcept
{
   return is_finite<N>(&v[0]);
}

//template <class T>
precise_complex_type fast_log(const precise_complex_type& z) noexcept;
/*{
   return precise_complex_type(log(abs(z)), arg(z));
}	*/

} // namespace flexiblesusy

#endif
