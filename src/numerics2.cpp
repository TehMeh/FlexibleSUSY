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

#include "numerics2.hpp"
#include <cmath>
#include <cstddef>

namespace flexiblesusy {

bool is_finite(const precise_real_type* v, long length)
{
   bool is_finite = true;

   for (long i = 0; i < length; i++)
      is_finite = is_finite && isfinite(v[i]);

   return is_finite;
}

bool is_equal(precise_real_type a, precise_real_type b, 
	precise_real_type prec) noexcept
{
   return is_zero(precise_real_type(a - b), prec);
}

//template <typename T>
bool is_equal(precise_complex_type a, precise_complex_type b,
              precise_real_type prec) noexcept
{
   return (is_equal(a.real(), b.real(), prec)
           && is_equal(a.imag(), b.imag(), prec));
}


precise_complex_type fast_log(const precise_complex_type& z) noexcept
{
   return precise_complex_type(log(abs(z)), arg(z));
}

} // namespace flexiblesusy
