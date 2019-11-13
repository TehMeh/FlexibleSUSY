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

#include "scan.hpp"
#include <cmath>
#include <cstddef>

/**
 * @file scan.cpp
 * @brief contains helper functions and classes for parameter scans
 */

namespace flexiblesusy {

/**
 * Returns vector with number_of_steps floating point values between
 * start and stop.  The endpoint is excluded.
 *
 * @param start smallest value
 * @param stop largest value (excluded)
 * @param number_of_steps number of values
 *
 * @return vector of floating point values
 */
std::vector<precise_real_type> float_range(precise_real_type start, precise_real_type stop,
                                std::size_t number_of_steps)
{
   const precise_real_type step_size = (stop - start) / number_of_steps;
   std::vector<precise_real_type> result(number_of_steps);

   for (std::size_t i = 0; i < number_of_steps; ++i) {
      const precise_real_type point = start + i * step_size;
      result[i] = point;
   }

   return result;
}

/**
 * Returns vector with number_of_steps floating point values between
 * start and stop.  The values are logarithmically distributed,
 * i.e. the difference of the logarithm of two consecutive points is
 * constant.  The endpoint is excluded.
 *
 * @param start smallest value
 * @param stop largest value (excluded)
 * @param number_of_steps number of values
 *
 * @return vector of floating point values
 */
std::vector<precise_real_type> float_range_log(precise_real_type start, precise_real_type stop,
                                    std::size_t number_of_steps)
{
   if (number_of_steps == 0)
      return std::vector<precise_real_type>();

   const precise_real_type step_size = (log(stop) - log(start)) / number_of_steps;
   std::vector<precise_real_type> result(number_of_steps);
   result[0] = start;

   for (std::size_t i = 1; i < number_of_steps; ++i) {
      const precise_real_type point = exp(step_size + log(result[i-1]));
      result[i] = point;
   }

   return result;
}

} // namespace flexiblesusy
