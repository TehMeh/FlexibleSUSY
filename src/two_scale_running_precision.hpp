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

#ifndef TWO_SCALE_RUNNING_PRECISION_H
#define TWO_SCALE_RUNNING_PRECISION_H

namespace flexiblesusy {

class Two_scale_running_precision {
public:
   virtual ~Two_scale_running_precision() = default;
   virtual precise_real_type get_precision(int) = 0;
};

class Two_scale_constant_precision : public Two_scale_running_precision {
public:
   explicit Two_scale_constant_precision(precise_real_type);
   virtual ~Two_scale_constant_precision();
   virtual precise_real_type get_precision(int);
private:
   precise_real_type precision;
};

class Two_scale_increasing_precision : public Two_scale_running_precision {
public:
   Two_scale_increasing_precision(precise_real_type, precise_real_type);
   virtual ~Two_scale_increasing_precision();
   virtual precise_real_type get_precision(int);
private:
   precise_real_type decreasing_factor;
   precise_real_type minimum_precision;
};

} // namespace flexiblesusy

#endif
