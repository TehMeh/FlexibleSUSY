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

#ifndef SM_THREELOOPHIGGS_H
#define SM_THREELOOPHIGGS_H

namespace flexiblesusy {
namespace sm_threeloophiggs {

/// SM Higgs self-energy 3-loop, only O(alpha_t alpha_s^2)
precise_real_type delta_mh_3loop_at_as_as_sm(
   precise_real_type scale, precise_real_type mt, precise_real_type yt, precise_real_type g3);

/// SM Higgs self-energy 3-loop, only O(alpha_t^2 alpha_s)
precise_real_type delta_mh_3loop_at_at_as_sm(
   precise_real_type scale, precise_real_type mt, precise_real_type yt, precise_real_type g3);

/// SM Higgs self-energy 3-loop, only O(alpha_t^3)
precise_real_type delta_mh_3loop_at_at_at_sm(
   precise_real_type scale, precise_real_type mt, precise_real_type yt, precise_real_type mh);

} // namespace sm_threeloophiggs
} // namespace flexiblesusy

#endif
