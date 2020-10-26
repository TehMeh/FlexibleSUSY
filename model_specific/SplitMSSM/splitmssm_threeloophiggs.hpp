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

#ifndef SPLIT_THREELOOPHIGGS_H
#define SPLIT_THREELOOPHIGGS_H

namespace flexiblesusy {
namespace splitmssm_threeloophiggs {

/// Higgs self-energy 3-loop, gluino contribution O(alpha_t alpha_s^2)
precise_real_type delta_mh_3loop_gluino_split(
   precise_real_type scale, precise_real_type mt, precise_real_type yt, precise_real_type g3, precise_real_type mg);

} // namespace splitmssm_threeloophiggs
} // namespace flexiblesusy

#endif
