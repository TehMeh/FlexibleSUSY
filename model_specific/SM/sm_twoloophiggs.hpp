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

#ifndef SM_TWOLOOPHIGGS_H
#define SM_TWOLOOPHIGGS_H

namespace flexiblesusy {
namespace sm_twoloophiggs {

/// SM Higgs 1-loop contribution
precise_real_type delta_mh_1loop_sm(
   precise_real_type p, precise_real_type scale, precise_real_type mt, precise_real_type yt,
   precise_real_type v, precise_real_type gY, precise_real_type g2, precise_real_type lambda);

/// SM Higgs 1-loop contribution, only O(alpha_t)
precise_real_type delta_mh_1loop_at_sm(
   precise_real_type p, precise_real_type scale, precise_real_type mt, precise_real_type yt);

/// SM Higgs self-energy 2-loop, only O(alpha_t alpha_s)
precise_real_type self_energy_higgs_2loop_at_as_sm(
   precise_real_type p2, precise_real_type scale, precise_real_type mt, precise_real_type yt, precise_real_type g3);

/// SM Higgs self-energy 2-loop, only O(alpha_b alpha_s)
precise_real_type self_energy_higgs_2loop_ab_as_sm(
   precise_real_type p2, precise_real_type scale, precise_real_type mb, precise_real_type yb, precise_real_type g3);

/// SM Higgs self-energy 2-loop, only O((alpha_b + alpha_t)^2)
precise_real_type self_energy_higgs_2loop_at_at_sm(
   precise_real_type p2, precise_real_type scale, precise_real_type mt, precise_real_type yt, precise_real_type mb);

/// SM Higgs self-energy 2-loop, only O(alpha_tau^2)
precise_real_type self_energy_higgs_2loop_atau_atau_sm(
   precise_real_type p2, precise_real_type scale, precise_real_type mtau, precise_real_type ytau);

/// SM Higgs tadpole 1-loop, only O(alpha_t)
precise_real_type tadpole_higgs_1loop_at_sm(
   precise_real_type scale, precise_real_type mt, precise_real_type yt);

/// SM Higgs tadpole 2-loop, only O(alpha_t alpha_s)
precise_real_type tadpole_higgs_2loop_at_as_sm(
   precise_real_type scale, precise_real_type mt, precise_real_type yt, precise_real_type g3);

/// SM Higgs tadpole 2-loop, only O(alpha_b alpha_s)
precise_real_type tadpole_higgs_2loop_ab_as_sm(
   precise_real_type scale, precise_real_type mb, precise_real_type yb, precise_real_type g3);

/// SM Higgs tadpole 2-loop, only O((alpha_b + alpha_t)^2)
precise_real_type tadpole_higgs_2loop_at_at_sm(
   precise_real_type scale, precise_real_type mt, precise_real_type yt, precise_real_type mb);

/// SM Higgs tadpole 2-loop, only O(alpha_tau^2)
precise_real_type tadpole_higgs_2loop_atau_atau_sm(
   precise_real_type scale, precise_real_type mtau, precise_real_type ytau);

/// SM Higgs 2-loop contribution, only O(alpha_t alpha_s)
precise_real_type delta_mh_2loop_at_as_sm(
   precise_real_type p2, precise_real_type scale, precise_real_type mt, precise_real_type yt, precise_real_type g3);

/// SM Higgs 2-loop contribution, only O(alpha_t alpha_s)
precise_real_type delta_mh_2loop_ab_as_sm(
   precise_real_type p2, precise_real_type scale, precise_real_type mb, precise_real_type yb, precise_real_type g3);

/// SM Higgs 2-loop contribution, only O((alpha_b + alpha_t)^2)
precise_real_type delta_mh_2loop_at_at_sm(
   precise_real_type p2, precise_real_type scale, precise_real_type mt, precise_real_type yt, precise_real_type mb);

/// SM Higgs 2-loop contribution, only O(alpha_tau^2)
precise_real_type delta_mh_2loop_atau_atau_sm(
   precise_real_type p2, precise_real_type scale, precise_real_type mtau, precise_real_type ytau);

/// SM Higgs 1-loop contribution from SUSYHD 1.0.2
precise_real_type delta_mh_1loop_sm_SUSYHD(
   precise_real_type vev, precise_real_type Mt, precise_real_type mh, precise_real_type MW, precise_real_type MZ, precise_real_type Q);

/// SM Higgs 2-loop contribution from SUSYHD 1.0.2
precise_real_type delta_mh_2loop_sm_SUSYHD(
   precise_real_type vev, precise_real_type Mt, precise_real_type Mh, precise_real_type g3);

} // namespace sm_twoloophiggs
} // namespace flexiblesusy

#endif
