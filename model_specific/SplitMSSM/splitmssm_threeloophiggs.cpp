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

#include "splitmssm_threeloophiggs.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {
namespace splitmssm_threeloophiggs {

/**
 * Split-SUSY Higgs self-energy 3-loop contribution from a gluino,
 * \f$O(\alpha_t \alpha_s^2)\f$.  Taken from arxiv:1312.5220,
 * Eq. (4.8).
 *
 * @note The result contains the 3-loop tadpole diagrams.  It is
 * therefore not 1-particle irreducible (1PI).
 *
 * @warning The result is in Landau gauge (\f$\xi = 0\f$).
 *
 * @param scale renormalization scale
 * @param mt MS-bar top mass
 * @param yt MS-bar Yukawa coupling
 * @param g3 MS-bar strong gauge coupling
 * @param mg MS-bar gluino mass
 *
 * @return real part of 3-loop gluino contribution to Higgs self-energy
 */
precise_real_type delta_mh_3loop_gluino_split(
   precise_real_type scale, precise_real_type mt, precise_real_type yt, precise_real_type g3, precise_real_type mg)
{
   const precise_real_type yt2 = Sqr(yt);
   const precise_real_type mt2 = Sqr(mt);
   const precise_real_type g34 = pow(g3, 4);
   const precise_real_type mg2 = Sqr(mg);
   const precise_real_type Q2 = Sqr(scale);
   const precise_real_type LogG = FiniteLog(mg2 / Q2);
   const precise_real_type LogG3 = pow(LogG, 3);

   const precise_real_type result = 64*g34*yt2*mt2*LogG3;

   return result * threeLoop;
}

} // namespace splitmssm_threeloophiggs
} // namespace flexiblesusy
