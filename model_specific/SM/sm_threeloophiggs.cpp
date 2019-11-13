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

#include "sm_threeloophiggs.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {
namespace sm_threeloophiggs {

/**
 * Standard Model Higgs self-energy 3-loop, \f$O(\alpha_t
 * \alpha_s^2)\f$.  Taken from arxiv:1407.4336, Eq. (3.3).
 *
 * @note The result contains the 3-loop tadpole diagrams.  It is
 * therefore not 1-particle irreducible (1PI).
 *
 * @param scale renormalization scale
 * @param mt MS-bar top mass
 * @param yt MS-bar Yukawa coupling
 * @param g3 MS-bar strong gauge coupling
 *
 * @return real part of 3-loop correction \f$O(\alpha_t \alpha_s^2)\f$
 */
precise_real_type delta_mh_3loop_at_as_as_sm(
   precise_real_type scale, precise_real_type mt, precise_real_type yt, precise_real_type g3)
{
   const precise_real_type yt2 = Sqr(yt);
   const precise_real_type mt2 = Sqr(mt);
   const precise_real_type g34 = Power4(g3);
   const precise_real_type Q2 = Sqr(scale);
   const precise_real_type LogT = FiniteLog(mt2 / Q2);
   const precise_real_type LogT2 = Sqr(LogT);
   const precise_real_type LogT3 = Power3(LogT);

   const precise_real_type result =
      g34*yt2*mt2*(248.1215180432007 + 839.1966169377614*LogT
                   + 160*LogT2 - 736*LogT3);

   return result * threeLoop;
}

/**
 * Standard Model Higgs self-energy 3-loop,
 * \f$O(\alpha_t^2\alpha_s)\f$.  Taken from arxiv:1407.4336, Eq. (3.3).
 *
 * @note The result contains the 3-loop tadpole diagrams.  It is
 * therefore not 1-particle irreducible (1PI).
 *
 * @param scale renormalization scale
 * @param mt MS-bar top mass
 * @param yt MS-bar Yukawa coupling
 * @param g3 MS-bar strong gauge coupling
 *
 * @return real part of 3-loop correction \f$O(\alpha_t^2\alpha_s) \f$
 */
precise_real_type delta_mh_3loop_at_at_as_sm(
   precise_real_type scale, precise_real_type mt, precise_real_type yt, precise_real_type g3)
{
   const precise_real_type yt4 = Power4(yt);
   const precise_real_type mt2 = Sqr(mt);
   const precise_real_type g32 = Sqr(g3);
   const precise_real_type Q2 = Sqr(scale);
   const precise_real_type LogT = FiniteLog(mt2 / Q2);
   const precise_real_type LogT2 = Sqr(LogT);
   const precise_real_type LogT3 = Power3(LogT);

   const precise_real_type result =
      g32*yt4*mt2*(2764.365124334015 + 1283.715638285500*LogT
                   - 360*LogT2 + 240*LogT3);

   return result * threeLoop;
}

/**
 * Standard Model Higgs self-energy 3-loop, \f$O(\alpha_t^3)\f$.
 * Taken from arxiv:1407.4336, Eq. (3.4).
 *
 * @note The result contains the 3-loop tadpole diagrams.  It is
 * therefore not 1-particle irreducible (1PI).
 *
 * @param scale renormalization scale
 * @param mt MS-bar top mass
 * @param yt MS-bar Yukawa coupling
 * @param mh MS-bar Higgs mass
 *
 * @return real part of 3-loop correction \f$O(\alpha_t^3) \f$
 */
precise_real_type delta_mh_3loop_at_at_at_sm(
   precise_real_type scale, precise_real_type mt, precise_real_type yt, precise_real_type mh)
{
   const precise_real_type yt6 = Power6(yt);
   const precise_real_type mt2 = Sqr(mt);
   const precise_real_type mh2 = Sqr(mh);
   const precise_real_type Q2 = Sqr(scale);
   const precise_real_type LogH = FiniteLog(mh2 / Q2);
   const precise_real_type LogT = FiniteLog(mt2 / Q2);
   const precise_real_type LogT2 = Sqr(LogT);
   const precise_real_type LogT3 = Power3(LogT);

   const precise_real_type result =
      yt6*mt2*(-3199.016554815089 + 36*LogH - 2653.510765697467*LogT
               + 756*LogH*LogT + 27.*0.5*LogT2 + 324*LogH*LogT2 - 225*LogT3);

   return result * threeLoop;
}

} // namespace sm_threeloophiggs
} // namespace flexiblesusy
