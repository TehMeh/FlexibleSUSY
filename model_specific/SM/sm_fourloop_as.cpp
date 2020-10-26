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

#include "sm_fourloop_as.hpp"
#include <cmath>
#include <ostream>

namespace flexiblesusy {
namespace sm_fourloop_as {

namespace {
   constexpr precise_real_type Pi = 3.1415926535897932384626433832795_p;
   template <typename T> constexpr T power2(T x) { return x*x; }
   template <typename T> constexpr T power3(T x) { return x*x*x; }
   template <typename T> constexpr T power4(T x) { return x*x*x*x; }
} // anonymous namespace

/**
 * 1-loop O(alpha_s) contributions to Delta alpha_s, Eq (23) of
 * [hep-ph/0004189]
 *
 * @param pars parameters
 *
 * @return 1-loop Delta alpha_s
 */
precise_real_type delta_alpha_s_1loop_as(const Parameters& pars)
{
   const precise_real_type as = pars.as;
   const precise_real_type mt2 = power2(pars.mt);
   const precise_real_type Q2 = power2(pars.Q);
   const precise_real_type L = log(Q2/mt2);

   return as / Pi * (1./6. * L);
}

/**
 * 2-loop O(alpha_s^2) contributions to Delta alpha_s, Eq (23) of
 * [hep-ph/0004189]
 *
 * @param pars parameters
 *
 * @return 2-loop Delta alpha_s
 */
precise_real_type delta_alpha_s_2loop_as_as(const Parameters& pars)
{
   const precise_real_type as = pars.as;
   const precise_real_type mt2 = power2(pars.mt);
   const precise_real_type Q2 = power2(pars.Q);
   const precise_real_type L = log(Q2/mt2);

   return power2(as / Pi) * (-11._p/72._p + 11._p/24_p*L + 1._p/36._p * power2(L));
}

/**
 * 3-loop O(alpha_s^3) contributions to Delta alpha_s, Eq (23) of
 * [hep-ph/0004189]
 *
 * @param pars parameters
 *
 * @return 3-loop Delta alpha_s
 */
precise_real_type delta_alpha_s_3loop_as_as_as(const Parameters& pars)
{
   const precise_real_type as = pars.as;
   const precise_real_type mt2 = power2(pars.mt);
   const precise_real_type Q2 = power2(pars.Q);
   const precise_real_type L = log(Q2/mt2);
   const precise_real_type L2 = power2(L);
   const precise_real_type L3 = power3(L);
   const precise_real_type nl = 5._p;
   const precise_real_type zeta3 = 1.202056903159594_p;

   return power3(as / Pi) * (
      - 564731._p/124416._p
      + 82043._p/27648._p * zeta3
      + 2645._p/1728._p * L
      + 167._p/576._p * L2
      + 1._p/216._p * L3
      + nl * (2633._p/31104._p - 67._p/576._p * L + 1._p/36._p * L2)
   );
}

/**
 * 4-loop O(alpha_s^4) contributions to Delta alpha_s, Eq (38) of
 * [hep-ph/0512060]
 *
 * @param pars parameters
 *
 * @return 4-loop Delta alpha_s
 */
precise_real_type delta_alpha_s_4loop_as_as_as_as(const Parameters& pars)
{
   const precise_real_type as = pars.as;
   const precise_real_type mt2 = power2(pars.mt);
   const precise_real_type Q2 = power2(pars.Q);
   const precise_real_type L = log(Q2/mt2);
   const precise_real_type L2 = power2(L);
   const precise_real_type L3 = power3(L);
   const precise_real_type L4 = power4(L);
   const precise_real_type nl = 5._p;
   const precise_real_type nl2 = power2(nl);
   const precise_real_type zeta3 = 1.202056903159594_p;
   const precise_real_type delta_MS_4 =
      5.170346990805882_p - 1.00993152453019_p * nl - 0.0219783748689228_p * nl2;

   return power4(as / Pi) * (
      + 121._p/1728._p
      - delta_MS_4
      - 11093717._p/746496._p * L
      + 3022001._p/165888._p * zeta3 * L
      + 1837._p/1152._p * L2
      + 2909._p/10368._p * L3
      + 1._p/1296._p * L4
      + nl * (
         + 141937._p/373248._p * L
         - 110779._p/82944._p * zeta3 * L
         + 277._p/10368._p * L2
         + 271._p/5184._p * L3
      )
      + nl2 * (
         - 6865._p/186624._p * L
         + 77._p/20736._p * L2
         - 1._p/324._p * L3
      )
   );
}

/**
 * Calculates alpha_s(SM(6)) from alpha_s(SM(5)) in QCD using
 * (inverted) Eq (14) from [hep-ph/0512060].
 *
 * @param pars parameters
 *
 * @return alpha_s(SM(6))
 */
precise_real_type calc_alpha_s(const Parameters& pars, int loops)
{
   const auto as = pars.as; // SM(5)
   const auto d1 = loops < 1 ? 0.0 : delta_alpha_s_1loop_as(pars);
   const auto d2 = loops < 2 ? 0.0 : delta_alpha_s_2loop_as_as(pars);
   const auto d3 = loops < 3 ? 0.0 : delta_alpha_s_3loop_as_as_as(pars);
   const auto d4 = loops < 4 ? 0.0 : delta_alpha_s_4loop_as_as_as_as(pars);

   return as * (1.0_p + d1 + d2 + d3 + d4);
}


/**
 * Calculates alpha_s(SM(6)) from alpha_s(SM(5)) in QCD using
 * the expressions given in [hep-ph/0512060].
 *
 * @param pars parameters
 *
 * @return alpha_s(SM(6))
 */
precise_real_type calc_alpha_s_alternative(const Parameters& pars, int loops)
{
   const auto as = pars.as; // SM(5)

   const auto d1 = loops < 1 ? 0.0 : delta_alpha_s_1loop_as(pars);
   const auto d2 = loops < 2 ? 0.0 : delta_alpha_s_2loop_as_as(pars);
   const auto d3 = loops < 3 ? 0.0 : delta_alpha_s_3loop_as_as_as(pars);
   const auto d4 = loops < 4 ? 0.0 : delta_alpha_s_4loop_as_as_as_as(pars);

   const auto del1 = d1;
   const auto del2 = d2 - power2(d1);
   const auto del3 = d3 + power3(d1) - 2.0_p*d1*d2;
   const auto del4 = d4 - 2.0_p*d1*d3 - power2(d2) + 3.0_p*power2(d1)*d2 - power4(d1);

   const auto delta = del1 + del2 + del3 + del4;

   return as / (1.0_p - delta);
}

std::ostream& operator<<(std::ostream& out, const Parameters& pars)
{
   out <<
      "Delta alpha_s(SM) parameters:\n"
      "alpha_s(SM(5)) = " <<  pars.as   << '\n' <<
      "mt             = " <<  pars.mt   << '\n' <<
      "Q              = " <<  pars.Q    << '\n';

   return out;
}

} // namespace sm_fourloop_as
} // namespace flexiblesusy
