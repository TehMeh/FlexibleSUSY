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

#include "sm_threeloop_as.hpp"
#include <cmath>
#include <ostream>

namespace flexiblesusy {
namespace sm_threeloop_as {

namespace {
   const precise_real_type Pi = 3.1415926535897932384626433832795_p;
   template <typename T> T power2(T x)  { return x*x; }
   template <typename T> T power3(T x)  { return x*x*x; }
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

   return power2(precise_real_type(as / Pi)) * (-11._p/72._p + 11._p/24_p*L + 1._p/36._p * power2(L));
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
   const precise_real_type nl = 5.;
   const precise_real_type zeta3 = 1.202056903159594_p;

   return power3<precise_real_type>(as / Pi) * (
      - 564731._p/124416._p
      + 82043._p/27648._p * zeta3
      + 2645._p/1728._p * L
      + 167._p/576._p * L2
      + 1._p/216._p * L3
      + nl * (2633._p/31104._p - 67._p/576._p * L + 1._p/36._p * L2)
   );
}

std::ostream& operator<<(std::ostream& out, const Parameters& pars)
{
   out <<
      "Delta alpha_s(SM) 2L parameters:\n"
      "alpha_s(SM(5)) = " <<  pars.as   << '\n' <<
      "mt             = " <<  pars.mt   << '\n' <<
      "Q              = " <<  pars.Q    << '\n';

   return out;
}

} // namespace sm_threeloop_as
} // namespace flexiblesusy
