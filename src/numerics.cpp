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

/** \file numerics.cpp
   - Project:     SOFTSUSY
   - Author:      Ben Allanach, Alexander Voigt
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

/// Comment if you want default softsusy behaviour
// #define USE_LOOPTOOLS

#include "numerics.h"
#include "numerics2.hpp"
#ifdef USE_LOOPTOOLS
#include <clooptools.h>
#endif
#include <algorithm>
#include <cmath>

namespace softsusy {

namespace {
/* constexpr --> const S.D. 
   constexpr is used to store numbers that are to big for available types? 
   with arbitrary precsiion this should not be a problem ?*/
const precise_real_type EPSTOL = 1.0e-11; ///< underflow accuracy
const precise_real_type TOL = 1e-4;

const precise_real_type dabs(precise_real_type a) noexcept { return a >= 0. ? a : -a; }
const precise_real_type sqr(precise_real_type a) noexcept { return a*a; }
const precise_real_type pow3(precise_real_type a) noexcept { return a*a*a; }
const precise_real_type pow6(precise_real_type a) noexcept { return a*a*a*a*a*a; }

const bool is_zero(precise_real_type m, precise_real_type tol) noexcept
{
   const precise_real_type am = dabs(m);
   const precise_real_type mtol = tol * am;

   if (mtol == 0.0 && am != 0.0 && tol != 0.0)
      return am <= tol;

   return am <= mtol;
}

const bool is_close(precise_real_type m1, precise_real_type m2, precise_real_type tol) noexcept
{
   const precise_real_type mmax = max(dabs(m1), dabs(m2));
   const precise_real_type mmin = min(dabs(m1), dabs(m2));
   const precise_real_type max_tol = tol * mmax;

   if (max_tol == 0.0 && mmax != 0.0 && tol != 0.0)
      return mmax - mmin <= tol;

   return mmax - mmin <= max_tol;
}

/// returns a/b if a/b is finite, otherwise returns numeric_limits::max()
template <typename T>
constexpr T divide_finite(T a, T b) noexcept {
   T result = a / b;
   if (!isfinite(result))
      result = std::numeric_limits<T>::max();
   return result;
}

// can be made constexpr in C++20
precise_real_type fB(const precise_complex_type& a) noexcept
{  
   precise_complex_type temp1,temp2,temp3;
   precise_real_type res;
   using flexiblesusy::fast_log;
   const precise_real_type x = a.real();

   if (fabs(x) < EPSTOL)
      return -1. - x + sqr(x) * 0.5;

   if (is_close(x, 1., EPSTOL))
      return -1.;

   temp1 = precise_complex_type(1.,0)-a;
   temp2 = precise_complex_type(1.0,0)-1.0 / a;
   temp3 = fast_log(temp1)-precise_complex_type(1.,0)-a*fast_log(temp2);
   res = real(temp3);
   return res;
}

} // anonymous namespace

precise_real_type a0(precise_real_type m, precise_real_type q) noexcept {
   //using std::fabs;
   //using std::log;
   const precise_real_type TOL = precise_real_type(1*pow(10,-4));
   if (fabs(m) < TOL) return 0.;
   return sqr(m) * (1.0 - 2. * log(fabs(m / q)));
}

precise_real_type ffn(precise_real_type p, precise_real_type m1, precise_real_type m2, precise_real_type q) noexcept {
   return a0(m1, q) - 2.0 * a0(m2, q) -
      (2.0 * sqr(p) + 2.0 * sqr(m1) - sqr(m2)) *
      b0(p, m1, m2, q);
}

precise_real_type gfn(precise_real_type p, precise_real_type m1, precise_real_type m2, precise_real_type q) noexcept {
   return (sqr(p) - sqr(m1) - sqr(m2)) * b0(p, m1, m2, q) - a0(m1, q)
      - a0(m2, q);
}

precise_real_type hfn(precise_real_type p, precise_real_type m1, precise_real_type m2, precise_real_type q) noexcept {
   return 4.0 * b22(p, m1, m2, q) + gfn(p, m1, m2, q);
}

precise_real_type b22bar(precise_real_type p, precise_real_type m1, precise_real_type m2, precise_real_type q) noexcept {
   return b22(p, m1, m2, q) - 0.25 * a0(m1, q) - 0.25 * a0(m2, q);
}

/**
 * Returns Re(B0(p,m1,m2,q)), from hep-ph/9606211
 */
precise_real_type b0(precise_real_type p, precise_real_type m1, precise_real_type m2, precise_real_type q) noexcept
{
   //using std::fabs;
  // using std::log;
  // using std::sqrt;

#ifdef USE_LOOPTOOLS
   setmudim(q*q);
   precise_real_type b0l = B0(p*p, m1*m1, m2*m2).real();
   //  return B0(p*p, m1*m1, m2*m2).real();
#endif

   // protect against infrared divergence
   if (is_zero(p, EPSTOL) && is_zero(m1, EPSTOL) && is_zero(m2, EPSTOL))
      return 0.0;

   const precise_real_type mMin = min(fabs(m1), fabs(m2));
   const precise_real_type mMax = max(fabs(m1), fabs(m2));

   const precise_real_type pSq = sqr(p), mMinSq = sqr(mMin), mMaxSq = sqr(mMax);
   /// Try to increase the accuracy of s
   const precise_real_type dmSq = mMaxSq - mMinSq;
   const precise_real_type s = pSq + dmSq;

   const precise_real_type pTest = divide_finite(pSq, mMaxSq);
   /// Decides level at which one switches to p=0 limit of calculations
   const precise_real_type pTolerance = precise_real_type(1*pow(10,-10));

   /// p is not 0
   if (pTest > pTolerance) {
      const precise_complex_type ieps(0.0, EPSTOL * mMaxSq);
      const precise_complex_type x = s + sqrt(sqr(s) - 4. * pSq * (mMaxSq - ieps));
      const precise_complex_type xPlus = x / (2. * pSq);
      const precise_complex_type xMinus = 2. * (mMaxSq - ieps) / x;

      return -2.0 * log(p / q) - fB(xPlus) - fB(xMinus);
   }

   if (is_close(m1, m2, EPSTOL)) {
      return - log(sqr(m1 / q));
   }

   const precise_real_type Mmax2 = mMaxSq, Mmin2 = mMinSq;

   if (Mmin2 < 1.e-30) {
      return 1.0 - log(Mmax2 / sqr(q));
   }

   return 1.0 - log(Mmax2 / sqr(q)) + Mmin2 * log(Mmax2 / Mmin2)
      / (Mmin2 - Mmax2);
}

/// Note that b1 is NOT symmetric in m1 <-> m2!!!
precise_real_type b1(precise_real_type p, precise_real_type m1, precise_real_type m2, precise_real_type q) noexcept
{
   //using std::fabs;
  //using std::log;

#ifdef USE_LOOPTOOLS
   setmudim(q*q);
   precise_real_type b1l = -B1(p*p, m1*m1, m2*m2).real();
   //    return b1l;
#endif

   // protect against infrared divergence
   if (is_zero(p, EPSTOL) && is_zero(m1, EPSTOL) && is_zero(m2, EPSTOL))
      return 0.0;

   const precise_real_type p2 = sqr(p), m12 = sqr(m1), m22 = sqr(m2), q2 = sqr(q);
   const precise_real_type pTest = divide_finite(p2, max(m12, m22));

   /// Decides level at which one switches to p=0 limit of calculations
   const precise_real_type pTolerance = precise_real_type(1*pow(10,-4));

   if (pTest > pTolerance) {
      return (a0(m2, q) - a0(m1, q) + (p2 + m12 - m22)
              * b0(p, m1, m2, q)) / (2.0 * p2);
   }

   if (fabs(m1) > 1.0e-15 && fabs(m2) > 1.0e-15) {
      const precise_real_type m14 = sqr(m12), m24 = sqr(m22);
      const precise_real_type m16 = m12*m14 , m26 = m22*m24;
      const precise_real_type m18 = sqr(m14), m28 = sqr(m24);
      const precise_real_type p4 = sqr(p2);

      if (fabs(m12 - m22) < pTolerance * max(m12, m22)) {
         return 0.08333333333333333*p2/m22
            + 0.008333333333333333*p4/m24
            + sqr(m12 - m22)*(0.041666666666666664/m24 +
                              0.016666666666666666*p2/m26 +
                              0.005357142857142856*p4/m28)
            + (m12 - m22)*(-0.16666666666666666/m22 -
                           0.03333333333333333*p2/m24 -
                           0.007142857142857142*p4/m26)
            - 0.5*log(m22/q2);
      }

      const precise_real_type l12 = log(m12/m22);

      return (3*m14 - 4*m12*m22 + m24 - 2*m14*l12)/(4.*sqr(m12 - m22))
         + (p2*(4*pow3(m12 - m22)*
                (2*m14 + 5*m12*m22 - m24) +
                (3*m18 + 44*m16*m22 - 36*m14*m24 - 12*m12*m26 + m28)*p2
                - 12*m14*m22*(2*sqr(m12 - m22) + (2*m12 + 3*m22)*p2)*l12))/
         (24.*pow6(m12 - m22)) - 0.5*log(m22/q2);
   }

   return (m12 > m22)
      ? -0.5*log(m12/q2) + 0.75
      : -0.5*log(m22/q2) + 0.25;
}

precise_real_type b22(precise_real_type p,  precise_real_type m1, precise_real_type m2, precise_real_type q) noexcept
{
   //using std::fabs;
   //using std::log;

#ifdef USE_LOOPTOOLS
   setmudim(q*q);
   precise_real_type b22l = B00(p*p, m1*m1, m2*m2).real();
#endif

   // protect against infrared divergence
   if (is_zero(p, EPSTOL) && is_zero(m1, EPSTOL) && is_zero(m2, EPSTOL))
      return 0.0;

   /// Decides level at which one switches to p=0 limit of calculations
   const precise_real_type p2 = sqr(p), m12 = sqr(m1), m22 = sqr(m2);
   const precise_real_type pTolerance = precise_real_type(1*pow(10,-10));

   if (p2 < pTolerance * max(m12, m22) ) {
      // m1 == m2 with good accuracy
      if (is_close(m1, m2, EPSTOL)) {
         return -m12 * log(sqr(m1 / q)) * 0.5 + m12 * 0.5;
      }
      // p == 0 limit
      if (fabs(m1) > EPSTOL && fabs(m2) > EPSTOL) {
         return 0.375 * (m12 + m22) - 0.25 *
            (sqr(m22) * log(sqr(m2 / q)) - sqr(m12) *
             log(sqr(m1 / q))) / (m22 - m12);
      }
      return (fabs(m1) < EPSTOL)
         ? 0.375 * m22 - 0.25 * m22 * log(sqr(m2 / q))
         : 0.375 * m12 - 0.25 * m12 * log(sqr(m1 / q));
   }

   const precise_real_type b0Save = b0(p, m1, m2, q);
   const precise_real_type a01 = a0(m1, q);
   const precise_real_type a02 = a0(m2, q);

   return 1.0 / 6.0 *
      (0.5 * (a01 + a02) + (m12 + m22 - 0.5 * p2)
       * b0Save + (m22 - m12) / (2.0 * p2) *
       (a02 - a01 - (m22 - m12) * b0Save) +
       m12 + m22 - p2 / 3.0);
}

precise_real_type d0(precise_real_type m1, precise_real_type m2, precise_real_type m3, precise_real_type m4) noexcept
{
   //using std::log;

   const precise_real_type m1sq = sqr(m1);
   const precise_real_type m2sq = sqr(m2);

   if (is_close(m1, m2, EPSTOL)) {
      const precise_real_type m3sq = sqr(m3);
      const precise_real_type m4sq = sqr(m4);

      if (is_zero(m2, EPSTOL)) {
         // d0 is undefined for m1 == m2 == 0
         return 0.;
      } else if (is_zero(m3, EPSTOL)) {
         return (-m2sq + m4sq - m2sq * log(m4sq/m2sq))/
            sqr(m2 * m2sq - m2 * m4sq);
      } else if (is_zero(m4, EPSTOL)) {
         return (-m2sq + m3sq - m2sq * log(m3sq/m2sq))/
            sqr(m2 * m2sq - m2 * m3sq);
      } else if (is_close(m2, m3, EPSTOL) && is_close(m2, m4, EPSTOL)) {
         return 1.0 / (6.0 * sqr(m2sq));
      } else if (is_close(m2, m3, EPSTOL)) {
         return (sqr(m2sq) - sqr(m4sq) + 2.0 * m4sq * m2sq * log(m4sq / m2sq)) /
            (2.0 * m2sq * sqr(m2sq - m4sq) * (m2sq - m4sq));
      } else if (is_close(m2, m4, EPSTOL)) {
         return (sqr(m2sq) - sqr(m3sq) + 2.0 * m3sq * m2sq * log(m3sq / m2sq)) /
            (2.0 * m2sq * sqr(m2sq - m3sq) * (m2sq - m3sq));
      } else if (is_close(m3, m4, EPSTOL)) {
         return -1.0 / sqr(m2sq - m3sq) *
            ((m2sq + m3sq) / (m2sq - m3sq) * log(m3sq / m2sq) + 2.0);
      }

      return
         (m4sq / sqr(m2sq - m4sq) * log(m4sq / m2sq) +
          m4sq / (m2sq * (m2sq - m4sq)) -
          m3sq / sqr(m2sq - m3sq) * log(m3sq / m2sq) -
          m3sq / (m2sq * (m2sq - m3sq))) / (m3sq - m4sq);
   }
   return (c0(m1, m3, m4) - c0(m2, m3, m4)) / (m1sq - m2sq);
}

precise_real_type d27(precise_real_type m1, precise_real_type m2, precise_real_type m3, precise_real_type m4) noexcept
{
   if (is_close(m1, m2, EPSTOL))
      m1 += TOL * 0.01;

   const precise_real_type m12 = sqr(m1), m22 = sqr(m2);

   return (m12 * c0(m1, m3, m4) - m22 * c0(m2, m3, m4))
      / (4.0 * (m12 - m22));
}

precise_real_type c0(precise_real_type m1, precise_real_type m2, precise_real_type m3) noexcept
{
   //using std::log;

#ifdef USE_LOOPTOOLS
   precise_real_type q = 100.;
   setmudim(q*q);
   precise_real_type psq = 0.;
   precise_real_type c0l = C0(psq, psq, psq, m1*m1, m2*m2, m3*m3).real();
#endif

   const precise_real_type m12 = sqr(m1), m22 = sqr(m2), m32 = sqr(m3);
   const bool m1_is_zero = is_zero(m1, EPSTOL);
   const bool m2_is_zero = is_zero(m2, EPSTOL);
   const bool m3_is_zero = is_zero(m3, EPSTOL);

   if ((m1_is_zero && m2_is_zero && m3_is_zero) ||
       (m2_is_zero && m3_is_zero) ||
       (m1_is_zero && m3_is_zero) ||
       (m1_is_zero && m2_is_zero)) {
      return 0.;
   }

   if (m1_is_zero) {
      if (is_close(m2,m3,EPSTOL))
         return -1./m22;
      return log(m32/m22)/(m22 - m32);
   }

   if (m2_is_zero) {
      if (is_close(m1,m3,EPSTOL))
         return -1./m12;
      return log(m32/m12)/(m12 - m32);
   }

   if (m3_is_zero) {
      if (is_close(m1,m2,EPSTOL))
         return -1./m12;
      return log(m22/m12)/(m12 - m22);
   }

   if (is_close(m2, m3, EPSTOL)) {
      if (is_close(m1, m2, EPSTOL))
         return - 0.5 / m22;
      return m12 / sqr(m12-m22) * log(m22/m12) + 1.0 / (m12 - m22);
   }

   if (is_close(m1, m2, EPSTOL)) {
      return - (1.0 + m32 / (m22-m32) * log(m32/m22)) / (m22-m32);
   }

   if (is_close(m1, m3, EPSTOL)) {
      return - (1.0 + m22 / (m32-m22) * log(m22/m32)) / (m32-m22);
   }

   return (m22 / (m12 - m22) * log(m22 / m12) -
           m32 / (m12 - m32) * log(m32 / m12)) / (m22 - m32);
}

/**
 * Derivative of B0(p^2, m1^2, m2^2, Q^2) w.r.t. p^2.
 *
 * @note Implemented only in the p^2 = 0 limit.
 *
 * @param p2 squared momentum
 * @param m2a squared mass
 * @param m2b squared mass
 *
 * @return derivative of B0 w.r.t. p^2
 */
precise_real_type d1_b0(precise_real_type /* p2 */, precise_real_type m2a, precise_real_type m2b) noexcept
{
   //using std::abs;

   const precise_real_type m4a = m2a * m2a;
   const precise_real_type m4b = m2b * m2b;

   if ((abs(m2a) < 0.0001) != (abs(m2b) < 0.0001)) {
      return (m4a - m4b) / (2. * pow3(m2a - m2b));
   } else if (abs(m2a) < 0.0001 && abs(m2b) < 0.0001) {
      return 0.;
   } else if (abs(m2b - m2a) < 0.001) {
      return 1./(6. * m2a) + (m2a - m2b)/(12.* m4a);
   }

   return (m4a - m4b + 2. * m2a * m2b * log(m2b/m2a))
      /(2. * pow3(m2a - m2b));
}

} // namespace softsusy
