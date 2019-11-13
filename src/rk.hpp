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

/**
 * @file rk.hpp
 * @brief Integration of ODEs by Runge-Kutta
 * @author Ben Allanach, Alexander Voigt
 *
 * The implementation of the Runge-Kutta routines have been derived
 * from SOFTSUSY [hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305].
 */

#include "precise.hpp"

#ifndef RK_H
#define RK_H

#include <algorithm>
#include <cmath>
#include <functional>

#include <Eigen/Core>

#include "logger.hpp"
#include "error.hpp"

namespace flexiblesusy {

namespace runge_kutta {

namespace {
/// Returns |a| with sign of b in front
inline precise_real_type sign(precise_real_type a, precise_real_type b) noexcept {
   return b >= 0 ? (precise_real_type)fabs(a) : (precise_real_type)(-fabs(a));
}
} // anonymous namespace

/// A single step of Runge Kutta (5th order), input:
/// y and dydx (derivative of y), x is independent variable. yout is value
/// after step. derivs is a user-supplied function
template <typename ArrayType, typename Derivs>
void rungeKuttaStep(const ArrayType& y, const ArrayType& dydx, precise_real_type x,
		    precise_real_type h, ArrayType& yout, ArrayType& yerr, Derivs derivs)
{
   const precise_real_type a2 = 0.2;
   const precise_real_type a3 = 0.3;
   const precise_real_type a4 = 0.6;
   const precise_real_type a5 = 1.0;
   const precise_real_type a6 = 0.875;
   const precise_real_type b21 = 0.2;
   const precise_real_type b31 = 3.0 / 40.0;
   const precise_real_type b32 = 9.0 / 40.0;
   const precise_real_type b41 = 0.3;
   const precise_real_type b42 = -0.9;
   const precise_real_type b43 = 1.2;
   const precise_real_type b51 = -11.0 / 54.0;
   const precise_real_type b52 = 2.5;
   const precise_real_type b53 = -70.0 / 27.0;
   const precise_real_type b54 = 35.0 / 27.0;
   const precise_real_type b61 = 1631.0 / 55296.0;
   const precise_real_type b62 = 175.0 / 512.0;
   const precise_real_type b63 = 575.0 / 13824.0;
   const precise_real_type b64 = 44275.0 / 110592.0;
   const precise_real_type b65 = 253.0 / 4096.0;
   const precise_real_type c1 = 37.0 / 378.0;
   const precise_real_type c3 = 250.0 / 621.0;
   const precise_real_type c4 = 125.0 / 594.0;
   const precise_real_type c6 = 512.0 / 1771.0;
   const precise_real_type dc5 = -277.00 / 14336.0;
   const precise_real_type dc1 = c1 - 2825.0 / 27648.0;
   const precise_real_type dc3 = c3 - 18575.0 / 48384.0;
   const precise_real_type dc4 = c4 - 13525.0 / 55296.0;
   const precise_real_type dc6 = c6 - 0.25;

   ArrayType ytemp = b21 * h * dydx + y;
   const ArrayType ak2 = derivs(x + a2 * h, ytemp);

   // Allowing piece-wise calculating of ytemp for speed reasons
   ytemp = y + h * (b31 * dydx + b32 * ak2);
   const ArrayType ak3 = derivs(x + a3 * h, ytemp);

   ytemp = y + h * (b41 * dydx + b42 * ak2 + b43 * ak3);
   const ArrayType ak4 = derivs(x+a4*h,ytemp);

   ytemp = y + h * (b51 * dydx + b52 * ak2 + b53 * ak3 + b54 * ak4);
   const ArrayType ak5 = derivs(x + a5 * h, ytemp);

   ytemp = y + h * (b61 * dydx + b62 * ak2 + b63 * ak3 + b64 * ak4 + b65 * ak5);
   const ArrayType ak6 = derivs(x + a6 * h, ytemp);

   yout = y + h * (c1 * dydx + c3 * ak3 + c4 * ak4 + c6 * ak6);
   yerr = h * (dc1 * dydx + dc3 * ak3 + dc4 * ak4 + dc5 * ak5 + dc6 * ak6);
}

/// organises the variable step-size for Runge-Kutta evolution
template <typename ArrayType, typename Derivs>
precise_real_type odeStepper(ArrayType& y, const ArrayType& dydx, precise_real_type& x, precise_real_type htry,
                  precise_real_type eps, const ArrayType& yscal, Derivs derivs,
                  int& max_step_dir)
{
   const precise_real_type SAFETY = 0.9;
   const precise_real_type PGROW = -0.2;
   const precise_real_type PSHRNK = -0.25;
   const precise_real_type ERRCON = 1.89e-4;
   const int n = y.size();
   precise_real_type errmax;
   precise_real_type h = htry;
   ArrayType yerr(n);
   ArrayType ytemp(n);

   for (;;) {
      rungeKuttaStep(y, dydx, x, h, ytemp, yerr, derivs);
      errmax = (yerr / yscal).abs().maxCoeff(&max_step_dir);
      errmax  /= eps;
      if (!isfinite(errmax)) {
#ifdef ENABLE_VERBOSE
         ERROR("odeStepper: non-perturbative running at Q = "
               << exp(x) << " GeV of parameter y(" << max_step_dir
               << ") = " << y(max_step_dir) << ", dy(" << max_step_dir
               << ")/dx = " << dydx(max_step_dir));
#endif
         throw NonPerturbativeRunningError(exp(x), max_step_dir, y(max_step_dir));
      }
      if (errmax <= 1.0) {
         break;
      }
      const precise_real_type htemp = SAFETY * h * pow(errmax, PSHRNK);
      h = (h >= 0.0 ? max(htemp, 0.1 * h) : min(htemp, 0.1 * h));
      if (x + h == x) {
#ifdef ENABLE_VERBOSE
         ERROR("At Q = " << exp(x) << " GeV "
               "stepsize underflow in odeStepper in parameter y("
               << max_step_dir << ") = " << y(max_step_dir) << ", dy("
               << max_step_dir << ")/dx = " << dydx(max_step_dir));
#endif
         throw NonPerturbativeRunningError(exp(x), max_step_dir, y(max_step_dir));
      }
   }
   x += h;
   y = ytemp;

   return errmax > ERRCON ? (precise_real_type)(SAFETY * h * pow(errmax,PGROW)) : (precise_real_type)(5.0 * h);
}

/// Organises integration of 1st order system of ODEs
template <typename ArrayType, typename Derivs,
          typename Stepper = decltype(runge_kutta::odeStepper<ArrayType,Derivs>)>
void integrateOdes(ArrayType& ystart, precise_real_type from, precise_real_type to, precise_real_type eps,
                   precise_real_type h1, precise_real_type hmin, Derivs derivs,
                   Stepper rkqs = runge_kutta::odeStepper<ArrayType,Derivs>, int max_steps = 400)
{
   const int nvar = ystart.size();
   const precise_real_type TINY = 1.0e-16;
   precise_real_type x = from;
   precise_real_type h = sign(h1, to - from);
   ArrayType yscal(nvar);
   ArrayType y(ystart);
   ArrayType dydx;
   int max_step_dir;

   for (int nstp = 0; nstp < max_steps; ++nstp) {
      dydx = derivs(x, y);
      yscal = y.abs() + (dydx * h).abs() + TINY;
      if ((x + h - to) * (x + h - from) > 0.0) {
         h = to - x;
      }

      const precise_real_type hnext = rkqs(y, dydx, x, h, eps, yscal, derivs, max_step_dir);

      if ((x - to) * (to - from) >= 0.0) {
         ystart = y;
         return;
      }

      h = hnext;

      if (fabs(hnext) <= hmin) {
         break;
      }
   }

#ifdef ENABLE_VERBOSE
   ERROR("Bailed out of rk.cpp:too many steps in integrateOdes\n"
         "********** Q = " << exp(x) << " *********");
   ERROR("max step in direction of " << max_step_dir);
   for (int i = 0; i < nvar; i++)
      ERROR("y(" << i << ") = " << y(i) << " dydx(" << i <<
            ") = " << dydx(i));
#endif

   throw NonPerturbativeRunningError(exp(x), max_step_dir, y(max_step_dir));
}

} // namespace runge_kutta

} // namespace flexiblesusy

#endif // RK_H
