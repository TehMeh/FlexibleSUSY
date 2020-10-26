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

#include "mssm_twoloophiggs.hpp"
#include "mssm_twoloophiggs.h"
#include "config.h"
#include "dilog.hpp"
#include <cmath>
#include <limits>
#include <utility>

#ifdef ENABLE_THREADS
   #include <mutex>
   #define LOCK_MUTEX() std::lock_guard<std::mutex> lg(mtx_mssm)
   namespace flexiblesusy {namespace mssm_twoloophiggs {
      static std::mutex mtx_mssm; /// locks MSSM fortran functions
   } // namespace mssm_twoloophiggs
} // namespace flexiblesusy
#else
   #define LOCK_MUTEX()
#endif

namespace flexiblesusy {
namespace mssm_twoloophiggs {

namespace {

precise_real_type sqr(precise_real_type a) { return a * a; }
precise_real_type pow3(precise_real_type a) { return a * a * a; }
precise_real_type sqrtabs(precise_real_type a) { return sqrt(abs(a)); }
precise_real_type logabs(precise_real_type x) { return log(abs(x)); }

precise_real_type phi(precise_real_type x, precise_real_type y, precise_real_type z)
{
   //* using std::log;

   const precise_real_type u = x/z, v = y/z;
   const precise_real_type lambda = sqrtabs(sqr(1 - u - v) - 4*u*v);
   const precise_real_type xp = 0.5 * (1 + (u - v) - lambda);
   const precise_real_type xm = 0.5 * (1 - (u - v) - lambda);

   return 1./lambda * (2*logabs(xp)*logabs(xm) - logabs(u)*logabs(v) -
                       2*(dilog(xp) + dilog(xm)) + M_PI*M_PI/3.);
}

/// First derivative of phi[t,T,g] w.r.t. T
precise_real_type dphi_010(precise_real_type t, precise_real_type T, precise_real_type g)
{
   //* usingstd::fabs;
   //* usingstd::sqrt;
   //* usingstd::log;
   //* usingstd::pow;

  precise_real_type Pi2 = M_PI * M_PI;
   const precise_real_type g2 = sqr(g);
   const precise_real_type abbr = (-4*t*T)/g2 + sqr(1 - t/g - T/g);
   const precise_real_type rabbr = sqrtabs(abbr);

   return ((g + t - T)*(Pi2 - 6*dilog(precise_real_type((g - rabbr*g + t - T)/(2.*g))) -
      6*dilog(precise_real_type((g - rabbr*g - t + T)/(2.*g))) -
      3*logabs(t/g)*logabs(T/g) + 6*logabs((g - rabbr*g + t -
      T)/(2.*g))*logabs((g - rabbr*g - t + T)/(2.*g))) + (3*rabbr*g* (
      rabbr*g*((-1 + rabbr)*g + t - T)*logabs(t/g) +
      2*T*(-2*g*logabs(4.) + (g + rabbr*g + t - T)*logabs((g - rabbr*g
      + t - T)/g) + (g + rabbr*g + t - T)*logabs((g + rabbr*g + t -
      T)/g) + g*logabs((g - rabbr*g - t + T)/g) - rabbr*g*logabs((g -
      rabbr*g - t + T)/g) - t*logabs((g - rabbr*g - t + T)/g) +
      T*logabs((g - rabbr*g - t + T)/g) + g*logabs((g + rabbr*g - t +
      T)/g) - rabbr*g*logabs((g + rabbr*g - t + T)/g) - t*logabs((g +
      rabbr*g - t + T)/g) + T*logabs((g + rabbr*g - t + T)/g)) ) ) /
      (T*(g - rabbr*g - t + T)))/(3.*pow(fabs(abbr),1.5)*g2);
}

/// First derivative of phi[g,t,T] w.r.t. T
precise_real_type dphi_001(precise_real_type g, precise_real_type t, precise_real_type T)
{
   //* usingstd::sqrt;

   const precise_real_type Pi2 = 9.869604401089359_p;
   const precise_real_type T2 = sqr(T);
   const precise_real_type T3 = T2*T;
   const precise_real_type x = sqrt(sqr(1 - g/T - t/T) - 4*g*t/T2);
   const precise_real_type y = -(2*(g/T2 + t/T2)*(1 - g/T - t/T) + (8*g*t)/T3)/(2*x);
   const precise_real_type ym = -g/T + t/T;
   const precise_real_type yp = g/T - t/T;
   const precise_real_type lgT = logabs(g/T);
   const precise_real_type ltT = logabs(t/T);
   const precise_real_type lxmym = logabs(0.5_p*(1 - x + ym));
   const precise_real_type lxmyp = logabs(0.5_p*(1 - x + yp));
   const precise_real_type lxpym = logabs(0.5_p*(1 + x + ym));
   const precise_real_type lxpyp = logabs(0.5_p*(1 + x + yp));
   const precise_real_type li2xym = dilog(precise_real_type(0.5_p*(1 - x + ym)));
   const precise_real_type li2xyp = dilog(precise_real_type(0.5_p*(1 - x + yp)));

   return ((t*(t - T) - g*(2*t + T) + sqr(g))*
       (-6*li2xym - 6*li2xyp - 3*lgT*ltT + 6*lxmym*lxmyp + Pi2)
       + 3*T*(lgT*T + ltT*T +
              +2*lxmyp*(g - t + y*T2)/(1 - x + ym)
              +2*lxmym*(-g + t + y*T2)/(1 - x + yp)
              -2*lxpyp*(g - t + y*T2)/(-1 + x - ym)
              -2*lxpym*(-g + t + y*T2)/(-1 + x - yp)
              )*sqr(x))/
      (3.*pow3(T)*pow3(x));
}

precise_real_type calc_At(precise_real_type mt2, precise_real_type mst12, precise_real_type mst22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type mu, precise_real_type tanb)
{
   const precise_real_type s2t = 2*cxt*sxt;
   const precise_real_type Xt = (mst12 - mst22)*s2t/2./sqrtabs(mt2);
   const precise_real_type At = Xt - mu/tanb;

   return At;
}

/// limit st -> 0
Eigen::Matrix<precise_real_type, 2, 1> tadpole_higgs_2loop_at_as_mssm_st_0(
   precise_real_type mt2, precise_real_type mg, precise_real_type mst12, precise_real_type mst22,
   precise_real_type /* sxt */, precise_real_type /* cxt */, precise_real_type scale2,
   precise_real_type mu, precise_real_type tanb, precise_real_type vev2, precise_real_type gs)
{
   //* usingstd::atan;
   //* usingstd::sin;
   //* usingstd::cos;

   const precise_real_type gs2 = sqr(gs);
   const precise_real_type g = sqr(mg);
   const precise_real_type q = scale2;
   const precise_real_type q2 = sqr(q);
   const precise_real_type t = mt2;
   const precise_real_type T1 = mst12;
   const precise_real_type T2 = mst22;
   const precise_real_type g2 = sqr(g);
   const precise_real_type v = sqrtabs(vev2);
   const precise_real_type beta = atan(tanb);
   const precise_real_type v2 = v * sin(beta);
   const precise_real_type v1 = v * cos(beta);
   const precise_real_type ltg = logabs(precise_real_type(t/g));
   const precise_real_type ltq = logabs(precise_real_type(t/g));
   const precise_real_type lgq = logabs(precise_real_type(g/q));
   const precise_real_type lT1q = logabs(precise_real_type(T1/q));
   const precise_real_type lT2q = logabs(precise_real_type(T2/q));
   const precise_real_type lgtq2 = logabs(precise_real_type(g*t/q2));
   const precise_real_type del1 = g2 + sqr(t) + sqr(T1) - 2*(g*t + g*T1 + t*T1);
   const precise_real_type del2 = g2 + sqr(t) + sqr(T2) - 2*(g*t + g*T2 + t*T2);

   const precise_real_type t1 =
      (16*mg*mu*(T1*T2*(g*(-lT1q + lT2q)*ltg + lT1q*ltg*t - lT2q*ltg*t + 5*T1 +
                        (-4 + lgtq2)*lT1q*T1 - lgq*ltq*T1 - 5*T2 + 4*lT2q*T2 -
                        lgtq2*lT2q*T2 + lgq*ltq*T2) + del1*T2*phi(g,t,T1) -
                 del2*T1*phi(g,t,T2)))/(T1*(T1 - T2)*T2*tanb*sqr(v1));

   const precise_real_type t2 =
      (16*(-(T2*(del1*mg*mu - (g2 + g*(t - T1))*(T1 - T2)*tanb)*phi(g,t,T1)) +
           T1*(mg*mu*T2*(g*(lT1q - lT2q)*ltg + lT2q*ltg*t - 5*T1 + lgq*ltq*T1 -
                         lT1q*(ltg*t + (-4 + lgtq2)*T1) +
                         (5 + (-4 + lgtq2)*lT2q - lgq*ltq)*T2) +
               T2*(-T1 + T2)*(g*(2 - 2*lgq*ltq + lT1q*(-2 + lgq + ltq) +
                                 lT2q*(-2 + lgq + ltq)) +
                              2*(-5 + lT1q*(-1 + ltq) + lT2q*(-1 + ltq) - 3*(-2 + ltq)*ltq)*
                              t + (-4 + lT1q)*lT1q*T1 + (-4 + lT2q)*lT2q*T2 + 5*(T1 + T2))*
               tanb + (del2*mg*mu + (g2 + g*(t - T2))*(T1 - T2)*tanb)*phi(g,t,T2)
              )))/(T1*(T1 - T2)*T2*tanb*sqr(v2));

   Eigen::Matrix<precise_real_type, 2, 1> result;
   result << t1, t2;

   const precise_real_type k2 = 0.00004010149318236068; // 1/(4 Pi)^4
   const precise_real_type pref = k2*mt2*gs2;

   return -result*pref;
}

/// limit st -> 0 and mst1 -> mst2
Eigen::Matrix<precise_real_type, 2, 1> tadpole_higgs_2loop_at_as_mssm_st_0_mst1_eq_mst2(
   precise_real_type mt2, precise_real_type mg, precise_real_type mst12, precise_real_type /* mst22 */,
   precise_real_type /* sxt */, precise_real_type /* cxt */, precise_real_type scale2,
   precise_real_type mu, precise_real_type tanb, precise_real_type vev2, precise_real_type gs)
{
   //* usingstd::atan;
   //* usingstd::sin;
   //* usingstd::cos;

   const precise_real_type gs2 = sqr(gs);
   const precise_real_type q = scale2;
   const precise_real_type q2 = sqr(q);
   const precise_real_type g = sqr(mg);
   const precise_real_type g2 = sqr(g);
   const precise_real_type t = mt2;
   const precise_real_type T = mst12;
   const precise_real_type Tsqr = sqr(T);
   const precise_real_type ltg = logabs(t/g);
   const precise_real_type lTq = logabs(T/q);
   const precise_real_type ltq = logabs(t/q);
   const precise_real_type lgq = logabs(g/q);
   const precise_real_type lgtq2 = logabs(g*t/q2);
   const precise_real_type v = sqrtabs(vev2);
   const precise_real_type beta = atan(tanb);
   const precise_real_type v2 = v * sin(beta);
   const precise_real_type v1 = v * cos(beta);

   const precise_real_type t1 =
      (16*mg*mu*(-((g2 - 2*g*t + sqr(t) - Tsqr)*phi(g,t,T)) +
                 T*(ltg*(-g + t) + (1 + lgtq2 - lgq*ltq - 4*lTq + lgtq2*lTq)*T +
                    (g2 - 2*t*T - 2*g*(t + T) + sqr(t) + Tsqr)*
                    dphi_001(g,t,T))))/(tanb*Tsqr*sqr(v1));

   const precise_real_type t2 =
      (-16*(-((g2*(mg*mu + 2*T*tanb) + mg*mu*(sqr(t) - Tsqr) -
               2*g*(mg*mu*t - t*T*tanb + tanb*Tsqr))*phi(g,t,T)) +
            T*(ltg*mg*mu*(-g + t) +
               T*((1 + lgtq2 - lgq*ltq - 4*lTq + lgtq2*lTq)*mg*mu +
                  2*tanb*(g*(1 + (-2 + ltq)*lTq + lgq*(-ltq + lTq)) +
                          t*(-5 - 2*lTq + 2*ltq*(3 + lTq) - 3*sqr(ltq)) +
                          T*(5 - 4*lTq + sqr(lTq)))) +
               mg*mu*(g2 - 2*t*T - 2*g*(t + T) + sqr(t) + Tsqr)*
               dphi_001(g,t,T))))/(tanb*Tsqr*sqr(v2));

   Eigen::Matrix<precise_real_type, 2, 1> result;
   result << t1, t2;

   const precise_real_type k2 = 0.00004010149318236068; // 1/(4 Pi)^4
   const precise_real_type pref = k2*mt2*gs2;

   return -result*pref;
}

/// Pietro Slavich implementation
Eigen::Matrix<precise_real_type, 2, 1> tadpole_higgs_2loop_at_as_mssm_general(
   precise_real_type mt2, precise_real_type mg, precise_real_type mst12, precise_real_type mst22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type scale2,
   precise_real_type mu, precise_real_type tanb, precise_real_type vev2, precise_real_type gs)
{
   Eigen::Matrix<precise_real_type, 2, 1> result;

   ewsb2loop_(&mt2, &mg, &mst12, &mst22, &sxt, &cxt, &scale2,
              &mu, &tanb, &vev2, &gs, &result(0), &result(1));

   // workaround for intel or Eigen bug causing unexpected behaviour
   // of result.allFinite()
   if (!isfinite(result(0)) || !isfinite(result(1)))
       result.setZero();

   return -result;
}

/// limit st -> 0 and mst1 -> mst2
Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_at_as_mssm_with_tadpoles_st_0_mst1_eq_mst2(
   precise_real_type mt2, precise_real_type mg, precise_real_type mst12, precise_real_type /* mst22 */,
   precise_real_type /* sxt */, precise_real_type /* cxt */, precise_real_type scale2, precise_real_type mu,
   precise_real_type tanb, precise_real_type vev2, precise_real_type gs, int /* scheme */)
{
   //* usingstd::atan;
   //* usingstd::sin;

   const precise_real_type gs2 = sqr(gs);
   const precise_real_type g = sqr(mg);
   const precise_real_type g2 = sqr(g);
   const precise_real_type q = scale2;
   const precise_real_type t = mt2;
   const precise_real_type T = mst12;
   const precise_real_type t2 = sqr(t);
   const precise_real_type t3 = t2*t;
   const precise_real_type Tsqr = sqr(T);
   const precise_real_type Tcub = Tsqr*T;
   const precise_real_type ltg = logabs(t/g);
   const precise_real_type lTg = logabs(T/g);
   const precise_real_type lTq = logabs(T/q);
   const precise_real_type ltq = logabs(t/q);
   const precise_real_type lgq = logabs(g/q);
   const precise_real_type lT2t2 = logabs(Tsqr/t2);
   const precise_real_type del = g2 + t2 + Tsqr - 2*(g*t + g*T + t*T);
   const precise_real_type sb = sin(atan(tanb));
   const precise_real_type ht2 = 2./vev2*mt2/sqr(sb);

   Eigen::Matrix<precise_real_type, 2, 2> result;

   result(0,0) = 0.;
   result(0,1) = (32*mg*mu*(-1 + ltq - T*dphi_010(t,T,g)))/T;
   result(1,0) = result(0,1);
   result(1,1) =
      (8*(8*del*mg*mu - 8*del*ltq*mg*mu + 8*del*g*tanb - 8*del*g*lgq*tanb +
          8*del*t*tanb - 8*del*lgq*t*tanb - 8*del*T*tanb +
          5*del*lT2t2*T*tanb + 8*del*ltg*T*tanb + 4*g2*ltg*T*tanb -
          4*del*lTg*T*tanb + 16*g2*lTg*T*tanb + 6*del*ltq*T*tanb +
          2*del*lTq*T*tanb + 40*g*ltg*t*T*tanb + 8*g*ltg*t2*tanb +
          12*ltg*T*t2*tanb - 8*ltg*t3*tanb - 4*ltg*tanb*Tcub +
          8*g*(g + t - T)*T*tanb*phi(t,T,g) + del*T*tanb*sqr(lT2t2) +
          8*del*T*tanb*sqr(ltq) - 8*del*T*tanb*sqr(lTq) +
          8*del*mg*mu*T*dphi_010(t,T,g)))/(del*T*tanb);

   const precise_real_type k2 = 0.00004010149318236068; // 1/(4 Pi)^4
   const precise_real_type pref = k2*ht2*mt2*gs2;

   return -result*pref;
}

/// Pietro Slavich implementation
Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_at_as_mssm_with_tadpoles_general(
   precise_real_type mt2, precise_real_type mg, precise_real_type mst12, precise_real_type mst22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type scale2, precise_real_type mu,
   precise_real_type tanb, precise_real_type vev2, precise_real_type gs, int scheme)
{
   Eigen::Matrix<precise_real_type, 2, 2> result;

   dszhiggs_(&mt2, &mg, &mst12, &mst22, &sxt, &cxt, &scale2, &mu,
             &tanb, &vev2, &gs, &scheme,
             &result(0,0), &result(1,1), &result(0,1));

   result(1,0) = result(0,1);

   return -result;
}

precise_real_type self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles_mst1_eq_mst2(
   precise_real_type mt2, precise_real_type mg, precise_real_type mst12, precise_real_type mst22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type scale2, precise_real_type mu,
   precise_real_type tanb, precise_real_type vev2, precise_real_type gs)
{
   //* usingstd::atan;
   //* usingstd::log;
   //* usingstd::sin;

   precise_real_type Pi2 = M_PI * M_PI;
   const precise_real_type g = sqr(mg);
   const precise_real_type g2 = sqr(g);
   const precise_real_type q = scale2;
   const precise_real_type q2 = sqr(scale2);
   const precise_real_type t = mt2;
   const precise_real_type T = mst12;
   const precise_real_type sb = sin(atan(tanb));
   const precise_real_type ht2 = 2./vev2*mt2/sqr(sb);
   const precise_real_type At = calc_At(mt2, mst12, mst22, sxt, cxt, mu, tanb);

   const precise_real_type result = (-2*(g*(2*At*g + 2*At*t - At*T + mg*T + mg*(g
      - t)*logabs(g/t) - At*T*logabs(g/q)*logabs(t/q) -
      mg*T*logabs(g/q)*logabs(t/q) - 4*mg*T*logabs(T/q) -
      2*At*T*sqr(logabs(T/q)) + logabs((g*t)/q2)*(-(At*(g + t - T)) +
      mg*T + (At + mg)*T*logabs(T/q))) - 2*(At + mg)*(g + t -
      T)*T*phi(t,T,g) + T*(At*(g2 + sqr(t - T) - 2*g*T) + mg*(g2 +
      sqr(t - T) - 2*g*(t + T)))*dphi_010(t,T,g)))/ (g*T);

   const precise_real_type pref = 4*sqr(gs)/sqr(16*Pi2) * ht2*mu*(1./tanb + tanb);

   return -pref * result;
}

precise_real_type self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles_general(
   precise_real_type mt2, precise_real_type mg, precise_real_type mst12, precise_real_type mst22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type scale2, precise_real_type mu,
   precise_real_type tanb, precise_real_type vev2, precise_real_type gs)
{
   precise_real_type result;

   dszodd_(&mt2, &mg, &mst12, &mst22, &sxt, &cxt, &scale2, &mu,
           &tanb, &vev2, &gs, &result);

   return -result;
}

} // anonymous namespace

Eigen::Matrix<precise_real_type, 2, 1> tadpole_higgs_2loop_at_as_mssm(
   precise_real_type mt2, precise_real_type mg, precise_real_type mst12, precise_real_type mst22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type scale2,
   precise_real_type mu, precise_real_type tanb, precise_real_type vev2, precise_real_type gs)
{
   if (abs(sxt) < 1e-8) {
      if (abs((mst12 - mst22)/mst12) < 1e-6)
         return tadpole_higgs_2loop_at_as_mssm_st_0_mst1_eq_mst2(
            mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);

      return tadpole_higgs_2loop_at_as_mssm_st_0(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);
   }

   return tadpole_higgs_2loop_at_as_mssm_general(
      mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);
}

Eigen::Matrix<precise_real_type, 2, 1> tadpole_higgs_2loop_at_at_mssm(
   precise_real_type mt2, precise_real_type mb2, precise_real_type mA2, precise_real_type mst12,
   precise_real_type mst22, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type sxb, precise_real_type cxb,
   precise_real_type scale2, precise_real_type mu, precise_real_type tanb, precise_real_type vev2)
{
   Eigen::Matrix<precise_real_type, 2, 1> result;

   {
      LOCK_MUTEX();

      ddstad_(&mt2, &mb2, &mA2, &mst12, &mst22, &msb12, &msb22,
              &sxt, &cxt, &sxb, &cxb, &scale2, &mu, &tanb, &vev2,
              &result(0), &result(1));
   }

   // workaround for intel or Eigen bug causing unexpected behaviour
   // of result.allFinite()
   if (!isfinite(result(0)) || !isfinite(result(1)))
       result.setZero();

   return -result;
}

Eigen::Matrix<precise_real_type, 2, 1> tadpole_higgs_2loop_ab_as_mssm(
   precise_real_type mb2, precise_real_type mg, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxb, precise_real_type cxb, precise_real_type scale2,
   precise_real_type mu, precise_real_type cotb, precise_real_type vev2, precise_real_type gs)
{
   Eigen::Matrix<precise_real_type, 2, 1> result(tadpole_higgs_2loop_at_as_mssm(
      mb2, mg, msb12, msb22, sxb, cxb, scale2,
      mu, cotb, vev2, gs));

   std::swap(result(0), result(1));

   return result;
}

Eigen::Matrix<precise_real_type, 2, 1> tadpole_higgs_2loop_atau_atau_mssm(
   precise_real_type mtau2, precise_real_type mA2, precise_real_type msv2, precise_real_type mstau12,
   precise_real_type mstau22, precise_real_type sintau, precise_real_type costau, precise_real_type scale2,
   precise_real_type mu, precise_real_type tanb, precise_real_type vev2)
{
   Eigen::Matrix<precise_real_type, 2, 1> result;

   tausqtad_(&mtau2, &mA2, &msv2, &mstau12, &mstau22, &sintau,
             &costau, &scale2, &mu, &tanb, &vev2, &result(0), &result(1));

   // workaround for intel or Eigen bug causing unexpected behaviour
   // of result.allFinite()
   if (!isfinite(result(0)) || !isfinite(result(1)))
       result.setZero();

   return -result;
}

Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_at_as_mssm_with_tadpoles(
   precise_real_type mt2, precise_real_type mg, precise_real_type mst12, precise_real_type mst22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type scale2, precise_real_type mu,
   precise_real_type tanb, precise_real_type vev2, precise_real_type gs, int scheme)
{
   if (abs((mst12 - mst22)/mst12) < 1e-8)
      return self_energy_higgs_2loop_at_as_mssm_with_tadpoles_st_0_mst1_eq_mst2(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs, scheme);

   return self_energy_higgs_2loop_at_as_mssm_with_tadpoles_general(
      mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs, scheme);
}

Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_at_at_mssm_with_tadpoles(
   precise_real_type mt2, precise_real_type mb2, precise_real_type mA2, precise_real_type mst12,
   precise_real_type mst22, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type sxb, precise_real_type cxb,
   precise_real_type scale2, precise_real_type mu, precise_real_type tanb, precise_real_type vev2)
{
   Eigen::Matrix<precise_real_type, 2, 2> result;

   {
      LOCK_MUTEX();

      ddshiggs_(&mt2, &mb2, &mA2, &mst12, &mst22, &msb12, &msb22,
                &sxt, &cxt, &sxb, &cxb, &scale2, &mu, &tanb, &vev2,
                &result(0,0), &result(0,1), &result(1,1));
   }

   result(1,0) = result(0,1);

   return -result;
}

Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_ab_as_mssm_with_tadpoles(
   precise_real_type mb2, precise_real_type mg, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxb, precise_real_type cxb, precise_real_type scale2, precise_real_type mu,
   precise_real_type cotb, precise_real_type vev2, precise_real_type gs, int scheme)
{
   Eigen::Matrix<precise_real_type, 2, 2> result(self_energy_higgs_2loop_at_as_mssm_with_tadpoles(
      mb2, mg, msb12, msb22, sxb, cxb, scale2, mu,
      cotb, vev2, gs, scheme));

   std::swap(result(0,0), result(1,1));

   return result;
}

Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_atau_atau_mssm_with_tadpoles(
   precise_real_type mtau2, precise_real_type mA2, precise_real_type msv2, precise_real_type mstau12,
   precise_real_type mstau22, precise_real_type sintau, precise_real_type costau, precise_real_type scale2,
   precise_real_type mu, precise_real_type tanb, precise_real_type vev2, int scheme)
{
   Eigen::Matrix<precise_real_type, 2, 2> result;

   tausqhiggs_(&mtau2, &mA2, &msv2, &mstau12, &mstau22, &sintau,
               &costau, &scale2, &mu, &tanb, &vev2, &scheme,
               &result(0,0), &result(1,1), &result(0,1));

   result(1,0) = result(0,1);

   return -result;
}

precise_real_type self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
   precise_real_type mt2, precise_real_type mg, precise_real_type mst12, precise_real_type mst22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type scale2, precise_real_type mu,
   precise_real_type tanb, precise_real_type vev2, precise_real_type gs)
{
   if (abs((mst12 - mst22)/mst12) < 1e-8) {
      const precise_real_type At = calc_At(mt2, mst12, mst22, sxt, cxt, mu, tanb);

      // if At = 0 => mu = 0 => dMA(2L) = 0
      if (abs(At) < std::numeric_limits<precise_real_type>::epsilon())
         return 0.;

      return self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles_mst1_eq_mst2(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);
   }

   return self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles_general(
      mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);
}

precise_real_type self_energy_pseudoscalar_2loop_at_at_mssm_with_tadpoles(
   precise_real_type mt2, precise_real_type mb2, precise_real_type mA2, precise_real_type mst12,
   precise_real_type mst22, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type sxb, precise_real_type cxb,
   precise_real_type scale2, precise_real_type mu, precise_real_type tanb, precise_real_type vev2)
{
   precise_real_type result;

   {
      LOCK_MUTEX();

      ddsodd_(&mt2, &mb2, &mA2, &mst12, &mst22, &msb12, &msb22,
              &sxt, &cxt, &sxb, &cxb, &scale2, &mu, &tanb, &vev2, &result);
   }

   return -result;
}

precise_real_type self_energy_pseudoscalar_2loop_ab_as_mssm_with_tadpoles(
   precise_real_type mb2, precise_real_type mg, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxb, precise_real_type cxb, precise_real_type scale2, precise_real_type mu,
   precise_real_type cotb, precise_real_type vev2, precise_real_type gs)
{
   return self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
      mb2, mg, msb12, msb22, sxb, cxb, scale2, mu,
      cotb, vev2, gs);
}

precise_real_type self_energy_pseudoscalar_2loop_atau_atau_mssm_with_tadpoles(
   precise_real_type mtau2, precise_real_type mA2, precise_real_type msv2, precise_real_type mstau12,
   precise_real_type mstau22, precise_real_type sintau, precise_real_type costau, precise_real_type scale2,
   precise_real_type mu, precise_real_type tanb, precise_real_type vev2)
{
   precise_real_type result;

   tausqodd_(&mtau2, &mA2, &msv2, &mstau12, &mstau22, &sintau,
             &costau, &scale2, &mu, &tanb, &vev2, &result);

   return -result;
}

// self-energies without tadpoles

Eigen::Matrix<precise_real_type, 2, 2> rotate_scalar(
   precise_real_type self_energy, precise_real_type tanb)
{
   const precise_real_type tanb2 = sqr(tanb);
   const precise_real_type sinb = tanb / sqrtabs(1. + tanb2);
   const precise_real_type cosb = 1. / sqrtabs(1. + tanb2);

   Eigen::Matrix<precise_real_type, 2, 2> result;

   result(0,0) = self_energy * sqr(sinb);
   result(0,1) = - self_energy * sinb * cosb;
   result(1,0) = result(0,1);
   result(1,1) = self_energy * sqr(cosb);

   return result;
}

Eigen::Matrix<precise_real_type, 2, 2> subtract_mssm_tadpoles_scalar(
   precise_real_type self_energy, const Eigen::Matrix<precise_real_type, 2, 1>& tadpoles,
   precise_real_type tanb)
{
   return rotate_scalar(self_energy, tanb) + Eigen::Matrix<precise_real_type, 2, 2>(tadpoles.asDiagonal());
}

Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_at_as_mssm(
   precise_real_type mt2, precise_real_type mg, precise_real_type mst12, precise_real_type mst22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type scale2, precise_real_type mu,
   precise_real_type tanb, precise_real_type vev2, precise_real_type gs, int scheme)
{
   const Eigen::Matrix<precise_real_type, 2, 2> result =
      self_energy_higgs_2loop_at_as_mssm_with_tadpoles(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu,
         tanb, vev2, gs, scheme);

   const precise_real_type dMA = self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
      mt2, mg, mst12, mst22, sxt, cxt, scale2, mu,
      tanb, vev2, gs);

   const Eigen::Matrix<precise_real_type, 2, 1> tadpoles =
      tadpole_higgs_2loop_at_as_mssm(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);

   const Eigen::Matrix<precise_real_type, 2, 2> tM =
      subtract_mssm_tadpoles_scalar(dMA, tadpoles, tanb);

   return result + tM;
}

Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_at_at_mssm(
   precise_real_type mt2, precise_real_type mb2, precise_real_type mA2, precise_real_type mst12,
   precise_real_type mst22, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type sxb, precise_real_type cxb,
   precise_real_type scale2, precise_real_type mu, precise_real_type tanb, precise_real_type vev2)
{
   const Eigen::Matrix<precise_real_type, 2, 2> result =
      self_energy_higgs_2loop_at_at_mssm_with_tadpoles(
         mt2, mb2, mA2, mst12, mst22, msb12, msb22,
         sxt, cxt, sxb, cxb, scale2, mu, tanb, vev2);

   const precise_real_type dMA = self_energy_pseudoscalar_2loop_at_at_mssm_with_tadpoles(
      mt2, mb2, mA2, mst12, mst22, msb12, msb22,
      sxt, cxt, sxb, cxb, scale2, mu, tanb, vev2);

   const Eigen::Matrix<precise_real_type, 2, 1> tadpoles =
      tadpole_higgs_2loop_at_at_mssm(
         mt2, mb2, mA2, mst12, mst22, msb12, msb22,
         sxt, cxt, sxb, cxb, scale2, mu, tanb, vev2);

   const Eigen::Matrix<precise_real_type, 2, 2> tM =
      subtract_mssm_tadpoles_scalar(dMA, tadpoles, tanb);

   return result + tM;
}

Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_ab_as_mssm(
   precise_real_type mb2, precise_real_type mg, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxb, precise_real_type cxb, precise_real_type scale2, precise_real_type mu,
   precise_real_type cotb, precise_real_type vev2, precise_real_type gs, int scheme)
{
   const Eigen::Matrix<precise_real_type, 2, 2> result =
      self_energy_higgs_2loop_ab_as_mssm_with_tadpoles(
         mb2, mg, msb12, msb22, sxb, cxb, scale2, mu,
         cotb, vev2, gs, scheme);

   const precise_real_type dMA = self_energy_pseudoscalar_2loop_ab_as_mssm_with_tadpoles(
      mb2, mg, msb12, msb22, sxb, cxb, scale2, mu, cotb, vev2, gs);

   const Eigen::Matrix<precise_real_type, 2, 1> tadpoles =
      tadpole_higgs_2loop_ab_as_mssm(
         mb2, mg, msb12, msb22, sxb, cxb, scale2, mu, cotb, vev2, gs);

   const Eigen::Matrix<precise_real_type, 2, 2> tM =
      subtract_mssm_tadpoles_scalar(dMA, tadpoles, 1./cotb);

   return result + tM;
}

Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_atau_atau_mssm(
   precise_real_type mtau2, precise_real_type mA2, precise_real_type msv2, precise_real_type mstau12,
   precise_real_type mstau22, precise_real_type sintau, precise_real_type costau, precise_real_type scale2,
   precise_real_type mu, precise_real_type tanb, precise_real_type vev2, int scheme)
{
   const Eigen::Matrix<precise_real_type, 2, 2> result =
      self_energy_higgs_2loop_atau_atau_mssm_with_tadpoles(
         mtau2, mA2, msv2, mstau12, mstau22, sintau, costau, scale2,
         mu, tanb, vev2, scheme);

   const precise_real_type dMA = self_energy_pseudoscalar_2loop_atau_atau_mssm_with_tadpoles(
      mtau2, mA2, msv2, mstau12, mstau22, sintau, costau, scale2,
      mu, tanb, vev2);

   const Eigen::Matrix<precise_real_type, 2, 1> tadpoles =
      tadpole_higgs_2loop_atau_atau_mssm(
         mtau2, mA2, msv2, mstau12, mstau22, sintau, costau, scale2,
         mu, tanb, vev2);

   const Eigen::Matrix<precise_real_type, 2, 2> tM =
      subtract_mssm_tadpoles_scalar(dMA, tadpoles, tanb);

   return result + tM;
}

Eigen::Matrix<precise_real_type, 2, 2> rotate_pseudoscalar(
   precise_real_type self_energy, precise_real_type tanb)
{
   const precise_real_type tanb2 = sqr(tanb);
   const precise_real_type sinb = tanb / sqrtabs(1. + tanb2);
   const precise_real_type cosb = 1. / sqrtabs(1. + tanb2);

   Eigen::Matrix<precise_real_type, 2, 2> result;

   // see hep-ph/0105096 Eq. (9)
   result(0,0) = self_energy * sqr(sinb);
   result(0,1) = self_energy * sinb * cosb;
   result(1,0) = result(0,1);
   result(1,1) = self_energy * sqr(cosb);

   return result;
}

Eigen::Matrix<precise_real_type, 2, 2> subtract_mssm_tadpoles_pseudoscalar(
   precise_real_type self_energy, const Eigen::Matrix<precise_real_type, 2, 1>& tadpoles,
   precise_real_type tanb)
{
   return rotate_pseudoscalar(self_energy, tanb) + Eigen::Matrix<precise_real_type, 2, 2>(tadpoles.asDiagonal());
}

Eigen::Matrix<precise_real_type, 2, 2> self_energy_pseudoscalar_2loop_at_as_mssm(
   precise_real_type mt2, precise_real_type mg, precise_real_type mst12, precise_real_type mst22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type scale2, precise_real_type mu,
   precise_real_type tanb, precise_real_type vev2, precise_real_type gs)
{
   const precise_real_type se = self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
      mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);

   const Eigen::Matrix<precise_real_type, 2, 1> tadpoles =
      tadpole_higgs_2loop_at_as_mssm(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);

   return subtract_mssm_tadpoles_pseudoscalar(se, tadpoles, tanb);
}

Eigen::Matrix<precise_real_type, 2, 2> self_energy_pseudoscalar_2loop_at_at_mssm(
   precise_real_type mt2, precise_real_type mb2, precise_real_type mA2, precise_real_type mst12,
   precise_real_type mst22, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type sxb, precise_real_type cxb,
   precise_real_type scale2, precise_real_type mu, precise_real_type tanb, precise_real_type vev2)
{
   const precise_real_type se = self_energy_pseudoscalar_2loop_at_at_mssm_with_tadpoles(
      mt2, mb2, mA2, mst12, mst22, msb12, msb22,
      sxt, cxt, sxb, cxb, scale2, mu, tanb, vev2);

   const Eigen::Matrix<precise_real_type, 2, 1> tadpoles =
      tadpole_higgs_2loop_at_at_mssm(
         mt2, mb2, mA2, mst12, mst22, msb12, msb22,
         sxt, cxt, sxb, cxb, scale2, mu, tanb, vev2);

   return subtract_mssm_tadpoles_pseudoscalar(se, tadpoles, tanb);
}

Eigen::Matrix<precise_real_type, 2, 2> self_energy_pseudoscalar_2loop_ab_as_mssm(
   precise_real_type mb2, precise_real_type mg, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxb, precise_real_type cxb, precise_real_type scale2, precise_real_type mu,
   precise_real_type cotb, precise_real_type vev2, precise_real_type gs)
{
   const precise_real_type se = self_energy_pseudoscalar_2loop_ab_as_mssm_with_tadpoles(
      mb2, mg, msb12, msb22, sxb, cxb, scale2, mu, cotb, vev2, gs);

   const Eigen::Matrix<precise_real_type, 2, 1> tadpoles =
      tadpole_higgs_2loop_ab_as_mssm(
         mb2, mg, msb12, msb22, sxb, cxb, scale2, mu, cotb, vev2, gs);

   return subtract_mssm_tadpoles_pseudoscalar(se, tadpoles, 1./cotb);
}

Eigen::Matrix<precise_real_type, 2, 2> self_energy_pseudoscalar_2loop_atau_atau_mssm(
   precise_real_type mtau2, precise_real_type mA2, precise_real_type msv2, precise_real_type mstau12,
   precise_real_type mstau22, precise_real_type sintau, precise_real_type costau, precise_real_type scale2,
   precise_real_type mu, precise_real_type tanb, precise_real_type vev2)
{
   const precise_real_type se = self_energy_pseudoscalar_2loop_atau_atau_mssm_with_tadpoles(
      mtau2, mA2, msv2, mstau12, mstau22, sintau, costau, scale2,
      mu, tanb, vev2);

   const Eigen::Matrix<precise_real_type, 2, 1> tadpoles =
      tadpole_higgs_2loop_atau_atau_mssm(
         mtau2, mA2, msv2, mstau12, mstau22, sintau, costau, scale2,
         mu, tanb, vev2);

   return subtract_mssm_tadpoles_pseudoscalar(se, tadpoles, tanb);
}

} // namespace mssm_twoloophiggs
} // namespace flexiblesusy
