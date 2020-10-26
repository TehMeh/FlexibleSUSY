// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "dilog.hpp"
#include <cmath>
#include <limits>

namespace flexiblesusy {

namespace {
   template <typename T>
   T sqr(T x) noexcept { return x*x; }
} // namespace

/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param x real argument
 * @note Implementation translated by R.Brun from CERNLIB DILOG function C332
 * @return \f$\mathrm{Li}_2(z)\f$
 */
/*double dilog(double x) noexcept {
   const double PI = M_PI;
   const double HF  = 0.5;
   const double PI2 = PI*PI;
   const double PI3 = PI2/3;
   const double PI6 = PI2/6;
   const double PI12 = PI2/12;
   const double C[20] = {0.42996693560813697, 0.40975987533077105,
     -0.01858843665014592, 0.00145751084062268,-0.00014304184442340,
      0.00001588415541880,-0.00000190784959387, 0.00000024195180854,
     -0.00000003193341274, 0.00000000434545063,-0.00000000060578480,
      0.00000000008612098,-0.00000000001244332, 0.00000000000182256,
     -0.00000000000027007, 0.00000000000004042,-0.00000000000000610,
      0.00000000000000093,-0.00000000000000014, 0.00000000000000002};

   double T,H,Y,S,A,ALFA,B1,B2,B0;

   if (x == 1) {
       H = PI6;
   } else if (x == -1) {
       H = -PI12;
   } else {
       T = -x;
       if (T <= -2) {
           Y = -1/(1+T);
           S = 1;
           B1= log(-T);
           B2= log(1+1/T);
           A = -PI3+HF*(B1*B1-B2*B2);
       } else if (T < -1) {
           Y = -1-T;
           S = -1;
           A = log(-T);
           A = -PI6+A*(A+log(1+1/T));
       } else if (T <= -0.5) {
           Y = -(1+T)/T;
           S = 1;
           A = log(-T);
           A = -PI6+A*(-HF*A+log(1+T));
       } else if (T < 0) {
           Y = -T/(1+T);
           S = -1;
           B1= log(1+T);
           A = HF*B1*B1;
       } else if (T <= 1) {
           Y = T;
           S = 1;
           A = 0;
       } else {
           Y = 1/T;
           S = -1;
           B1= log(T);
           A = PI6+HF*B1*B1;
       }
       H    = Y+Y-1;
       ALFA = H+H;
       B1   = 0;
       B2   = 0;
       for (int i=19;i>=0;i--){
          B0 = C[i] + ALFA*B1-B2;
          B2 = B1;
          B1 = B0;
       }
       H = -(S*(B0-H*B2)+A);
    }
    return H;
}*/

/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param x real argument
 * @note Implementation translated by R.Brun from CERNLIB DILOG function C332
 * @return \f$\mathrm{Li}_2(z)\f$
 *
 * @note not full long double precision yet!
 */
/*long double dilog(long double x) noexcept {
   return dilog(static_cast<double>(x));
}*/

//S.D. precise types
//real
precise_real_type dilog(precise_real_type x) noexcept {
   const precise_real_type PI = M_PI;
   const precise_real_type HF  = 0.5_p;
   const precise_real_type PI2 = PI*PI;
   const precise_real_type PI3 = PI2/3;
   const precise_real_type PI6 = PI2/6;
   const precise_real_type PI12 = PI2/12;
   const precise_real_type C[20] = {0.42996693560813697_p, 0.40975987533077105_p,
     -0.01858843665014592_p, 0.00145751084062268_p,-0.00014304184442340_p,
      0.00001588415541880_p,-0.00000190784959387_p, 0.00000024195180854_p,
     -0.00000003193341274_p, 0.00000000434545063_p,-0.00000000060578480_p,
      0.00000000008612098_p,-0.00000000001244332_p, 0.00000000000182256_p,
     -0.00000000000027007_p, 0.00000000000004042_p,-0.00000000000000610_p,
      0.00000000000000093_p,-0.00000000000000014_p, 0.00000000000000002_p};

   precise_real_type T,H,Y,S,A,ALFA,B1,B2,B0;

   if (x == 1) {
       H = PI6;
   } else if (x == -1) {
       H = -PI12;
   } else {
       T = -x;
       if (T <= -2) {
           Y = -1/(1+T);
           S = 1;
           B1= log(-T);
           B2= log(1+1/T);
           A = -PI3+HF*(B1*B1-B2*B2);
       } else if (T < -1) {
           Y = -1-T;
           S = -1;
           A = log(-T);
           A = -PI6+A*(A+log(1+1/T));
       } else if (T <= -0.5_p) {
           Y = -(1+T)/T;
           S = 1;
           A = log(-T);
           A = -PI6+A*(-HF*A+log(1+T));
       } else if (T < 0) {
           Y = -T/(1+T);
           S = -1;
           B1= log(1+T);
           A = HF*B1*B1;
       } else if (T <= 1) {
           Y = T;
           S = 1;
           A = 0;
       } else {
           Y = 1/T;
           S = -1;
           B1= log(T);
           A = PI6+HF*B1*B1;
       }
       H    = Y+Y-1;
       ALFA = H+H;
       B1   = 0;
       B2   = 0;
       for (int i=19;i>=0;i--){
          B0 = C[i] + ALFA*B1-B2;
          B2 = B1;
          B1 = B0;
       }
       H = -(S*(B0-H*B2)+A);
    }
    return H;
}
// complex

precise_complex_type dilog(const precise_complex_type& z) noexcept
{
   const precise_real_type PI = 3.141592653589793238_p;
   static const int N = 12;

   // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
   // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 11}]
   const precise_real_type bf[N] = {
      - 1._p/4._p                  , 1._p/36._p                 ,
      - 1._p/36.e2_p               , 1._p/21168.e1_p            ,
      - 1._p/108864.e2_p           , 1._p/52690176.e1_p         ,
      - 4.064761645144225527e-11_p, 8.921691020456452555e-13_p,
      - 1.993929586072107569e-14_p, 4.518980029619918192e-16_p,
      - 1.035651761218124701e-17_p, 2.395218621026186746e-19_p
   };

   const precise_real_type rz = real(z);
   const precise_real_type iz = imag(z);
   const precise_real_type az = abs(z);

   // special cases
   if (iz == 0._p) {
      if (rz <= 1._p)
         return precise_complex_type(dilog(rz), 0._p);
      if (rz > 1._p)
         return precise_complex_type(dilog(rz), -PI*log(rz));
   } else if (az < std::numeric_limits<precise_real_type>::epsilon()) {
      return z;
   }

   precise_complex_type cy, cz;
   int jsgn, ipi12;

   // transformation to |z|<1, Re(z)<=0.5
   if (rz <= 0.5_p) {
      if (az > 1._p) {
         cy = -0.5_p * sqr(precise_complex_type(log(-z)));
         cz = -log(1._p - 1._p / z);
         jsgn = -1;
         ipi12 = -2;
      } else { // (az <= 1.)
         cy = 0;
         cz = -log(1._p - z);
         jsgn = 1;
         ipi12 = 0;
      }
   } else { // rz > 0.5
      if (az <= sqrt(2*rz)) {
         cz = -log(z);
         cy = cz * log(1._p - z);
         jsgn = -1;
         ipi12 = 2;
      } else { // (az > sqrt(2*rz))
         cy = -0.5_p * sqr(precise_complex_type(log(-z)));
         cz = -log(1._p - 1._p / z);
         jsgn = -1;
         ipi12 = -2;
      }
   }

   // the dilogarithm
   const precise_complex_type cz2(sqr(cz));
   const precise_complex_type sum =
      cz +
      cz2 * (bf[0] +
      cz  * (bf[1] +
      cz2 * (bf[2] +
      cz2 * (bf[3] +
      cz2 * (bf[4] +
      cz2 * (bf[5] +
      cz2 * (bf[6] +
      cz2 * (bf[7] +
      cz2 * (bf[8] +
      cz2 * (bf[9] +
      cz2 * (bf[10] +
      cz2 * (bf[11]))))))))))));

   return static_cast<precise_real_type>(jsgn) * sum + cy + ipi12 * PI * PI / 12._p;
} // S.D.


/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param z complex argument
 * @note Implementation translated from SPheno to C++
 * @return \f$\mathrm{Li}_2(z)\f$
 */
/*std::complex<long double> dilog(const std::complex<long double>& z) noexcept
{
   const long double PI = 3.141592653589793238l;
   static const int N = 12;

   // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
   // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 11}]
   const long double bf[N] = {
      - 1.l/4.l                  , 1.l/36.l                 ,
      - 1.l/36.e2l               , 1.l/21168.e1l            ,
      - 1.l/108864.e2l           , 1.l/52690176.e1l         ,
      - 4.064761645144225527e-11l, 8.921691020456452555e-13l,
      - 1.993929586072107569e-14l, 4.518980029619918192e-16l,
      - 1.035651761218124701e-17l, 2.395218621026186746e-19l
   };

   const long double rz = std::real(z);
   const long double iz = std::imag(z);
   const long double az = std::abs(z);

   // special cases
   if (iz == 0.l) {
      if (rz <= 1.l)
         return std::complex<long double>(dilog(rz), 0.l);
      if (rz > 1.l)
         return std::complex<long double>(dilog(rz), -PI*std::log(rz));
   } else if (az < std::numeric_limits<long double>::epsilon()) {
      return z;
   }

   std::complex<long double> cy, cz;
   int jsgn, ipi12;

   // transformation to |z|<1, Re(z)<=0.5
   if (rz <= 0.5l) {
      if (az > 1.l) {
         cy = -0.5l * sqr(std::log(-z));
         cz = -std::log(1.l - 1.l / z);
         jsgn = -1;
         ipi12 = -2;
      } else { // (az <= 1.)
         cy = 0;
         cz = -std::log(1.l - z);
         jsgn = 1;
         ipi12 = 0;
      }
   } else { // rz > 0.5
      if (az <= std::sqrt(2*rz)) {
         cz = -std::log(z);
         cy = cz * std::log(1.l - z);
         jsgn = -1;
         ipi12 = 2;
      } else { // (az > sqrt(2*rz))
         cy = -0.5l * sqr(std::log(-z));
         cz = -std::log(1.l - 1.l / z);
         jsgn = -1;
         ipi12 = -2;
      }
   }

   // the dilogarithm
   const std::complex<long double> cz2(sqr(cz));
   const std::complex<long double> sum =
      cz +
      cz2 * (bf[0] +
      cz  * (bf[1] +
      cz2 * (bf[2] +
      cz2 * (bf[3] +
      cz2 * (bf[4] +
      cz2 * (bf[5] +
      cz2 * (bf[6] +
      cz2 * (bf[7] +
      cz2 * (bf[8] +
      cz2 * (bf[9] +
      cz2 * (bf[10] +
      cz2 * (bf[11]))))))))))));

   return static_cast<long double>(jsgn) * sum + cy + ipi12 * PI * PI / 12.l;
}
*/
/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param z complex argument
 * @note Implementation translated from SPheno to C++
 * @return \f$\mathrm{Li}_2(z)\f$
 */
/*std::complex<double> dilog(const std::complex<double>& z) noexcept
{
   const double PI = 3.141592653589793;
   static const int N = 10;

   // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
   // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 9}]
   const double bf[N] = {
      - 1./4.,
      + 1./36.,
      - 1./3600.,
      + 1./211680.,
      - 1./10886400.,
      + 1./526901760.,
      - 4.064761645144226e-11,
      + 8.921691020456453e-13,
      - 1.993929586072108e-14,
      + 4.518980029619918e-16
   };

   const double rz = std::real(z);
   const double iz = std::imag(z);
   const double az = std::abs(z);

   // special cases
   if (iz == 0.) {
      if (rz <= 1.)
         return {dilog(rz), 0.};
      else // (rz > 1.)
         return {dilog(rz), -PI*std::log(rz)};
   } else if (az < std::numeric_limits<double>::epsilon()) {
      return z;
   }

   std::complex<double> cy, cz;
   int jsgn, ipi12;

   // transformation to |z|<1, Re(z)<=0.5
   if (rz <= 0.5) {
      if (az > 1.) {
         cy = -0.5 * sqr(std::log(-z));
         cz = -std::log(1. - 1. / z);
         jsgn = -1;
         ipi12 = -2;
      } else { // (az <= 1.)
         cy = 0;
         cz = -std::log(1. - z);
         jsgn = 1;
         ipi12 = 0;
      }
   } else { // rz > 0.5
      if (az <= std::sqrt(2*rz)) {
         cz = -std::log(z);
         cy = cz * std::log(1. - z);
         jsgn = -1;
         ipi12 = 2;
      } else { // (az > sqrt(2*rz))
         cy = -0.5 * sqr(std::log(-z));
         cz = -std::log(1. - 1. / z);
         jsgn = -1;
         ipi12 = -2;
      }
   }

   // the dilogarithm
   const std::complex<double> cz2(sqr(cz));
   const std::complex<double> sum =
      cz +
      cz2 * (bf[0] +
      cz  * (bf[1] +
      cz2 * (bf[2] +
      cz2 * (bf[3] +
      cz2 * (bf[4] +
      cz2 * (bf[5] +
      cz2 * (bf[6] +
      cz2 * (bf[7] +
      cz2 * (bf[8] +
      cz2 * (bf[9]))))))))));

   return double(jsgn) * sum + cy + ipi12 * PI * PI / 12.;
}*/

/**
 * @brief Clausen function \f$\mathrm{Cl}_2(\theta)\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_2(\theta)\f$
 */
/*double clausen_2(double x) noexcept
{
   using std::exp;
   const double PI = 3.141592653589793;
   const std::complex<double> i(0.,1.);

   while (x >= 2*PI)
      x -= 2*PI;

   while (x < 0.)
      x += 2*PI;

   if (std::abs(x) < std::numeric_limits<double>::epsilon() ||
       std::abs(x - PI) < std::numeric_limits<double>::epsilon() ||
       std::abs(x - 2*PI) < std::numeric_limits<double>::epsilon())
      return 0.;

   return std::imag(dilog(exp(i*x)));
}*/

precise_real_type clausen_2(precise_real_type x) noexcept
{
   //using std::exp;
   const precise_real_type PI = 3.141592653589793;
   const precise_complex_type i(0.,1.);

   while (x >= 2*PI)
      x -= 2*PI;

   while (x < 0.)
      x += 2*PI;

   if (abs(x) < std::numeric_limits<precise_real_type>::epsilon() ||
       abs(x - PI) < std::numeric_limits<precise_real_type>::epsilon() ||
       abs(x - 2*PI) < std::numeric_limits<precise_real_type>::epsilon())
      return 0.;

   return imag(dilog(precise_complex_type(exp(i*x))));
}


} // namespace flexiblesusy
