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

#include "wrappers.hpp"
#include "dilog.hpp"
#include "numerics2.hpp"

namespace flexiblesusy {

/*double AbsSqr(double z) noexcept
{
   return z * z;
}

precise_real_type AbsSqr(precise_real_type z) noexcept
{
   return z * z;
}

double AbsSqr(const std::complex<double>& z) noexcept
{
   return std::norm(z);
}

precise_real_type AbsSqr(const precise_complex_type& z) noexcept
{
   return norm(z);
}*/
/*
double AbsSqrt(double x) noexcept
{
   return std::sqrt(std::fabs(x));
}*/

precise_real_type AbsSqrt(precise_real_type x) noexcept
{
   return sqrt(fabs(x));
}

/*double ArcTan(double a) noexcept
{
   return std::atan(a);
}
*/
precise_real_type ArcTan(precise_real_type a) noexcept
{
   return atan(a);
}

/*double ArcSin(double a) noexcept
{
   return std::asin(a);
}*/

precise_real_type ArcSin(precise_real_type a) noexcept
{
   return asin(a);
}

/*double ArcCos(double a) noexcept
{
   return std::acos(a);
}*/

precise_real_type ArcCos(precise_real_type a) noexcept
{
   return acos(a);
}

/*double Arg(const std::complex<double>& z) noexcept
{
   return std::arg(z);
}*/

precise_real_type Arg(const precise_complex_type& z) noexcept
{
   return arg(z);
}
/*
double Conj(double a) noexcept
{
   return a;
}

precise_real_type Conj(const precise_real_type a) noexcept
{
   return a;
}

std::complex<double> Conj(const std::complex<double>& a) noexcept
{
   return std::conj(a);
}

precise_complex_type Conj(const precise_complex_type& a) noexcept
{
   return conj(a);
}*/

/*double Tan(double a) noexcept
{
   return std::tan(a);
}*/

precise_real_type Tan(precise_real_type a) noexcept
{
   return tan(a);
}

/*double Cot(double a) noexcept
{
   return 1./Tan(a);
}*/

precise_real_type Cot(precise_real_type a) noexcept
{
   return 1./Tan(a);
}

/*double Cos(double x) noexcept
{
   return std::cos(x);
}
*/
precise_real_type Cos(precise_real_type x) noexcept
{
   return cos(x);
}

/*double Sin(double x) noexcept
{
   return std::sin(x);
}*/

precise_real_type Sin(precise_real_type x) noexcept
{
   return sin(x);
}

/*double Sec(double x) noexcept
{
   return 1./Cos(x);
}
*/
precise_real_type Sec(precise_real_type x) noexcept
{
   return 1./Cos(x);
}

/*double Csc(double x) noexcept
{
   return 1./Sin(x);
}*/

precise_real_type Csc(precise_real_type x) noexcept
{
   return 1./Sin(x);
}

int Delta(int i, int j) noexcept
{
   return i == j;
}

bool IsClose(double a, double b, double eps) noexcept
{
   return std::abs(a - b) < eps;
}

bool IsClose(precise_real_type a, precise_real_type b, precise_real_type eps) noexcept
{
   return abs(a - b) < eps;
}

bool IsCloseRel(double a, double b, double eps) noexcept
{
   if (IsClose(a, b, std::numeric_limits<double>::epsilon()))
      return true;

   if (std::abs(a) < std::numeric_limits<double>::epsilon())
      return IsClose(a, b, eps);

   return std::abs((a - b)/a) < eps;
}

bool IsCloseRel(precise_real_type a, precise_real_type b, precise_real_type eps) noexcept
{
   if (IsClose(a, b, std::numeric_limits<precise_real_type>::epsilon()))
      return true;

   if (abs(a) < std::numeric_limits<precise_real_type>::epsilon())
      return IsClose(a, b, eps);

   return abs((a - b)/a) < eps;
}


bool IsFinite(double x) noexcept
{
   return std::isfinite(x);
}

bool IsFinite(precise_real_type x) noexcept
{
   return isfinite(x);
}

bool IsFinite(const std::complex<double>& x) noexcept
{
   return std::isfinite(x.real()) && std::isfinite(x.imag());
}

bool IsFinite(const precise_complex_type& x) noexcept
{
   return isfinite(x.real()) && isfinite(x.imag());
}

int KroneckerDelta(int i, int j) noexcept
{
   return i == j;
}

/*double Log(double a) noexcept
{
   return std::log(a);
}*/

precise_real_type Log(precise_real_type a) noexcept
{
   return log(a);
}

/*std::complex<double> ComplexLog(double a) noexcept
{
   return fast_log(std::complex<double>(a,0.));
}*/

precise_complex_type ComplexLog(precise_real_type a) noexcept
{
   return fast_log(precise_complex_type(a,0.));
}

/*std::complex<double> ComplexLog(const std::complex<double>& z) noexcept
{
   return fast_log(z);
}*/


/*precise_complex_type fast_log(const precise_complex_type& z) noexcept //Defined elsewhere, numerics.hpp or numerics2.hpp
{
   return precise_complex_type(log(abs(z)),arg(z));
}*/

precise_complex_type ComplexLog(const precise_complex_type& z) noexcept
{
   return fast_log(z);
}

double FiniteLog(double a) noexcept
{
   const double l = std::log(a);
   return std::isfinite(l) ? l : 0.;
}

precise_real_type FiniteLog(precise_real_type a) noexcept
{
   const precise_real_type l = log(a);
   return isfinite(l) ? l : 0.;
}

/*double MaxAbsValue(double x) noexcept
{
   return Abs(x);
}*/

precise_real_type MaxAbsValue(precise_real_type x) noexcept
{
   return Abs(x);
}

/*double MaxAbsValue(const std::complex<double>& x) noexcept
{
   return Abs(x);
}*/

precise_real_type MaxAbsValue(const precise_complex_type& x) noexcept
{
   return Abs(x);
}

double MaxRelDiff(double a, double b)
{
   const double sTin = fabs(a);
   const double sTout = fabs(b);
   const double maxx = std::max(sTin, sTout);
   const double underflow = 1.0e-20;

   if (maxx < underflow)
      return 0.0;

   return std::abs((a - b) / maxx);
}

precise_real_type MaxRelDiff(precise_real_type a, precise_real_type b)
{
   const precise_real_type sTin = fabs(a);
   const precise_real_type sTout = fabs(b);
   const precise_real_type maxx = max(sTin, sTout);
   const precise_real_type underflow = 1.0e-20;

   if (maxx < underflow)
      return 0.0;

   return abs((a - b) / maxx);
}

double MaxRelDiff(const std::complex<double>& a, const std::complex<double>& b)
{
   const double sTin = std::abs(a);
   const double sTout = std::abs(b);
   const double maxx = std::max(sTin, sTout);
   const double underflow = 1.0e-20;

   if (maxx < underflow)
      return 0.0;

   return std::abs(a - b) / maxx;
}

precise_real_type MaxRelDiff(const precise_complex_type a, const precise_complex_type& b)
{
   const precise_real_type sTin = abs(a);
   const precise_real_type sTout = abs(b);
   const precise_real_type maxx = max(sTin, sTout);
   const precise_real_type underflow = 1.0e-20;

   if (maxx < underflow)
      return 0.0;

   return abs(a - b) / maxx;
}

/*double PolyLog(int n, double z)
{
   if (n == 2)
      return dilog(z);
   throw SetupError("PolyLog(n!=2) not implemented");
}*/

precise_real_type PolyLog(int n, precise_real_type z)
{
   if (n == 2)
      return dilog(z);
   throw SetupError("PolyLog(n!=2) not implemented");
}

/*std::complex<double> PolyLog(int n, const std::complex<double>& z)
{
   if (n == 2)
      return dilog(z);
   throw SetupError("PolyLog(n!=2) not implemented");
}*/

precise_complex_type PolyLog(int n, const precise_complex_type& z)
{
   if (n == 2)
      return dilog(z);
   throw SetupError("PolyLog(n!=2) not implemented");
}

/*double Re(double x) noexcept
{
   return x;
}
*/
/*precise_real_type Re(precise_real_type x) noexcept
{
   return x;
}

double Re(const std::complex<double>& x) noexcept
{
   return std::real(x);
}

precise_real_type Re(const precise_complex_type& x) noexcept
{
   return real(x);
}*/

int Round(double a) noexcept
{
   return static_cast<int>(a >= 0. ? a + 0.5 : a - 0.5);
}

int Round(precise_real_type a) noexcept
{
   return static_cast<int>(a >= 0. ? (precise_real_type)(a + 0.5) : (precise_real_type)(a - 0.5));
}

/*double Im(double) noexcept
{
   return 0.;
}

precise_real_type Im(precise_real_type) noexcept
{
   return 0.;
}

double Im(const std::complex<double>& x) noexcept
{
   return std::imag(x);
}

precise_real_type Im(const precise_complex_type& x) noexcept
{
   return imag(x);
}*/

int Sign(double x) noexcept
{
   return (x >= 0.0 ? 1 : -1);
}

int Sign(precise_real_type x) noexcept
{
   return (x >= 0.0 ? 1 : -1);
}

int Sign(int x) noexcept
{
   return (x >= 0 ? 1 : -1);
}

/*double SignedAbsSqrt(double a) noexcept
{
   return Sign(a) * AbsSqrt(a);
}*/

precise_real_type SignedAbsSqrt(precise_real_type a) noexcept
{
   return Sign(a) * AbsSqrt(a);
}

precise_real_type Sqrt(precise_real_type a) noexcept
{
   return sqrt(a);
}


double Total(double a) noexcept
{
   return a;
}

precise_real_type Total(precise_real_type a) noexcept
{
   return a;
}

std::complex<double> Total(const std::complex<double>& a) noexcept
{
   return a;
}

precise_complex_type Total(const precise_complex_type& a) noexcept
{
   return a;
}

/// unit vector of length N into direction i
Eigen::VectorXdp UnitVector(int N, int i)
{
   Eigen::VectorXdp v = Eigen::VectorXdp::Zero(N);
   v(i) = 1;

   return v;
}

/// unit matrix projector of size MxN into direction i, j
Eigen::MatrixXdp MatrixProjector(int M, int N, int i, int j)
{
   Eigen::MatrixXdp m = Eigen::MatrixXdp::Zero(M,N);
   m(i,j) = 1;

   return m;
}

double ZeroSqrt(double x) noexcept
{
   return (x > 0.0 ? std::sqrt(x) : 0.0);
}

precise_real_type ZeroSqrt(precise_real_type x) noexcept
{
   return (x > 0.0 ? sqrt(x) : (precise_real_type)0.0);
}

} // namespace flexiblesusy
