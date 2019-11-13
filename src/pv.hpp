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

#ifndef pv_hpp
#define pv_hpp

#include <complex>
#include "cextensions.hpp"
#include "config.h"

// override macros defined in config.h when building test suite
#ifdef TEST_PV_FFLITE
#  undef  ENABLE_LOOPTOOLS
#  define ENABLE_FFLITE 1
#elif defined(TEST_PV_LOOPTOOLS)
#  undef  ENABLE_FFLITE
#  define ENABLE_LOOPTOOLS 1
#elif defined(TEST_PV_SOFTSUSY)
#  undef  ENABLE_FFLITE
#  undef  ENABLE_LOOPTOOLS
#endif

namespace flexiblesusy {

namespace passarino_veltman {

// functions in LoopTools depend on global variables (such as delta,
// mu, lambda, minmass as well as the internal cache) and therefore
// are not thread-safe
#ifdef ENABLE_LOOPTOOLS
#  define PVATTR noexcept ATTR(pure)
#else
#  define PVATTR noexcept ATTR(const)
#endif

#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)

precise_complex_type A0 (precise_real_type m2, precise_real_type scl2) PVATTR;
precise_complex_type B0 (precise_real_type p2, precise_real_type m2a, precise_real_type m2b, precise_real_type scl2)PVATTR;
precise_complex_type B1 (precise_real_type p2, precise_real_type m2a, precise_real_type m2b, precise_real_type scl2)PVATTR;
precise_complex_type B00(precise_real_type p2, precise_real_type m2a, precise_real_type m2b, precise_real_type scl2)PVATTR;

precise_complex_type A0 (precise_complex_type m2, precise_real_type scl2) PVATTR;
precise_complex_type B0 (precise_complex_type p2, precise_complex_type m2a,
			 precise_complex_type m2b, precise_real_type scl2) PVATTR;
precise_complex_type B1 (precise_complex_type p2, precise_complex_type m2a,
			 precise_complex_type m2b, precise_real_type scl2) PVATTR;
precise_complex_type B00(precise_complex_type p2, precise_complex_type m2a,
			 precise_complex_type m2b, precise_real_type scl2) PVATTR;

template<class T> precise_complex_type B22(T, T, T, precise_real_type) PVATTR;
template<class T> precise_complex_type F0 (T, T, T, precise_real_type) PVATTR;
template<class T> precise_complex_type G0 (T, T, T, precise_real_type) PVATTR;
template<class T> precise_complex_type H0 (T, T, T, precise_real_type) PVATTR;

// CHECK: are the following correct complexifications of B22, H0, F0, G0?

template<class T>
precise_complex_type B22(T p2, T m2a, T m2b, precise_real_type scl2) noexcept
{
    return B00(p2, m2a, m2b, scl2) - A0(m2a, scl2)/4.0 - A0(m2b, scl2)/4.0;
}

template<class T>
precise_complex_type F0(T p2, T m2a, T m2b, precise_real_type scl2) noexcept
{
    return A0(m2a, scl2) - 2.0*A0(m2b, scl2)
	   - (2.0*p2 + 2.0*m2a - m2b) * B0(p2, m2a, m2b, scl2);
}

template<class T>
precise_complex_type G0(T p2, T m2a, T m2b, precise_real_type scl2) noexcept
{
    return (p2 - m2a - m2b) * B0(p2, m2a, m2b, scl2)
	   - A0(m2a, scl2) - A0(m2b, scl2);
}

template<class T>
precise_complex_type H0(T p2, T m2a, T m2b, precise_real_type scl2) noexcept
{
    return 4.0*B00(p2, m2a, m2b, scl2) + G0(p2, m2a, m2b, scl2);
}

/// Derivative of B0(p^2,m1^2,m2^2,Q^2) w.r.t. p^2
precise_complex_type D1B0(precise_real_type p2, precise_real_type m2a, precise_real_type m2b) PVATTR;
precise_complex_type D1B0(precise_complex_type p2, precise_complex_type m2a,
                          precise_complex_type m2b) PVATTR;

/// Derivative of F0(p^2,m1^2,m2^2,Q^2) w.r.t. p^2
template<class T> precise_complex_type D1F0(T, T, T, precise_real_type) PVATTR;
template<class T> precise_complex_type D1F0(T p2, T m2a, T m2b, precise_real_type scl2) noexcept
{
    return - (2.0*p2 + 2.0*m2a - m2b) * D1B0(p2, m2a, m2b)
       - 2.0 * B0(p2, m2a, m2b, scl2);
}

/// Derivative of G0(p^2,m1^2,m2^2,Q^2) w.r.t. p^2
template<class T> precise_complex_type D1G0(T, T, T, precise_real_type) PVATTR;
template<class T> precise_complex_type D1G0(T p2, T m2a, T m2b, precise_real_type scl2) noexcept
{
    return (p2 - m2a - m2b) * D1B0(p2, m2a, m2b)
       + B0(p2, m2a, m2b, scl2);
}

#endif

// the following are mainly for interfacing with loop function
// implementations from softsusy since they come only with precise_real_type
// return type.  If LoopTools or FF is in use, they reduce simply to
// A0(m2, scl2).real(), etc.
precise_real_type ReA0 (precise_real_type m2, precise_real_type scl2) PVATTR;
precise_real_type ReB0 (precise_real_type p2, precise_real_type m2a, precise_real_type m2b, precise_real_type scl2) PVATTR;
precise_real_type ReB1 (precise_real_type p2, precise_real_type m2a, precise_real_type m2b, precise_real_type scl2) PVATTR;
precise_real_type ReB00(precise_real_type p2, precise_real_type m2a, precise_real_type m2b, precise_real_type scl2) PVATTR;
precise_real_type ReB22(precise_real_type p2, precise_real_type m2a, precise_real_type m2b, precise_real_type scl2) PVATTR;
precise_real_type ReH0 (precise_real_type p2, precise_real_type m2a, precise_real_type m2b, precise_real_type scl2) PVATTR;
precise_real_type ReF0 (precise_real_type p2, precise_real_type m2a, precise_real_type m2b, precise_real_type scl2) PVATTR;
precise_real_type ReG0 (precise_real_type p2, precise_real_type m2a, precise_real_type m2b, precise_real_type scl2) PVATTR;

/// Real part of derivative of B0(p^2,m1^2,m2^2,Q^2) w.r.t. p^2
precise_real_type ReD1B0(precise_real_type p2, precise_real_type m2a, precise_real_type m2b) PVATTR;
/// Real part of derivative of F0(p^2,m1^2,m2^2,Q^2) w.r.t. p^2
precise_real_type ReD1F0(precise_real_type p2, precise_real_type m2a, precise_real_type m2b, precise_real_type scl2) PVATTR;
/// Real part of derivative of G0(p^2,m1^2,m2^2,Q^2) w.r.t. p^2
precise_real_type ReD1G0(precise_real_type p2, precise_real_type m2a, precise_real_type m2b, precise_real_type scl2) PVATTR;

} // namespace passarino_veltman

} // namespace flexiblesusy

#endif // pv_hpp
