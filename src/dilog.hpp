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

#include "precise.hpp"

#ifndef DILOG_H
#define DILOG_H

#include <complex>

#define DILOGATTR noexcept

namespace flexiblesusy {

/// real dilogarithm
/*double dilog(double) DILOGATTR;

/// real dilogarithm
long double dilog(long double) DILOGATTR;

/// complex dilogarithm
std::complex<double> dilog(const std::complex<double>&) DILOGATTR;

/// complex dilogarithm
std::complex<long double> dilog(const std::complex<long double>&) DILOGATTR;*/

//S.D. precise types

precise_real_type dilog(const precise_real_type) DILOGATTR;

precise_complex_type dilog(const precise_complex_type&) DILOGATTR;

/// Clausen function Cl_2(x)
/*double clausen_2(double) DILOGATTR;*/

precise_real_type clausen_2(precise_real_type) DILOGATTR;

} // namespace flexiblesusy

#undef DILOGATTR

#endif
