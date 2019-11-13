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

#ifndef THRESHOLD_LOOP_FUNCTIONS_H
#define THRESHOLD_LOOP_FUNCTIONS_H

#define TCF(n) threshold_loop_functions::F ## n
#define TCf(n) threshold_loop_functions::f ## n
#define TCfth(n) threshold_loop_functions::fth ## n
#define TCf0 threshold_loop_functions::f
#define TCg0 threshold_loop_functions::g
#define TCIabc threshold_loop_functions::Iabc
#define TCB0 threshold_loop_functions::B0
#define TCDB0 threshold_loop_functions::DB0
#define TCC0 threshold_loop_functions::C0
#define TCD0 threshold_loop_functions::D0
#define TCD2t threshold_loop_functions::D2t
#define TCD4t threshold_loop_functions::D4t
#define TCW threshold_loop_functions::W
#define TDelta threshold_loop_functions::delta_xyz
#define TPhi threshold_loop_functions::phi_xyz
#define TCD1F(n) threshold_loop_functions::D1F ## n
#define TCD2F(n) threshold_loop_functions::D2F ## n
#define TCD1f(n) threshold_loop_functions::D1f ## n
#define TCD10f(n) threshold_loop_functions::D10f ## n
#define TCD01f(n) threshold_loop_functions::D01f ## n
#define TCD1f0 threshold_loop_functions::D1f
#define TCD1g0 threshold_loop_functions::D1g

#include "cextensions.hpp"

namespace flexiblesusy {
namespace threshold_loop_functions {

// loop functions from arXiv:1407.4081

#define TCFATTR noexcept ATTR(const)

precise_real_type F1(precise_real_type) TCFATTR;
precise_real_type F2(precise_real_type) TCFATTR;
precise_real_type F3(precise_real_type) TCFATTR;
precise_real_type F4(precise_real_type) TCFATTR;
precise_real_type F5(precise_real_type) TCFATTR;
precise_real_type F6(precise_real_type) TCFATTR;
precise_real_type F7(precise_real_type) TCFATTR;
precise_real_type F8(precise_real_type, precise_real_type) TCFATTR;
precise_real_type F9(precise_real_type, precise_real_type) TCFATTR;

precise_real_type f(precise_real_type) TCFATTR;
precise_real_type g(precise_real_type) TCFATTR;

precise_real_type f1(precise_real_type) TCFATTR;
precise_real_type f2(precise_real_type) TCFATTR;
precise_real_type f3(precise_real_type) TCFATTR;
precise_real_type f4(precise_real_type) TCFATTR;
precise_real_type f5(precise_real_type, precise_real_type) TCFATTR;
precise_real_type f6(precise_real_type, precise_real_type) TCFATTR;
precise_real_type f7(precise_real_type, precise_real_type) TCFATTR;
precise_real_type f8(precise_real_type, precise_real_type) TCFATTR;

// 2-loop threshold function fth[y] from MhEFT-1.1
precise_real_type fth1(precise_real_type) TCFATTR;
precise_real_type fth2(precise_real_type) TCFATTR;
precise_real_type fth3(precise_real_type) TCFATTR;

// first derivatives

precise_real_type D1F1(precise_real_type) TCFATTR;
precise_real_type D1F2(precise_real_type) TCFATTR;
precise_real_type D1F3(precise_real_type) TCFATTR;
precise_real_type D1F4(precise_real_type) TCFATTR;
precise_real_type D1F5(precise_real_type) TCFATTR;
precise_real_type D1F6(precise_real_type) TCFATTR;
precise_real_type D1F7(precise_real_type) TCFATTR;
precise_real_type D1f(precise_real_type) TCFATTR;
precise_real_type D1g(precise_real_type) TCFATTR;
precise_real_type D1f1(precise_real_type) TCFATTR;
precise_real_type D1f2(precise_real_type) TCFATTR;
precise_real_type D1f3(precise_real_type) TCFATTR;
precise_real_type D1f4(precise_real_type) TCFATTR;
precise_real_type D10f5(precise_real_type, precise_real_type) TCFATTR;
precise_real_type D01f5(precise_real_type, precise_real_type) TCFATTR;
precise_real_type D10f6(precise_real_type, precise_real_type) TCFATTR;
precise_real_type D01f6(precise_real_type, precise_real_type) TCFATTR;
precise_real_type D10f7(precise_real_type, precise_real_type) TCFATTR;
precise_real_type D01f7(precise_real_type, precise_real_type) TCFATTR;
precise_real_type D10f8(precise_real_type, precise_real_type) TCFATTR;
precise_real_type D01f8(precise_real_type, precise_real_type) TCFATTR;

// second derivatives

precise_real_type D2F1(precise_real_type) TCFATTR;
precise_real_type D2F2(precise_real_type) TCFATTR;
precise_real_type D2F3(precise_real_type) TCFATTR;
precise_real_type D2F4(precise_real_type) TCFATTR;
precise_real_type D2F5(precise_real_type) TCFATTR;
precise_real_type D2F6(precise_real_type) TCFATTR;
precise_real_type D2F7(precise_real_type) TCFATTR;

/// \f$I_{abc}(a,b,c)\f$ (arguments are interpreted as unsquared)
precise_real_type Iabc(precise_real_type, precise_real_type, precise_real_type) TCFATTR;

/// \f$Delta_{xyz}(x,y,z)\f$ (arguments are interpreted as squared masses)
precise_real_type delta_xyz(precise_real_type, precise_real_type, precise_real_type) TCFATTR;

/// \f$phi_{xyz}(x,y,z)\f$ (arguments are interpreted as squared masses)
precise_real_type phi_xyz(precise_real_type, precise_real_type, precise_real_type) TCFATTR;

/// \f$B_0(p=0,m_1,m_2,Q)\f$ (arguments are interpreted as unsquared)
precise_real_type B0(precise_real_type, precise_real_type, precise_real_type) TCFATTR;

/// \f$B_0'(p=0,m_1,m_2)\f$ (arguments are interpreted as unsquared)
precise_real_type DB0(precise_real_type, precise_real_type) TCFATTR;

/// \f$C_0(p=0,m_1,m_2,m_3)\f$ (arguments are interpreted as unsquared)
precise_real_type C0(precise_real_type, precise_real_type, precise_real_type) TCFATTR;

/// \f$D_0(p=0,m_1,m_2,m_3,m_4)\f$ (arguments are interpreted as unsquared)
precise_real_type D0(precise_real_type, precise_real_type, precise_real_type, precise_real_type) TCFATTR;

/// \f$\tilde{D}_2(m_1,m_2,m_3,m_4)\f$ (arguments are interpreted as unsquared)
precise_real_type D2t(precise_real_type, precise_real_type, precise_real_type, precise_real_type) TCFATTR;

/// \f$\tilde{D}_4(m_1,m_2,m_3,m_4,Q)\f$ (arguments are interpreted as unsquared)
precise_real_type D4t(precise_real_type, precise_real_type, precise_real_type, precise_real_type, precise_real_type) TCFATTR;

/// \f$Q(m_1,m_2,Q)\f$ (arguments are interpreted as unsquared)
precise_real_type W(precise_real_type, precise_real_type, precise_real_type) TCFATTR;

} // namespace threshold_loop_functions
} // namespace flexiblesusy

#endif
