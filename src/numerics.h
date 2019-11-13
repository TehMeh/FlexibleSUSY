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

/** \file numerics.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach, Alexander Voigt
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief loop functions
*/
#include "precise.hpp"

#ifndef NUMERICS_H
#define NUMERICS_H

namespace softsusy {

precise_real_type a0(precise_real_type m, precise_real_type q) noexcept;
precise_real_type b0(precise_real_type p, precise_real_type m1, precise_real_type m2, precise_real_type q) noexcept;
precise_real_type b1(precise_real_type p, precise_real_type m1, precise_real_type m2, precise_real_type q) noexcept;
precise_real_type b22(precise_real_type p,  precise_real_type m1, precise_real_type m2, precise_real_type q) noexcept;
precise_real_type c0(precise_real_type m1, precise_real_type m2, precise_real_type m3) noexcept;
precise_real_type d27(precise_real_type m1, precise_real_type m2, precise_real_type m3, precise_real_type m4) noexcept;
precise_real_type d0(precise_real_type m1, precise_real_type m2, precise_real_type m3, precise_real_type m4) noexcept;
precise_real_type ffn(precise_real_type p, precise_real_type m1, precise_real_type m2, precise_real_type q) noexcept;
precise_real_type gfn(precise_real_type p, precise_real_type m1, precise_real_type m2, precise_real_type q) noexcept;
precise_real_type hfn(precise_real_type p, precise_real_type m1, precise_real_type m2, precise_real_type q) noexcept;
precise_real_type b22bar(precise_real_type p, precise_real_type m1, precise_real_type m2, precise_real_type q) noexcept;

precise_real_type d1_b0(precise_real_type p, precise_real_type m1, precise_real_type m2) noexcept;

} // namespace softsusy

#endif
