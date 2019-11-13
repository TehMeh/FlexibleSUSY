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

#ifndef EFFECTIVE_COUPLINGS_H
#define EFFECTIVE_COUPLINGS_H

#include <complex>
#include <map>

namespace flexiblesusy {

namespace effective_couplings {

// loop functions of hep-ph/9504378, Eq. (53)
precise_complex_type scaling_function(precise_real_type tau);
precise_complex_type AS0(precise_real_type tau);
precise_complex_type AS12(precise_real_type tau);
precise_complex_type AS1(precise_real_type tau);
precise_complex_type AP12(precise_real_type tau);

// QCD corrections to the diphoton decay width
// for scalars and pseudoscalars, at the scale
// m_decay / 2
precise_complex_type scalar_diphoton_fermion_loop(
   precise_real_type m_decay, precise_real_type m_loop);
precise_complex_type pseudoscalar_diphoton_fermion_loop(
   precise_real_type m_decay, precise_real_type m_loop);

precise_complex_type linear_interpolation(
   precise_real_type x, const std::map<precise_real_type,precise_complex_type >& data);
precise_complex_type quadratic_interpolation(
   precise_real_type x, const std::map<precise_real_type,precise_complex_type >& data);

std::map<precise_real_type,precise_complex_type > get_scalar_fermion_loop_data();
std::map<precise_real_type,precise_complex_type > get_pseudoscalar_fermion_loop_data();

} // namespace effective_couplings

} // namespace flexiblesusy

#endif
