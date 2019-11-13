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

#ifndef SFERMIONS_H
#define SFERMIONS_H

#include <boost/array.hpp>
#include <Eigen/Core>

namespace flexiblesusy {
namespace sfermions {

enum Sparticles {
   up = 0,
   down = 1,
   neutrino = 2,
   electron = 3,
   NUMBER_OF_MSSM_SPARTICLES
};

extern const boost::array<precise_real_type, NUMBER_OF_MSSM_SPARTICLES> Isospin;
extern const boost::array<precise_real_type, NUMBER_OF_MSSM_SPARTICLES> Hypercharge_left;
extern const boost::array<precise_real_type, NUMBER_OF_MSSM_SPARTICLES> Hypercharge_right;

/**
 * parameters needed to fill 2 x 2 sfermion mass matrix
 */
struct Mass_data {
   precise_real_type ml2{};      ///< soft mass of left-handed sfermion
   precise_real_type mr2{};      ///< soft mass of right-handed sfermion
   precise_real_type yf{};       ///< Yukawa coupling
   precise_real_type vd{}, vu{}; ///< Higgs VEVs
   precise_real_type gY{}, g2{}; ///< gauge couplings (not GUT normalized)
   precise_real_type Tyf{};      ///< trilinear coupling
   precise_real_type mu{};       ///< Superpotential parameter
   precise_real_type T3{};       ///< weak isospin
   precise_real_type Yl{};       ///< Hypercharge of left-handed sfermion
   precise_real_type Yr{};       ///< Hypercharge of right-handed sfermion
};

precise_real_type diagonalize_sfermions_2x2(const Mass_data&,
                                 Eigen::Array<precise_real_type,2,1>&);

} // namespace sfermions
} // namespace flexiblesusy

#endif
