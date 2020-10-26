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

#ifndef MSSM_TWOLOOPHIGGS_H
#define MSSM_TWOLOOPHIGGS_H

#include "precise.hpp"
#include <Eigen/Core>

/**
 * @file mssm_twoloophiggs.hpp
 * @brief function declarations for 2-loop MSSM Higgs self-energies
 *        and tadpoles
 *
 * Notation:
 *
 * mt2    : squared DR-bar top mass in the MSSM
 * mb2    : squared DR-bar bottom mass in the MSSM
 * mtau2  : squared DR-bar tau mass in the MSSM
 * mg     : DR-bar gluino mass in the MSSM
 * mA2    : squared DR-bar CP-odd Higgs mass in the MSSM
 * mst12  : squared DR-bar lightest stop mass
 * mst22  : squared DR-bar heaviest stop mass
 * msb12  : squared DR-bar lightest sbottom mass
 * msb22  : squared DR-bar heaviest sbottom mass
 * mstau12: squared DR-bar lightest stau mass
 * mstau22: squared DR-bar heaviest stau mass
 * msv2   : squared DR-bar tau sneutrino mass
 *
 * sxt    : sine of DR-bar stop mixing angle in the MSSM
 * cxt    : cosine of DR-bar stop mixing angle in the MSSM
 * sxb    : sine of DR-bar sbottom mixing angle in the MSSM
 * cxb    : cosine of DR-bar sbottom mixing angle in the MSSM
 * sintau : sine of DR-bar stau mixing angle in the MSSM
 * costau : cosine of DR-bar stau mixing angle in the MSSM
 *
 * gs     : DR-bar strong gauge coupling g3 in the MSSM
 * mu     : DR-bar mu-parameter in the MSSM (arXiv:0907.4682)
 * tanb   : DR-bar tan(beta) = vu/vd in the MSSM
 * cotb   : DR-bar 1/tan(beta) in the MSSM
 * vev2   : squared DR-bar vev^2 = (vu^2 + vd^2) in the MSSM
 *
 * scheme : DR-bar scheme (0) or on-shell scheme (1)
 */

namespace flexiblesusy {
namespace mssm_twoloophiggs {

// tadpoles

Eigen::Matrix<precise_real_type, 2, 1> tadpole_higgs_2loop_at_as_mssm(
   precise_real_type mt2, precise_real_type mg, precise_real_type mst12, precise_real_type mst22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type scale2,
   precise_real_type mu, precise_real_type tanb, precise_real_type vev2, precise_real_type gs);

Eigen::Matrix<precise_real_type, 2, 1> tadpole_higgs_2loop_at_at_mssm(
   precise_real_type mt2, precise_real_type mb2, precise_real_type mA2, precise_real_type mst12,
   precise_real_type mst22, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type sxb, precise_real_type cxb,
   precise_real_type scale2, precise_real_type mu, precise_real_type tanb, precise_real_type vev2);

Eigen::Matrix<precise_real_type, 2, 1> tadpole_higgs_2loop_ab_as_mssm(
   precise_real_type mb2, precise_real_type mg, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxb, precise_real_type cxb, precise_real_type scale2,
   precise_real_type mu, precise_real_type cotb, precise_real_type vev2, precise_real_type gs);

Eigen::Matrix<precise_real_type, 2, 1> tadpole_higgs_2loop_atau_atau_mssm(
   precise_real_type mtau2, precise_real_type mA2, precise_real_type msv2, precise_real_type mstau12,
   precise_real_type mstau22, precise_real_type sintau, precise_real_type costau, precise_real_type scale2,
   precise_real_type mu, precise_real_type tanb, precise_real_type vev2);

// self-energies

Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_at_as_mssm(
   precise_real_type mt2, precise_real_type mg, precise_real_type mst12, precise_real_type mst22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type scale2, precise_real_type mu,
   precise_real_type tanb, precise_real_type vev2, precise_real_type gs, int scheme = 0);

Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_at_at_mssm(
   precise_real_type mt2, precise_real_type mb2, precise_real_type mA2, precise_real_type mst12,
   precise_real_type mst22, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type sxb, precise_real_type cxb,
   precise_real_type scale2, precise_real_type mu, precise_real_type tanb, precise_real_type vev2);

Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_ab_as_mssm(
   precise_real_type mb2, precise_real_type mg, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxb, precise_real_type cxb, precise_real_type scale2, precise_real_type mu,
   precise_real_type cotb, precise_real_type vev2, precise_real_type gs, int scheme = 0);

Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_atau_atau_mssm(
   precise_real_type mtau2, precise_real_type mA2, precise_real_type msv2, precise_real_type mstau12,
   precise_real_type mstau22, precise_real_type sintau, precise_real_type costau, precise_real_type scale2,
   precise_real_type mu, precise_real_type tanb, precise_real_type vev2, int scheme = 0);


Eigen::Matrix<precise_real_type, 2, 2> self_energy_pseudoscalar_2loop_at_as_mssm(
   precise_real_type mt2, precise_real_type mg, precise_real_type mst12, precise_real_type mst22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type scale2, precise_real_type mu,
   precise_real_type tanb, precise_real_type vev2, precise_real_type gs);

Eigen::Matrix<precise_real_type, 2, 2> self_energy_pseudoscalar_2loop_at_at_mssm(
   precise_real_type mt2, precise_real_type mb2, precise_real_type mA2, precise_real_type mst12,
   precise_real_type mst22, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type sxb, precise_real_type cxb,
   precise_real_type scale2, precise_real_type mu, precise_real_type tanb, precise_real_type vev2);

Eigen::Matrix<precise_real_type, 2, 2> self_energy_pseudoscalar_2loop_ab_as_mssm(
   precise_real_type mb2, precise_real_type mg, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxb, precise_real_type cxb, precise_real_type scale2, precise_real_type mu,
   precise_real_type cotb, precise_real_type vev2, precise_real_type gs);

Eigen::Matrix<precise_real_type, 2, 2> self_energy_pseudoscalar_2loop_atau_atau_mssm(
   precise_real_type mtau2, precise_real_type mA2, precise_real_type msv2, precise_real_type mstau12,
   precise_real_type mstau22, precise_real_type sintau, precise_real_type costau, precise_real_type scale2,
   precise_real_type mu, precise_real_type tanb, precise_real_type vev2);

// self-energies with tadpoles added

Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_at_as_mssm_with_tadpoles(
   precise_real_type mt2, precise_real_type mg, precise_real_type mst12, precise_real_type mst22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type scale2, precise_real_type mu,
   precise_real_type tanb, precise_real_type vev2, precise_real_type gs, int scheme = 0);

Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_at_at_mssm_with_tadpoles(
   precise_real_type mt2, precise_real_type mb2, precise_real_type mA2, precise_real_type mst12,
   precise_real_type mst22, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type sxb, precise_real_type cxb,
   precise_real_type scale2, precise_real_type mu, precise_real_type tanb, precise_real_type vev2);

Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_ab_as_mssm_with_tadpoles(
   precise_real_type mb2, precise_real_type mg, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxb, precise_real_type cxb, precise_real_type scale2, precise_real_type mu,
   precise_real_type cotb, precise_real_type vev2, precise_real_type gs, int scheme = 0);

Eigen::Matrix<precise_real_type, 2, 2> self_energy_higgs_2loop_atau_atau_mssm_with_tadpoles(
   precise_real_type mtau2, precise_real_type mA2, precise_real_type msv2, precise_real_type mstau12,
   precise_real_type mstau22, precise_real_type sintau, precise_real_type costau, precise_real_type scale2,
   precise_real_type mu, precise_real_type tanb, precise_real_type vev2, int scheme = 0);


precise_real_type self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
   precise_real_type mt2, precise_real_type mg, precise_real_type mst12, precise_real_type mst22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type scale2, precise_real_type mu,
   precise_real_type tanb, precise_real_type vev2, precise_real_type gs);

precise_real_type self_energy_pseudoscalar_2loop_at_at_mssm_with_tadpoles(
   precise_real_type mt2, precise_real_type mb2, precise_real_type mA2, precise_real_type mst12,
   precise_real_type mst22, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxt, precise_real_type cxt, precise_real_type sxb, precise_real_type cxb,
   precise_real_type scale2, precise_real_type mu, precise_real_type tanb, precise_real_type vev2);

precise_real_type self_energy_pseudoscalar_2loop_ab_as_mssm_with_tadpoles(
   precise_real_type mb2, precise_real_type mg, precise_real_type msb12, precise_real_type msb22,
   precise_real_type sxb, precise_real_type cxb, precise_real_type scale2, precise_real_type mu,
   precise_real_type cotb, precise_real_type vev2, precise_real_type gs);

precise_real_type self_energy_pseudoscalar_2loop_atau_atau_mssm_with_tadpoles(
   precise_real_type mtau2, precise_real_type mA2, precise_real_type msv2, precise_real_type mstau12,
   precise_real_type mstau22, precise_real_type sintau, precise_real_type costau, precise_real_type scale2,
   precise_real_type mu, precise_real_type tanb, precise_real_type vev2);

} // namespace mssm_twoloophiggs
} // namespace flexiblesusy

#endif
