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

#ifndef SplitMSSM_FullMSSM_THRESHOLDS_H
#define SplitMSSM_FullMSSM_THRESHOLDS_H

#include <Eigen/Core>

/**
 * @file splitmssm_thresholds.hpp
 *
 * Contains function declarations for threshold corrections of the
 * MSSM to the Standard Model or the Split-MSSM from arXiv:1407.4081.
 *
 * Example for matching of the MSSM to the Standard Model:
@code
   Parameters pars;
   // fill parameters ...

   precise_real_type lambda = lambda_tree_level(pars);

   if (loopLevel > 0) {
      lambda +=
         + delta_lambda_1loop_reg(pars)
         + delta_lambda_1loop_phi(pars)
         + delta_lambda_1loop_chi_1(pars)
         + delta_lambda_1loop_chi_2(pars);
   }

   if (loopLevel > 1)
      lambda += delta_lambda_2loop_phi_HSS(pars);
@endcode
 *
 * Example for matching of the MSSM to the Split-MSSM:
@code
   Parameters pars;
   // fill parameters ...
   // Note: g1 and g2 are defined in the Split-MSSM

   precise_real_type lambda_split = lambda_tree_level(pars);
   precise_real_type gYu = gYu_tree_level(pars);
   precise_real_type gYd = gYd_tree_level(pars);
   precise_real_type g2u = g2u_tree_level(pars);
   precise_real_type g2d = g2d_tree_level(pars);

   if (loopLevel > 0) {
      lambda_split +=
         + delta_lambda_1loop_reg(pars)
         + delta_lambda_1loop_phi(pars);
      gYu += delta_gYu_1loop(pars);
      gYd += delta_gYd_1loop(pars);
      g2u += delta_g2u_1loop(pars);
      g2d += delta_g2d_1loop(pars);
   }

   if (loopLevel > 1)
      lambda_split += delta_lambda_2loop_phi(pars);
@endcode
 *
 * Example for matching of the Split-MSSM to the Standard Model:
@code
   Parameters pars;
   // fill parameters ...
   // Note: parameters are defined in the Split-MSSM

   precise_real_type lambda_SM = lambda_split;

   if (loopLevel > 0) {
      lambda_SM += delta_lambda_1loop_chi_1(
         scale, mu, lambda_split, gYu, gYd, g2u, g2d, m1, m2);
   }
@endcode
 */

namespace flexiblesusy {
namespace splitmssm_thresholds {

/**
 * @class Parameters
 * @brief Parameters for MSSM threshold corrections to the SM or Split-MSSM
 *
 * Contains the running MS-bar parameters of the EFT (either the SM or
 * the Split-MSSM) and the runnign DR-bar parameters of the MSSM.  See
 * arXiv:1407.4081 for the parameter definition.
 */
struct Parameters {
   precise_real_type g1{0.}, g2{0.}, g3{0.}; ///< MS-bar gauge couplings in the EFT (GUT normalized)
   precise_real_type gt{0.};          ///< MS-bar top Yukawa coupling of the SM or Split-MSSM
   precise_real_type At{0.};          ///< DR-bar trilinear coupling for the stops in the MSSM
   precise_real_type mu{0.};          ///< bilinear Higgsino coupling
   precise_real_type mA{0.};          ///< mass of the heavy Higgs precise_real_typett
   precise_real_type m1{0.};          ///< bino mass parameter
   precise_real_type m2{0.};          ///< wino mass parameter
   precise_real_type tan_beta{0.};    ///< mixing angle of the heavy Higgs precise_real_typett in the MSSM
   precise_real_type scale{0.};       ///< renormalization scale
   Eigen::Matrix<precise_real_type,3,3> mq2{Eigen::Matrix<precise_real_type,3,3>::Zero()}; ///< DR-bar squared soft-breaking left-handed squark mass parameters in the MSSM
   Eigen::Matrix<precise_real_type,3,3> mu2{Eigen::Matrix<precise_real_type,3,3>::Zero()}; ///< DR-bar squared soft-breaking right-handed up-squark mass parameters in the MSSM
   Eigen::Matrix<precise_real_type,3,3> md2{Eigen::Matrix<precise_real_type,3,3>::Zero()}; ///< DR-bar squared soft-breaking right-handed down-squark mass parameters in the MSSM
   Eigen::Matrix<precise_real_type,3,3> ml2{Eigen::Matrix<precise_real_type,3,3>::Zero()}; ///< DR-bar squared soft-breaking left-handed selectron mass parameters in the MSSM
   Eigen::Matrix<precise_real_type,3,3> me2{Eigen::Matrix<precise_real_type,3,3>::Zero()}; ///< DR-bar squared soft-breaking right-handed selectron mass parameters in the MSSM
};

std::ostream& operator<<(std::ostream&, const Parameters&);

precise_real_type lambda_tree_level(const Parameters&);
precise_real_type gYu_tree_level(const Parameters&);
precise_real_type gYd_tree_level(const Parameters&);
precise_real_type g2u_tree_level(const Parameters&);
precise_real_type g2d_tree_level(const Parameters&);

precise_real_type delta_lambda_1loop_reg(const Parameters&);
precise_real_type delta_lambda_1loop_phi(const Parameters&);
precise_real_type delta_lambda_1loop_chi_1(const Parameters&);
precise_real_type delta_lambda_1loop_chi_1(
   precise_real_type scale, precise_real_type mu, precise_real_type lambda, precise_real_type gYu, precise_real_type gYd,
   precise_real_type g2u, precise_real_type g2d, precise_real_type m1, precise_real_type m2);
precise_real_type delta_lambda_1loop_chi_2(const Parameters&);
precise_real_type delta_lambda_1loop_chi_2(
   precise_real_type scale, precise_real_type mu, precise_real_type m2, precise_real_type g1, precise_real_type g2, precise_real_type tan_beta);
precise_real_type delta_lambda_2loop_phi(const Parameters&);
precise_real_type delta_lambda_2loop_phi_HSS(const Parameters&);
precise_real_type delta_gYu_1loop(const Parameters&);
precise_real_type delta_gYd_1loop(const Parameters&);
precise_real_type delta_g2u_1loop(const Parameters&);
precise_real_type delta_g2d_1loop(const Parameters&);
precise_real_type delta_gt_1loop_chi(
   precise_real_type scale, precise_real_type mu, precise_real_type gYu, precise_real_type gYd,
   precise_real_type g2u, precise_real_type g2d, precise_real_type m1, precise_real_type m2);
precise_real_type delta_m2_1loop_chi(const Parameters&);

} // namespace splitmssm_thresholds
} // namespace flexiblesusy

#endif
