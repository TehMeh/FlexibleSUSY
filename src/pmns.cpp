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

#include "pmns.hpp"
#include "ew_input.hpp"
#include "numerics2.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

void PMNS_parameters::reset_to_diagonal()
{
   theta_12 = 0.;
   theta_13 = 0.;
   theta_23 = 0.;
   delta    = 0.;
   alpha_1  = 0.;
   alpha_2  = 0.;
}

void PMNS_parameters::reset_to_observation()
{
   theta_12 = Electroweak_constants::PMNS_THETA12;
   theta_13 = Electroweak_constants::PMNS_THETA13;
   theta_23 = Electroweak_constants::PMNS_THETA23;
   delta    = Electroweak_constants::PMNS_DELTA;
   alpha_1  = Electroweak_constants::PMNS_ALPHA1;
   alpha_2  = Electroweak_constants::PMNS_ALPHA2;
}

Eigen::Matrix<precise_real_type,3,3> PMNS_parameters::get_real_pmns() const
{
   const precise_complex_type eID(polar(1.0, delta));
   const precise_real_type s12 = Sin(theta_12);
   const precise_real_type s13 = Sin(theta_13);
   const precise_real_type s23 = Sin(theta_23);
   const precise_real_type c12 = Cos(theta_12);
   const precise_real_type c13 = Cos(theta_13);
   const precise_real_type c23 = Cos(theta_23);

   // set phase factor e^(i delta) to +1 or -1 depending on the sign
   // of its real part
   const int pf = Sign(Re(eID));

   Eigen::Matrix<precise_real_type,3,3> pmns_matrix;
   pmns_matrix(0, 0) = c12 * c13;
   pmns_matrix(0, 1) = s12 * c13;
   pmns_matrix(0, 2) = pf * s13;
   pmns_matrix(1, 0) = -s12 * c23 - pf * c12 * s23 * s13;
   pmns_matrix(1, 1) = c12 * c23 - pf * s12 * s23 * s13;
   pmns_matrix(1, 2) = s23 * c13;
   pmns_matrix(2, 0) = s12 * s23 - pf * c12 * c23 * s13;
   pmns_matrix(2, 1) = -c12 * s23 - pf * s12 * c23 * s13;
   pmns_matrix(2, 2) = c23 * c13;

   return pmns_matrix;
}

Eigen::Matrix<precise_complex_type,3,3> PMNS_parameters::get_complex_pmns() const
{
   const precise_complex_type eID(polar(1.0, delta));
   const precise_complex_type eIAlpha1(polar(1.0, 0.5*alpha_1));
   const precise_complex_type eIAlpha2(polar(1.0, 0.5*alpha_2));
   const precise_real_type s12 = Sin(theta_12);
   const precise_real_type s13 = Sin(theta_13);
   const precise_real_type s23 = Sin(theta_23);
   const precise_real_type c12 = Cos(theta_12);
   const precise_real_type c13 = Cos(theta_13);
   const precise_real_type c23 = Cos(theta_23);

   Eigen::Matrix<precise_complex_type,3,3> pmns_matrix;
   pmns_matrix(0, 0) = c12 * c13 * eIAlpha1;
   pmns_matrix(0, 1) = s12 * c13 * eIAlpha2;
   pmns_matrix(0, 2) = s13 / eID;
   pmns_matrix(1, 0) = (-s12 * c23 - c12 * s23 * s13 * eID) * eIAlpha1;
   pmns_matrix(1, 1) = (c12 * c23 - s12 * s23 * s13 * eID) * eIAlpha2;
   pmns_matrix(1, 2) = s23 * c13;
   pmns_matrix(2, 0) = (s12 * s23 - c12 * c23 * s13 * eID) * eIAlpha1;
   pmns_matrix(2, 1) = (-c12 * s23 - s12 * c23 * s13 * eID) * eIAlpha2;
   pmns_matrix(2, 2) = c23 * c13;

   return pmns_matrix;
}

void PMNS_parameters::to_pdg_convention(Eigen::Matrix<precise_real_type,3,3>& Vv,
                                        Eigen::Matrix<precise_real_type,3,3>& Ve,
                                        Eigen::Matrix<precise_real_type,3,3>& Ue)
{
   Eigen::Matrix<precise_real_type,3,3> pmns(Ve*Vv.adjoint());
   to_pdg_convention(pmns, Vv, Ve, Ue);
}

void PMNS_parameters::to_pdg_convention(Eigen::Matrix<precise_real_type,3,3>& pmns,
                                        Eigen::Matrix<precise_real_type,3,3>& Vv,
                                        Eigen::Matrix<precise_real_type,3,3>& Ve,
                                        Eigen::Matrix<precise_real_type,3,3>& Ue)
{
   Eigen::Matrix<precise_real_type,3,3> signs_E(Eigen::Matrix<precise_real_type,3,3>::Identity());
   Eigen::Matrix<precise_real_type,3,3> signs_V(Eigen::Matrix<precise_real_type,3,3>::Identity());

   // make 33 element positive
   if (pmns(2, 2) < 0.) {
      signs_E(2, 2) = -1.;
      for (int j = 0; j < 3; ++j) {
         pmns(2, j) *= -1.;
      }
   }

   // make 23 element positive
   if (pmns(1, 2) < 0.) {
      signs_V(2, 2) = -1;
      signs_E(2, 2) *= -1;
      for (int j = 0; j < 3; ++j) {
         pmns(2, j) *= -1;
         pmns(j, 2) *= -1;
      }
   }

   Ve = signs_E * Ve;
   Ue = signs_E * Ue;
   Vv = signs_V * Vv;
}

void PMNS_parameters::to_pdg_convention(
   Eigen::Matrix<precise_complex_type,3,3>& Vv,
   Eigen::Matrix<precise_complex_type,3,3>& Ve,
   Eigen::Matrix<precise_complex_type,3,3>& Ue)
{
   Eigen::Matrix<precise_complex_type,3,3> pmns(Ve*Vv.adjoint());
   to_pdg_convention(pmns, Vv, Ve, Ue);
}

namespace {
//ATTENTION S.D.
precise_complex_type phase(const precise_complex_type& z)
{
   precise_real_type r = abs(z);
   return r == 0 ? precise_complex_type(1) : precise_complex_type(z/r);
}

void calc_phase_factors(
   const Eigen::Matrix<precise_complex_type,3,3>& pmns,
   const precise_complex_type& p,
   precise_complex_type& o,
   Eigen::DiagonalMatrix<precise_complex_type,3>& l)
{  
   o = conj(phase(p * pmns(0,2)));
   l.diagonal().bottomRightCorner<2,1>() = (o * pmns.bottomRightCorner<2,1>()).
      unaryExpr(std::ptr_fun(phase)).conjugate();
}

/// restrict sin or cos to interval [-1,1]
precise_real_type sanitize_hypot(precise_real_type sc)
{
   if (sc < -1.) sc = -1.;
   if (sc > 1.) sc = 1.;
   return sc;
}

} // anonymous namespace

void PMNS_parameters::to_pdg_convention(
   Eigen::Matrix<precise_complex_type,3,3>& pmns,
   Eigen::Matrix<precise_complex_type,3,3>& Vv,
   Eigen::Matrix<precise_complex_type,3,3>& Ve,
   Eigen::Matrix<precise_complex_type,3,3>& Ue)
{
   precise_complex_type o;
   Eigen::DiagonalMatrix<precise_complex_type,3> l(1,1,1);

   const precise_real_type s13 = sanitize_hypot(abs(pmns(0,2)));
   const precise_real_type c13_sq = 1. - Sqr(s13);
   if (is_zero(c13_sq)) {
      o = conj(phase(pmns(0,2)));
      const auto rel_phase = sqrt(phase(pmns(1,0) * pmns(2,1)));
      const auto p = conj(rel_phase * rel_phase)
         * phase(pmns(1,0) * pmns(1,1));
      l.diagonal()[1] = conj(o * rel_phase);
      l.diagonal()[2] = conj(o * rel_phase) * p;
   } else {
      const precise_real_type c13 = sqrt(c13_sq);
      const precise_real_type s12 = sanitize_hypot(abs(pmns(0,1)) / c13);
      const precise_real_type c12 = sqrt(1 - Sqr(s12));
      const precise_real_type s23 = sanitize_hypot(abs(pmns(1,2)) / c13);
      const precise_real_type c23 = sqrt(1 - Sqr(s23));
      const precise_real_type jcp = imag(pmns(1,2) * conj(pmns(0,2))
                                   * pmns(0,1) * conj(pmns(1,1)));
      const precise_real_type side1 = s12*s23;
      const precise_real_type side2 = c12*c23*s13;
      const precise_real_type cosdelta = sanitize_hypot(
         is_zero(jcp) ? precise_real_type(1) // delta is removable
         : (Sqr(side1)+Sqr(side2)-norm(pmns(2,0))) / (2*side1*side2));
      const precise_real_type sindelta = sqrt(1 - Sqr(cosdelta));
      const precise_complex_type p(cosdelta, sindelta); // Exp[I delta]
      calc_phase_factors(pmns, p, o, l);

      const auto a1 = phase(o*pmns(0,0));
      const auto a2 = phase(o*pmns(0,1));

      Eigen::Matrix<precise_complex_type,2,2> pmnsBL{
         (o * l * pmns).bottomLeftCorner<2,2>()};
      Eigen::Array<precise_complex_type,2,2> maj_phases;
      maj_phases << a1, a2, a1, a2;
      const Eigen::Array<precise_real_type,2,2> imagBL{
         (maj_phases.conjugate() * pmnsBL.array()).imag()};
      if (!((imagBL <= 0).all() || (imagBL >= 0).all())) {
         calc_phase_factors(pmns, conj(p), o, l);
      }
   }

   Ve.transpose() *= l * o;
   Ue.transpose() *= (l * o).diagonal().conjugate().asDiagonal();
   pmns = Ve * Vv.adjoint();
}

} // namespace flexiblesusy
