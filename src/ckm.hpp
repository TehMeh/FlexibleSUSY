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

#ifndef CKM_H
#define CKM_H

//#include <Eigen/Core>

namespace flexiblesusy {

struct CKM_parameters {
   void reset_to_diagonal();
   void reset_to_observation();

   void set_from_wolfenstein(precise_real_type, precise_real_type, precise_real_type, precise_real_type);
   void get_wolfenstein(precise_real_type&, precise_real_type&, precise_real_type&, precise_real_type&);

   Eigen::Matrix<precise_real_type,3,3> get_real_ckm() const;
   Eigen::Matrix<precise_complex_type,3,3> get_complex_ckm() const;

   static void to_pdg_convention(Eigen::Matrix<precise_real_type,3,3>&,
                                 Eigen::Matrix<precise_real_type,3,3>&,
                                 Eigen::Matrix<precise_real_type,3,3>&,
                                 Eigen::Matrix<precise_real_type,3,3>&,
                                 Eigen::Matrix<precise_real_type,3,3>&);
   static void to_pdg_convention(Eigen::Matrix<precise_real_type,3,3>&,
                                 Eigen::Matrix<precise_real_type,3,3>&,
                                 Eigen::Matrix<precise_real_type,3,3>&,
                                 Eigen::Matrix<precise_real_type,3,3>&);
   static void to_pdg_convention(Eigen::Matrix<precise_complex_type,3,3>&,
                                 Eigen::Matrix<precise_complex_type,3,3>&,
                                 Eigen::Matrix<precise_complex_type,3,3>&,
                                 Eigen::Matrix<precise_complex_type,3,3>&,
                                 Eigen::Matrix<precise_complex_type,3,3>&);
   static void to_pdg_convention(Eigen::Matrix<precise_complex_type,3,3>&,
                                 Eigen::Matrix<precise_complex_type,3,3>&,
                                 Eigen::Matrix<precise_complex_type,3,3>&,
                                 Eigen::Matrix<precise_complex_type,3,3>&);

   precise_real_type theta_12{0.}, theta_13{0.}, theta_23{0.}, delta{0.};
};

} // namespace flexiblesusy

#endif
