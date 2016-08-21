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

#ifndef GSL_UTILS_H
#define GSL_UTILS_H

#include <gsl/gsl_vector.h>
#include <Eigen/Core>

namespace flexiblesusy {

class GSL_vector;

/// Returns true if GSL_vector contains only finite elements, false otherwise
bool is_finite(const GSL_vector&);
/// Returns true if GSL vector contains only finite elements, false otherwise
bool is_finite(const gsl_vector*);
Eigen::ArrayXd to_eigen_array(const gsl_vector*);
Eigen::VectorXd to_eigen_vector(const gsl_vector*);
Eigen::VectorXd to_eigen_vector(const GSL_vector&);
GSL_vector to_GSL_vector(const Eigen::VectorXd&);
GSL_vector to_GSL_vector(const gsl_vector*);
void copy(const Eigen::ArrayXd&, gsl_vector*);

}

#endif
