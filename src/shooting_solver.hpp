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

#ifndef SHOOTING_SOLVER_H
#define SHOOTING_SOLVER_H

#include "composite_root_finder.hpp"
#include <cstddef>
#include <vector>
#include <Eigen/Core>

namespace flexiblesusy {

/**
 * @class Shooting_solver
 * @brief Solves the BVP using the shooting method
 */
template<std::size_t N>
class Shooting_solver {
public:
   using Vec_t = Eigen::Matrix<double, N, 1>;
   using Fun_t = std::function<Vec_t(const Vec_t&)>;

   /// returns solution
   Vec_t get_solution() const { return solution; }
   /// set precision goal
   void set_precision(double p) { precision = p; }
   /// set maximum number of iterations
   void set_max_iterations(std::size_t n) { max_it = n; }
   /// clear all internal data
   void reset();
   /// solves the boundary value problem with initial guess
   void solve(const Fun_t&, const Vec_t&);

private:
   Vec_t solution{Vec_t::Zero()};
   double precision{1e-3};
   std::size_t max_it{100};
};

template<std::size_t N>
void Shooting_solver<N>::reset()
{
   *this = Shooting_solver<N>();
}

template<std::size_t N>
void Shooting_solver<N>::solve(const Fun_t& fun, const Vec_t& v0)
{
   Composite_root_solver<N> crs;
   crs.set_precision(precision);
   crs.set_max_iterations(max_it);
   crs.solve(fun, v0);

   solution = crs.get_solution();
}

} // namespace flexiblesusy

#endif
