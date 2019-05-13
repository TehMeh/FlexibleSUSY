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

#include "rg_flow.hpp"
#include "root_finder.hpp"
#include <cstddef>
#include <vector>
#include <Eigen/Core>

namespace flexiblesusy {

/**
 * @class Shooting_solver
 * @brief Boundary value problem solver for shooting method.  This
 * solver encapsulates the root_finder.
 */
template<std::size_t N>
class Shooting_solver {
public:
   using Vec_t = Eigen::Matrix<double, N, 1>;
   using Fun_t = std::function<Vec_t(const Vec_t&)>;

   /// returns solution
   Vec_t get_solution() const { return point; }
   /// set precision goal
   void set_precision(double p) { precision = p; }
   /// set maximum number of iterations
   void set_max_iterations(std::size_t n) { max_it = n; }
   /// clear all internal data
   void reset();
   /// solves the boundary value problem with initial guess
   void solve(const Fun_t&, const Vec_t&);

private:
   Vec_t point{Vec_t::Zero()};
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
   std::vector<Root_finder<N>> rfs = {
      Root_finder<N>(fun, max_it, precision, Root_finder<N>::GSLHybridS),
      Root_finder<N>(fun, max_it, precision, Root_finder<N>::GSLHybrid),
      Root_finder<N>(fun, max_it, precision, Root_finder<N>::GSLBroyden)
   };

   auto err = !GSL_SUCCESS;

   for (auto& rf: rfs) {
      err = rf.solve(v0);

      if (err == GSL_SUCCESS) {
         point = rf.get_solution();
         break;
      }
   }

   if (err != GSL_SUCCESS)
      throw NoConvergenceError(max_it);
}

} // namespace flexiblesusy

#endif
