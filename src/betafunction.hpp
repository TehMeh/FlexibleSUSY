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

/**
 * @file betafunction.hpp
 * @brief contains class Beta_function
 */

#include "precise.hpp"

#ifndef BETAFUNCTION_H
#define BETAFUNCTION_H

#include "basic_rk_integrator.hpp"

#include <Eigen/Core>

namespace flexiblesusy {

/**
 * @class Beta_function
 * @brief beta function interface
 *
 * Beta_function is the abstract base class for the beta functions of
 * the parameter classes of the generated models.  It defines the
 * basic RG running interface.  The run() and run_to() functions use
 * the Runge-Kutta algorithm to integrate the RGEs up to a given
 * scale.
 */
class Beta_function {
public:
   using Derivs = std::function<Eigen::ArrayXdp(precise_real_type, const Eigen::ArrayXdp&)>;
   using ODE_integrator = std::function<void(precise_real_type, precise_real_type, Eigen::ArrayXdp&, Derivs, precise_real_type)>;

   Beta_function() = default;
   Beta_function(const Beta_function&) = default;
   Beta_function(Beta_function&&) = default;
   virtual ~Beta_function() = default;
   Beta_function& operator=(const Beta_function&) = default;
   Beta_function& operator=(Beta_function&&) = default;

   void set_scale(precise_real_type s) { scale = s; }
   void set_number_of_parameters(int pars) { num_pars = pars; }
   void set_loops(precise_real_type l) { loops = (int)l; }
   void set_thresholds(precise_real_type t) { thresholds = (int)t; }
   void set_zero_threshold(precise_real_type t) { zero_threshold = t; }
   void set_integrator(const ODE_integrator& i) { integrator = i; }

   precise_real_type get_scale() const { return scale; }
   int get_number_of_parameters() const { return num_pars; }
   int get_loops() const { return loops; }
   int get_thresholds() const { return thresholds; }
   precise_real_type get_zero_threshold() const { return zero_threshold; }

   void reset();

   virtual Eigen::ArrayXdp get() const = 0;
   virtual void set(const Eigen::ArrayXdp&) = 0;
   virtual Eigen::ArrayXdp beta() const = 0;

   virtual void run(precise_real_type, precise_real_type, precise_real_type eps = -1.0);
   virtual void run_to(precise_real_type, precise_real_type eps = -1.0);

protected:
   void call_rk(precise_real_type, precise_real_type, Eigen::ArrayXdp&, Derivs, precise_real_type eps = -1.0);

private:
   int num_pars{0};              ///< number of parameters
   int loops{0};                 ///< to what loop order does the RG evolution run
   int thresholds{0};            ///< threshold correction loop order
   precise_real_type scale{0.};             ///< current renormalization scale
   precise_real_type tolerance{1.e-4};      ///< running tolerance
   precise_real_type min_tolerance{1.e-11}; ///< minimum tolerance allowed
   precise_real_type zero_threshold{1.e-11};///< threshold for treating values as zero
   ODE_integrator integrator{
      runge_kutta::Basic_rk_integrator<Eigen::ArrayXdp>()};

   Eigen::ArrayXdp derivatives(precise_real_type, const Eigen::ArrayXdp&);
   precise_real_type get_tolerance(precise_real_type eps);
};

} // namespace flexiblesusy

#endif
