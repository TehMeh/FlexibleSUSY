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

#ifndef lattice_numerical_constraint_hpp
#define lattice_numerical_constraint_hpp


#include <gsl/gsl_deriv.h>

#include "lattice_foreign_constraint.hpp"

namespace flexiblesusy {

struct NumericalConstraintCommon {
    static const std::vector<size_t> empty_vector;
    static constexpr Real default_epsilon = 1e-8;
};

class NumericalConstraint :
    public NumericalConstraintCommon,
    public ForeignConstraint {
public:
    NumericalConstraint(std::vector<size_t> dependence = empty_vector,
			Real epsilon = default_epsilon);
    void init(RGFlow<Lattice> *flow, size_t theory, size_t site);
    void operator()();
    using ForeignConstraint::init;
protected:
    virtual Real c(const Real *x) const = 0;
private:
    std::vector<double> x_local;
    std::vector<size_t> nonzeros;
    std::vector<bool> depends_on;
    gsl_function F_gsl;
    static double c_wrap(double xj, void *params);
    size_t j;
    Real deriv_epsilon;
};

class AnyNumericalConstraint : public NumericalConstraint {
public:
    AnyNumericalConstraint
    (std::function<Real(const AnyNumericalConstraint *, const Real *x)> fxn,
     std::vector<size_t> dependence = empty_vector,
     Real epsilon = default_epsilon) :
	NumericalConstraint(dependence, epsilon), fxn_(fxn)
	{}
protected:
    Real c(const Real *x) const { return fxn_(this, x); }
private:
    std::function<Real(const AnyNumericalConstraint *, const Real *x)> fxn_;
};

class NumericalMatching :
    public NumericalConstraintCommon,
    public ForeignMatching {
public:
    NumericalMatching(std::vector<size_t> depL = empty_vector,
		      std::vector<size_t> depH = empty_vector,
		      Real epsilon = default_epsilon);
    void init(RGFlow<Lattice> *flow, size_t lower_theory);
    void operator()();
    using ForeignMatching::init;
protected:
    virtual Real c(const Real *w, const Real *x) const = 0;
private:
    std::vector<double> w_local;
    std::vector<double> x_local;
    std::vector<size_t> nonzerosL, nonzerosH;
    std::vector<bool> depends_on;
    gsl_function F_gsl;
    static double c_wrap(double wxj, void *params);
    size_t j;
    Real &wx_local(size_t j)
    { return j < w_local.size() ? w_local[j] : x_local[j - w_local.size()]; }
    Real deriv_epsilon;
};

class AnyNumericalMatching : public NumericalMatching {
public:
    AnyNumericalMatching
    (std::function<Real(const AnyNumericalMatching *,
			const Real *w, const Real *x)> fxn,
     std::vector<size_t> depL = empty_vector,
     std::vector<size_t> depH = empty_vector,
     Real epsilon = default_epsilon) :
	NumericalMatching(depL, depH, epsilon), fxn_(fxn)
	{}
protected:
    Real c(const Real *w, const Real *x) const { return fxn_(this, w, x); }
private:
    std::function<Real(const AnyNumericalMatching *,
		       const Real *w, const Real *x)> fxn_;
};

}

#endif // lattice_numerical_constraint_hpp
