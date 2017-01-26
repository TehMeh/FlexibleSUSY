#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMSemiAnalytic_ewsb

#include <boost/test/unit_test.hpp>

#include "test_CMSSMSemiAnalytic.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_ewsb_solver.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_solutions.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_CMSSMSemiAnalytic_ewsb_tree_level_solution )
{
   const double precision = 1.0e-5;
   CMSSMSemiAnalytic_input_parameters input;
   CMSSMSemiAnalytic_mass_eigenstates model(input);
   setup_CMSSMSemiAnalytic(model, input);

   const double high_scale = 2.e16;
   model.run_to(high_scale);

   Boundary_values values;
   setup_high_scale_CMSSMSemiAnalytic(model, values);

   CMSSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);
   solns.set_AzeroBasis(values.Azero);
   solns.set_m12Basis(values.m12);
   solns.set_m0SqBasis(values.m0Sq);
   solns.set_BMu0Basis(values.BMu0);
   solns.set_MuBasis(values.Mu);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);

   CMSSMSemiAnalytic_ewsb_solver<Semi_analytic> ewsb_solver;
   ewsb_solver.set_loop_order(0);
   ewsb_solver.set_semi_analytic_solutions(&solns);
   model.set_ewsb_solver(&ewsb_solver);

   std::cout << "initially, m0Sq = " << model.get_m0Sq() << '\n';
   model.set_ewsb_iteration_precision(precision);
   const int error = model.solve_ewsb_tree_level();
   std::cout << "finally, m0Sq = " << model.get_m0Sq() << '\n';
   BOOST_CHECK_EQUAL(error, 0);

   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(model.get_ewsb_eq_hh_2(), precision);
}
