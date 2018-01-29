
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_higgs_loop_corrections

#include <boost/test/unit_test.hpp>

#include "test_SM.hpp"
#include "wrappers.hpp"
#include "pv.hpp"
#include "SM_two_scale_model.hpp"
#include "standard_model.hpp"
#include "sm_twoloophiggs.hpp"

using namespace flexiblesusy;
using namespace passarino_veltman;
using namespace flexiblesusy::sm_twoloophiggs;

BOOST_AUTO_TEST_CASE( test_SM_1loop_alpha_t )
{
   SM_input_parameters input;
   input.LambdaIN = 0.25;
   SM<Two_scale> m;
   setup_SM_const(m, input);

   m.calculate_DRbar_masses();

   const double mt = m.get_MFu(2);
   const double mt2 = Sqr(mt);
   const double yt = m.get_Yu(2,2);
   const double p = 125.;
   const double p2 = Sqr(p);
   const double scale = mt;
   const double scale2 = Sqr(scale);
   const double v = m.get_v();

   const double se_smh = delta_mh_1loop_at_sm(
      p, scale, mt, yt);

   // top loop for p = 0, Drees p. 8
   const double se_fs_eff =
      -6 * Sqr(yt) * oneOver16PiSqr * (
         ReA0(mt2, scale2)
         + 2*mt2*ReB0(p2, mt2, mt2, scale2));

   const double se_fs =
      se_fs_eff
      + 3.*Sqr(yt)*p2*ReB0(p2,mt2,mt2,scale2) * oneOver16PiSqr;

   // tadpole
   const double t_fs =
      6 * oneOver16PiSqr * (2*yt/Sqrt(2)) * mt / v * ReA0(mt2, scale2);

   BOOST_CHECK_CLOSE_FRACTION(se_smh, se_fs + t_fs, 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_SM_1loop_full )
{
   SM_input_parameters input;
   input.LambdaIN = 0.25;
   SM<Two_scale> m;
   setup_SM_const(m, input);

   m.set_Ye(ZEROMATRIX(3,3));
   m.set_Yd(ZEROMATRIX(3,3));
   m.set_Yu(0,0,0);
   m.set_Yu(1,1,0);
   m.solve_ewsb_tree_level();
   m.calculate_DRbar_masses();

   const double mt = m.get_MFu(2);
   const double mh = m.get_Mhh();
   const double mt2 = Sqr(mt);
   const double yt = m.get_Yu(2,2);
   const double p = mh; // this is assumed in the diagram with a Higgs loop
   const double p2 = Sqr(p);
   const double scale = mt;
   const double scale2 = Sqr(scale);
   const double v = m.get_v();
   const double gY = m.get_g1() * SM_info::normalization_g1;
   const double g2 = m.get_g2() * SM_info::normalization_g2;
   const double lambda = m.get_Lambdax()/2;

   const double se_smh = delta_mh_1loop_sm(
      p, scale, mt, yt, v, gY, g2, lambda);

   m.set_scale(scale);

   const double se_fs = Re(m.self_energy_hh_1loop(p));
   const double t_fs  = -Re(m.tadpole_hh_1loop() / v);

   // BOOST_CHECK_CLOSE_FRACTION(se_smh, se_fs + t_fs, 1.0e-10);
}

namespace {

std::pair<SM_mass_eigenstates, standard_model::Standard_model> make_point(int loops)
{
   SM_mass_eigenstates m1;
   standard_model::Standard_model m2;

   const double Q = 91.;
   const double vev = 246.;
   const double g1 = 0.4;
   const double g2 = 0.5;
   const double g3 = 1.1;

   m1.set_scale(Q);
   m1.set_v(vev);
   m1.set_mu2(vev*vev);
   m1.set_g1(g1);
   m1.set_g2(g2);
   m1.set_g3(g3);
   m1.set_Yu(0, 0, 0.001   * Sqrt(2.) / vev);
   m1.set_Yu(1, 1, 0.010   * Sqrt(2.) / vev);
   m1.set_Yu(2, 2, 165.0   * Sqrt(2.) / vev);
   m1.set_Yd(0, 0, 0.001   * Sqrt(2.) / vev);
   m1.set_Yd(1, 1, 0.010   * Sqrt(2.) / vev);
   m1.set_Yd(2, 2, 2.9     * Sqrt(2.) / vev);
   m1.set_Ye(0, 0, 0.001   * Sqrt(2.) / vev);
   m1.set_Ye(1, 1, 0.010   * Sqrt(2.) / vev);
   m1.set_Ye(2, 2, 1.77699 * Sqrt(2.) / vev);
   m1.set_Lambdax(0.2);

   m2.set_scale(Q);
   m2.set(m1.get());

   m1.do_calculate_sm_pole_masses(true);
   m1.set_loops(loops);
   m1.set_pole_mass_loop_order(loops);
   m1.set_ewsb_loop_order(loops);
   m1.solve_ewsb_tree_level();
   m1.calculate_DRbar_masses();
   m1.solve_ewsb();
   m1.calculate_pole_masses();

   m2.set_loops(loops);
   m2.set_pole_mass_loop_order(loops);
   m2.set_ewsb_loop_order(loops);
   m2.solve_ewsb_tree_level();
   m2.calculate_DRbar_masses();
   m2.solve_ewsb();
   m2.calculate_pole_masses();

   return std::make_pair(m1, m2);
}

void compare_Mh(int loops)
{
   const auto m = make_point(loops);

   double Mh_1 = 0., Mh_2 = 0.;

   if (loops > 0) {
      Mh_1 = m.first.get_physical().Mhh;
      Mh_2 = m.second.get_physical().Mhh;
   } else {
      Mh_1 = m.first.get_Mhh();
      Mh_2 = m.second.get_Mhh();
   }

   BOOST_REQUIRE(!m.first.get_problems().have_problem());
   BOOST_REQUIRE(!m.second.get_problems().have_problem());

   BOOST_TEST_MESSAGE("Mh(" << loops << "L) = " << Mh_1 << " = " << Mh_2);
   BOOST_CHECK_GT(Mh_1, 0.);
   BOOST_CHECK_CLOSE_FRACTION(Mh_1, Mh_2, 1e-10);
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE( test_SM_nloop_literature )
{
   const double eps = 1e-10;

   compare_Mh(0);
   compare_Mh(1);
   compare_Mh(2);
   compare_Mh(3);
   compare_Mh(4);
}
