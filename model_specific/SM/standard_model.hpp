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
 * @file standard_model.hpp
 *
 * @brief contains class for Standard model running and self-energies
 *
 */

#include "precise.hpp"

#ifndef STANDARD_MODEL_H
#define STANDARD_MODEL_H

#include "betafunction.hpp"
#include "standard_model_physical.hpp"
#include "loop_corrections.hpp"
#include "threshold_corrections.hpp"
#include "error.hpp"
#include "problems.hpp"
#include "physical_input.hpp"

#include <array>
#include <iosfwd>
#include <string>

#include <Eigen/Core>

namespace softsusy {
   class QedQcd;
} // namespace softsusy

namespace flexiblesusy {

class EWSB_solver;

namespace standard_model_info {

   enum Particles : int {VG, Hp, Fv, Ah, hh, VP, VZ, Fd, Fu, Fe, VWp,
      NUMBER_OF_PARTICLES};

   enum Parameters : int {g1, g2, g3, Lambdax, Yu0_0, Yu0_1, Yu0_2, Yu1_0,
      Yu1_1, Yu1_2, Yu2_0, Yu2_1, Yu2_2, Yd0_0, Yd0_1, Yd0_2, Yd1_0, Yd1_1, Yd1_2
      , Yd2_0, Yd2_1, Yd2_2, Ye0_0, Ye0_1, Ye0_2, Ye1_0, Ye1_1, Ye1_2, Ye2_0,
      Ye2_1, Ye2_2, mu2, v, NUMBER_OF_PARAMETERS};

   enum Mixings : int {ReVd00, ImVd00, ReVd01, ImVd01, ReVd02, ImVd02,
      ReVd10, ImVd10, ReVd11, ImVd11, ReVd12, ImVd12, ReVd20, ImVd20, ReVd21,
      ImVd21, ReVd22, ImVd22, ReUd00, ImUd00, ReUd01, ImUd01, ReUd02, ImUd02,
      ReUd10, ImUd10, ReUd11, ImUd11, ReUd12, ImUd12, ReUd20, ImUd20, ReUd21,
      ImUd21, ReUd22, ImUd22, ReVu00, ImVu00, ReVu01, ImVu01, ReVu02, ImVu02,
      ReVu10, ImVu10, ReVu11, ImVu11, ReVu12, ImVu12, ReVu20, ImVu20, ReVu21,
      ImVu21, ReVu22, ImVu22, ReUu00, ImUu00, ReUu01, ImUu01, ReUu02, ImUu02,
      ReUu10, ImUu10, ReUu11, ImUu11, ReUu12, ImUu12, ReUu20, ImUu20, ReUu21,
      ImUu21, ReUu22, ImUu22, ReVe00, ImVe00, ReVe01, ImVe01, ReVe02, ImVe02,
      ReVe10, ImVe10, ReVe11, ImVe11, ReVe12, ImVe12, ReVe20, ImVe20, ReVe21,
      ImVe21, ReVe22, ImVe22, ReUe00, ImUe00, ReUe01, ImUe01, ReUe02, ImUe02,
      ReUe10, ImUe10, ReUe11, ImUe11, ReUe12, ImUe12, ReUe20, ImUe20, ReUe21,
      ImUe21, ReUe22, ImUe22, NUMBER_OF_MIXINGS};

   extern const precise_real_type normalization_g1;
   extern const precise_real_type normalization_g2;
   extern const precise_real_type normalization_g3;

   extern const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities;
   extern const std::array<std::string, NUMBER_OF_PARTICLES> particle_names;
   extern const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names;
   extern const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names;
   extern const std::array<std::string, NUMBER_OF_MIXINGS> particle_mixing_names;
   extern const std::string model_name;
   constexpr bool is_low_energy_model = false;
   constexpr bool is_supersymmetric_model = false;

   class Standard_model_particle_names : public Names {
   public:
      virtual ~Standard_model_particle_names() = default;
      virtual const std::string& get(int index) const override {
         return particle_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARAMETERS;
      }
   };

   class Standard_model_parameter_names : public Names {
   public:
      virtual ~Standard_model_parameter_names() = default;
      virtual const std::string& get(int index) const override {
         return parameter_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARTICLES;
      }
   };

   const Standard_model_particle_names  particle_names_getter{};
   const Standard_model_parameter_names parameter_names_getter{};

} // namespace standard_model_info

namespace standard_model {

template <class T>
class StandardModel;

/**
 * @class Standard_model
 * @brief model class with routines for SM running and self-energies
 */
class Standard_model : public Beta_function {
public:

   Standard_model();
   Standard_model(precise_real_type scale_, precise_real_type loops_, precise_real_type thresholds_
   , precise_real_type g1_, precise_real_type g2_, precise_real_type g3_, precise_real_type Lambdax_, const Eigen::Matrix<
   precise_real_type,3,3>& Yu_, const Eigen::Matrix<precise_real_type,3,3>& Yd_, const Eigen::Matrix<
   precise_real_type,3,3>& Ye_, precise_real_type mu2_, precise_real_type v_);
   Standard_model(const Standard_model&) = default;
   Standard_model(Standard_model&&) = default;

   virtual ~Standard_model() = default;

   Standard_model& operator=(const Standard_model&) = default;
   Standard_model& operator=(Standard_model&&) = default;

   /// number of EWSB equations
   static const int number_of_ewsb_equations = 1;

   void calculate_DRbar_masses();
   void calculate_pole_masses();
   void check_pole_masses_for_tachyons();
   void do_force_output(bool);
   bool do_force_output() const;
   void set_ewsb_iteration_precision(precise_real_type);
   void set_ewsb_loop_order(int);
   void set_loop_corrections(const Loop_corrections&);
   const Loop_corrections& get_loop_corrections() const;
   void set_threshold_corrections(const Threshold_corrections&);
   const Threshold_corrections& get_threshold_corrections() const;
   void set_pole_mass_loop_order(int);
   int get_pole_mass_loop_order() const;
   void set_physical(const Standard_model_physical&);
   precise_real_type get_ewsb_iteration_precision() const;
   precise_real_type get_ewsb_loop_order() const;
   const Standard_model_physical& get_physical() const;
   Standard_model_physical& get_physical();
   const Problems& get_problems() const;
   Problems& get_problems();
   int solve_ewsb_tree_level();
   int solve_ewsb_one_loop();
   int solve_ewsb();            ///< solve EWSB at ewsb_loop_order level

   virtual Eigen::ArrayXdp beta() const override;
   virtual Eigen::ArrayXdp get() const override;
   void print(std::ostream& out = std::cerr) const;
   virtual void set(const Eigen::ArrayXdp&) override;

   Standard_model calc_beta() const;
   Standard_model calc_beta(int) const;
   void clear();
   void clear_running_parameters();
   void clear_DRbar_parameters();
   void clear_problems();

   void calculate_spectrum();
   std::string name() const;
   virtual void run_to(precise_real_type scale, precise_real_type eps = -1.0) override;
   void set_precision(precise_real_type);
   precise_real_type get_precision() const;

   void set_physical_input(const Physical_input& input_) { input = input_; }
   const Physical_input& get_physical_input() const { return input; }
   Physical_input& get_physical_input() { return input; }

   void initialise_from_input(const softsusy::QedQcd&);

   void set_g1(precise_real_type g1_) { g1 = g1_; }
   void set_g2(precise_real_type g2_) { g2 = g2_; }
   void set_g3(precise_real_type g3_) { g3 = g3_; }
   void set_Lambdax(precise_real_type Lambdax_) { Lambdax = Lambdax_; }
   void set_Yu(const Eigen::Matrix<precise_real_type,3,3>& Yu_) { Yu = Yu_; }
   void set_Yu(int i, int k, precise_real_type value) { Yu(i,k) = value; }
   void set_Yd(const Eigen::Matrix<precise_real_type,3,3>& Yd_) { Yd = Yd_; }
   void set_Yd(int i, int k, precise_real_type value) { Yd(i,k) = value; }
   void set_Ye(const Eigen::Matrix<precise_real_type,3,3>& Ye_) { Ye = Ye_; }
   void set_Ye(int i, int k, precise_real_type value) { Ye(i,k) = value; }
   void set_mu2(precise_real_type mu2_) { mu2 = mu2_; }
   void set_v(precise_real_type v_) { v = v_; }

   precise_real_type get_g1() const { return g1; }
   precise_real_type get_g2() const { return g2; }
   precise_real_type get_g3() const { return g3; }
   precise_real_type get_Lambdax() const { return Lambdax; }
   const Eigen::Matrix<precise_real_type,3,3>& get_Yu() const { return Yu; }
   precise_real_type get_Yu(int i, int k) const { return Yu(i,k); }
   const Eigen::Matrix<precise_real_type,3,3>& get_Yd() const { return Yd; }
   precise_real_type get_Yd(int i, int k) const { return Yd(i,k); }
   const Eigen::Matrix<precise_real_type,3,3>& get_Ye() const { return Ye; }
   precise_real_type get_Ye(int i, int k) const { return Ye(i,k); }
   precise_real_type get_mu2() const { return mu2; }
   precise_real_type get_v() const { return v; }

   precise_real_type get_MVG() const { return MVG; }
   precise_real_type get_MHp() const { return MHp; }
   const Eigen::Array<precise_real_type,3,1>& get_MFv() const { return MFv; }
   precise_real_type get_MFv(int i) const { return MFv(i); }
   precise_real_type get_MAh() const { return MAh; }
   precise_real_type get_Mhh() const { return Mhh; }
   precise_real_type get_MVP() const { return MVP; }
   precise_real_type get_MVZ() const { return MVZ; }
   const Eigen::Array<precise_real_type,3,1>& get_MFd() const { return MFd; }
   precise_real_type get_MFd(int i) const { return MFd(i); }
   const Eigen::Array<precise_real_type,3,1>& get_MFu() const { return MFu; }
   precise_real_type get_MFu(int i) const { return MFu(i); }
   const Eigen::Array<precise_real_type,3,1>& get_MFe() const { return MFe; }
   precise_real_type get_MFe(int i) const { return MFe(i); }
   precise_real_type get_MVWp() const { return MVWp; }
   const Eigen::Array<precise_real_type,2,1>& get_MVPVZ() const { return MVPVZ; }
   precise_real_type get_MVPVZ(int i) const { return MVPVZ(i); }

   const Eigen::Matrix<precise_complex_type,3,3>& get_Vd() const { return Vd; }
   const precise_complex_type& get_Vd(int i, int k) const { return Vd(i,k); }
   const Eigen::Matrix<precise_complex_type,3,3>& get_Ud() const { return Ud; }
   const precise_complex_type& get_Ud(int i, int k) const { return Ud(i,k); }
   const Eigen::Matrix<precise_complex_type,3,3>& get_Vu() const { return Vu; }
   const precise_complex_type& get_Vu(int i, int k) const { return Vu(i,k); }
   const Eigen::Matrix<precise_complex_type,3,3>& get_Uu() const { return Uu; }
   const precise_complex_type& get_Uu(int i, int k) const { return Uu(i,k); }
   const Eigen::Matrix<precise_complex_type,3,3>& get_Ve() const { return Ve; }
   const precise_complex_type& get_Ve(int i, int k) const { return Ve(i,k); }
   const Eigen::Matrix<precise_complex_type,3,3>& get_Ue() const { return Ue; }
   const precise_complex_type& get_Ue(int i, int k) const { return Ue(i,k); }
   const Eigen::Matrix<precise_real_type,2,2>& get_ZZ() const { return ZZ; }
   precise_real_type get_ZZ(int i, int k) const { return ZZ(i,k); }

   precise_real_type get_mass_matrix_VG() const;
   void calculate_MVG();
   precise_real_type get_mass_matrix_Hp() const;
   void calculate_MHp();
   Eigen::Matrix<precise_real_type,3,3> get_mass_matrix_Fv() const;
   void calculate_MFv();
   precise_real_type get_mass_matrix_Ah() const;
   void calculate_MAh();
   precise_real_type get_mass_matrix_hh() const;
   void calculate_Mhh();
   precise_real_type get_mass_matrix_VP() const;
   void calculate_MVP();
   precise_real_type get_mass_matrix_VZ() const;
   void calculate_MVZ();
   Eigen::Matrix<precise_real_type,3,3> get_mass_matrix_Fd() const;
   void calculate_MFd();
   Eigen::Matrix<precise_real_type,3,3> get_mass_matrix_Fu() const;
   void calculate_MFu();
   Eigen::Matrix<precise_real_type,3,3> get_mass_matrix_Fe() const;
   void calculate_MFe();
   precise_real_type get_mass_matrix_VWp() const;
   void calculate_MVWp();
   Eigen::Matrix<precise_real_type,2,2> get_mass_matrix_VPVZ() const;
   void calculate_MVPVZ();

   precise_real_type get_ewsb_eq_hh_1() const;

   precise_real_type CpconjHpHphh() const;
   precise_real_type CpconjHpVWpVP() const;
   precise_real_type CpconjHpVZVWp() const;
   precise_real_type CpHpgWpCbargZ() const;
   precise_real_type CpconjHpbargWpCgZ() const;
   precise_real_type CpHpgZbargWp() const;
   precise_real_type CpconjHpbargZgWp() const;
   precise_real_type CpHpconjHpAhAh() const;
   precise_real_type CpHpconjHphhhh() const;
   precise_real_type CpHpconjHpconjHpHp() const;
   precise_complex_type CpconjHpVWpAh() const;
   precise_real_type CpconjHpVWphh() const;
   precise_real_type CpconjHpVPHp() const;
   precise_real_type CpconjHpVZHp() const;
   precise_real_type CpHpconjHpconjVWpVWp() const;
   precise_complex_type CpHpconjHpVZVZ() const;
   precise_complex_type CpconjHpbarFdFuPR(int gI1, int gI2) const;
   precise_complex_type CpconjHpbarFdFuPL(int gI1, int gI2) const;
   precise_real_type CpconjHpbarFeFvPR(int , int ) const;
   precise_complex_type CpconjHpbarFeFvPL(int gI1, int gI2) const;
   precise_real_type CpAhhhAh() const;
   precise_complex_type CpAhbargWpgWp() const;
   precise_complex_type CpAhbargWpCgWpC() const;
   precise_real_type CpAhAhAhAh() const;
   precise_real_type CpAhAhhhhh() const;
   precise_real_type CpAhAhconjHpHp() const;
   precise_complex_type CpAhVZhh() const;
   precise_complex_type CpAhconjVWpHp() const;
   precise_real_type CpAhAhconjVWpVWp() const;
   precise_complex_type CpAhAhVZVZ() const;
   precise_complex_type CpAhbarFdFdPR(int gI1, int gI2) const;
   precise_complex_type CpAhbarFdFdPL(int gI1, int gI2) const;
   precise_complex_type CpAhbarFeFePR(int gI1, int gI2) const;
   precise_complex_type CpAhbarFeFePL(int gI1, int gI2) const;
   precise_complex_type CpAhbarFuFuPR(int gI1, int gI2) const;
   precise_complex_type CpAhbarFuFuPL(int gI1, int gI2) const;
   precise_real_type CphhAhAh() const;
   precise_real_type Cphhhhhh() const;
   precise_real_type CphhVZVZ() const;
   precise_real_type CphhbargWpgWp() const;
   precise_real_type CphhbargWpCgWpC() const;
   precise_real_type CphhbargZgZ() const;
   precise_real_type CphhconjHpHp() const;
   precise_real_type CphhconjVWpVWp() const;
   precise_real_type CphhhhAhAh() const;
   precise_real_type Cphhhhhhhh() const;
   precise_real_type CphhhhconjHpHp() const;
   precise_complex_type CphhVZAh() const;
   precise_real_type CphhconjVWpHp() const;
   precise_real_type CphhhhconjVWpVWp() const;
   precise_complex_type CphhhhVZVZ() const;
   precise_complex_type CphhbarFdFdPR(int gI1, int gI2) const;
   precise_complex_type CphhbarFdFdPL(int gI1, int gI2) const;
   precise_complex_type CphhbarFeFePR(int gI1, int gI2) const;
   precise_complex_type CphhbarFeFePL(int gI1, int gI2) const;
   precise_complex_type CphhbarFuFuPR(int gI1, int gI2) const;
   precise_complex_type CphhbarFuFuPL(int gI1, int gI2) const;
   precise_complex_type CpVZhhAh() const;
   precise_real_type CpVZVZhh() const;
   precise_real_type CpVZbargWpgWp() const;
   precise_real_type CpVZbargWpCgWpC() const;
   precise_real_type CpVZconjHpHp() const;
   precise_real_type CpVZconjVWpHp() const;
   precise_complex_type CpVZVZAhAh() const;
   precise_complex_type CpVZVZhhhh() const;
   precise_complex_type CpVZVZconjHpHp() const;
   precise_real_type CpVZconjVWpVWp() const;
   precise_real_type CpVZbarFdFdPL(int gI1, int gI2) const;
   precise_real_type CpVZbarFdFdPR(int gI1, int gI2) const;
   precise_real_type CpVZbarFeFePL(int gI1, int gI2) const;
   precise_real_type CpVZbarFeFePR(int gI1, int gI2) const;
   precise_real_type CpVZbarFuFuPL(int gI1, int gI2) const;
   precise_real_type CpVZbarFuFuPR(int gI1, int gI2) const;
   precise_real_type CpVZbarFvFvPL(int gI1, int gI2) const;
   precise_real_type CpVZbarFvFvPR(int , int ) const;
   precise_real_type CpVZVZconjVWpVWp1() const;
   precise_real_type CpVZVZconjVWpVWp2() const;
   precise_real_type CpVZVZconjVWpVWp3() const;
   precise_complex_type CpconjVWpHpAh() const;
   precise_real_type CpconjVWpHphh() const;
   precise_real_type CpconjVWpVPHp() const;
   precise_real_type CpconjVWpVWphh() const;
   precise_real_type CpconjVWpVZHp() const;
   precise_real_type CpconjVWpbargPgWp() const;
   precise_real_type CpconjVWpbargWpCgP() const;
   precise_real_type CpconjVWpbargWpCgZ() const;
   precise_real_type CpconjVWpbargZgWp() const;
   precise_real_type CpVWpconjVWpAhAh() const;
   precise_real_type CpVWpconjVWphhhh() const;
   precise_real_type CpVWpconjVWpconjHpHp() const;
   precise_real_type CpconjVWpVWpVP() const;
   precise_real_type CpconjVWpVZVWp() const;
   precise_complex_type CpconjVWpbarFdFuPL(int gI1, int gI2) const;
   precise_real_type CpconjVWpbarFdFuPR(int , int ) const;
   precise_complex_type CpconjVWpbarFeFvPL(int gI1, int gI2) const;
   precise_real_type CpconjVWpbarFeFvPR(int , int ) const;
   precise_real_type CpVWpconjVWpVPVP1() const;
   precise_real_type CpVWpconjVWpVPVP2() const;
   precise_real_type CpVWpconjVWpVPVP3() const;
   precise_real_type CpVWpconjVWpVZVZ1() const;
   precise_real_type CpVWpconjVWpVZVZ2() const;
   precise_real_type CpVWpconjVWpVZVZ3() const;
   precise_real_type CpVWpconjVWpconjVWpVWp1() const;
   precise_real_type CpVWpconjVWpconjVWpVWp2() const;
   precise_real_type CpVWpconjVWpconjVWpVWp3() const;
   precise_complex_type CpbarUFdFdAhPL(int gO2, int gI1) const;
   precise_complex_type CpbarUFdFdAhPR(int gO1, int gI1) const;
   precise_complex_type CpbarUFdhhFdPL(int gO2, int gI2) const;
   precise_complex_type CpbarUFdhhFdPR(int gO1, int gI2) const;
   precise_complex_type CpbarUFdVGFdPR(int gO2, int gI2) const;
   precise_complex_type CpbarUFdVGFdPL(int gO1, int gI2) const;
   precise_complex_type CpbarUFdVPFdPR(int gO2, int gI2) const;
   precise_complex_type CpbarUFdVPFdPL(int gO1, int gI2) const;
   precise_complex_type CpbarUFdVZFdPR(int gO2, int gI2) const;
   precise_complex_type CpbarUFdVZFdPL(int gO1, int gI2) const;
   precise_complex_type CpbarUFdconjHpFuPL(int gO2, int gI2) const;
   precise_complex_type CpbarUFdconjHpFuPR(int gO1, int gI2) const;
   precise_real_type CpbarUFdconjVWpFuPR(int , int ) const;
   precise_complex_type CpbarUFdconjVWpFuPL(int gO1, int gI2) const;
   precise_complex_type CpbarUFuFuAhPL(int gO2, int gI1) const;
   precise_complex_type CpbarUFuFuAhPR(int gO1, int gI1) const;
   precise_complex_type CpbarUFuhhFuPL(int gO2, int gI2) const;
   precise_complex_type CpbarUFuhhFuPR(int gO1, int gI2) const;
   precise_complex_type CpbarUFuHpFdPL(int gO2, int gI2) const;
   precise_complex_type CpbarUFuHpFdPR(int gO1, int gI2) const;
   precise_complex_type CpbarUFuVGFuPR(int gO2, int gI2) const;
   precise_complex_type CpbarUFuVGFuPL(int gO1, int gI2) const;
   precise_complex_type CpbarUFuVPFuPR(int gO2, int gI2) const;
   precise_complex_type CpbarUFuVPFuPL(int gO1, int gI2) const;
   precise_real_type CpbarUFuVWpFdPR(int , int ) const;
   precise_complex_type CpbarUFuVWpFdPL(int gO1, int gI2) const;
   precise_complex_type CpbarUFuVZFuPR(int gO2, int gI2) const;
   precise_complex_type CpbarUFuVZFuPL(int gO1, int gI2) const;
   precise_complex_type CpbarUFeFeAhPL(int gO2, int gI1) const;
   precise_complex_type CpbarUFeFeAhPR(int gO1, int gI1) const;
   precise_complex_type CpbarUFehhFePL(int gO2, int gI2) const;
   precise_complex_type CpbarUFehhFePR(int gO1, int gI2) const;
   precise_complex_type CpbarUFeVPFePR(int gO2, int gI2) const;
   precise_complex_type CpbarUFeVPFePL(int gO1, int gI2) const;
   precise_complex_type CpbarUFeVZFePR(int gO2, int gI2) const;
   precise_complex_type CpbarUFeVZFePL(int gO1, int gI2) const;
   precise_complex_type CpbarUFeconjHpFvPL(int gO2, int gI2) const;
   precise_real_type CpbarUFeconjHpFvPR(int , int ) const;
   precise_real_type CpbarUFeconjVWpFvPR(int , int ) const;
   precise_real_type CpbarUFeconjVWpFvPL(int gO1, int gI2) const;
   precise_complex_type CpbarFdFdAhPL(int gO2, int gI1) const;
   precise_complex_type CpbarFdFdAhPR(int gO1, int gI1) const;
   precise_complex_type CpbarFdhhFdPL(int gO2, int gI2) const;
   precise_complex_type CpbarFdhhFdPR(int gO1, int gI2) const;
   precise_real_type CpbarFdVZFdPR(int gO2, int gI2) const;
   precise_real_type CpbarFdVZFdPL(int gO1, int gI2) const;
   precise_complex_type CpbarFdconjHpFuPL(int gO2, int gI2) const;
   precise_complex_type CpbarFdconjHpFuPR(int gO1, int gI2) const;
   precise_real_type CpbarFdconjVWpFuPR(int , int ) const;
   precise_complex_type CpbarFdconjVWpFuPL(int gO1, int gI2) const;
   precise_complex_type CpbarFeFeAhPL(int gO2, int gI1) const;
   precise_complex_type CpbarFeFeAhPR(int gO1, int gI1) const;
   precise_complex_type CpbarFehhFePL(int gO2, int gI2) const;
   precise_complex_type CpbarFehhFePR(int gO1, int gI2) const;
   precise_real_type CpbarFeVZFePR(int gO2, int gI2) const;
   precise_real_type CpbarFeVZFePL(int gO1, int gI2) const;
   precise_complex_type CpbarFeconjHpFvPL(int gO2, int gI2) const;
   precise_real_type CpbarFeconjHpFvPR(int , int ) const;
   precise_real_type CpbarFeconjVWpFvPR(int , int ) const;
   precise_complex_type CpbarFeconjVWpFvPL(int gO1, int gI2) const;
   precise_complex_type CpbarFuFuAhPL(int gO2, int gI1) const;
   precise_complex_type CpbarFuFuAhPR(int gO1, int gI1) const;
   precise_complex_type CpbarFuhhFuPL(int gO2, int gI2) const;
   precise_complex_type CpbarFuhhFuPR(int gO1, int gI2) const;
   precise_complex_type CpbarFuHpFdPL(int gO2, int gI2) const;
   precise_complex_type CpbarFuHpFdPR(int gO1, int gI2) const;
   precise_real_type CpbarFuVPFuPR(int gO2, int gI2) const;
   precise_real_type CpbarFuVPFuPL(int gO1, int gI2) const;
   precise_real_type CpbarFuVWpFdPR(int , int ) const;
   precise_complex_type CpbarFuVWpFdPL(int gO1, int gI2) const;
   precise_real_type CpbarFuVZFuPR(int gO2, int gI2) const;
   precise_real_type CpbarFuVZFuPL(int gO1, int gI2) const;
   precise_complex_type self_energy_Hp_1loop(precise_real_type p ) const;
   precise_complex_type self_energy_Ah_1loop(precise_real_type p ) const;
   precise_complex_type self_energy_hh_1loop(precise_real_type p ) const;
   precise_complex_type self_energy_VZ_1loop(precise_real_type p ) const;
   precise_complex_type self_energy_VWp_1loop(precise_real_type p ) const;
   precise_complex_type self_energy_Fd_1loop_1(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fd_1loop_PR(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fd_1loop_PL(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fu_1loop_1(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fu_1loop_PR(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fu_1loop_PL(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fe_1loop_1(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fe_1loop_PR(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fe_1loop_PL(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fd_1loop_1_heavy_rotated(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fd_1loop_PR_heavy_rotated(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fd_1loop_PL_heavy_rotated(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fe_1loop_1_heavy_rotated(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fe_1loop_PR_heavy_rotated(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fe_1loop_PL_heavy_rotated(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fu_1loop_1_heavy_rotated(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fu_1loop_PR_heavy_rotated(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fu_1loop_PL_heavy_rotated(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fu_1loop_1_heavy(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fu_1loop_PR_heavy(precise_real_type p , int gO1, int gO2) const;
   precise_complex_type self_energy_Fu_1loop_PL_heavy(precise_real_type p , int gO1, int gO2) const;

   Eigen::Matrix<precise_complex_type,3,3> self_energy_Fd_1loop_1(precise_real_type p) const;
   Eigen::Matrix<precise_complex_type,3,3> self_energy_Fd_1loop_PR(precise_real_type p) const;
   Eigen::Matrix<precise_complex_type,3,3> self_energy_Fd_1loop_PL(precise_real_type p) const;
   Eigen::Matrix<precise_complex_type,3,3> self_energy_Fu_1loop_1(precise_real_type p) const;
   Eigen::Matrix<precise_complex_type,3,3> self_energy_Fu_1loop_PR(precise_real_type p) const;
   Eigen::Matrix<precise_complex_type,3,3> self_energy_Fu_1loop_PL(precise_real_type p) const;
   Eigen::Matrix<precise_complex_type,3,3> self_energy_Fe_1loop_1(precise_real_type p) const;
   Eigen::Matrix<precise_complex_type,3,3> self_energy_Fe_1loop_PR(precise_real_type p) const;
   Eigen::Matrix<precise_complex_type,3,3> self_energy_Fe_1loop_PL(precise_real_type p) const;

   precise_complex_type tadpole_hh_1loop() const;

   /// calculates the tadpoles at current loop order
   Eigen::Matrix<precise_real_type, number_of_ewsb_equations, 1> tadpole_equations() const;

   /// calculates Higgs 2-loop self-energy
   precise_real_type self_energy_hh_2loop(precise_real_type p) const;
   /// calculates Higgs 3-loop self-energy
   precise_real_type self_energy_hh_3loop() const;
   /// calculates Higgs 4-loop self-energy
   precise_real_type self_energy_hh_4loop() const;

   void calculate_MVG_pole();
   void calculate_MFv_pole();
   void calculate_Mhh_pole();
   void calculate_MVP_pole();
   void calculate_MVZ_pole();
   void calculate_MFd_pole();
   void calculate_MFu_pole();
   void calculate_MFe_pole();
   void calculate_MVWp_pole();
   precise_real_type calculate_MVWp_pole(precise_real_type);
   precise_real_type calculate_MVZ_pole(precise_real_type);

   precise_real_type calculate_MFv_DRbar(precise_real_type, int) const;
   precise_real_type calculate_MFe_DRbar(precise_real_type, int) const;
   precise_real_type calculate_MFu_DRbar(precise_real_type, int) const;
   precise_real_type calculate_MFd_DRbar(precise_real_type, int) const;
   precise_real_type calculate_MVP_DRbar(precise_real_type);
   precise_real_type calculate_MVZ_DRbar(precise_real_type);
   precise_real_type calculate_MVWp_DRbar(precise_real_type);
   precise_real_type calculate_Mhh_DRbar(precise_real_type);

   precise_real_type ThetaW() const;

   precise_real_type calculate_delta_alpha_em(precise_real_type alphaEm) const;
   precise_real_type calculate_delta_alpha_s(precise_real_type alphaS) const;
   void calculate_Lambdax_DRbar();
   precise_real_type calculate_theta_w(const softsusy::QedQcd&, precise_real_type alpha_em_drbar);
   void calculate_Yu_DRbar(const softsusy::QedQcd&);
   void calculate_Yd_DRbar(const softsusy::QedQcd&);
   void calculate_Ye_DRbar(const softsusy::QedQcd&);
   precise_real_type recalculate_mw_pole(precise_real_type);
   precise_real_type max_rel_diff(const Standard_model& old) const;

protected:

   // Running parameters
   precise_real_type g1{};
   precise_real_type g2{};
   precise_real_type g3{};
   precise_real_type Lambdax{};
   Eigen::Matrix<precise_real_type,3,3> Yu{Eigen::Matrix<precise_real_type,3,3>::Zero()};
   Eigen::Matrix<precise_real_type,3,3> Yd{Eigen::Matrix<precise_real_type,3,3>::Zero()};
   Eigen::Matrix<precise_real_type,3,3> Ye{Eigen::Matrix<precise_real_type,3,3>::Zero()};
   precise_real_type mu2{};
   precise_real_type v{};

private:

   static const int numberOfParameters = 33;

   struct Beta_traces {
      precise_real_type traceYdAdjYd{};
      precise_real_type traceYeAdjYe{};
      precise_real_type traceYuAdjYu{};
      precise_real_type traceYdAdjYdYdAdjYd{};
      precise_real_type traceYeAdjYeYeAdjYe{};
      precise_real_type traceYuAdjYuYuAdjYu{};
      precise_real_type traceYdAdjYuYuAdjYd{};
      precise_real_type traceYdAdjYdYdAdjYdYdAdjYd{};
      precise_real_type traceYdAdjYdYdAdjYuYuAdjYd{};
      precise_real_type traceYdAdjYuYuAdjYdYdAdjYd{};
      precise_real_type traceYdAdjYuYuAdjYuYuAdjYd{};
      precise_real_type traceYeAdjYeYeAdjYeYeAdjYe{};
      precise_real_type traceYuAdjYuYuAdjYuYuAdjYu{};
   };
   void calc_beta_traces(Beta_traces&) const;

   precise_real_type calc_beta_g1_one_loop(const Beta_traces&) const;
   precise_real_type calc_beta_g1_two_loop(const Beta_traces&) const;
   precise_real_type calc_beta_g1_three_loop(const Beta_traces&) const;
   precise_real_type calc_beta_g2_one_loop(const Beta_traces&) const;
   precise_real_type calc_beta_g2_two_loop(const Beta_traces&) const;
   precise_real_type calc_beta_g2_three_loop(const Beta_traces&) const;
   precise_real_type calc_beta_g3_one_loop(const Beta_traces&) const;
   precise_real_type calc_beta_g3_two_loop(const Beta_traces&) const;
   precise_real_type calc_beta_g3_three_loop(const Beta_traces&) const;
   precise_real_type calc_beta_g3_four_loop(const Beta_traces&) const;
   precise_real_type calc_beta_g3_five_loop(const Beta_traces&) const;
   precise_real_type calc_beta_Lambdax_one_loop(const Beta_traces&) const;
   precise_real_type calc_beta_Lambdax_two_loop(const Beta_traces&) const;
   precise_real_type calc_beta_Lambdax_three_loop(const Beta_traces&) const;
   precise_real_type calc_beta_Lambdax_four_loop(const Beta_traces&) const;
   Eigen::Matrix<precise_real_type,3,3> calc_beta_Yu_one_loop(const Beta_traces&) const;
   Eigen::Matrix<precise_real_type,3,3> calc_beta_Yu_two_loop(const Beta_traces&) const;
   Eigen::Matrix<precise_real_type,3,3> calc_beta_Yu_three_loop(const Beta_traces&) const;
   Eigen::Matrix<precise_real_type,3,3> calc_beta_Yu_four_loop(const Beta_traces&) const;
   Eigen::Matrix<precise_real_type,3,3> calc_beta_Yd_one_loop(const Beta_traces&) const;
   Eigen::Matrix<precise_real_type,3,3> calc_beta_Yd_two_loop(const Beta_traces&) const;
   Eigen::Matrix<precise_real_type,3,3> calc_beta_Yd_three_loop(const Beta_traces&) const;
   Eigen::Matrix<precise_real_type,3,3> calc_beta_Ye_one_loop(const Beta_traces&) const;
   Eigen::Matrix<precise_real_type,3,3> calc_beta_Ye_two_loop(const Beta_traces&) const;
   Eigen::Matrix<precise_real_type,3,3> calc_beta_Ye_three_loop(const Beta_traces&) const;
   precise_real_type calc_beta_mu2_one_loop(const Beta_traces&) const;
   precise_real_type calc_beta_mu2_two_loop(const Beta_traces&) const;
   precise_real_type calc_beta_mu2_three_loop(const Beta_traces&) const;
   precise_real_type calc_beta_v_one_loop(const Beta_traces&) const;
   precise_real_type calc_beta_v_two_loop(const Beta_traces&) const;
   precise_real_type calc_beta_v_three_loop(const Beta_traces&) const;

   using EWSB_vector_t = Eigen::Matrix<precise_real_type,number_of_ewsb_equations,1>;

   class EEWSBStepFailed : public Error {
   public:
      EEWSBStepFailed() : Error("Could not perform EWSB step") {}
      virtual ~EEWSBStepFailed() = default;
   };

   int ewsb_loop_order{4};
   int pole_mass_loop_order{4};
   bool force_output{false};      ///< switch to force output of pole masses
   precise_real_type precision{1e-4};        ///< RG running precision
   precise_real_type ewsb_iteration_precision{1e-5};
   Standard_model_physical physical{}; ///< contains the pole masses and mixings
   Problems problems{standard_model_info::model_name,
                     &standard_model_info::particle_names_getter,
                     &standard_model_info::parameter_names_getter};
   Loop_corrections loop_corrections{}; ///< used loop pole mass corrections
   Threshold_corrections threshold_corrections{}; ///< used low-energy threshold corrections
   Physical_input input{};

   int get_number_of_ewsb_iterations() const;
   int get_number_of_mass_iterations() const;
   int solve_ewsb_iteratively();
   int solve_ewsb_iteratively(int);
   int solve_ewsb_iteratively_with(EWSB_solver*, const Eigen::Matrix<precise_real_type, number_of_ewsb_equations, 1>&);
   int solve_ewsb_tree_level_custom();
   EWSB_vector_t ewsb_initial_guess();
   EWSB_vector_t ewsb_step() const;
   void copy_DRbar_masses_to_pole_masses();

   void initial_guess_for_parameters(const softsusy::QedQcd&);
   bool check_convergence(const Standard_model& old) const;

   // Passarino-Veltman loop functions
   precise_real_type A0(precise_real_type) const;
   precise_real_type B0(precise_real_type, precise_real_type, precise_real_type) const;
   precise_real_type B1(precise_real_type, precise_real_type, precise_real_type) const;
   precise_real_type B00(precise_real_type, precise_real_type, precise_real_type) const;
   precise_real_type B22(precise_real_type, precise_real_type, precise_real_type) const;
   precise_real_type H0(precise_real_type, precise_real_type, precise_real_type) const;
   precise_real_type F0(precise_real_type, precise_real_type, precise_real_type) const;
   precise_real_type G0(precise_real_type, precise_real_type, precise_real_type) const;

   // DR-bar masses
   precise_real_type MVG{};
   precise_real_type MHp{};
   Eigen::Array<precise_real_type,3,1> MFv{Eigen::Array<precise_real_type,3,1>::Zero()};
   precise_real_type MAh{};
   precise_real_type Mhh{};
   precise_real_type MVP{};
   precise_real_type MVZ{};
   Eigen::Array<precise_real_type,3,1> MFd{Eigen::Array<precise_real_type,3,1>::Zero()};
   Eigen::Array<precise_real_type,3,1> MFu{Eigen::Array<precise_real_type,3,1>::Zero()};
   Eigen::Array<precise_real_type,3,1> MFe{Eigen::Array<precise_real_type,3,1>::Zero()};
   precise_real_type MVWp{};
   Eigen::Array<precise_real_type,2,1> MVPVZ{Eigen::Array<precise_real_type,2,1>::Zero()};

   // DR-bar mixing matrices
   Eigen::Matrix<precise_complex_type,3,3> Vd{Eigen::Matrix<precise_complex_type,3,3>::Zero()};
   Eigen::Matrix<precise_complex_type,3,3> Ud{Eigen::Matrix<precise_complex_type,3,3>::Zero()};
   Eigen::Matrix<precise_complex_type,3,3> Vu{Eigen::Matrix<precise_complex_type,3,3>::Zero()};
   Eigen::Matrix<precise_complex_type,3,3> Uu{Eigen::Matrix<precise_complex_type,3,3>::Zero()};
   Eigen::Matrix<precise_complex_type,3,3> Ve{Eigen::Matrix<precise_complex_type,3,3>::Zero()};
   Eigen::Matrix<precise_complex_type,3,3> Ue{Eigen::Matrix<precise_complex_type,3,3>::Zero()};
   Eigen::Matrix<precise_real_type,2,2> ZZ{Eigen::Matrix<precise_real_type,2,2>::Zero()};


};

std::ostream& operator<<(std::ostream&, const Standard_model&);

} // namespace standard_model

} // namespace flexiblesusy

#endif
