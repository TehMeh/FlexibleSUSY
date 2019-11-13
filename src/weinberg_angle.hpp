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

#ifndef WEINBERG_ANGLE_H
#define WEINBERG_ANGLE_H

#include <Eigen/Core>

namespace flexiblesusy {

namespace weinberg_angle {

/**
 * @class Weinberg_angle
 * @brief Class to calculate the DR-bar weak mixing angle
 */
class Weinberg_angle {
public:
   /**
    * @class Data
    * @brief Model parameters necessary for calculating weak mixing angle
    *
    * @attention The W and Z self-energies are assumed to be
    * calculated using the top quark pole mass, instead of the top
    * quark DR-bar mass.
    */
   struct Data {
      Data();

      precise_real_type scale;                  ///< renormalization scale
      precise_real_type alpha_em_drbar;         ///< alpha_em(MZ, DR-bar, SUSY)
      precise_real_type fermi_contant;          ///< Fermi constant
      precise_real_type self_energy_z_at_mz;    ///< self-energy Z at p = MZ, mt = mt_pole
      precise_real_type self_energy_w_at_0;     ///< self-energy W at p = 0, mt = mt_pole
      precise_real_type self_energy_w_at_mw;    ///< self-energy W at p = MW, mt = mt_pole
      precise_real_type mw_pole;                ///< W pole mass
      precise_real_type mz_pole;                ///< Z pole mass
      precise_real_type mt_pole;                ///< top quark pole mass
      precise_real_type mh_drbar;               ///< lightest CP-even Higgs DR-bar mass
      precise_real_type hmix_12;                ///< CP-even Higgs mixing Cos(alpha)
      precise_real_type msel_drbar;             ///< left-handed selectron DR-bar mass
      precise_real_type msmul_drbar;            ///< left-handed smuon DR-bar mass
      precise_real_type msve_drbar;             ///< electron-sneutrino DR-bar mass
      precise_real_type msvm_drbar;             ///< muon-sneutrino DR-bar mass
      Eigen::ArrayXdp mn_drbar;       ///< Neutralino DR-bar mass
      Eigen::ArrayXdp mc_drbar;       ///< Chargino DR-bar mass
      Eigen::MatrixXcdp zn;           ///< Neutralino mixing matrix
      Eigen::MatrixXcdp um;           ///< Chargino mixing matrix
      Eigen::MatrixXcdp up;           ///< Chargino mixing matrix
      precise_real_type gY;                     ///< U(1)_Y gauge coupling
      precise_real_type g2;                     ///< SU(2)_L gauge coupling
      precise_real_type g3;                     ///< SU(3)_c gauge coupling
      precise_real_type tan_beta;               ///< tan(beta) = vu / vd
   };

   struct Self_energy_data {
      Self_energy_data();
      precise_real_type scale;                  ///< renormalization scale
      precise_real_type mt_pole;                ///< top quark pole mass
      precise_real_type mt_drbar;               ///< top quark DR-bar mass
      precise_real_type mb_drbar;               ///< bottom quark DR-bar mass
      precise_real_type gY;                     ///< U(1)_Y gauge coupling
      precise_real_type g2;                     ///< SU(2)_L gauge coupling
   };

   Weinberg_angle();

   void enable_susy_contributions(); ///< enable susy contributions
   void disable_susy_contributions(); ///< disable susy contributions

   void set_data(const Data&);       ///< set data necessary for the calculation
   void set_number_of_iterations(int); ///< maximum number of iterations
   void set_number_of_loops(int);    ///< set number of loops
   void set_precision_goal(precise_real_type);  ///< set precision goal
   precise_real_type get_rho_hat() const;       ///< returns the rho parameter
   precise_real_type get_sin_theta() const;     ///< returns sin(theta_w)

   /// calculates the sinus of the Weinberg angle
   int calculate(precise_real_type rho_start = 1.0, precise_real_type sin_start = 0.48);

   static precise_real_type replace_mtop_in_self_energy_z(precise_real_type, precise_real_type, const Self_energy_data&);
   static precise_real_type replace_mtop_in_self_energy_w(precise_real_type, precise_real_type, const Self_energy_data&);

private:
   int number_of_iterations; ///< maximum number of iterations
   int number_of_loops;      ///< number of loops
   precise_real_type precision_goal;         ///< precision goal
   precise_real_type rho_hat;                ///< output rho-hat parameter
   precise_real_type sin_theta;              ///< output sin(theta)
   Data data;
   bool susy_contributions;       ///< model type

   static precise_real_type calculate_delta_r(precise_real_type, precise_real_type, const Data&, bool add_susy_contributions = true, int number_of_loops = 2);
   static precise_real_type calculate_delta_rho(precise_real_type, precise_real_type, const Data&, bool add_susy_contributions = true, int number_of_loops = 2);
   static precise_real_type calculate_delta_vb(precise_real_type, precise_real_type, const Data&, bool add_susy_contributions = true);
   static precise_real_type calculate_delta_vb_sm(precise_real_type, precise_real_type, const Data&);
   static precise_real_type calculate_delta_vb_susy(precise_real_type, const Data&);
   static precise_real_type rho_2(precise_real_type);

   static precise_real_type calculate_self_energy_z_top(precise_real_type, precise_real_type, const Self_energy_data&);
   static precise_real_type calculate_self_energy_w_top(precise_real_type, precise_real_type, const Self_energy_data&);
};

} // namespace weinberg_angle

} // namespace flexiblesusy

#endif
