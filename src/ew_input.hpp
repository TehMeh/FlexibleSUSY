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

#ifndef ELECTROWEAK_INPUT_H
#define ELECTROWEAK_INPUT_H

#include <cmath>

namespace flexiblesusy {

/**
 * @namespace Electroweak_constants
 *
 * Contains Standard Model parameters at the electroweak scale
 */
namespace Electroweak_constants {
   namespace {
      const precise_real_type vev = 246.22;
      const precise_real_type root2 = sqrt(2.0);
      const precise_real_type mtoprun = 165;
      const precise_real_type mbrun = 2.9;
      const precise_real_type mtau = 1.77699;
      const precise_real_type yt = mtoprun * root2 / vev;
      const precise_real_type yb = mbrun * root2 / vev;
      const precise_real_type ytau = mtau * root2 / vev;
      const precise_real_type MZ = 91.1876;
      const precise_real_type Error_MZ = 0.0021; ///< uncertainty on MZ from PDG
      const precise_real_type MW = 80.385;
      const precise_real_type MH = 125.09; ///< Higgs mass from PDG (CMS and ATLAS combination)
      const precise_real_type Error_MH = 0.24; ///< uncertainty on MH from PDG - 0.11 (sys) and 0.21 stat combined in quadrature.
      const precise_real_type MUP = 2.4e-3; ///< default running quark mass from PDG
      const precise_real_type MDOWN = 4.75e-3; ///< default running quark mass from PDG
      const precise_real_type MSTRANGE = 0.104; ///< default running quark mass from PDG
      const precise_real_type MCHARM = 1.27; ///< default running quark mass from PDG
      const precise_real_type MBOTTOM = 4.18; ///< default running quark mass from PDG
      const precise_real_type MTOP = 165.0; ///< default running quark mass from PDG
      const precise_real_type MELECTRON = 5.10998902e-4; ///< default pole lepton mass from PDG
      const precise_real_type MMUON = 1.056583715e-1; ///< default pole lepton mass from PDG
      const precise_real_type MTAU = 1.77699; ///< default pole lepton mass from PDG
      const precise_real_type PMTOP = 173.34; ///< default pole mass from CDF/D0 Run II 1207.1069
      const precise_real_type PMBOTTOM = 4.9; ///< default pole mass from PDG
      const precise_real_type aem = 1.0 / 127.916; // at MZ
      const precise_real_type sinThetaW2 = 0.23122;
      const precise_real_type sinThetaW = sqrt(sinThetaW2);
      const precise_real_type cosThetaW2 = 1 - sinThetaW2;
      const precise_real_type cosThetaW = sqrt(cosThetaW2);
      const precise_real_type alpha1 = 5.0 * aem / (3.0 * (1.0 - sinThetaW2));
      const precise_real_type alpha2 = aem / sinThetaW2;
      const precise_real_type alpha3 = 0.1184; // at MZ from PDG
      const precise_real_type e  = sqrt(4.0 * M_PI * aem);
      const precise_real_type g1 = sqrt(4.0 * M_PI * alpha1);
      const precise_real_type g2 = sqrt(4.0 * M_PI * alpha2);
      const precise_real_type g3 = sqrt(4.0 * M_PI * alpha3);
      const precise_real_type gYSM = 3.57232027E-01;     ///< gY MS-bar in the SM at Q = MZ
      const precise_real_type g1SM = sqrt(5./3.) * gYSM; ///< g1 MS-bar in the SM at Q = MZ
      const precise_real_type g2SM = 6.51103848E-01;     ///< g2 MS-bar in the SM at Q = MZ
      const precise_real_type g3SM = 1.21087245E+00;     ///< g3 MS-bar in the SM at Q = MZ
      const precise_real_type CKM_THETA12 = 0.229206; ///< From Vus/Vud in global CKM fit, PDG
      const precise_real_type CKM_THETA13 = 0.003960; ///< From Vub in global CKM fit, PDG
      const precise_real_type CKM_THETA23 = 0.042223; ///< From Vcb/Vtb in global CKM fit, PDG
      const precise_real_type CKM_DELTA   = 0.;
      const precise_real_type PMNS_THETA12 = 0.5 * asin(sqrt(0.846));
      const precise_real_type PMNS_THETA13 = 0.5 * asin(sqrt(0.093));
      const precise_real_type PMNS_THETA23 = 0.5 * asin(sqrt(0.999));
      const precise_real_type PMNS_DELTA   = 0.;
      const precise_real_type PMNS_ALPHA1  = 0.;
      const precise_real_type PMNS_ALPHA2  = 0.;
      const precise_real_type gfermi = 1.1663787e-5; ///< Fermi constant G_F
      const precise_real_type yeSM = 2.85784368E-06; ///< Ye(1,1) MS-bar in the SM at Q = MZ
      const precise_real_type ymSM = 5.90911374E-04; ///< Ye(2,2) MS-bar in the SM at Q = MZ
      const precise_real_type ylSM = 9.95869693E-03; ///< Ye(3,3) MS-bar in the SM at Q = MZ
      const precise_real_type yuSM = 7.89527379E-06; ///< Yu(1,1) MS-bar in the SM at Q = MZ
      const precise_real_type ycSM = 3.60854291E-03; ///< Yu(2,2) MS-bar in the SM at Q = MZ
      const precise_real_type ytSM = 9.76017610E-01; ///< Yu(3,3) MS-bar in the SM at Q = MZ
      const precise_real_type ydSM = 1.56989573E-05; ///< Yd(1,1) MS-bar in the SM at Q = MZ
      const precise_real_type ysSM = 3.43724539E-04; ///< Yd(2,2) MS-bar in the SM at Q = MZ
      const precise_real_type ybSM = 1.64406299E-02; ///< Yd(3,3) MS-bar in the SM at Q = MZ
      const precise_real_type mu2SM = 7.67488232E+03; ///< mu^2 MS-bar in the SM at Q = MZ
      const precise_real_type lamSM = 2.79613357E-01; ///< lambda MS-bar in the SM at Q = MZ
      const precise_real_type vSM = 2.48997424E+02;   ///< VEV MS-bar in the SM at Q = MZ
   } // namespace
} // namespace Electroweak_constants

} // namespace flexiblesusy

#endif
