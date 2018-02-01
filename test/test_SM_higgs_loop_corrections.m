(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

Needs["TestSuite`", "TestSuite.m"];

t = (gt v / Sqrt[2])^2;
q = Q^2;
y = gt^2;
g = g3^2
L = Log[t/q];

DMh20L = \[Lambda] v^2;
DMh21L = -12 y t L;
DMh22L = 32 y g t (3 L^2 + L);
DMh23L = g^2 y t (248.1215180432007 + 839.1966169377614 L + 160 L^2 - 736 L^3);
DMh24L = 22080*g^3 y t (4.33441 + 1.16581 L - 3.74561 L^2 - 0.626087 L^3 + L^4);

Mh24L = DMh20L + h DMh21L + h^2 DMh22L + h^3 DMh23L + h^4 DMh24L;

betaRules = {
    g1 -> 0,
    g2 -> 0,
    g3 -> g3[Q],
    g\[Tau] -> 0,
    gb -> 0,
    gt -> gt[Q],
    \[Lambda] -> \[Lambda][Q],
    v -> v[Q],
    m2 -> m2[Q]
};

simp = {
    g3[Q] -> g3,
    gt[Q] -> gt,
    \[Lambda][Q] -> 0,
    v[Q] -> v,
    m2[Q] -> m2,
    v -> Sqrt[2] mt / gt
};

ass = { mt > 0, Q > 0, g3 > 0, gt > 0, v > 0 };

(* beta functions *)
SMdir = FileNameJoin[{Directory[], "meta", "SM"}];
betag3 = Get[FileNameJoin[{SMdir, "beta_g3.m"}]] /. betaRules;
betagt = Get[FileNameJoin[{SMdir, "beta_gt.m"}]] /. betaRules;
betala = Get[FileNameJoin[{SMdir, "beta_lambda.m"}]] /. betaRules;
betav  = Get[FileNameJoin[{SMdir, "beta_v.m"}]] /. betaRules;
betam2 = Get[FileNameJoin[{SMdir, "beta_m2.m"}]] /. betaRules;

(* define derivative of g3 *)
Derivative[1][g3][Q] = (
    h^1 betag3[[1]] +
    h^2 betag3[[2]] +
    h^3 betag3[[3]]
)/Q;

(* define derivative of gt *)
Derivative[1][gt][Q] = (
    h^1 betagt[[1]] +
    h^2 betagt[[2]] +
    h^3 betagt[[3]]
)/Q;

(* define derivative of lambda *)
Derivative[1][\[Lambda]][Q] = (
    h^1 betala[[1]] +
    h^2 betala[[2]] +
    h^3 betala[[3]] +
    h^4 betala[[4]]
)/Q;

(* define derivative of VEV *)
Derivative[1][v][Q] = (
    h^1 betav[[1]] +
    h^2 betav[[2]]
)/Q;

(* define derivative of m^2 *)
Derivative[1][m2][Q] = (
    h^1 betam2[[1]] +
    h^2 betam2[[2]]
)/Q;

deriv = D[Mh24L /. betaRules, Q] //. simp;

deriv = Collect[
    Normal[Series[deriv, {h, 0, 4}]],
    {h, gt}, Simplify[#, ass]&
];

(* ignore O(at^n) terms with n > 1 *)
deriv = deriv /. { gt^n_ :> 0 /; n > 2 };

Print["Testing 0L renormalization group invariance of Mh^2 in the SM ..."];
TestEquality[Coefficient[deriv, h, 0], 0];
Print["Testing 1L renormalization group invariance of Mh^2 in the SM ..."];
TestEquality[Coefficient[deriv, h, 1], 0];
Print["Testing 2L renormalization group invariance of Mh^2 in the SM ..."];
TestEquality[Coefficient[deriv, h, 2], 0];
Print["Testing 3L renormalization group invariance of Mh^2 in the SM ..."];
TestEquality[Coefficient[deriv, h, 3], 0];
Print["Testing 4L renormalization group invariance of Mh^2 in the SM ..."];
TestEquality[Coefficient[deriv, h, 4], 0];

PrintTestSummary[];