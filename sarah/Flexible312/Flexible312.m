(* ::Package:: *)

Off[General::spell]

Model`Name = "Flexible312";
Model`NameLaTeX ="Grimus-Neufeld+2HDM";



(* Gauge superfields*)

Gauge[[1]]={B, U[1], hypercharge, g1, False};
Gauge[[2]]={WB, SU[2], left, g2, True};
Gauge[[3]]={G, SU[3], color, g3, False};



(*Chiral Superfields*)

FermionFields[[1]] = {q, 3, {uL, dL}, 1/6, 2, 3};
FermionFields[[2]] = {l, 3, {vL, eL}, -1/2, 2, 1};
FermionFields[[3]] = {d, 3, conj[dR], 1/3, 1, -3};
FermionFields[[4]] = {u, 3, conj[uR], -2/3, 1, -3};
FermionFields[[5]] = {e, 3, conj[eR], 1, 1, 1};
FermionFields[[6]] = {N0, 1, conj[vR], 0, 1, 1};

ScalarFields[[1]] = {H1, 1, {H1p, H10}, 1/2, 2, 1};
ScalarFields[[2]] = {H2, 1, {H2p, H20}, 1/2, 2, 1};


(*DEFINITIONS*)

NameOfStates={GaugeES, EWSB};
(* Before EWSB *)
DEFINITION[GaugeES][Additional] = {
	{LagHC, {AddHC->True}},
	{LagNoHC, {AddHC->False}}
};


LagNoHC = -(M112 conj[H1].H1 + M222 conj[H2].H2 + Z1111 conj[H1].H1.conj[H1].H1+ \
			Z2222 conj[H2].H2.conj[H2].H2 + Z2211 conj[H2].H2.conj[H1].H1 + \
			Z2112 conj[H2].H1.conj[H1].H2);





LagHC =-( M12 conj[H1].H2 + Z2121 conj[H2].H1.conj[H2].H1 +\
		 Ye1 conj[H1].e.l + Ye2 conj[H2].e.l +\
		Z1112 conj[H1].H1.conj[H1].H2 + Z2212 conj[H2].H2.conj[H1].H2 + \
		Yv1 H1.N0.l + Yv2 H2.N0.l + Yd1 conj[H1].d.q + Yd2 conj[H2].d.q +\
		Yu1 H1.u.q + Yu2 H2.u.q + MR/2 N0.N0);
		
(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] = 
{
	{{VB, VWB[3]}, {VP, VZ}, ZZ}, (* B-boson, W3 boson, photon, Z boson, ZZ-mixing matrix*)
	{{VWB[1], VWB[2]}, {VWm, conj[VWm]}, ZW} (* W1 boson, W2 boson,W^- -boson, W^-* -boson, ZW-mixing matrix *)
};

 DEFINITION[GaugeES][Phases]= 
{    {H2p, Exp[I eta]},
     {H20, Exp[I eta]}
   };
(* VEVs *)

DEFINITION[EWSB][VEVs]=
{
	{H10, {v, 1/Sqrt[2]}, {sigma1, I/Sqrt[2]}, {phi1, 1/Sqrt[2]}},
	{H20, {0,0}, {sigma2, I/Sqrt[2]}, {phi2, 1/Sqrt[2]}}
};

DEFINITION[EWSB][MatterSector]=
{  
      {{phi1, phi2, sigma1, sigma2}, {hh, ZH}},
      {{conj[H1p],conj[H2p]},{Hm,ZP}},
	  {{{dL}, {conj[dR]}}, {{DL,Vd}, {DR,Ud}}},
      {{{uL}, {conj[uR]}}, {{UL,Vu}, {UR,Uu}}},
      {{{eL}, {conj[eR]}}, {{EL,Ve}, {ER,Ue}}},
   
      {{vL, conj[vR]}, {Ksim, Vpmns}}
};


(* Dirac spinors *)

DEFINITION[EWSB][DiracSpinors]=
{
	Fd->{DL, conj[DR]},
	Fe->{EL, conj[ER]},
	Fu->{UL, conj[UR]},
	Fv->{Ksim, conj[Ksim]}
};

DEFINITION[GaugeES][DiracSpinors]=
{  
	Fd1->{dL, 0},
	Fd2->{0, dR},
	Fu1->{uL, 0},
	Fu2->{0, uR},
	Fe1->{eL, 0},
	Fe2->{0, eR},
	Fv1->{vL, 0},
	Fv2->{0, vR}
};


