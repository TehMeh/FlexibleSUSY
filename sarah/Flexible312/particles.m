(* ::Package:: *)

ParticleDefinitions[GaugeES] = { 

      {H10,  {    PDG -> {0},
                 Width -> 0, 
                 Mass -> Automatic,
                 FeynArtsNr -> 1,
                 LaTeX -> "H^0",
                 OutputName -> "H0" }},                         
      
      
      {H1p,  {             PDG -> {0},
                 Width -> 0, 
                 Mass -> Automatic,
                 FeynArtsNr -> 2,
                 LaTeX -> "H^+",
                 OutputName -> "Hp" }},  
                 
    {H20, {LaTeX -> "H_2^0",
           OutputName -> "H20",
           FeynArtsNr -> 301}},
   {H2p, {LaTeX -> "H_2^+",
          OutputName -> "H2p",
          FeynArtsNr -> 302}},
                   
                    
      {VB,   { Description -> "B-Boson"}},                                                   
      {VG,   { Description -> "Gluon"}},          
      {VWB,  { Description -> "W-Bosons"}},
	
      {gB,   { Description -> "B-Boson Ghost"}},                                                   
      {gG,   { Description -> "Gluon Ghost" }},          
      {gWB,  { Description -> "W-Boson Ghost"}},
      
      {Fd1,   { Description -> "Left-Down-Quarks",
                     OutputName -> "Fd1",
                     FeynArtsNr->401}},   
      {Fd2,   { Description -> "Right-Down-Quarks",
                     OutputName -> "Fd2",
                     FeynArtsNr->402}},
      {Fe1,   { Description -> "Left-Electrons",
                     OutputName -> "Fe1",
                     FeynArtsNr->403}},   
      {Fe2,   { Description -> "Right-Electrons",
                     OutputName -> "Fe2",
                     FeynArtsNr->404}}, 
      {Fu1,   { Description -> "Left-Up-Quarks",
                     OutputName -> "Fu1",
                     FeynArtsNr->405}},   
      {Fu2,   { Description -> "Right-Up-Quarks",
                     OutputName -> "Fu2",
                     FeynArtsNr->406}}, 
      {Fv1,   { Description -> "Left-Neutrinos",
                     OutputName -> "Fv1",
                     FeynArtsNr->407}},   
      {Fv2,   { Description -> "Right-Neutrinos",
                     OutputName -> "Fv2",
                     FeynArtsNr->408}}
      
};


 ParticleDefinitions[EWSB] = {
   
(*      {hh ,  { Description -> "Higgs"}}, 
      {Ah ,  { Description -> "Pseudo-Scalar Higgs"}}, *)
      
      {hh ,  { Description -> "Higgs",
               PDG -> {0,25, 35,36},
               PDG.IX->{0,100000001,100000002,100000003} }}, 
               
      {Hm,  { Description -> "Charged Higgs"}},
      
      {VP,   { Description -> "Photon"}}, 
      {VZ,   { Description -> "Z-Boson", Goldstone -> hh[{1}]}}, 
      {VG,   { Description -> "Gluon" }},          
      {VWm,  { Description -> "W-Boson",
               Goldstone -> Hm[{1}] }},
     
      {gP,   { Description -> "Photon Ghost"}},                                                   
      {gWm,  { Description -> "Negative W-Boson Ghost"}}, 
      {gWmC, { Description -> "Positive W-Boson Ghost" }}, 
      {gZ,   { Description -> "Z-Boson Ghost" }},
      {gG,   { Description -> "Gluon Ghost" }}, 
   
   
      {Fd,   { Description -> "Down-Quarks"}},   
      {Fu,   { Description -> "Up-Quarks"}},   
      {Fe,   { Description -> "Leptons" }},
      {Fv,   { Description -> "Neutrinos", 
               PDG->{12,14,16,18}}}           
 };
 
 


 WeylFermionAndIndermediate =
 {
    {H,      {   PDG -> {0},
                 Width -> 0, 
                 Mass -> Automatic,
                 LaTeX -> "H",
                 OutputName -> "H" }},
   

   {sigma1, {LaTeX -> "\\sigma_1"}},
   {sigma2, {LaTeX -> "\\sigma_2"}},

   {phi1, {LaTeX -> "\\phi_1"}},
   {phi2, {LaTeX -> "\\phi_2"}},

   {dR,     {LaTeX -> "d_R" }},
   {eR,     {LaTeX -> "e_R" }},
   {uR,     {LaTeX -> "u_R" }},
   {eL,     {LaTeX -> "e_L" }},
   {dL,     {LaTeX -> "d_L" }},
   {uL,     {LaTeX -> "u_L" }},
   {vL,     {LaTeX -> "\\nu_L" }}, 
   {vR,     {LaTeX -> "\\nu_R"}},
   
   {DR,     {LaTeX -> "D_R" }},
   {ER,     {LaTeX -> "E_R" }},
   {UR,     {LaTeX -> "U_R" }},
   {EL,     {LaTeX -> "E_L" }},
   {DL,     {LaTeX -> "D_L" }},
   {UL,     {LaTeX -> "U_L" }},
   {Ksim,   {LaTeX -> "\\Xi_m"}}
 };
