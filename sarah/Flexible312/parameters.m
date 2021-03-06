(* ::Package:: *)

ParameterDefinitions =
{
	{g1, {Description -> "Hypercharge-Coupling"}},
	{g2, {Description -> "Left-Coupling"}},
	{g3, {Description -> "Strong-Coupling"}},
	
	{M112, {LaTeX -> "m^2_{1}",
			OutputName -> M112,
			Real->True,
               LesHouches -> {HMIX,20}}},
	{M222, {LaTeX -> "m^2_{2}",
			OutputName -> M222,
			Real->True,
               LesHouches -> {HMIX,21}}},
	{M12, {LaTeX -> "M_{12}^2",
			OutputName -> M12,
               LesHouches -> {HMIX,22}}},
	{MR, {LaTeX -> "m_{R}",
			OutputName -> MR,
			LesHouches->{HMIX, 23}}},
	{Z1111, {LaTeX -> "Z_{1111}",
			OutputName -> Z1111,
			Real->True,
               LesHouches -> {HMIX,31}}},
	{Z2222, {LaTeX -> "Z_{2222}",
			OutputName -> Z2222,
			Real->True,
               LesHouches -> {HMIX,32}}},
	{Z2211, {LaTeX -> "Z_{2211}",
			OutputName -> Z2211,
			Real->True,
               LesHouches -> {HMIX,33}}},
	{Z2112, {LaTeX -> "Z_{2112}",
			OutputName -> Z2112,
			Real->True,
               LesHouches -> {HMIX,34}}},
	{Z2121, {LaTeX -> "Z_{2121}",
			OutputName -> Z2121,
               LesHouches -> {HMIX,35}}},
	{Z2212, {LaTeX -> "Z_{2212}",
			OutputName -> Z2212,
               LesHouches -> {HMIX,36}}},
	{Z1112, {LaTeX -> "Z_{1112}",
			OutputName -> Z1112,
               LesHouches -> {HMIX,37}}},
               
	{Ye1, {Description -> "Lepton-Yukawa-Coupling",
			OutputName -> Ye1,
			DependenceNum->None,
			LesHouches->Ye1,
			LaTeX->"Y_{E1}"}},
    {Yd1, {Description -> "Down-Yukawa-Coupling",
			OutputName -> Yd1,
			DependenceNum->None,
			LesHouches->Yd1,
			LaTeX->"Y_{D1}"}},
	{Yu1, {Description -> "Up-Yukawa-Coupling",
			OutputName -> Yu1,
			DependenceNum->None,
			LesHouches->Yu1,
			LaTeX->"Y_{U1}"}},
    {Yu2, {Description -> "Up-Yukawa-Coupling-Dub2",
           OutputName -> Yu2,
           LesHouches-> Yu2,
           LaTeX-> "Y_{U2}"}},
    {Yd2, {Description -> "Down-Yukawa-Coupling-Dub2",
           OutputName -> Yd2,
           LesHouches-> Yd2,
           LaTeX->"U_{D2}"}},
    {Ye2, {Description -> "Lepton-Yukawa-Coupling-Dub2",
           OutputName -> Ye2,
           LesHouches-> Ye2,
           LaTeX->"Y_{E1}"}},
    {Yv1, {Description -> "N0-Yukawa-Coupling-Dub1",
           OutputName -> Yv1,
           LesHouches-> Yv1,
           LaTeX->"Y_{N1}"}},
    {Yv2, {Description -> "N0-Yukawa-Coupling-Dub2",
           OutputName -> Yv2,
           LesHouches-> Yv2,
           LaTeX->"Y_{N2}"}},
			
	{ThetaW,    { Description -> "Weinberg-Angle"}}, 
	{ZZ, {Description ->   "Photon-Z Mixing Matrix"}}, (*Is this enough to get the matrix from the "global(?)" parameter file?*)
	{ZW, {Description -> "W Mixing Matrix" }},
	
(*    {v1,        { Description -> "Down-VEV", LaTeX -> "v_1"}}, 
    {v2,        { Description -> "Up-VEV", LaTeX -> "v_2"}},       
    {v,         { Description -> "EW-VEV"}},*)
    {v,        { Description -> "EW-VEV",
                 DependenceNum -> Sqrt[4*Mass[VWp]^2/(g2^2)],
                 DependenceSPheno -> None  }},
    (*{vTemp,     { Description -> "VEV-STUB",
                 OutputName -> vStub, Value -> 0,
                 LesHouches->{HMIX, 104}  }},*)
                 
	{Vu,        {Description ->"Left-Up-Mixing-Matrix"}},
	{Vd,        {Description ->"Left-Down-Mixing-Matrix"}},
	{Uu,        {Description ->"Right-Up-Mixing-Matrix"}},
	{Ud,        {Description ->"Right-Down-Mixing-Matrix"}}, 
	{Ve,        {Description ->"Left-Lepton-Mixing-Matrix"}},
	{Ue,        {Description ->"Right-Lepton-Mixing-Matrix"}},
	{Vpmns,     {Description -> "Neutrino-Mixing-Matrix",
	             LaTeX->"V_{PMNS}",
	             Value->None,
	             LesHouches->VPMNS,
	             OutputName->Vpmns}},	
	
		{\[Beta],   { Description -> "Pseudo Scalar mixing angle"  }},             
		{\[Alpha],  { Description -> "Scalar mixing angle" }},  
	
	    (*{ZH,        { Description->"Scalar-Mixing-Matrix"}},
		{ZA,        { Description->"Pseudo-Scalar-Mixing-Matrix"}},
		{ZP,        { Description->"Charged-Mixing-Matrix"}},*)
		

		 {ZH,        { Description->"Scalar-Mixing-Matrix",
                       Dependence -> None,
                       DependenceNum -> None,
                       DependenceOptional -> None,
                       Real -> True }},

		{ZP,        { Description->"Charged-Mixing-Matrix",
              Real -> False,
              Dependence -> None,
              DependenceOptional -> None,
              DependenceNum -> None}},  
		(*^^^ from THDM-II ^^^*)
		
{		AlphaS,    {Description -> "Alpha Strong"}},	
		{e,         { Description -> "electric charge"}}, 

		{Gf,        { Description -> "Fermi's constant"}},
		{aEWinv,    { Description -> "inverse weak coupling constant at mZ"}},
		
{eta,       { Real -> True, 
              OutputName ->"eta",
              LaTeX -> "\\eta",
              LesHouches->{HMIX,500} }}(*,

{Cos[eta]+I Sin[eta], { OutputName ->Cos[eta]+I Sin[eta],
              LaTeX -> "\\cos(\\eta)+i\\sin(\\eta)",
              LesHouches->{HMIX,601} }}*)
};
