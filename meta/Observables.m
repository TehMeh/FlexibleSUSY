BeginPackage["Observables`", {"FlexibleSUSY`", "SARAH`", "BetaFunction`", "Parameters`", "TreeMasses`", "Utils`", "CConversion`", "TextFormatting`"}];

(* observables *)
Begin["FlexibleSUSYObservable`"];
FSObservables = { aMuonGM2Calc, aMuonGM2CalcUncertainty,
                  CpHiggsPhotonPhoton, CpHiggsGluonGluon,
                  CpPseudoScalarPhotonPhoton, CpPseudoScalarGluonGluon };
End[];

GetRequestedObservables::usage="";
CountNumberOfObservables::usage="";
CreateObservablesDefinitions::usage="";
CreateObservablesInitialization::usage="";
CreateSetAndDisplayObservablesFunctions::usage="";
CreateClearObservablesFunction::usage="";
CalculateObservables::usage="";

Begin["`Private`"];

GetRequestedObservables[blocks_] :=
    Module[{observables, dim},
           observables = DeleteDuplicates[Cases[blocks, a_?(MemberQ[FlexibleSUSYObservable`FSObservables,#]&) :> a, {0, Infinity}]];
           If[MemberQ[observables, FlexibleSUSYObservable`CpHiggsPhotonPhoton] ||
              MemberQ[observables, FlexibleSUSYObservable`CpHiggsGluonGluon],
              dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson]
              If[FreeQ[TreeMasses`GetParticles[], SARAH`HiggsBoson] ||
                 TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson] == 0,
                 Print["Warning: no physical Higgs boson found."];
                 Print["         Effective couplings for Higgs boson will not"];
                 Print["         be calculated."];
                 observables = DeleteCases[observables, a_ /; (a === FlexibleSUSYObservable`CpHiggsPhotonPhoton ||
                                                               a === FlexibleSUSYObservable`CpHiggsGluonGluon)];
                ];
             ];
           If[MemberQ[observables, FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton] ||
              MemberQ[observables, FlexibleSUSYObservable`CpPseudoScalarGluonGluon],
              If[FreeQ[TreeMasses`GetParticles[], SARAH`PseudoScalar] ||
                 TreeMasses`GetDimensionWithoutGoldstones[SARAH`PseudoScalar] == 0,
                 Print["Warning: no physical pseudoscalar boson found."];
                 Print["         Effective couplings for pseudoscalar boson will not"];
                 Print["         be calculated."];
                 observables = DeleteCases[observables, a_ /; (a === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton ||
                                                               a === FlexibleSUSYObservable`CpPseudoScalarGluonGluon)]; 
                ];
             ];
           observables
          ];

GetObservableName[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc] := "a_muon_gm2calc";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty] := "a_muon_gm2calc_uncertainty";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`CpHiggsPhotonPhoton] := "eff_cp_higgs_photon_photon";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`CpHiggsGluonGluon] := "eff_cp_higgs_gluon_gluon";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton] := "eff_cp_pseudoscalar_photon_photon";
GetObservableName[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarGluonGluon] := "eff_cp_pseudoscalar_gluon_gluon";

GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc] := "a_muon = (g-2)/2 of the muon (calculated with GM2Calc)";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty] := "uncertainty of (g-2)/2 of the muon (calculated with GM2Calc)";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`CpHiggsPhotonPhoton] := "effective H-Photon-Photon coupling";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`CpHiggsGluonGluon] := "effective H-Gluon-Gluon coupling";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton] := "effective A-Photon-Photon coupling";
GetObservableDescription[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarGluonGluon] := "effective A-Gluon-Gluon coupling";

GetObservableType[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc] := CConversion`ScalarType[CConversion`realScalarCType];
GetObservableType[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty] := CConversion`ScalarType[CConversion`realScalarCType];

GetObservableType[obs_ /; obs === FlexibleSUSYObservable`CpHiggsPhotonPhoton] :=
    Module[{dim, type},
           dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson];
           If[dim == 1,
              type = CConversion`ScalarType[CConversion`complexScalarCType],
              type = CConversion`ArrayType[CConversion`complexScalarCType, dim]
             ];
           type
          ];

GetObservableType[obs_ /; obs === FlexibleSUSYObservable`CpHiggsGluonGluon] :=
    Module[{dim, type},
           dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`HiggsBoson];
           If[dim == 1,
              type = CConversion`ScalarType[CConversion`complexScalarCType],
              type = CConversion`ArrayType[CConversion`complexScalarCType, dim]
             ];
           type
          ];

GetObservableType[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton] :=
    Module[{dim, type},
           dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`PseudoScalar];
           If[dim == 1,
              type = CConversion`ScalarType[CConversion`complexScalarCType],
              type = CConversion`ArrayType[CConversion`complexScalarCType, dim]
             ];
           type
          ];

GetObservableType[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarGluonGluon] :=
    Module[{dim, type},
           dim = TreeMasses`GetDimensionWithoutGoldstones[SARAH`PseudoScalar];
           If[dim == 1,
              type = CConversion`ScalarType[CConversion`complexScalarCType],
              type = CConversion`ArrayType[CConversion`complexScalarCType, dim]
             ];
           type
          ];

CountNumberOfObservables[observables_List] :=
    Module[{i, number = 0},
           For[i = 1, i <= Length[observables], i++,
               If[MemberQ[FlexibleSUSYObservable`FSObservables, observables[[i]]],
                  number += BetaFunction`CountNumberOfParameters[GetObservableType[observables[[i]]]];,
                  Print["Warning: ignoring invalid observable ", observables[[i]]];
                 ];
              ];
           number
          ];

CreateObservablesDefinitions[observables_List] :=
    Module[{i, type, name, description, definitions = ""},
           For[i = 1, i <= Length[observables], i++,
               If[MemberQ[FlexibleSUSYObservable`FSObservables, observables[[i]]],
                  name = GetObservableName[observables[[i]]];
                  description = GetObservableDescription[observables[[i]]];
                  type = CConversion`CreateCType[GetObservableType[observables[[i]]]];
                  definitions = definitions <> type <> " " <> name <> "; ///< " <> description <> "\n";,
                  Print["Warning: ignoring invalid observable ", observables[[i]]];
                 ];
              ];
           definitions
          ];

CreateObservablesInitialization[observables_List] :=
    Module[{i, name, type, init = ""},
           For[i = 1, i <= Length[observables], i++,
               If[MemberQ[FlexibleSUSYObservable`FSObservables, observables[[i]]],
                  name = GetObservableName[observables[[i]]];
                  type = GetObservableType[observables[[i]]];
                  If[init == "",
                     init = ": " <> CConversion`CreateDefaultConstructor[name, type] <> "\n";,
                     init = init <> ", " <> CConversion`CreateDefaultConstructor[name, type] <> "\n";
                    ];,
                  Print["Warning: ignoring invalid observable ", observables[[i]]];
                 ];
              ];
           init
          ];

CreateSetAndDisplayObservablesFunctions[observables_List] :=
    Module[{numObservables, i, name, type, paramCount = 0, nAssignments, assignment,
            display = "", displayNames = "", set = ""},
           numObservables = CountNumberOfObservables[observables];
           If[numObservables != 0,
              display = "Eigen::ArrayXd vec(" <> FlexibleSUSY`FSModelName
                        <> "_observables::NUMBER_OF_OBSERVABLES);\n\n";
              displayNames = "std::vector<std::string> names("
                             <> FlexibleSUSY`FSModelName
                             <> "_observables::NUMBER_OF_OBSERVABLES);\n\n";
              set = "assert(vec.rows() == " <> FlexibleSUSY`FSModelName
                    <> "_observables::NUMBER_OF_OBSERVABLES);\n\n";
              For[i = 1, i <= Length[observables], i++,
                  If[MemberQ[FlexibleSUSYObservable`FSObservables, observables[[i]]],
                     name = GetObservableName[observables[[i]]];
                     type = GetObservableType[observables[[i]]];
                     {assignment, nAssignments} = Parameters`CreateSetAssignment[name, paramCount, type, "vec"];
                     set = set <> assignment;
                     {assignment, nAssignments} = Parameters`CreateDisplayAssignment[name, paramCount, type, "vec"];
                     display = display <> assignment;
                     {assignment, nAssignments} = Parameters`CreateStdVectorNamesAssignment[name, paramCount, type];
                     displayNames = displayNames <> assignment;
                     paramCount += nAssignments;,
                     Print["Warning: ignoring invalid observable ", observables[[i]]];
                    ];
                 ];,
               display = "Eigen::ArrayXd vec(1);\n\nvec(0) = 0.;\n";
               set = "";
               displayNames = "std::vector<std::string> names(1);\n\n"
                              <> "names[0] = \"no observables defined\";\n";
             ];
           {display, displayNames, set}
          ];

CreateClearObservablesFunction[observables_List] :=
    Module[{i, name, type, result = ""},
           For[i = 1, i <= Length[observables], i++,
               If[MemberQ[FlexibleSUSYObservable`FSObservables, observables[[i]]],
                  name = GetObservableName[observables[[i]]];
                  type = GetObservableType[observables[[i]]];
                  result = result <> CConversion`SetToDefault[name, type];,
                  Print["Warning: ignoring invalid observable ", observables[[i]]];
                 ];
              ];
           result
          ];

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2Calc, structName_String] :=
    structName <> ".AMUGM2CALC = gm2calc_calculate_amu(gm2calc_data);";

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`aMuonGM2CalcUncertainty, structName_String] :=
    structName <> ".AMUGM2CALCUNCERTAINTY = gm2calc_calculate_amu_uncertainty(gm2calc_data);";

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`CpHiggsPhotonPhoton, structName_String] :=
    Module[{i, type, dim, result = ""},
           type = GetObservableType[obs];
           If[MatchQ[type, CConversion`ArrayType[_,_]],
              dim = type /. CConversion`ArrayType[_, d_] -> d;
              For[i = 1, i <= dim, i++,
                  result = result <> structName <> ".EFFCPHIGGSPHOTONPHOTON("
                           <> ToString[i-1] <> ") = effective_couplings.eff_Cp"
                           <> CConversion`ToValidCSymbolString[SARAH`HiggsBoson]
                           <> CConversion`ToValidCSymbolString[SARAH`VectorP]
                           <> CConversion`ToValidCSymbolString[SARAH`VectorP] <> "("
                           <> ToString[i-1] <> If[i != dim, ");\n", ");"];
                 ];,
              dim == 1;
              result = structName <> ".EFFCPHIGGSPHOTONPHOTON = effective_couplings.eff_Cp"
               <> CConversion`ToValidCSymbolString[SARAH`HiggsBoson]
               <> CConversion`ToValidCSymbolString[SARAH`VectorP]
               <> CConversion`ToValidCSymbolString[SARAH`VectorP] <> "();"
             ];
           result
          ];

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`CpHiggsGluonGluon, structName_String] :=
    Module[{i, type, dim, result = ""},
           type = GetObservableType[obs];
           If[MatchQ[type, CConversion`ArrayType[_,_]],
              dim = type /. CConversion`ArrayType[_, d_] -> d;
              For[i = 1, i <= dim, i++,
                  result = result <> structName <> ".EFFCPHIGGSGLUONGLUON("
                           <> ToString[i-1] <> ") = effective_couplings.eff_Cp"
                           <> CConversion`ToValidCSymbolString[SARAH`HiggsBoson]
                           <> CConversion`ToValidCSymbolString[SARAH`VectorG]
                           <> CConversion`ToValidCSymbolString[SARAH`VectorG] <> "("
                           <> ToString[i-1] <> If[i != dim, ");\n", ");"];
                 ];,
              dim == 1;
              result = structName <> ".EFFCPHIGGSGLUONGLUON = effective_couplings.eff_Cp"
               <> CConversion`ToValidCSymbolString[SARAH`HiggsBoson]
               <> CConversion`ToValidCSymbolString[SARAH`VectorG]
               <> CConversion`ToValidCSymbolString[SARAH`VectorG] <> "();"
             ];
           result
          ];

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton, structName_String] :=
    Module[{i, type, dim, result = ""},
           type = GetObservableType[obs];
           If[MatchQ[type, CConversion`ArrayType[_,_]],
              dim = type /. CConversion`ArrayType[_, d_] -> d;
              For[i = 1, i <= dim, i++,
                  result = result <> structName <> ".EFFCPPSEUDOSCALARPHOTONPHOTON("
                           <> ToString[i-1] <> ") = effective_couplings.eff_Cp"
                           <> CConversion`ToValidCSymbolString[SARAH`PseudoScalar]
                           <> CConversion`ToValidCSymbolString[SARAH`VectorP]
                           <> CConversion`ToValidCSymbolString[SARAH`VectorP] <> "("
                           <> ToString[i-1] <> If[i != dim, ");\n", ");"];
                 ];,
              dim == 1;
              result = structName <> ".EFFCPPSEUDOSCALARPHOTONPHOTON = effective_couplings.eff_Cp"
               <> CConversion`ToValidCSymbolString[SARAH`PseudoScalar]
               <> CConversion`ToValidCSymbolString[SARAH`VectorP]
               <> CConversion`ToValidCSymbolString[SARAH`VectorP] <> "();"
             ];
           result
          ];

CalculateObservable[obs_ /; obs === FlexibleSUSYObservable`CpPseudoScalarGluonGluon, structName_String] :=
    Module[{i, type, dim, result = ""},
           type = GetObservableType[obs];
           If[MatchQ[type, CConversion`ArrayType[_,_]],
              dim = type /. CConversion`ArrayType[_, d_] -> d;
              For[i = 1, i <= dim, i++,
                  result = result <> structName <> ".EFFCPPSEUDOSCALARGLUONGLUON("
                           <> ToString[i-1] <> ") = effective_couplings.eff_Cp"
                           <> CConversion`ToValidCSymbolString[SARAH`PseudoScalar]
                           <> CConversion`ToValidCSymbolString[SARAH`VectorG]
                           <> CConversion`ToValidCSymbolString[SARAH`VectorG] <> "("
                           <> ToString[i-1] <> If[i != dim, ");\n", ");"];
                 ];,
              dim == 1;
              result = structName <> ".EFFCPPSEUDOSCALARGLUONGLUON = effective_couplings.eff_Cp"
               <> CConversion`ToValidCSymbolString[SARAH`PseudoScalar]
               <> CConversion`ToValidCSymbolString[SARAH`VectorG]
               <> CConversion`ToValidCSymbolString[SARAH`VectorG] <> "();"
             ];
           result
          ];

FillGM2CalcInterfaceData[struct_String] :=
    Module[{filling, mwStr,
            w, pseudoscalar, smuon, muonsneutrino, chargino, neutralino,
            mu, m1, m2, m3, mq2, mu2, md2, ml2, me2, tu, td, te, yu, yd, ye},
           w             = Parameters`GetParticleFromDescription["W-Boson"];
           pseudoscalar  = Parameters`GetParticleFromDescription["Pseudo-Scalar Higgs"];
           smuon         = Parameters`GetParticleFromDescription["Smuon"];
           muonsneutrino = Parameters`GetParticleFromDescription["Muon Sneutrino"];
           chargino      = Parameters`GetParticleFromDescription["Charginos"];
           neutralino    = Parameters`GetParticleFromDescription["Neutralinos"];
           mu            = Parameters`GetParameterFromDescription["Mu-parameter"];
           m1            = Parameters`GetParameterFromDescription["Bino Mass parameter"];
           m2            = Parameters`GetParameterFromDescription["Wino Mass parameter"];
           m3            = Parameters`GetParameterFromDescription["Gluino Mass parameter"];
           mq2           = Parameters`GetParameterFromDescription["Softbreaking left Squark Mass"];
           mu2           = Parameters`GetParameterFromDescription["Softbreaking right Up-Squark Mass"];
           md2           = Parameters`GetParameterFromDescription["Softbreaking right Down-Squark Mass"];
           ml2           = Parameters`GetParameterFromDescription["Softbreaking left Slepton Mass"];
           me2           = Parameters`GetParameterFromDescription["Softbreaking right Slepton Mass"];
           tu            = Parameters`GetParameterFromDescription["Trilinear-Up-Coupling"];
           td            = Parameters`GetParameterFromDescription["Trilinear-Down-Coupling"];
           te            = Parameters`GetParameterFromDescription["Trilinear-Lepton-Coupling"];
           yu            = Parameters`GetParameterFromDescription["Up-Yukawa-Coupling"];
           yd            = Parameters`GetParameterFromDescription["Down-Yukawa-Coupling"];
           ye            = Parameters`GetParameterFromDescription["Lepton-Yukawa-Coupling"];
           mwStr         = "MODEL.get_physical()." <> CConversion`RValueToCFormString[FlexibleSUSY`M[w]];
           filling = \
           struct <> ".alpha_s_MZ = ALPHA_S_MZ;\n" <>
           struct <> ".MZ    = MZPole;\n" <>
           "if (!is_zero(" <> mwStr <> "))\n" <>
              TextFormatting`IndentText[struct <> ".MW = " <> mwStr <> ";"] <> "\n" <>
           "else if (!is_zero(MWPole))\n" <>
              TextFormatting`IndentText[struct <> ".MW = MWPole;"] <> "\n" <>
           struct <> ".mb_mb = MBMB;\n" <>
           struct <> ".MT    = MTPole;\n" <>
           struct <> ".MTau  = MTauPole;\n" <>
           struct <> ".MM    = MMPole;\n" <>
           struct <> ".MA0   = MODEL.get_physical()." <>
           CConversion`RValueToCFormString[FlexibleSUSY`M[pseudoscalar][1]] <> ";\n" <>
           struct <> ".MSvm  = MODEL.get_physical()." <>
           CConversion`RValueToCFormString[FlexibleSUSY`M[muonsneutrino]] <> ";\n" <>
           struct <> ".MSm   = MODEL.get_physical()." <>
           CConversion`RValueToCFormString[FlexibleSUSY`M[smuon]] <> ";\n" <>
           struct <> ".MCha  = MODEL.get_physical()." <>
           CConversion`RValueToCFormString[FlexibleSUSY`M[chargino]] <> ";\n" <>
           struct <> ".MChi  = MODEL.get_physical()." <>
           CConversion`RValueToCFormString[FlexibleSUSY`M[neutralino]] <> ";\n" <>
           struct <> ".scale = MODEL.get_scale();\n" <>
           struct <> ".TB    = MODEL.get_" <> CConversion`RValueToCFormString[SARAH`VEVSM2] <> "() / " <>
                              "MODEL.get_" <> CConversion`RValueToCFormString[SARAH`VEVSM1] <> "();\n" <>
           struct <> ".Mu    = MODEL.get_" <> CConversion`RValueToCFormString[mu] <> "();\n" <>
           struct <> ".M1    = MODEL.get_" <> CConversion`RValueToCFormString[m1] <> "();\n" <>
           struct <> ".M2    = MODEL.get_" <> CConversion`RValueToCFormString[m2] <> "();\n" <>
           struct <> ".M3    = MODEL.get_" <> CConversion`RValueToCFormString[m3] <> "();\n" <>
           struct <> ".mq2   = MODEL.get_" <> CConversion`RValueToCFormString[mq2] <> "();\n" <>
           struct <> ".mu2   = MODEL.get_" <> CConversion`RValueToCFormString[mu2] <> "();\n" <>
           struct <> ".md2   = MODEL.get_" <> CConversion`RValueToCFormString[md2] <> "();\n" <>
           struct <> ".ml2   = MODEL.get_" <> CConversion`RValueToCFormString[ml2] <> "();\n" <>
           struct <> ".me2   = MODEL.get_" <> CConversion`RValueToCFormString[me2] <> "();\n" <>
           struct <> ".Au    = div_save(MODEL.get_" <> CConversion`RValueToCFormString[tu] <>
                               "(), MODEL.get_" <> CConversion`RValueToCFormString[yu] <> "());\n" <>
           struct <> ".Ad    = div_save(MODEL.get_" <> CConversion`RValueToCFormString[td] <>
                               "(), MODEL.get_" <> CConversion`RValueToCFormString[yd] <> "());\n" <>
           struct <> ".Ae    = div_save(MODEL.get_" <> CConversion`RValueToCFormString[te] <>
                               "(), MODEL.get_" <> CConversion`RValueToCFormString[ye] <> "());\n";
           "GM2Calc_data " <> struct <> ";\n" <> filling
          ];

FillEffectiveCouplingsInterfaceData[struct_String] :=
    FlexibleSUSY`FSModelName <> "_effective_couplings " <> struct <> "(model, qedqcd);\n";

FillInterfaceData[{}] := "";

FillInterfaceData[obs_List] :=
    Module[{filled = ""},
           If[MemberQ[obs,FlexibleSUSYObservable`aMuonGM2Calc] ||
              MemberQ[obs,FlexibleSUSYObservable`aMuonGM2CalcUncertainty],
              filled = filled <> FillGM2CalcInterfaceData["gm2calc_data"];
             ];
           If[MemberQ[obs,FlexibleSUSYObservable`CpHiggsPhotonPhoton]         ||
              MemberQ[obs,FlexibleSUSYObservable`CpHiggsGluonGluon]           ||
              MemberQ[obs, FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton] ||
              MemberQ[obs, FlexibleSUSYObservable`CpPseudoScalarGluonGluon],
              filled = filled <> FillEffectiveCouplingsInterfaceData["effective_couplings"];
             ];
           filled
          ];

CalculateObservables[something_, structName_String] :=
    Module[{observables},
           observables = Cases[something, a_?(MemberQ[FlexibleSUSYObservable`FSObservables,#]&) :> a, {0, Infinity}];
           FillInterfaceData[observables] <> "\n" <>
           Utils`StringJoinWithSeparator[CalculateObservable[#,structName]& /@ observables, "\n"]
          ];

End[];

EndPackage[];
