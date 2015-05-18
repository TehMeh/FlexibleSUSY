
BeginPackage["GMuonMinus2`", {"SARAH`", "TextFormatting`", "TreeMasses`", "LoopMasses`", "Vertices`"}];

Clear[CreateCalculationFunctions];
Clear[CreateCalculation];

GetPhysicalName::usage="The name of the c++ data type that stores the resulting magnetic moment";
CreatePhysicalDefinition::usage="Returns the c++ code with declaration of the magnetic moment variable";

CreateParticles::usage="Returns the c++ code that contains all particle classes";
CreateDiagrams::usage="Returns the c++ code that contains all relevant diagram classes";
CreateVertexFunctions::usage="Returns the c++ code that contains all relevant vertex functions";

CreateCalculation::usage="Returns the c++ code that performs the actual calculation the magnetic moment";
CreateThreadedCalculation::usage="Same as above, threaded version";

CreateDefinitions::usage="Returns the c++ that contains all function definitions"

NPointFunctions::usage="Returns a list of all n point functions that are needed. Actually it is a list of fake functions to extract vertex functions...";

Begin["`Private`"];

CXXNamespace = (FlexibleSUSY`FSModelName <> "_GMuonMinus2");

(************* Begin public interface *******************)

GetPhysicalName[] := "GMuonMinus2";
CreatePhysicalDefinition[] := "double " <> GMuonMinus2`GetPhysicalName[] <> ";";

(* Create c++ classes for all particles *)
CreateParticles[] := Module[{particles, muonIndex, code},
                            (* Get a list of all particles *)
                            particles = TreeMasses`GetParticles[];
                            muonIndex = GetMuonIndex[];
                            
                            code = ("// Particles (SARAH-style)\n" <>
                                    "struct Particle {};\n\n" <>
                                    
                                    StringJoin @ Riffle[("struct " <> ParticleToString[#] <>
                                                         ": public Particle { " <>
                                                         "static const int numberOfGenerations = " <>
                                                         ToString @ TreeMasses`GetDimension[#] <>
                                                         "; };" &) /@ particles, "\n"] <> "\n\n" <>
                                    "// Special particle families\n" <>
                                    "using Photon = " <> ParticleToString @ GetPhoton[] <> ";\n" <>
                                    "using MuonFamily = " <> ParticleToString @ GetMuonFamily[] <> ";\n\n" <>
                                    
                                    "static const unsigned int muonIndex( void )\n" <>
                                    "{ unsigned int muonIndex" <>
                                    If[muonIndex =!= Null, " = " <> ToString[muonIndex-1], ""] <>
                                    "; return muonIndex; }\n" <>
                                    "static const double muonCharge( const EvaluationContext &context )\n" <>
                                    "{ return context.model." <>
                                    NameOfCouplingFunction[{GetPhoton[], GetMuonFamily[], SARAH`bar[GetMuonFamily[]]}] <> "PL" <>
                                    If[muonIndex =!= Null,
                                       "( " <> ToString[muonIndex-1] <> ", " <> ToString[muonIndex-1] <> " )",
                                       "()"] <> "; }\n\n" <>
                                    
                                    "// AntiParticles (SARAH-style)\n\n" <>
                                    "template<class P> struct bar { using type = bar<P>; };\n" <>
                                    "template<class P> struct conj { using type = conj<P>; };\n\n" <>
                                    
                                    "template<class P> struct bar<bar<P>> { using type = P; };\n" <>
                                    "template<class P> struct conj<conj<P>> { using type = P; };\n\n" <>
                                    
                                    "// Particles that are their own antiparticles\n" <>
                                    StringJoin @ Riffle[("template<> struct " <>
                                                         If[IsScalar[#] || IsVector[#],
                                                            "conj<" <> ParticleToString[#] <> ">",
                                                            "bar<" <> ParticleToString[#] <> ">"] <>
                                                         " { using type = " <> ParticleToString[#] <> "; };"
                                                         &) /@ Select[particles, (# == AntiParticle[#] &)],
                                                        "\n"]
                                    );
                            
                            Return[code];
                            ];

CreateDiagrams[] := Module[{diagramTypes, diagramTypeHeads, code},
                           diagrams = contributingFeynmanDiagramTypes;
                           diagramHeads = DeleteDuplicates @ (Head /@ diagrams);
                           
                           code = "// The different diagram types that contribute to the muon magnetic moment\n";
                           code = (code <>
                                   StringJoin @ Riffle[("template<unsigned int> class " <> SymbolName[#] <> ";" &)
                                                       /@ diagramHeads, "\n"] <>
                                   "\n\n");
                           
                           code = (code <> "// Indexed diagram types\n" <>
                                   StringJoin @ Riffle[("template<> class " <> SymbolName[Head[#]] <>
                                                        "<" <> ToString @ #[[1]] <> "> {};" &)
                                                       /@ diagrams, "\n"]);
                           
                           code = (code <> "\n\n" <>
                                   "template<class ...Args> struct DiagramEvaluator;\n\n" <>
                                   StringJoin @ Riffle[CreateDiagramEvaluatorClass /@ contributingFeynmanDiagramTypes, "\n\n"]);
                           
                           Return[code];
                           ];

CreateVertexFunctions[] := CreateVertices[][[1]];

CreateDiagramEvaluatorClass[type_OneLoopDiagram] := ("template<class PhotonEmitter, class ExchangeParticle>\n" <>
                                                     "struct DiagramEvaluator<OneLoopDiagram<" <>
                                                     ToString @ type[[1]] <>
                                                     ">, PhotonEmitter, ExchangeParticle>\n" <>
                                                     "{ static double value( const EvaluationContext &context ); };");

calculationCode = Null;
CreateCalculation[] := Module[{code},
                              (* If we have been here before return the old result *)
                              If[calculationCode =!= Null, Return[calculationCode]];
                              
                              code = "/********** GMuonMinus2.m generated calculation code **********/\n\n";
                              
                              (* Generate code that simply adds up all contributions *)
                              code = (code <>
                                      "EvaluationContext context{ model };\n" <>
                                      "double val = 0.0;\n\n" <>
                                      StringJoin @ Riffle[("val += " <> # <> "::value( context );" &) /@ ConcreteDiagramEvaluators[],
                                                          "\n"] <> "\n\n" <>
                                      "return val;"
                                      );
                              
                              calculationCode = code;
                              Return[code];
                              ];

CreateThreadedCalculation[] := CreateCalculation[];

CreateDefinitions[] := (CreateEvaluationContextSpecializations[] <> "\n\n" <>
                        CreateVertices[][[2]]);

NPointFunctions[] := (If[nPointFunctions === Null, CreateVertices[]]; Return[nPointFunctions]);


(**************** End public interface *****************)

CreateEvaluationContextSpecializations[] :=
Module[{particles, code},
       particles = TreeMasses`GetParticles[];
       particles = Select[particles, (! TreeMasses`IsGhost[#] &)];
       
       code = (StringJoin @
               Riffle[("template<> double EvaluationContext::mass<" <> ToString[#] <> ">( " <>
                       If[TreeMasses`GetDimension[#] === 1, "void", "int index"] <> " ) const\n" <>
                       "{ return model.get_M" <> ParticleToString[#] <>
                       If[TreeMasses`GetDimension[#] === 1, "()", "( index )"] <> "; }"
                       &) /@ particles, "\n\n"]);
       
       Return[code];
       ];

(************************ Begin helper routines *******************************)

GetMuonFamily[] := If[TreeMasses`GetDimension[SARAH`Electron] =!= 1,
                        SARAH`Electron,
                        Cases[SARAH`ParticleDefinitions[FlexibleSUSY`FSEigenstates],
                              {p_, {Description -> "Muon", ___}} -> p, 1][[1]]
                        ];
GetMuonIndex[] := If[TreeMasses`GetDimension[SARAH`Electron] =!= 1, 2, Null];

GetPhoton[] := SARAH`Photon;

IsLorentzIndex[index_] := StringMatchQ[ToString @ index, "lt" ~~ __];

StripLorentzIndices[p_Symbol] := p;
StripLorentzIndices[SARAH`bar[p_]] := SARAH`bar[StripLorentzIndices[p]];
StripLorentzIndices[Susyno`LieGroups`conj[p_]] := Susyno`LieGroups`conj[StripLorentzIndices[p]];
StripLorentzIndices[p_] := Module[{remainingIndices},
                                  remainingIndices = Select[p[[1]], (!IsLorentzIndex[#] &)];
                                  If[Length[remainingIndices] === 0, Head[p],
                                     Head[p][remainingIndices]]
                                  ];
SetAttributes[StripLorentzIndices, {Listable}];

(* Takes a SARAH particle and returns its antiparticle *)
AntiParticle[SARAH`bar[p_]] := p;
AntiParticle[Susyno`LieGroups`conj[p_]] := p;
AntiParticle[p_] := Module[{pNoIndices = Vertices`StripFieldIndices[p]},
                           If[IsScalar[pNoIndices] || IsVector[pNoIndices],
                              Susyno`LieGroups`conj[p],
                              SARAH`bar[p]]];
SetAttributes[AntiParticle, {Listable}];

ParticleToString[p_] := SymbolName[p];
ParticleToString[SARAH`bar[p_]] := "bar<" <> SymbolName[p] <> ">::type";
ParticleToString[Susyno`LieGroups`conj[p_]] := CXXNamespace <> "::conj<" <> SymbolName[p] <> ">::type";

ChangeIndexNumber[index_,number_Integer] := Symbol @ StringReplace[index,
                       Shortest[name__] ~~ NumberString :> name ~~ ToString[number]];

IndexReplacementRulesForNewParticleIndexNumber[particle_, number_Integer] :=
    ((# -> ChangeIndexNumber[ToString[#], number] &) /@ Vertices`FieldIndexList[particle]);

IndexReplacementRulesForParticleReordering[particles_List, ordering_] :=
Module[{indices = Table[i, {i, Length[particles]}], index, fieldIndexList},
       Flatten[(IndexReplacementRulesForNewParticleIndexNumber[particles[[ordering[[#]]]], #] &) /@ indices]];

OrderParticles[particles_List, ordering_] := Module[{indexRules},
       indexRules = IndexReplacementRulesForParticleReordering[particles, ordering];
       Return[particles[[ordering]] /. indexRules];
                                                    ];

OrderVertex[vertex_, ordering_] := Module[{indexRules, particles, expr, newVertex},
                                          indexRules = IndexReplacementRulesForParticleReordering[vertex[[1]], ordering];
                                          
                                          particles = vertex[[1]][[ordering]];
                                          expr = vertex[[2;;]];
                                          
                                          newVertex = (Join[{particles}, expr] /. indexRules);
                                          Return[newVertex];
                                          ];

(* MemoizingVertex[] works just like SARAH`Vertex[], but it caches the results *)
(* MemoizingVertex[] only works when __no__ indices are specified!!! *)
(* Use of memoization gives 30% speedup for the MSSM! *)
memoizedVertices = {};
MemoizingVertex[particles_List, options : OptionsPattern[SARAH`Vertex]] :=
    Module[{memo, ordering, orderedParticles},
           (* First we sort the particles *)
           ordering = Ordering[particles];
           orderedParticles = particles[[ordering]];
           
           memo = Select[memoizedVertices, MatchesMemoizedVertex[orderedParticles], 1];
           If[memo =!= {}, memo = memo[[1]],
              (* Create a new entry *)
              memo = SARAH`Vertex[orderedParticles, options];
              AppendTo[memoizedVertices, memo];];
           
           (* Now return the particles to their original order *)
           memo = OrderVertex[memo, Ordering[ordering]];
           Return[memo]];

MatchesMemoizedVertex[particles_List][vertex_] := MatchQ[particles, Vertices`StripFieldIndices /@ vertex[[1]]];

(* Argument is a SARAH Vertex[] result *)
IsNonZeroVertex[v_] := MemberQ[v[[2 ;;]][[All, 1]], Except[0]];

(********************** End helper routines **************************)

(* The different diagram types that should be taken into consideration *)
(* They need to be called DIAGRAMTYPENAME[_Integer]! See CreateDiagramClasses[] below. *)
(* And do not try to cause an integer overflow in the c++ converted code... *)
contributingFeynmanDiagramTypes = {
    OneLoopDiagram[3], (* Photon is emitted by a fermion, exchange particle is a scalar *)
    OneLoopDiagram[4]  (* Photon is emitted by a scalar, exchange particle is a fermion *)
};

(* This is just a convenient way to help ContributingDiagramsOfType[] *)
OneLoopDiagram[3][fermions_, scalars_, vectors_] := {fermions, scalars};
OneLoopDiagram[4][fermions_, scalars_, vectors_] := {scalars, fermions};

(* Find all diagrams of the type type_, testing all corresponding combinations of particles *)
(* For now only LoopDiagrams 3 through 4 (see further above) are supported.
 If you are adding more diagram types, you should probably make a new overload
 of ContributingDiagramsOfType[] instead of extending this one. *)
(* IMPORTANT: Return value should be a list of Diagram[DIAGRAMTYPENAME[_Integer], Particles___]
 This is important for the c++ conversion that assumes every argument after the type
 is a particle and uses ParticleToString for conversion *)

ContributingDiagramsOfType[type : (OneLoopDiagram[3] | OneLoopDiagram[4]), fermions_, scalars_, vectors_] :=
    Module[{photonEmitters, exchangeParticles, photonVertices, muonVertices, test},
           (* Get the photon emitter and the exchange particle categories corresponding to the
            diagram type *)
           {photonEmitters, exchangeParticles} = type[fermions, scalars, vectors];
           
           (* For every potential photon emitter, check whether it actually can emit a photon *)
           photonVertices = (MemoizingVertex[{GetPhoton[], #, AntiParticle[#]},
                              SARAH`UseDependences -> True,
                              SARAH`Eigenstates -> FlexibleSUSY`FSEigenstates] &) /@ photonEmitters;
           photonVertices = Select[photonVertices, IsNonZeroVertex];
           
           (* From SARAH's more or less cryptically formatted result extract the particles' names *)
           photonEmitters = (Vertices`StripFieldIndices /@ photonVertices)[[All, 1]][[All, 2]];
           
           (* Since we do not know anything about the particles, we have to include their
            corresponding antiparticles as well *)
           photonEmitters = Union[photonEmitters, AntiParticle[photonEmitters]];
           exchangeParticles = Union[exchangeParticles, AntiParticle[exchangeParticles]];
           
           (* Now we check which of the photon emitting particles actually can interact with
            a muon in the way we want. *)
           muonVertices = Outer[(MemoizingVertex[{GetMuonFamily[], #1, #2},
                                                 SARAH`UseDependences -> True,
                                                 SARAH`Eigenstates -> FlexibleSUSY`FSEigenstates] &),
                                photonEmitters, exchangeParticles];
           muonVertices = Select[#, IsNonZeroVertex] & /@ muonVertices;
           muonVertices = Cases[muonVertices, Except[{}]];
           
           (* We return the antiparticles of the particles we just found to *)
           (* This is just a convention and nothing serious. The returned
            particles are the decay products of the muon:
            i.e. if muon ---> p1 + p2
            Then p1 and p2 are returned *)
           test = (Diagram[type, AntiParticle[#[[1]]], AntiParticle[#[[2]]]] &) /@
                (Vertices`StripFieldIndices /@ # &) /@ Flatten[muonVertices[[All, All, 1, 2 ;;]], 1]
           ];

createdVertices = Null;
nPointFunctions = Null;
CreateVertices[] := Module[{contributingDiagrams, vertices,
                            vertexClassesPrototypes, vertexClassesDefinitions},
                           If[createdVertices =!= Null, Return[createdVertices]];
                           
                           contributingDiagrams = ContributingDiagrams[];
                           
                           vertices = Flatten[VerticesForDiagram /@ contributingDiagrams, 1];
                           vertices = DeleteDuplicates[vertices];
                           
                           {vertexClassesPrototypes, vertexClassesDefinitions} = Transpose @ (CreateVertexFunction /@ vertices);
                           vertexClassesPrototypes = Cases[vertexClassesPrototypes, Except[""]];
                           vertexClassesDefinitions = Cases[vertexClassesDefinitions, Except[""]];
                           
                           AppendTo[vertices, {GetPhoton[], GetMuonFamily[], SARAH`bar[GetMuonFamily[]]}];
                           vertices = (OrderParticles[#, Ordering[(Vertices`StripFieldIndices /@ #)]] &) /@ vertices;
                           vertices = DeleteDuplicates[vertices,
                                                       (Vertices`StripFieldIndices[#1] === Vertices`StripFieldIndices[#2] &)];
                           
                           nPointFunctions = Flatten[(Module[{couplings, particles, temp},
                                                             couplings = {ReplacePart[#, 0 -> SARAH`Cp]};
                                                             
                                                             particles = Vertices`StripFieldIndices /@ #;
                                                             (* FIXME: Diagram types 5-6 not supported *)
                                                             If[MemberQ[TreeMasses`IsFermion /@ particles, True],
                                                                couplings = {couplings[[1]][SARAH`PL],
                                                                    couplings[[1]][SARAH`PR]}];
                                                             
                                                             (Null[Null, #] &) /@ couplings
                                                             ] &) /@ vertices];
                           
                           createdVertices = {vertexClassesPrototypes, vertexClassesDefinitions};
                           createdVertices = (StringJoin @ Riffle[#, "\n\n"] &) /@ createdVertices;
                           
                           Return[createdVertices];
                           ];

VerticesForDiagram[Diagram[loopDiagram_OneLoopDiagram, photonEmitter_, exchangeParticle_]] :=
    Module[{photonVertex, muonVertex},
           photonVertex = MemoizingVertex[{GetPhoton[], photonEmitter, AntiParticle[photonEmitter]}];
           muonVertex = MemoizingVertex[{GetMuonFamily[], AntiParticle[photonEmitter], AntiParticle[exchangeParticle]}];
           
           photonVertex = StripLorentzIndices @ photonVertex[[1]];
           muonVertex = StripLorentzIndices @ muonVertex[[1]];

           Return[{photonVertex, muonVertex}];
           ];

vertexFunctions = {};
CreateVertexFunction[indexedParticles_List] :=
    (Module[{prototypes, definitions = "", ordering, particles, orderedParticles,
             orderedIndexedParticles, addSpacing = True},
            particles = Vertices`StripFieldIndices /@ indexedParticles;
            If[MemberQ[vertexFunctions, particles], Return[{"",""}]];
            
            ordering = Ordering[particles];
            orderedParticles = particles[[ordering]];
            orderedIndexedParticles = OrderParticles[indexedParticles, ordering];
            
            If[MemberQ[vertexFunctions, orderedParticles] === True,
               (* There is already an entry *)
               prototypes = "";
               addSpacing = False,
               (* There is no entry yet, create it *)
               {prototypes, definitions} = CreateOrderedVertexFunction[orderedIndexedParticles];
               AppendTo[vertexFunctions, orderedParticles];
               ];
            
            If[ordering === Table[i, {i, 1, Length[ordering]}],
               Return[{prototypes, definitions}]];
            
            orderedVertexFunction = ("VertexFunction<" <>
                                     StringJoin @ Riffle[ParticleToString /@ orderedParticles, ", "] <>
                                     ">");
            
            prototypes = (prototypes <> If[addSpacing, "\n\n", ""] <>
                          "template<> struct VertexFunctionData<" <>
                          StringJoin @ Riffle[ParticleToString /@ particles, ", "] <>
                          ">\n" <>
                          "{\n" <>
                          IndentText @
                          ("static const bool is_permutation = true;\n" <>
                           "using orig_type = " <> orderedVertexFunction <> ";\n" <>
                           "using particlePermutation = permutation<" <>
                           StringJoin @ Riffle[ToString /@ (Ordering[ordering] - 1), ", "] <>
                           ">;\n"
                           ) <>
                          "};");
            
            AppendTo[vertexFunctions, particles];
            Return[{prototypes, definitions}];
            ];);

(* ParsedVertex structure:
    ParsedVertex[
        {numP1Indices, numP2Indices, ...},
        {{minIndex1, minIndex2, ...}, {maxIndex1+1, maxIndex2+1, ...}},
        VertexClassName,
        VertexFunctionBody
    ]
 
 Getters are available!
 *)

NameOfCouplingFunction[particles_List] :=
    ((* FIXME: Not upwards compatible if naming conventions change *)
     "Cp" <> StringReplace[StringJoin @ (ParticleToString /@ Sort[particles]), "<" | ">" | (CXXNamespace <> "::") | "::type" -> ""]);

ParseVertex[indexedParticles_List] :=
    Module[{particles, numberOfIndices, indexParameters,
        parsedVertex, vertexClassName, vertexFunctionBody,
        sarahParticles, particleInfo, indexBounds},
           numberOfIndices = ((Length @ FieldIndexList[#] &) /@ indexedParticles);
           particles = Vertices`StripFieldIndices /@ indexedParticles;
           
           indexParameters = StringJoin @ Riffle[Table["indices[" <> ToString[i] <> "]",
                                                       {i, 0, Total[numberOfIndices] - 1}],
                                                 ", "];
           If[indexParameters =!= "", indexParameters = " " <> indexParameters <> " "];
           
           (* FIXME: Diagram types 5-6 not supported *)
           If[MemberQ[TreeMasses`IsFermion /@ particles, True],
              vertexClassName = "LeftAndRightComponentedVertex";
              
              vertexFunctionBody = ("std::complex<double> left = context.model." <>
                                    NameOfCouplingFunction[particles] <> "PL" <>
                                    "(" <> indexParameters <> ");\n" <>
                                    "std::complex<double> right = context.model." <>
                                    NameOfCouplingFunction[particles] <> "PR" <>
                                    "(" <> indexParameters <> ");\n\n" <>
                                    
                                    "return vertex_type( left, right );"),
              vertexClassName = "SingleComponentedVertex";
              vertexFunctionBody = ("return vertex_type( context.model." <>
                                    NameOfCouplingFunction[particles] <>
                                    "(" <> indexParameters <> ") );")];
           
           sarahParticles = SARAH`getParticleName /@ particles;
           particleInfo = Flatten[(Cases[SARAH`Particles[FlexibleSUSY`FSEigenstates], {#, ___}] &) /@
                                  sarahParticles, 1];
           
           (* INFO: I do not think this ever occurs... *)
           particleInfo = DeleteCases[particleInfo, {SARAH`generation, 1}, {3}];
           particleInfo = DeleteCases[particleInfo, {SARAH`lorentz, _}, {3}];
           
           indexBounds = (With[{particleIndex = #},
                               (If[#[[1]] === SARAH`generation,
                                   {particleInfo[[particleIndex, 2]]-1, particleInfo[[particleIndex, 3]]},
                                   {1, #[[2]]}]
                                &) /@ particleInfo[[particleIndex, 5]]]
                          &) /@ Table[i, {i, Length[particles]}];
           indexBounds = Cases[Flatten[indexBounds, 1], Except[{}]];
           
           If[indexBounds === {},
              indexBounds = {{},{}},
              indexBounds = Transpose @ indexBounds];
           
           parsedVertex = ParsedVertex[numberOfIndices,
                                       indexBounds,
                                       vertexClassName,
                                       vertexFunctionBody];
           
           Return[parsedVertex];
           ];

NumberOfIndices[parsedVertex_ParsedVertex] := Total[parsedVertex[[1]]];
NumberOfIndices[parsedVertex_ParsedVertex, pIndex_Integer] := parsedVertex[[1, pIndex]];

IndexBounds[parsedVertex_ParsedVertex] := parsedVertex[[2]];

VertexClassName[parsedVertex_ParsedVertex] := parsedVertex[[3]];
VertexFunctionBody[parsedVertex_ParsedVertex] := parsedVertex[[4]];

CreateOrderedVertexFunction[orderedIndexedParticles_List] :=
    (Module[{prototype, definition, orderedParticles, dataClassName, functionClassName,
        parsedVertex, particleIndexStartF, particleIndexStart, indexBounds},
            orderedParticles = Vertices`StripFieldIndices /@ orderedIndexedParticles;
            parsedVertex = ParseVertex[orderedIndexedParticles];
            dataClassName = "VertexFunctionData<" <> StringJoin @ Riffle[ParticleToString /@ orderedParticles, ", "] <> ">";
            functionClassName = "VertexFunction<" <> StringJoin @ Riffle[ParticleToString /@ orderedParticles, ", "] <> ">";
            
            prototype = ("template<> struct " <> dataClassName <> "\n" <>
                         "{\n" <>
                         IndentText @
                         ("static const bool is_permutation = false;\n\n" <>
                          
                          "using indices_type = std::array<int, " <> ToString @ NumberOfIndices[parsedVertex] <> ">;\n" <>
                          "using index_bounds = IndexBounds<" <> ToString @ NumberOfIndices[parsedVertex] <> ">;\n" <>
                          "using vertex_type = " <> VertexClassName[parsedVertex] <> ";\n\n" <>
                          
                          "static const unsigned int particleIndexStart[" <> ToString[Length[orderedParticles]+1] <> "];\n" <>
                          "static const index_bounds indexB;\n"
                          ) <>
                         "};");
            
            particleIndexStartF[1] = 0;
            particleIndexStartF[pIndex_] := particleIndexStartF[pIndex-1] + NumberOfIndices[parsedVertex, pIndex-1];
            particleIndexStartF[Length[orderedParticles]+1] = NumberOfIndices[parsedVertex];
            
            particleIndexStart = Table[particleIndexStartF[i], {i, 1, Length[orderedParticles] + 1}];
            
            indexBounds = IndexBounds[parsedVertex];
            
            prototype = (prototype <> "\n" <>
                         "const unsigned int " <> dataClassName <> "::particleIndexStart[" <> ToString[Length[orderedParticles]+1]
                         <> "] = { " <> StringJoin @ Riffle[ToString /@ particleIndexStart, ", "] <> " };");
            
            If[NumberOfIndices[parsedVertex] =!= 0,
               prototype = (prototype <> "\n" <>
                            "const " <> dataClassName <> "::index_bounds " <> dataClassName <> "::indexB = { " <>
                            "{ " <> StringJoin @ Riffle[ToString /@ indexBounds[[1]], ", "] <> " }, " <>
                            "{ " <> StringJoin @ Riffle[ToString /@ indexBounds[[2]], ", "] <> " } };"
                            );];
            definition = ("template<> template<> " <> functionClassName <> "::vertex_type\n" <>
                          functionClassName <> "::vertex( const indices_type &indices, const EvaluationContext &context )\n" <>
                          "{\n" <>
                          IndentText @ VertexFunctionBody[parsedVertex] <> "\n" <>
                          "}");
            
            Return[{prototype, definition}];
            ];);

(* Find all contributing diagrams *)
cachedContributingDiagrams = Null;
ContributingDiagrams[] := Module[{particles, fermions, scalars, vectors},
                                 If[cachedContributingDiagrams =!= Null, Return[cachedContributingDiagrams]];
                                 
                                 particles = TreeMasses`GetParticles[];
                                 fermions = Select[particles, TreeMasses`IsFermion];
                                 scalars = Select[particles, TreeMasses`IsScalar];
                                 vectors = Select[particles, TreeMasses`IsVector];
                                 
                                 cachedContributingDiagrams = Flatten[(ContributingDiagramsOfType[#, fermions, scalars, vectors] &)
                                                                      /@ contributingFeynmanDiagramTypes
                                                                      , 1];
                                 Return[cachedContributingDiagrams];
                                 ];

(* Returns a list of all concrete diagram evaluators (as strings)
 e.g. "DiagramEvaluator<OneLoopDiagram<1>, Fe, VP>"
 that need to be invoked in our calculation *)
ConcreteDiagramEvaluators[] := (("DiagramEvaluator<" <> SymbolName @ Head @ #[[1]] <> "<" <>
                                 ToString @ #[[1,1]] <> ">, " <>
                                 StringJoin @ (Riffle[ParticleToString /@ ReplacePart[#[[2;;]], 0 -> List], ", "]) <>
                                 ">" &)
                                /@ ContributingDiagrams[]);

End[];

EndPackage[];