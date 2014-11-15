
FSModelName = "MSSMNoFV";

(* CMSSM input parameters *)

MINPAR = {
    {3, TanBeta},
    {4, Sign[\[Mu]]}
};

EXTPAR = { {0, Qin},
           {1, M1},
           {2, M2},
           {3, M3},
           {11, AtIN},
           {12, AbIN},
           {13, AtauIN},
           {14, AcIN},
           {15, AsIN},
           {16, AmuonIN},
           {17, AuIN},
           {18, AdIN},
           {19, AeIN},
           {21, mHd2IN},
           {22, mHu2IN},
           {31, ml11IN},
           {32, ml22IN},
           {33, ml33IN},
           {34, me11IN},
           {35, me22IN},
           {36, me33IN},
           {41, mq11IN},
           {42, mq22IN},
           {43, mq33IN},
           {44, mu11IN},
           {45, mu22IN},
           {46, mu33IN},
           {47, md11IN},
           {48, md22IN},
           {49, md33IN} };

EWSBOutputParameters = { B[\[Mu]], \[Mu] };

SUSYScale = Sqrt[M[St[1]] M[St[2]]];

SUSYScaleFirstGuess = SM[MZ];

SUSYScaleInput = {};

HighScale == Qin;

HighScaleFirstGuess = 2.0 10^16;

HighScaleInput = {
   {T[Ye][1,1], AeIN*Ye[1,1]},
   {T[Ye][2,2], AmuonIN*Ye[2,2]},
   {T[Ye][3,3], AtauIN*Ye[3,3]},
   {T[Yd][1,1], AdIN*Yd[1,1]},
   {T[Yd][2,2], AsIN*Yd[2,2]},
   {T[Yd][3,3], AbIN*Yd[3,3]},
   {T[Yu][1,1], AuIN*Yu[1,1]},
   {T[Yu][2,2], AcIN*Yu[2,2]},
   {T[Yu][3,3], AtIN*Yu[3,3]},
   {mHd2, mHd2IN},
   {mHu2, mHu2IN},
   {mq2[1,1],mq11IN^2},
   {mq2[2,2],mq22IN^2},
   {mq2[3,3],mq33IN^2},
   {ml2[1,1],ml11IN^2},
   {ml2[2,2],ml22IN^2},
   {ml2[3,3],ml33IN^2},
   {md2[1,1],md11IN^2},
   {md2[2,2],md22IN^2},
   {md2[3,3],md33IN^2},
   {mu2[1,1],mu11IN^2},
   {mu2[2,2],mu22IN^2},
   {mu2[3,3],mu33IN^2},
   {me2[1,1],me11IN^2},
   {me2[2,2],me22IN^2},
   {me2[3,3],me33IN^2},
   {MassB, M1},
   {MassWB,M2},
   {MassG,M3}
};

LowScale = SM[MZ];

LowScaleFirstGuess = SM[MZ];

LowScaleInput = {
   {Yu, Automatic},
   {Yd, Automatic},
   {Ye, Automatic},
   {vd, 2 MZDRbar / Sqrt[GUTNormalization[g1]^2 g1^2 + g2^2] Cos[ArcTan[TanBeta]]},
   {vu, 2 MZDRbar / Sqrt[GUTNormalization[g1]^2 g1^2 + g2^2] Sin[ArcTan[TanBeta]]}
};

InitialGuessAtLowScale = {
   {vd, SM[vev] Cos[ArcTan[TanBeta]]},
   {vu, SM[vev] Sin[ArcTan[TanBeta]]},
   {Yu, Automatic},
   {Yd, Automatic},
   {Ye, Automatic}
};

InitialGuessAtHighScale = {
   {\[Mu]   , 1.0},
   {B[\[Mu]], 0.0}
};

UseHiggs2LoopMSSM = True;
EffectiveMu = \[Mu];

ExtraSLHAOutputBlocks = {
   {ALPHA, {{ArcSin[Pole[ZH[2,2]]]}}},
   {HMIX , {{1, \[Mu]},
            {2, vu / vd},
            {3, Sqrt[vu^2 + vd^2]},
            {4, M[Ah[2]]^2},
            {101, B[\[Mu]]},
            {102, vd},
            {103, vu} } },
   {Au,    {{1, 1, T[Yu][1,1] / Yu[1,1]},
            {2, 2, T[Yu][2,2] / Yu[2,2]},
            {3, 3, T[Yu][3,3] / Yu[3,3]} } },
   {Ad,    {{1, 1, T[Yd][1,1] / Yd[1,1]},
            {2, 2, T[Yd][2,2] / Yd[2,2]},
            {3, 3, T[Yd][3,3] / Yd[3,3]} } },
   {Ae,    {{1, 1, T[Ye][1,1] / Ye[1,1]},
            {2, 2, T[Ye][2,2] / Ye[2,2]},
            {3, 3, T[Ye][3,3] / Ye[3,3]} } },
   {MSOFT, {{1, MassB},
            {2, MassWB},
            {3, MassG},
            {21, mHd2},
            {22, mHu2},
            {31, Sqrt[ml2[1,1]]},
            {32, Sqrt[ml2[2,2]]},
            {33, Sqrt[ml2[3,3]]},
            {34, Sqrt[me2[1,1]]},
            {35, Sqrt[me2[2,2]]},
            {36, Sqrt[me2[3,3]]},
            {41, Sqrt[mq2[1,1]]},
            {42, Sqrt[mq2[2,2]]},
            {43, Sqrt[mq2[3,3]]},
            {44, Sqrt[mu2[1,1]]},
            {45, Sqrt[mu2[2,2]]},
            {46, Sqrt[mu2[3,3]]},
            {47, Sqrt[md2[1,1]]},
            {48, Sqrt[md2[2,2]]},
            {49, Sqrt[md2[3,3]]} } },
   {MASS,  {{1000021, Pole[M[Glu]]    },
            {1000012, Pole[M[SveL]]   },
            {1000014, Pole[M[SvmL]]   },
            {1000016, Pole[M[SvtL]]   },
            {1000024, Pole[M[Cha[1]]] },
            {1000037, Pole[M[Cha[2]]] },
            {     25, Pole[M[hh[1]]]  },
            {     35, Pole[M[hh[2]]]  },
            {     37, Pole[M[Hpm[2]]] },
            {     36, Pole[M[Ah[2]]]  },
            {1000001, Sum[conj[Pole[ZD[i,1]]] Pole[ZD[i,1]] Pole[M[Sd[i]]], {i,1,2}] }, (* ~d_L *)
            {2000001, Sum[conj[Pole[ZD[i,2]]] Pole[ZD[i,2]] Pole[M[Sd[i]]], {i,1,2}] }, (* ~d_R *)
            {1000003, Sum[conj[Pole[ZS[i,1]]] Pole[ZS[i,1]] Pole[M[Ss[i]]], {i,1,2}] }, (* ~s_L *)
            {2000003, Sum[conj[Pole[ZS[i,2]]] Pole[ZS[i,2]] Pole[M[Ss[i]]], {i,1,2}] }, (* ~s_R *)
            {1000005, Pole[M[Sb[1]]]  },
            {2000005, Pole[M[Sb[2]]]  },
            {1000011, Sum[conj[Pole[ZE[i,1]]] Pole[ZE[i,1]] Pole[M[Se[i]]], {i,1,2}] }, (* ~e_L *)
            {2000011, Sum[conj[Pole[ZE[i,2]]] Pole[ZE[i,2]] Pole[M[Se[i]]], {i,1,2}] }, (* ~e_R *)
            {1000013, Sum[conj[Pole[ZM[i,1]]] Pole[ZM[i,1]] Pole[M[Sm[i]]], {i,1,2}] }, (* ~m_L *)
            {2000013, Sum[conj[Pole[ZM[i,2]]] Pole[ZM[i,2]] Pole[M[Sm[i]]], {i,1,2}] }, (* ~m_R *)
            {1000015, Pole[M[Stau[1]]]},
            {2000015, Pole[M[Stau[2]]]},
            {1000002, Sum[conj[Pole[ZU[i,1]]] Pole[ZU[i,1]] Pole[M[Su[i]]], {i,1,2}] }, (* ~u_L *)
            {2000002, Sum[conj[Pole[ZU[i,2]]] Pole[ZU[i,2]] Pole[M[Su[i]]], {i,1,2}] }, (* ~u_R *)
            {1000004, Sum[conj[Pole[ZC[i,1]]] Pole[ZC[i,1]] Pole[M[Sc[i]]], {i,1,2}] }, (* ~c_L *)
            {2000004, Sum[conj[Pole[ZC[i,2]]] Pole[ZC[i,2]] Pole[M[Sc[i]]], {i,1,2}] }, (* ~c_R *)
            {1000006, Pole[M[St[1]]]  },
            {2000006, Pole[M[St[2]]]  },
            {1000022, Pole[M[Chi[1]]] },
            {1000023, Pole[M[Chi[2]]] },
            {1000025, Pole[M[Chi[3]]] },
            {1000035, Pole[M[Chi[4]]] },
            {     24, SM[MW]          } } }
};

TopQuark    = Ft;
BottomQuark = Fb;
Neutrino    = Fvt;
Electron    = Ftau;

SMParticles = {
    VectorP, VectorZ, VectorG, VectorW,
    Fu, Fd, Fc, Fs, Ft, Fb, Fe, Fm, Ftau, Fve, Fvm, Fvt
};

SA`Casimir[Ss, SARAH`color] = 4/3;
SA`Casimir[Sc, SARAH`color] = 4/3;
SA`Casimir[Sb, SARAH`color] = 4/3;
SA`Casimir[St, SARAH`color] = 4/3;
