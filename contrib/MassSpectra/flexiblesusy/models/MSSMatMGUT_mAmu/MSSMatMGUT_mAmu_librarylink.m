Print["================================"];
Print["FlexibleSUSY 2.0.1"];
Print["MSSMatMGUT_mAmu"];
Print["http://flexiblesusy.hepforge.org"];
Print["================================"];

libMSSMatMGUT_mAmu = FileNameJoin[{Directory[], "models", "MSSMatMGUT_mAmu", "MSSMatMGUT_mAmu_librarylink.so"}];

FSMSSMatMGUT_mAmuGetSettings = LibraryFunctionLoad[libMSSMatMGUT_mAmu, "FSMSSMatMGUT_mAmuGetSettings", LinkObject, LinkObject];
FSMSSMatMGUT_mAmuGetSMInputParameters = LibraryFunctionLoad[libMSSMatMGUT_mAmu, "FSMSSMatMGUT_mAmuGetSMInputParameters", LinkObject, LinkObject];
FSMSSMatMGUT_mAmuGetInputParameters = LibraryFunctionLoad[libMSSMatMGUT_mAmu, "FSMSSMatMGUT_mAmuGetInputParameters", LinkObject, LinkObject];
FSMSSMatMGUT_mAmuGetProblems = LibraryFunctionLoad[libMSSMatMGUT_mAmu, "FSMSSMatMGUT_mAmuGetProblems", LinkObject, LinkObject];
FSMSSMatMGUT_mAmuGetWarnings = LibraryFunctionLoad[libMSSMatMGUT_mAmu, "FSMSSMatMGUT_mAmuGetWarnings", LinkObject, LinkObject];
FSMSSMatMGUT_mAmuToSLHA = LibraryFunctionLoad[libMSSMatMGUT_mAmu, "FSMSSMatMGUT_mAmuToSLHA", LinkObject, LinkObject];

FSMSSMatMGUT_mAmuOpenHandleLib = LibraryFunctionLoad[libMSSMatMGUT_mAmu, "FSMSSMatMGUT_mAmuOpenHandle", {{Real,1}}, Integer];
FSMSSMatMGUT_mAmuCloseHandle = LibraryFunctionLoad[libMSSMatMGUT_mAmu, "FSMSSMatMGUT_mAmuCloseHandle", {Integer}, Void];

FSMSSMatMGUT_mAmuSetLib = LibraryFunctionLoad[libMSSMatMGUT_mAmu, "FSMSSMatMGUT_mAmuSet", {Integer, {Real,1}}, Void];

FSMSSMatMGUT_mAmuCalculateSpectrum = LibraryFunctionLoad[libMSSMatMGUT_mAmu, "FSMSSMatMGUT_mAmuCalculateSpectrum", LinkObject, LinkObject];
FSMSSMatMGUT_mAmuCalculateObservables = LibraryFunctionLoad[libMSSMatMGUT_mAmu, "FSMSSMatMGUT_mAmuCalculateObservables", LinkObject, LinkObject];

FSMSSMatMGUT_mAmuCalculateSpectrum::error = "`1`";
FSMSSMatMGUT_mAmuCalculateSpectrum::warning = "`1`";

FSMSSMatMGUT_mAmuCalculateObservables::error = "`1`";
FSMSSMatMGUT_mAmuCalculateObservables::warning = "`1`";

FSMSSMatMGUT_mAmu::info = "`1`";
FSMSSMatMGUT_mAmu::nonum = "Error: `1` is not a numeric input value!";
FSMSSMatMGUT_mAmuMessage[s_] := Message[FSMSSMatMGUT_mAmu::info, s];

FSMSSMatMGUT_mAmuCheckIsNumeric[a_?NumericQ] := a;
FSMSSMatMGUT_mAmuCheckIsNumeric[a_] := (Message[FSMSSMatMGUT_mAmu::nonum, a]; Abort[]);

fsDefaultSettings = {
      precisionGoal -> 1.*^-4,           (* FlexibleSUSY[0] *)
      maxIterations -> 0,                (* FlexibleSUSY[1] *)
      solver -> 1,     (* FlexibleSUSY[2] *)
      calculateStandardModelMasses -> 0, (* FlexibleSUSY[3] *)
      poleMassLoopOrder -> 2,            (* FlexibleSUSY[4] *)
      ewsbLoopOrder -> 2,                (* FlexibleSUSY[5] *)
      betaFunctionLoopOrder -> 3,        (* FlexibleSUSY[6] *)
      thresholdCorrectionsLoopOrder -> 2,(* FlexibleSUSY[7] *)
      higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
      higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
      higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
      higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
      forceOutput -> 0,                  (* FlexibleSUSY[12] *)
      topPoleQCDCorrections -> 1,        (* FlexibleSUSY[13] *)
      betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
      forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
      poleMassScale -> 0,                (* FlexibleSUSY[17] *)
      eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
      eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
      eftMatchingLoopOrderUp -> 2,       (* FlexibleSUSY[20] *)
      eftMatchingLoopOrderDown -> 1,     (* FlexibleSUSY[21] *)
      eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
      calculateBSMMasses -> 1,           (* FlexibleSUSY[23] *)
      thresholdCorrections -> 123111321, (* FlexibleSUSY[24] *)
      higgs3loopCorrectionRenScheme -> 0,(* FlexibleSUSY[25] *)
      higgs3loopCorrectionAtAsAs -> 1,   (* FlexibleSUSY[26] *)
      higgs3loopCorrectionAbAsAs -> 1,   (* FlexibleSUSY[27] *)
      higgs3loopCorrectionAtAtAs -> 1,   (* FlexibleSUSY[28] *)
      higgs3loopCorrectionAtAtAt -> 1,   (* FlexibleSUSY[29] *)
      parameterOutputScale -> 0          (* MODSEL[12] *)
};

fsDefaultSMParameters = {
    alphaEmMZ -> 1/127.916, (* SMINPUTS[1] *)
    GF -> 1.1663787*^-5,    (* SMINPUTS[2] *)
    alphaSMZ -> 0.1184,     (* SMINPUTS[3] *)
    MZ -> 91.1876,          (* SMINPUTS[4] *)
    mbmb -> 4.18,           (* SMINPUTS[5] *)
    Mt -> 173.34,           (* SMINPUTS[6] *)
    Mtau -> 1.777,          (* SMINPUTS[7] *)
    Mv3 -> 0,               (* SMINPUTS[8] *)
    MW -> 80.385,           (* SMINPUTS[9] *)
    Me -> 0.000510998902,   (* SMINPUTS[11] *)
    Mv1 -> 0,               (* SMINPUTS[12] *)
    Mm -> 0.1056583715,     (* SMINPUTS[13] *)
    Mv2 -> 0,               (* SMINPUTS[14] *)
    md2GeV -> 0.00475,      (* SMINPUTS[21] *)
    mu2GeV -> 0.0024,       (* SMINPUTS[22] *)
    ms2GeV -> 0.104,        (* SMINPUTS[23] *)
    mcmc -> 1.27,           (* SMINPUTS[24] *)
    CKMTheta12 -> 0,
    CKMTheta13 -> 0,
    CKMTheta23 -> 0,
    CKMDelta -> 0,
    PMNSTheta12 -> 0,
    PMNSTheta13 -> 0,
    PMNSTheta23 -> 0,
    PMNSDelta -> 0,
    PMNSAlpha1 -> 0,
    PMNSAlpha2 -> 0,
    alphaEm0 -> 1/137.035999074,
    Mh -> 125.09
};

fsMSSMatMGUT_mAmuDefaultInputParameters = {
   TanBeta -> 0,
   SignMu -> 0,
   MuInput -> 0,
   mA2Input -> 0,
   Aeij -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   Adij -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   Auij -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   mq2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   ml2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   md2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   mu2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   me2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   MassBInput -> 0,
   MassWBInput -> 0,
   MassGInput -> 0
};

Options[FSMSSMatMGUT_mAmuOpenHandle] = {
    Sequence @@ fsDefaultSettings,
    Sequence @@ fsDefaultSMParameters,
    Sequence @@ fsMSSMatMGUT_mAmuDefaultInputParameters
};

FSMSSMatMGUT_mAmuOpenHandle[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSMSSMatMGUT_mAmuOpenHandle[a, Sequence @@ s, r];

FSMSSMatMGUT_mAmuOpenHandle[OptionsPattern[]] :=
    FSMSSMatMGUT_mAmuOpenHandleLib[
        FSMSSMatMGUT_mAmuCheckIsNumeric /@ {
            (* spectrum generator settings *)
            OptionValue[precisionGoal],
            OptionValue[maxIterations],
            OptionValue[solver],
            OptionValue[calculateStandardModelMasses],
            OptionValue[poleMassLoopOrder],
            OptionValue[ewsbLoopOrder],
            OptionValue[betaFunctionLoopOrder],
            OptionValue[thresholdCorrectionsLoopOrder],
            OptionValue[higgs2loopCorrectionAtAs],
            OptionValue[higgs2loopCorrectionAbAs],
            OptionValue[higgs2loopCorrectionAtAt],
            OptionValue[higgs2loopCorrectionAtauAtau],
            OptionValue[forceOutput],
            OptionValue[topPoleQCDCorrections],
            OptionValue[betaZeroThreshold],
            OptionValue[forcePositiveMasses],
            OptionValue[poleMassScale],
            OptionValue[eftPoleMassScale],
            OptionValue[eftMatchingScale],
            OptionValue[eftMatchingLoopOrderUp],
            OptionValue[eftMatchingLoopOrderDown],
            OptionValue[eftHiggsIndex],
            OptionValue[calculateBSMMasses],
            OptionValue[thresholdCorrections],
            OptionValue[higgs3loopCorrectionRenScheme],
            OptionValue[higgs3loopCorrectionAtAsAs],
            OptionValue[higgs3loopCorrectionAbAsAs],
            OptionValue[higgs3loopCorrectionAtAtAs],
            OptionValue[higgs3loopCorrectionAtAtAt],
            OptionValue[parameterOutputScale],

            (* Standard Model input parameters *)
            OptionValue[alphaEmMZ],
            OptionValue[GF],
            OptionValue[alphaSMZ],
            OptionValue[MZ],
            OptionValue[mbmb],
            OptionValue[Mt],
            OptionValue[Mtau],
            OptionValue[Mv3],
            OptionValue[MW],
            OptionValue[Me],
            OptionValue[Mv1],
            OptionValue[Mm],
            OptionValue[Mv2],
            OptionValue[md2GeV],
            OptionValue[mu2GeV],
            OptionValue[ms2GeV],
            OptionValue[mcmc],
            OptionValue[CKMTheta12],
            OptionValue[CKMTheta13],
            OptionValue[CKMTheta23],
            OptionValue[CKMDelta],
            OptionValue[PMNSTheta12],
            OptionValue[PMNSTheta13],
            OptionValue[PMNSTheta23],
            OptionValue[PMNSDelta],
            OptionValue[PMNSAlpha1],
            OptionValue[PMNSAlpha2],
            OptionValue[alphaEm0],
            OptionValue[Mh]

            (* MSSMatMGUT_mAmu input parameters *)
            ,
            OptionValue[TanBeta],
            OptionValue[SignMu],
            OptionValue[MuInput],
            OptionValue[mA2Input],
            OptionValue[Aeij][[1,1]],
            OptionValue[Aeij][[1,2]],
            OptionValue[Aeij][[1,3]],
            OptionValue[Aeij][[2,1]],
            OptionValue[Aeij][[2,2]],
            OptionValue[Aeij][[2,3]],
            OptionValue[Aeij][[3,1]],
            OptionValue[Aeij][[3,2]],
            OptionValue[Aeij][[3,3]],
            OptionValue[Adij][[1,1]],
            OptionValue[Adij][[1,2]],
            OptionValue[Adij][[1,3]],
            OptionValue[Adij][[2,1]],
            OptionValue[Adij][[2,2]],
            OptionValue[Adij][[2,3]],
            OptionValue[Adij][[3,1]],
            OptionValue[Adij][[3,2]],
            OptionValue[Adij][[3,3]],
            OptionValue[Auij][[1,1]],
            OptionValue[Auij][[1,2]],
            OptionValue[Auij][[1,3]],
            OptionValue[Auij][[2,1]],
            OptionValue[Auij][[2,2]],
            OptionValue[Auij][[2,3]],
            OptionValue[Auij][[3,1]],
            OptionValue[Auij][[3,2]],
            OptionValue[Auij][[3,3]],
            OptionValue[mq2Input][[1,1]],
            OptionValue[mq2Input][[1,2]],
            OptionValue[mq2Input][[1,3]],
            OptionValue[mq2Input][[2,1]],
            OptionValue[mq2Input][[2,2]],
            OptionValue[mq2Input][[2,3]],
            OptionValue[mq2Input][[3,1]],
            OptionValue[mq2Input][[3,2]],
            OptionValue[mq2Input][[3,3]],
            OptionValue[ml2Input][[1,1]],
            OptionValue[ml2Input][[1,2]],
            OptionValue[ml2Input][[1,3]],
            OptionValue[ml2Input][[2,1]],
            OptionValue[ml2Input][[2,2]],
            OptionValue[ml2Input][[2,3]],
            OptionValue[ml2Input][[3,1]],
            OptionValue[ml2Input][[3,2]],
            OptionValue[ml2Input][[3,3]],
            OptionValue[md2Input][[1,1]],
            OptionValue[md2Input][[1,2]],
            OptionValue[md2Input][[1,3]],
            OptionValue[md2Input][[2,1]],
            OptionValue[md2Input][[2,2]],
            OptionValue[md2Input][[2,3]],
            OptionValue[md2Input][[3,1]],
            OptionValue[md2Input][[3,2]],
            OptionValue[md2Input][[3,3]],
            OptionValue[mu2Input][[1,1]],
            OptionValue[mu2Input][[1,2]],
            OptionValue[mu2Input][[1,3]],
            OptionValue[mu2Input][[2,1]],
            OptionValue[mu2Input][[2,2]],
            OptionValue[mu2Input][[2,3]],
            OptionValue[mu2Input][[3,1]],
            OptionValue[mu2Input][[3,2]],
            OptionValue[mu2Input][[3,3]],
            OptionValue[me2Input][[1,1]],
            OptionValue[me2Input][[1,2]],
            OptionValue[me2Input][[1,3]],
            OptionValue[me2Input][[2,1]],
            OptionValue[me2Input][[2,2]],
            OptionValue[me2Input][[2,3]],
            OptionValue[me2Input][[3,1]],
            OptionValue[me2Input][[3,2]],
            OptionValue[me2Input][[3,3]],
            OptionValue[MassBInput],
            OptionValue[MassWBInput],
            OptionValue[MassGInput]
        }
];

Options[FSMSSMatMGUT_mAmuSet] = Options[FSMSSMatMGUT_mAmuOpenHandle];

FSMSSMatMGUT_mAmuSet[handle_Integer, a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSMSSMatMGUT_mAmuSet[handle, a, Sequence @@ s, r];

FSMSSMatMGUT_mAmuSet[handle_Integer, p:OptionsPattern[]] :=
    FSMSSMatMGUT_mAmuSetLib[
        handle,
        ReleaseHold[Hold[FSMSSMatMGUT_mAmuCheckIsNumeric /@ {
            (* spectrum generator settings *)
            OptionValue[precisionGoal],
            OptionValue[maxIterations],
            OptionValue[solver],
            OptionValue[calculateStandardModelMasses],
            OptionValue[poleMassLoopOrder],
            OptionValue[ewsbLoopOrder],
            OptionValue[betaFunctionLoopOrder],
            OptionValue[thresholdCorrectionsLoopOrder],
            OptionValue[higgs2loopCorrectionAtAs],
            OptionValue[higgs2loopCorrectionAbAs],
            OptionValue[higgs2loopCorrectionAtAt],
            OptionValue[higgs2loopCorrectionAtauAtau],
            OptionValue[forceOutput],
            OptionValue[topPoleQCDCorrections],
            OptionValue[betaZeroThreshold],
            OptionValue[forcePositiveMasses],
            OptionValue[poleMassScale],
            OptionValue[eftPoleMassScale],
            OptionValue[eftMatchingScale],
            OptionValue[eftMatchingLoopOrderUp],
            OptionValue[eftMatchingLoopOrderDown],
            OptionValue[eftHiggsIndex],
            OptionValue[calculateBSMMasses],
            OptionValue[thresholdCorrections],
            OptionValue[higgs3loopCorrectionRenScheme],
            OptionValue[higgs3loopCorrectionAtAsAs],
            OptionValue[higgs3loopCorrectionAbAsAs],
            OptionValue[higgs3loopCorrectionAtAtAs],
            OptionValue[higgs3loopCorrectionAtAtAt],
            OptionValue[parameterOutputScale],

            (* Standard Model input parameters *)
            OptionValue[alphaEmMZ],
            OptionValue[GF],
            OptionValue[alphaSMZ],
            OptionValue[MZ],
            OptionValue[mbmb],
            OptionValue[Mt],
            OptionValue[Mtau],
            OptionValue[Mv3],
            OptionValue[MW],
            OptionValue[Me],
            OptionValue[Mv1],
            OptionValue[Mm],
            OptionValue[Mv2],
            OptionValue[md2GeV],
            OptionValue[mu2GeV],
            OptionValue[ms2GeV],
            OptionValue[mcmc],
            OptionValue[CKMTheta12],
            OptionValue[CKMTheta13],
            OptionValue[CKMTheta23],
            OptionValue[CKMDelta],
            OptionValue[PMNSTheta12],
            OptionValue[PMNSTheta13],
            OptionValue[PMNSTheta23],
            OptionValue[PMNSDelta],
            OptionValue[PMNSAlpha1],
            OptionValue[PMNSAlpha2],
            OptionValue[alphaEm0],
            OptionValue[Mh]

            (* MSSMatMGUT_mAmu input parameters *)
            ,
            OptionValue[TanBeta],
            OptionValue[SignMu],
            OptionValue[MuInput],
            OptionValue[mA2Input],
            OptionValue[Aeij][[1,1]],
            OptionValue[Aeij][[1,2]],
            OptionValue[Aeij][[1,3]],
            OptionValue[Aeij][[2,1]],
            OptionValue[Aeij][[2,2]],
            OptionValue[Aeij][[2,3]],
            OptionValue[Aeij][[3,1]],
            OptionValue[Aeij][[3,2]],
            OptionValue[Aeij][[3,3]],
            OptionValue[Adij][[1,1]],
            OptionValue[Adij][[1,2]],
            OptionValue[Adij][[1,3]],
            OptionValue[Adij][[2,1]],
            OptionValue[Adij][[2,2]],
            OptionValue[Adij][[2,3]],
            OptionValue[Adij][[3,1]],
            OptionValue[Adij][[3,2]],
            OptionValue[Adij][[3,3]],
            OptionValue[Auij][[1,1]],
            OptionValue[Auij][[1,2]],
            OptionValue[Auij][[1,3]],
            OptionValue[Auij][[2,1]],
            OptionValue[Auij][[2,2]],
            OptionValue[Auij][[2,3]],
            OptionValue[Auij][[3,1]],
            OptionValue[Auij][[3,2]],
            OptionValue[Auij][[3,3]],
            OptionValue[mq2Input][[1,1]],
            OptionValue[mq2Input][[1,2]],
            OptionValue[mq2Input][[1,3]],
            OptionValue[mq2Input][[2,1]],
            OptionValue[mq2Input][[2,2]],
            OptionValue[mq2Input][[2,3]],
            OptionValue[mq2Input][[3,1]],
            OptionValue[mq2Input][[3,2]],
            OptionValue[mq2Input][[3,3]],
            OptionValue[ml2Input][[1,1]],
            OptionValue[ml2Input][[1,2]],
            OptionValue[ml2Input][[1,3]],
            OptionValue[ml2Input][[2,1]],
            OptionValue[ml2Input][[2,2]],
            OptionValue[ml2Input][[2,3]],
            OptionValue[ml2Input][[3,1]],
            OptionValue[ml2Input][[3,2]],
            OptionValue[ml2Input][[3,3]],
            OptionValue[md2Input][[1,1]],
            OptionValue[md2Input][[1,2]],
            OptionValue[md2Input][[1,3]],
            OptionValue[md2Input][[2,1]],
            OptionValue[md2Input][[2,2]],
            OptionValue[md2Input][[2,3]],
            OptionValue[md2Input][[3,1]],
            OptionValue[md2Input][[3,2]],
            OptionValue[md2Input][[3,3]],
            OptionValue[mu2Input][[1,1]],
            OptionValue[mu2Input][[1,2]],
            OptionValue[mu2Input][[1,3]],
            OptionValue[mu2Input][[2,1]],
            OptionValue[mu2Input][[2,2]],
            OptionValue[mu2Input][[2,3]],
            OptionValue[mu2Input][[3,1]],
            OptionValue[mu2Input][[3,2]],
            OptionValue[mu2Input][[3,3]],
            OptionValue[me2Input][[1,1]],
            OptionValue[me2Input][[1,2]],
            OptionValue[me2Input][[1,3]],
            OptionValue[me2Input][[2,1]],
            OptionValue[me2Input][[2,2]],
            OptionValue[me2Input][[2,3]],
            OptionValue[me2Input][[3,1]],
            OptionValue[me2Input][[3,2]],
            OptionValue[me2Input][[3,3]],
            OptionValue[MassBInput],
            OptionValue[MassWBInput],
            OptionValue[MassGInput]
        }] /. HoldPattern[OptionValue[param_]] :> param /.
        { p } /.
        FSMSSMatMGUT_mAmuGetSettings[handle] /.
        FSMSSMatMGUT_mAmuGetSMInputParameters[handle] /.
        FSMSSMatMGUT_mAmuGetInputParameters[handle]]];
