(* ::Package:: *)

(* ==========================================================================
   Predictive scaling of topological tolerance classes
   --------------------------------------------------------------------------
   Reproduces Section 10 of the paper:
     - Eq. (2)-(3): fit of the mean drift mu(M) and of the standard
       deviation sigma(M), trained only on 10 <= M <= 100.
     - Fig. 11 : predicted vs. empirical mu(M) and sigma(M) up to M = 300.
     - Table 15: empirical vs. predicted tolerance bounds, 120 <= M <= 300.
     - Extreme test at M = 2000.

   Requirements: Mathematica 12 or later. No external packages.
   Runtime: about 3 minutes on a standard laptop. The M = 2000 test is the
   slowest single step.

   Note on randomness: no random seed is fixed, so each run draws a fresh set
   of noise realisations. The reported values are therefore reproduced within
   the statistical fluctuation of 1000 Monte Carlo trials.
   ========================================================================== *)

ClearAll["Global`*"];


(* --- 1. Parameters ------------------------------------------------------ *)
noiseLevel = 0.6;        (* noise bound s: noise is uniform in [-s, s]      *)
sigmaGauss = 0.8;        (* standard deviation of each Gaussian pulse       *)
stepRes = 0.1;           (* sampling resolution, in spatial units           *)
dist = 5;                (* distance d between two consecutive peaks        *)
minWidthSamples = 8;     (* W_min = 0.8 units = 8 samples at step 0.1       *)
amplitudeThreshold = 0.5;(* half-maximum level used for the filtration      *)
numTrials = 1000;        (* Monte Carlo trials per configuration            *)
targetAccuracy = 0.95;   (* target accuracy tau                             *)


(* --- 2. Signal generation and topological decoding ---------------------- *)

(* Clean pulse train: M Gaussian pulses of unit amplitude, spaced by d. *)
GenerateCleanSignal[M_Integer, maxSpace_, step_] :=
  Module[{grid, peakPos},
   grid = Range[0, maxSpace, step];
   peakPos = Range[1, M]*dist;
   Total[Table[Exp[-((grid - p)^2)/(2*sigmaGauss^2)], {p, peakPos}]]];

(* Band-limited noise: uniform values on the integer grid, then interpolated
   to the fine sampling grid. This is the channel model of Section 4. *)
GenerateNoise[maxSpace_, numPts_] :=
  ArrayResample[
   RandomReal[{-noiseLevel, noiseLevel}, maxSpace + 1], {numPts}];

(* Estimator of Eq. (1): count the connected components of the super-level
   set at threshold 0.5 that are at least W_min samples wide. *)
DecodeTDA[noisySignal_] :=
  Count[Split[UnitStep[noisySignal - amplitudeThreshold]],
   l : {1 ..} /; Length[l] >= minWidthSamples];

(* One full Monte Carlo run for a given M: returns the list of estimates. *)
RunTrials[M_Integer, trials_Integer] :=
  Module[{maxSpace, numPts, cleanSignal},
   maxSpace = M*dist + 5;
   numPts = Round[maxSpace/stepRes] + 1;
   cleanSignal = GenerateCleanSignal[M, maxSpace, stepRes];
   Table[DecodeTDA[cleanSignal + GenerateNoise[maxSpace, numPts]],
    {trials}]];

(* Data-driven window expansion (Algorithm 2 of Section 11.2). The window
   starts at [M, M] and grows towards the adjacent class with the higher
   empirical frequency, until the coverage reaches the target. *)
FindAsymmetricWindow[estimates_, targetAcc_, M_] :=
 Module[{freq, prob, wMin, wMax, currentP, minV, maxV, pLeft, pRight},
  freq = Counts[estimates];
  prob[k_] := N[Lookup[freq, k, 0]/Length[estimates]];
  wMin = M; wMax = M; currentP = prob[M];
  {minV, maxV} = {Min[Keys[freq]], Max[Keys[freq]]};
  While[currentP < targetAcc,
   If[wMin <= minV && wMax >= maxV, Break[]];
   pLeft = If[wMin > minV, prob[wMin - 1], -1.0];
   pRight = If[wMax < maxV, prob[wMax + 1], -1.0];
   If[pLeft >= pRight,
    wMin--; currentP += prob[wMin],
    wMax++; currentP += prob[wMax]]];
  {wMin, wMax, wMax - wMin + 1}];

FormatBounds[wMin_, wMax_] :=
  "[" <> ToString[wMin] <> ", " <> ToString[wMax] <> "]";


(* --- 3. Empirical data, M = 10 to 300 ----------------------------------- *)
Print["--- Generating experimental data (M = 10 to 300) ... ---"];

allM = Join[Range[10, 100, 10], Range[120, 300, 20]];
meanErrors = Association[];
stdDevs = Association[];
allEstimates = Association[];

Do[
  Module[{estimates, errors},
   estimates = RunTrials[M, numTrials];
   errors = estimates - M;
   allEstimates[M] = estimates;
   meanErrors[M] = Mean[errors];
   stdDevs[M] = StandardDeviation[errors]],
  {M, allM}];


(* --- 4. Calibration, trained only on M <= 100 --------------------------- *)
Print["\n--- Fitting the model (training set: M <= 100) ---"];

trainM = Range[10, 100, 10];
dataMeanTrain = Table[{M, meanErrors[M]}, {M, trainM}];
dataStdTrain = Table[{M, stdDevs[M]}, {M, trainM}];

(* Eq. (2): the mean bias of each pulse adds coherently, so mu grows as M. *)
fitMean = FindFit[dataMeanTrain, a*M, {a}, M];
aFit = a /. fitMean;

(* Eq. (3): the random part adds incoherently, so sigma grows as Sqrt[M]. *)
fitStd = FindFit[dataStdTrain, k*Sqrt[M], {k}, M];
kFit = k /. fitStd;

Print["Mean drift : mu(M) = ", aFit, " * M"];
Print["Std. dev.  : sigma(M) = ", kFit, " * Sqrt[M]"];

(* Eq. (4): predicted tolerance interval at the 95% target (z = 1.96). *)
zScore = 1.96;
PredictBounds[M_] :=
  Module[{center, halfWidth, lower, upper},
   center = M + aFit*M;
   halfWidth = zScore*kFit*Sqrt[M];
   lower = Round[center - halfWidth];
   upper = Round[center + halfWidth];
   {lower, upper, upper - lower + 1}];


(* --- 5. Figure 11: predicted vs. empirical mu(M) and sigma(M) ----------- *)
dataMeanAll = Table[{M, meanErrors[M]}, {M, allM}];
dataStdAll = Table[{M, stdDevs[M]}, {M, allM}];

(* The legend sits inside the frame, so no space is wasted beside the plot.
   Points and line need two separate legend objects. *)
plotLegend = Column[{
    PointLegend[{Red}, {"Empirical data"}, LegendMarkerSize -> 14],
    LineLegend[{Directive[Blue, Dashed, AbsoluteThickness[2]]},
     {"Prediction"}, LegendMarkerSize -> 30]},
   Spacings -> 0];

(* Shared formatting: closed frame, ticks and numbers on the bottom and left
   sides only, no cartesian axes crossing the plot. *)
frameFormat = Sequence[
   Frame -> True,
   Axes -> False,
   FrameStyle -> Directive[Black, AbsoluteThickness[1.2]],
   FrameTicks -> {{Automatic, None}, {Automatic, None}},
   LabelStyle -> {FontSize -> 14, Black},
   ImageSize -> 430];

plotMean = Show[
   ListPlot[dataMeanAll, PlotStyle -> Directive[Red, PointSize[0.018]]],
   Plot[aFit*M, {M, 10, 300},
    PlotStyle -> Directive[Blue, Dashed, AbsoluteThickness[2]]],
   PlotRange -> {{0, 310}, {-2.3, 0.35}},
   (* Faint dashed line at zero: the drift is always negative. *)
   GridLines -> {{}, {{0, Directive[GrayLevel[0.65], Dashed]}}},
   FrameLabel -> {{Style["Mean error (drift)", 16], None},
                  {Style["Number of pulses M", 16], None}},
   Epilog -> Inset[plotLegend, Scaled[{0.97, 0.95}], {Right, Top}],
   frameFormat];

plotStd = Show[
   ListPlot[dataStdAll, PlotStyle -> Directive[Red, PointSize[0.018]]],
   Plot[kFit*Sqrt[M], {M, 10, 300},
    PlotStyle -> Directive[Blue, Dashed, AbsoluteThickness[2]]],
   PlotRange -> {{0, 310}, {0, 5}},
   FrameLabel -> {{Style["Standard deviation", 16], None},
                  {Style["Number of pulses M", 16], None}},
   Epilog -> Inset[plotLegend, Scaled[{0.05, 0.95}], {Left, Top}],
   frameFormat];

combinedPlot = Row[{plotMean, Spacer[25], plotStd}];
Print[combinedPlot];


(* --- 6. Table 15: empirical vs. predicted bounds, 120 <= M <= 300 ------- *)
Print["\n--- Building Table 15 (validation on unseen lengths) ---"];

validationM = Range[120, 300, 20];

validationTable = Table[
   Module[{eMin, eMax, eW, pMin, pMax, pW},
    {eMin, eMax, eW} = FindAsymmetricWindow[allEstimates[M], targetAccuracy, M];
    {pMin, pMax, pW} = PredictBounds[M];
    {M, FormatBounds[eMin, eMax], eW, FormatBounds[pMin, pMax], pW}],
   {M, validationM}];

Print[Grid[
   Prepend[validationTable,
    {"M", "Empirical Bounds", "Width", "Predicted Bounds", "Width"}],
   Frame -> All]];


(* --- 7. Extreme test at M = 2000 ---------------------------------------- *)
(* Twenty times the upper limit of the training range. The coefficients are
   kept fixed: nothing is refitted here. *)
Print["\n--- Extreme validation test (M = 2000) ---"];
PrintTemporary["Simulating 2000 pulses. This is the slowest step."];

extremeM = 2000;
estimatesExt = RunTrials[extremeM, numTrials];

{eMinExt, eMaxExt, eWidthExt} =
  FindAsymmetricWindow[estimatesExt, targetAccuracy, extremeM];
{pMinExt, pMaxExt, pWidthExt} = PredictBounds[extremeM];

Print["Empirical bounds : ", FormatBounds[eMinExt, eMaxExt],
  "  width: ", eWidthExt];
Print["Predicted bounds : ", FormatBounds[pMinExt, pMaxExt],
  "  width: ", pWidthExt];

Print["\n--- Done. ---"];


(* --- 8. Export ----------------------------------------------------------- *)
(* Both plots go into a single PDF. A dialog asks where to save it. *)
With[{dir = SystemDialogInput["Directory"]},
 If[StringQ[dir],
  Export[FileNameJoin[{dir, "PredictiveScaling.pdf"}], combinedPlot];
  Print["Figure exported."],
  Print["Export cancelled."]]];
