(* ::Package:: *)

(* ==========================================================================
   Predictive scaling of topological tolerance classes
   --------------------------------------------------------------------------
   Reproduces Section 10 of the paper:
     - Eq. (2)-(3): fit of the mean drift mu(L) and of the standard
       deviation sigma(L), trained only on 10 <= L <= 100.
     - Fig. 11 : predicted vs. empirical mu(L) and sigma(L) up to L = 300.
     - Table 15: empirical vs. predicted tolerance bounds, 120 <= L <= 300.
     - Extreme test at L = 2000.

   Requirements: Mathematica 12 or later. No external packages.
   Runtime: about 3 minutes on a standard laptop. The L = 2000 test is the
   slowest single step.

   Author: <name>
   License: MIT (see LICENSE)
   ========================================================================== *)

ClearAll["Global`*"];

(* --- 0. Note on randomness ---------------------------------------------- *)
(* No random seed is fixed, so each run draws a fresh set of noise
   realisations. The reported values are therefore reproduced within the
   statistical fluctuation of 1000 Monte Carlo trials: in our runs the fitted
   coefficients are stable to the digits quoted in the paper, and the
   tolerance bounds move by at most one pulse. *)


(* --- 1. Parameters ------------------------------------------------------ *)
noiseLevel = 0.6;        (* noise bound s: noise is uniform in [-s, s]      *)
sigmaGauss = 0.8;        (* standard deviation of each Gaussian pulse       *)
stepRes = 0.1;           (* sampling resolution, in spatial units           *)
dist = 5;                (* distance d between two consecutive peaks        *)
minWidthSamples = 8;     (* W_min = 0.8 units = 8 samples at step 0.1       *)
amplitudeThreshold = 0.5;(* half-maximum level (HML) used for the filtration*)
numTrials = 1000;        (* Monte Carlo trials per configuration            *)
targetAccuracy = 0.95;   (* target accuracy tau                             *)


(* --- 2. Signal generation and topological decoding ---------------------- *)

(* Clean pulse train: L Gaussian pulses of unit amplitude, spaced by d. *)
GenerateCleanSignal[L_Integer, maxSpace_, step_] :=
  Module[{grid, peakPos},
   grid = Range[0, maxSpace, step];
   peakPos = Range[1, L]*dist;
   Total[Table[Exp[-((grid - p)^2)/(2*sigmaGauss^2)], {p, peakPos}]]];

(* Band-limited noise: uniform values on the integer grid, then interpolated
   to the fine sampling grid. This is the channel model of Section 4. *)
GenerateNoise[maxSpace_, numPts_] :=
  ArrayResample[
   RandomReal[{-noiseLevel, noiseLevel}, maxSpace + 1], {numPts}];

(* Estimator L-hat of Eq. (1): count the connected components of the
   super-level set at threshold 0.5 that are at least W_min samples wide.
   This is Algorithm 1 of Section 11.1. *)
DecodeTDA[noisySignal_] :=
  Count[Split[UnitStep[noisySignal - amplitudeThreshold]],
   l : {1 ..} /; Length[l] >= minWidthSamples];

(* One full Monte Carlo run for a given L: returns the list of estimates. *)
RunTrials[L_Integer, trials_Integer] :=
  Module[{maxSpace, numPts, cleanSignal},
   maxSpace = L*dist + 5;
   numPts = Round[maxSpace/stepRes] + 1;
   cleanSignal = GenerateCleanSignal[L, maxSpace, stepRes];
   Table[DecodeTDA[cleanSignal + GenerateNoise[maxSpace, numPts]],
    {trials}]];

(* Data-driven window expansion (Algorithm 2 of Section 11.2).
   The window starts at [L, L] and grows towards the adjacent class with the
   higher empirical frequency, until the coverage reaches the target. *)
FindAsymmetricWindow[estimates_, targetAcc_, L_] :=
 Module[{freq, prob, wMin, wMax, currentP, minV, maxV, pLeft, pRight},
  freq = Counts[estimates];
  prob[k_] := N[Lookup[freq, k, 0]/Length[estimates]];
  wMin = L; wMax = L; currentP = prob[L];
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


(* --- 3. Empirical data, L = 10 to 300 ----------------------------------- *)
Print["--- Generating experimental data (L = 10 to 300) ... ---"];

allL = Join[Range[10, 100, 10], Range[120, 300, 20]];
meanErrors = Association[];
stdDevs = Association[];
allEstimates = Association[];

Do[
  Module[{estimates, errors},
   estimates = RunTrials[L, numTrials];
   errors = estimates - L;
   allEstimates[L] = estimates;
   meanErrors[L] = Mean[errors];
   stdDevs[L] = StandardDeviation[errors]],
  {L, allL}];


(* --- 4. Calibration, trained only on L <= 100 --------------------------- *)
Print["\n--- Fitting the model (training set: L <= 100) ---"];

trainL = Range[10, 100, 10];
dataMeanTrain = Table[{L, meanErrors[L]}, {L, trainL}];
dataStdTrain = Table[{L, stdDevs[L]}, {L, trainL}];

(* Eq. (2): the mean bias of each pulse adds coherently, so mu grows as L. *)
fitMean = FindFit[dataMeanTrain, m*L, {m}, L];
mFit = m /. fitMean;

(* Eq. (3): the random part adds incoherently, so sigma grows as Sqrt[L]. *)
fitStd = FindFit[dataStdTrain, k*Sqrt[L], {k}, L];
kFit = k /. fitStd;

Print["Mean drift : mu(L) = ", mFit, " * L"];
Print["Std. dev.  : sigma(L) = ", kFit, " * Sqrt[L]"];

(* Eq. (4): predicted tolerance interval at the 95% target (z = 1.96). *)
zScore = 1.96;
PredictBounds[L_] :=
  Module[{center, halfWidth, lower, upper},
   center = L + mFit*L;
   halfWidth = zScore*kFit*Sqrt[L];
   lower = Round[center - halfWidth];
   upper = Round[center + halfWidth];
   {lower, upper, upper - lower + 1}];


(* --- 5. Figure 11: predicted vs. empirical mu(L) and sigma(L) ----------- *)
dataMeanAll = Table[{L, meanErrors[L]}, {L, allL}];
dataStdAll = Table[{L, stdDevs[L]}, {L, allL}];

plotMean = Show[
   ListPlot[dataMeanAll, PlotStyle -> Directive[Red, PointSize[Large]],
    PlotLegends -> {"Empirical Data"}],
   Plot[mFit*L, {L, 10, 300}, PlotStyle -> Directive[Blue, Dashed, Thick],
    PlotLegends -> {"Prediction"}],
   AxesLabel -> {"Number of Pulses (L)", "Mean Error (Drift)"},
   PlotLabel -> "Predictive Power: Mean Error up to 300",
   ImageSize -> Medium];

plotStd = Show[
   ListPlot[dataStdAll, PlotStyle -> Directive[Red, PointSize[Large]],
    PlotLegends -> {"Empirical Data"}],
   Plot[kFit*Sqrt[L], {L, 10, 300}, PlotStyle -> Directive[Blue, Dashed, Thick],
    PlotLegends -> {"Prediction"}],
   AxesLabel -> {"Number of Pulses (L)", "Standard Deviation"},
   PlotLabel -> "Predictive Power: StdDev up to 300",
   ImageSize -> Medium];

combinedPlot = GraphicsRow[{plotMean, plotStd}, ImageSize -> 1000];
Print[combinedPlot];


(* --- 6. Table 15: empirical vs. predicted bounds, 120 <= L <= 300 ------- *)
Print["\n--- Building Table 15 (validation on unseen lengths) ---"];

validationL = Range[120, 300, 20];

validationTable = Table[
   Module[{eMin, eMax, eW, pMin, pMax, pW},
    {eMin, eMax, eW} = FindAsymmetricWindow[allEstimates[L], targetAccuracy, L];
    {pMin, pMax, pW} = PredictBounds[L];
    {L, FormatBounds[eMin, eMax], eW, FormatBounds[pMin, pMax], pW}],
   {L, validationL}];

Print[Grid[
   Prepend[validationTable,
    {"L", "Empirical Bounds", "Width", "Predicted Bounds", "Width"}],
   Frame -> All]];


(* --- 7. Extreme test at L = 2000 ---------------------------------------- *)
(* Twenty times the upper limit of the training range. The coefficients are
   kept fixed: nothing is refitted here. *)
Print["\n--- Extreme validation test (L = 2000) ---"];
PrintTemporary["Simulating 2000 pulses. This is the slowest step."];

extremeL = 2000;
estimatesExt = RunTrials[extremeL, numTrials];

{eMinExt, eMaxExt, eWidthExt} =
  FindAsymmetricWindow[estimatesExt, targetAccuracy, extremeL];
{pMinExt, pMaxExt, pWidthExt} = PredictBounds[extremeL];

Print["Empirical bounds : ", FormatBounds[eMinExt, eMaxExt],
  "  width: ", eWidthExt];
Print["Predicted bounds : ", FormatBounds[pMinExt, pMaxExt],
  "  width: ", pWidthExt];

Print["\n--- Done. ---"];
