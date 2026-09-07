(* ::Package:: *)

(* ==========================================================================
   PNM vs OOK: performance comparison in a high-noise channel
   --------------------------------------------------------------------------
   Reproduces Section 9 of the paper:
     - Table 10: per-bitstring accuracy of PNM and OOK at s = 0.6.
     - Table 11: summary metrics (average, best and worst case).
     - Table 12: single-bit reliability of the OOK receiver.

   Setting noiseLevel = 0.7 reproduces Tables 13 and 14 (extreme noise).

   Requirements: Mathematica 12 or later. No external packages.
   Runtime: a few minutes on a standard laptop.

   Note on randomness: no random seed is fixed, so each run draws a fresh set
   of noise realisations. The reported values are reproduced within the
   statistical fluctuation of the Monte Carlo sampling.
   ========================================================================== *)

ClearAll["Global`*"];


(* --- 1. Parameters ------------------------------------------------------ *)
numTrials = 1000;        (* Monte Carlo trials per configuration            *)
numTrials1Bit = 2000;    (* trials for the single-bit test of Table 12      *)
sigma = 0.8;             (* standard deviation of each Gaussian pulse       *)
noiseLevel = 0.6;        (* noise bound s: set to 0.7 for Tables 13-14      *)
stepRes = 0.1;           (* sampling resolution, in spatial units           *)
dist = 5;                (* distance d between two consecutive pulse slots  *)
minWidthSamples = 8;     (* W_min = 0.8 units = 8 samples at step 0.1       *)
amplitudeThreshold = 0.5;(* half-maximum level (HML)                        *)

silenceSlots = 5;        (* auto-gating: silent gap that ends the message   *)
silenceMaxSpace = silenceSlots*dist;   (* converted to spatial units *)

(* PNM configuration (Table 9): 8 contiguous bins of width 10. *)
gridL = 80;
binWidth = 10;
numClasses = 8;
startValue = 0;

(* OOK configuration: 3 pulse slots only. *)
gridLOOK = 3;


(* --- 2. Shared building blocks ------------------------------------------ *)

(* Clean pulse train: unit-amplitude Gaussian pulses at the given positions. *)
GenerateCleanSignal[peakPos_List, maxSpace_, numPts_] :=
  If[Length[peakPos] > 0,
   Total[Table[Exp[-((Range[0, maxSpace, stepRes] - p)^2)/(2*sigma^2)],
     {p, peakPos}]],
   ConstantArray[0, numPts]];

(* Band-limited noise: uniform values on the integer grid, interpolated to
   the fine sampling grid (channel model of Section 4). *)
GenerateNoise[maxSpace_, numPts_] :=
  ArrayResample[
   RandomReal[{-noiseLevel, noiseLevel}, maxSpace + 1], {numPts}];

(* PNM receiver: counts the valid connected components, and stops as soon as
   it has seen a long enough silence (auto-gating). *)
AnalyzeWithAutoGating[noisySignal_] :=
 Module[{runs, count = 0, silence = 0},
  runs = Split[UnitStep[noisySignal - amplitudeThreshold]];
  Do[
   If[First[run] == 1,
    (* a component: keep it only if it is wide enough *)
    If[Length[run] >= minWidthSamples, count++; silence = 0],
    (* a gap: accumulate its length and stop if the message has ended *)
    silence += Length[run]*stepRes;
    If[silence >= silenceMaxSpace, Break[]]],
   {run, runs}];
  count];

(* OOK receiver: looks for a valid pulse inside the time slot centred on
   tCentre, and returns 1 or 0. *)
DecodeOOKSlot[noisySignal_, tCentre_, numPts_] :=
 Module[{idxStart, idxEnd, window},
  idxStart = Clip[Round[(tCentre - 2.5)/stepRes], {1, numPts}];
  idxEnd = Clip[Round[(tCentre + 2.5)/stepRes], {1, numPts}];
  window = noisySignal[[idxStart ;; idxEnd]];
  If[Count[Split[UnitStep[window - amplitudeThreshold]],
     l : {1 ..} /; Length[l] >= minWidthSamples] >= 1, 1, 0]];

FormatBin[{lower_, upper_}] :=
  "[" <> ToString[lower] <> " - " <> ToString[upper] <> "]";


(* --- 3. PNM: accuracy per bin ------------------------------------------- *)
Print["--- Testing PNM (noise level s = ", noiseLevel, ") ---"];

maxSpacePNM = gridL*dist + 5;
numPtsPNM = Round[maxSpacePNM/stepRes] + 1;

(* Bins {0,9}, {10,19}, ... as in Table 9. *)
binsList =
  Table[{startValue + (i - 1)*binWidth, startValue + i*binWidth - 1},
   {i, numClasses}];

pnmResults = Table[
   Module[{lowerB, upperB, targetL, peakPos, clean, correct = 0, estimate},
    {lowerB, upperB} = bin;
    targetL = Ceiling[Mean[bin]];
    (* The clean signal depends only on targetL, so it is built once. *)
    peakPos = Table[5 + (j - 1)*dist, {j, targetL}];
    clean = GenerateCleanSignal[peakPos, maxSpacePNM, numPtsPNM];
    Do[
     estimate =
      AnalyzeWithAutoGating[clean + GenerateNoise[maxSpacePNM, numPtsPNM]];
     If[lowerB <= estimate <= upperB, correct++],
     {numTrials}];
    {targetL, {lowerB, upperB}, N[correct/numTrials]}],
   {bin, binsList}];


(* --- 4. OOK: accuracy per 3-bit string ---------------------------------- *)
Print["--- Testing OOK (noise level s = ", noiseLevel, ") ---"];

maxSpaceOOK = gridLOOK*dist + 5;
numPtsOOK = Round[maxSpaceOOK/stepRes] + 1;
testStrings = Tuples[{0, 1}, 3];

ookResults = Table[
   Module[{peakPos, clean, correct = 0, noisy, decoded},
    peakPos = Flatten[Table[If[str[[j]] == 1, {j*dist}, {}], {j, gridLOOK}]];
    clean = GenerateCleanSignal[peakPos, maxSpaceOOK, numPtsOOK];
    Do[
     noisy = clean + GenerateNoise[maxSpaceOOK, numPtsOOK];
     decoded = Table[DecodeOOKSlot[noisy, j*dist, numPtsOOK], {j, gridLOOK}];
     If[decoded === str, correct++],
     {numTrials}];
    N[correct/numTrials]],
   {str, testStrings}];


(* --- 5. Table 10: per-bitstring comparison ------------------------------ *)
(* The i-th bin of Table 9 encodes the i-th bitstring, so the two lists are
   compared row by row. *)
comparisonTable = Table[
   Module[{str, targetL, bin, accPNM, accOOK},
    str = StringJoin[ToString /@ testStrings[[i]]];
    {targetL, bin, accPNM} = pnmResults[[i]];
    accOOK = ookResults[[i]];
    {str, targetL, FormatBin[bin],
     100. accPNM, 100. accOOK, 100. (accPNM - accOOK)}],
   {i, numClasses}];

Print["\n=== Table 10: PNM vs OOK per bitstring (s = ", noiseLevel, ") ==="];
Print[Grid[
   Join[
    {{"Bitstring", "Target L", "Bin Range", "PNM (%)", "OOK (%)", "Gain (%)"}},
    comparisonTable,
    {{"Average", "", "",
      Mean[comparisonTable[[All, 4]]],
      Mean[comparisonTable[[All, 5]]],
      Mean[comparisonTable[[All, 6]]]}}],
   Frame -> All]];


(* --- 6. Table 11: summary metrics --------------------------------------- *)
accPNMall = comparisonTable[[All, 4]];
accOOKall = comparisonTable[[All, 5]];

Print["\n=== Table 11: summary (s = ", noiseLevel, ") ==="];
Print[Grid[
   {{"Metric", "PNM (Proposed)", "OOK (Traditional)"},
    {"Average Accuracy (%)", Mean[accPNMall], Mean[accOOKall]},
    {"Best Case Accuracy (%)", Max[accPNMall], Max[accOOKall]},
    {"Worst Case Accuracy (%)", Min[accPNMall], Min[accOOKall]},
    {"Spatial Units Used", gridL*dist, gridLOOK*dist}},
   Frame -> All]];


(* --- 7. Table 12: single-bit reliability of OOK -------------------------- *)
Print["\n--- Testing single-bit OOK reliability ---"];

maxSpace1Bit = dist + 5;
numPts1Bit = Round[maxSpace1Bit/stepRes] + 1;

results1Bit = Association @ Table[
    Module[{peakPos, clean, correct = 0, noisy, decoded},
     peakPos = If[bit == 1, {dist}, {}];
     clean = GenerateCleanSignal[peakPos, maxSpace1Bit, numPts1Bit];
     Do[
      noisy = clean + GenerateNoise[maxSpace1Bit, numPts1Bit];
      decoded = DecodeOOKSlot[noisy, dist, numPts1Bit];
      If[decoded == bit, correct++],
      {numTrials1Bit}];
     bit -> N[correct/numTrials1Bit]],
    {bit, {0, 1}}];

Print["\n=== Table 12: single-bit OOK accuracy (s = ", noiseLevel, ") ==="];
Print[Grid[
   {{"Target Bit", "OOK Accuracy (%)"},
    {"Bit 0", 100. results1Bit[0]},
    {"Bit 1", 100. results1Bit[1]},
    {"Average", 100. Mean[Values[results1Bit]]}},
   Frame -> All]];

Print["\n--- Done. ---"];
