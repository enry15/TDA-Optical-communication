(* ::Package:: *)

(* ==========================================================================
   Empirical determination of the optimal tolerance classes
   --------------------------------------------------------------------------
   Reproduces the tables of Section 5 of the paper, for one target accuracy
   at a time:
     - total number of required classes C,
     - specific class intervals [Lmin, Lmax],
   over a grid of signal lengths L and noise levels s.

   Both tables are printed as ready-to-paste LaTeX code. The caption and the
   label follow the value of targetAccuracy, so changing that single number
   produces the whole set of tables for another accuracy.

   Requirements: Mathematica 12 or later. No external packages.
   Runtime: the full 10 x 10 grid is the heavy part of the script. Shorten
   valL or valNoise for a quick test.

   Note on randomness: no random seed is fixed, so each run draws a fresh set
   of noise realisations. The reported values are reproduced within the
   statistical fluctuation of the Monte Carlo sampling.

   Author: <name>
   License: MIT (see LICENSE)
   ========================================================================== *)

ClearAll["Global`*"];


(* --- 1. Parameters ------------------------------------------------------ *)
targetAccuracy = 0.95;   (* tau: change this to get the 90, 85 or 80% tables *)
numTrials = 1000;        (* Monte Carlo trials per configuration             *)

sigmaVal = 0.8;          (* standard deviation of each Gaussian pulse        *)
stepRes = 0.1;           (* sampling resolution, in spatial units            *)
dist = 5;                (* distance d between two consecutive peaks         *)
minWidthSamples = 8;     (* W_min = 0.8 units = 8 samples at step 0.1        *)
amplitudeThreshold = 0.5;(* half-maximum level (HML)                         *)

valL = Range[10, 100, 10];
valNoise = Range[0.51, 0.60, 0.01];


(* --- 2. Signal generation and topological filter ------------------------ *)

(* Clean pulse train: L Gaussian pulses of unit amplitude, spaced by d. *)
GenerateCleanSignal[L_, maxSpace_] :=
  Quiet[
   Module[{xArray, peakPositions},
    xArray = Range[0, maxSpace, stepRes];
    peakPositions = Range[dist, L*dist, dist];
    Total[Table[Exp[-((xArray - pos)^2)/(2*sigmaVal^2)], {pos, peakPositions}]]],
   General::munfl];

(* Band-limited noise: uniform values on the integer grid, interpolated to
   the fine sampling grid. *)
GenerateNoise[noise_, maxSpace_, numPts_] :=
  ArrayResample[RandomReal[{-noise, noise}, maxSpace + 1], {numPts}];

(* Estimator of Eq. (1): connected components above the threshold that are
   at least W_min samples wide. *)
CountValidPulses[noisySignal_] :=
  Count[Split[UnitStep[noisySignal - amplitudeThreshold]],
   l : {1 ..} /; Length[l] >= minWidthSamples];


(* --- 3. Data-driven window expansion ------------------------------------ *)
(* The window starts at [L, L] and grows towards the adjacent class with the
   higher empirical frequency, until the coverage reaches the target. *)
FindOptimalWindow[Linput_, noiseInput_] :=
 Module[{maxSpace, numPts, clean, estimates, freq, prob,
   wMin, wMax, currentP, minV, maxV, pLeft, pRight},
  maxSpace = Linput*dist + 5;
  numPts = Round[maxSpace/stepRes] + 1;
  clean = GenerateCleanSignal[Linput, maxSpace];
  estimates =
   Table[CountValidPulses[clean + GenerateNoise[noiseInput, maxSpace, numPts]],
    {numTrials}];
  freq = Counts[estimates];
  prob[k_] := N[Lookup[freq, k, 0]/numTrials];
  wMin = Linput; wMax = Linput; currentP = prob[Linput];
  {minV, maxV} = {Min[Keys[freq]], Max[Keys[freq]]};
  While[currentP < targetAccuracy,
   If[wMin <= minV && wMax >= maxV, Break[]];
   pLeft = If[wMin > minV, prob[wMin - 1], -1.0];
   pRight = If[wMax < maxV, prob[wMax + 1], -1.0];
   If[pLeft >= pRight,
    wMin--; currentP += prob[wMin],
    wMax++; currentP += prob[wMax]]];
  {wMax - wMin + 1, "[" <> ToString[wMin] <> ", " <> ToString[wMax] <> "]"}];


(* --- 4. LaTeX exporter --------------------------------------------------- *)
GenerateLaTeXCode[matrix_, headers_, label_, captionText_] :=
 Module[{colFormat, headerStr, dataRows, finalCaption},
  colFormat = "c|" <> StringJoin[ConstantArray["c", Length[headers] - 1]];
  headerStr =
   StringRiffle[Map["\\textbf{" <> ToString[#] <> "}" &, headers], " & "] <> " \\\\";
  dataRows = Map[StringRiffle[ToString /@ #, " & "] <> " \\\\" &, matrix];
  finalCaption =
   captionText <> " (Target Accuracy: " <>
    ToString[Round[targetAccuracy*100]] <> "\\%).";
  StringJoin[
   "\\begin{table}[H]\n\\centering\n",
   "\\makebox[\\textwidth][c]{\n",
   "\\begin{tabular}{", colFormat, "}\n\\toprule\n",
   headerStr, "\n\\midrule\n",
   StringRiffle[dataRows, "\n"],
   "\n\\bottomrule\n\\end{tabular}}\n",
   "\\caption{", finalCaption, "}\n",
   "\\label{", label, "}\n\\end{table}\n"]];


(* --- 5. Execution -------------------------------------------------------- *)
Print["--- Computing tolerance classes at ",
  Round[targetAccuracy*100], "% target accuracy ---"];

rawMatrix = Table[
   PrintTemporary["Calculating L = ", L, " ..."];
   Table[FindOptimalWindow[L, noise], {noise, valNoise}],
   {L, valL}];

matrixClasses =
  Table[Prepend[rawMatrix[[i, All, 1]], valL[[i]]], {i, Length[valL]}];
matrixIntervals =
  Table[Prepend[rawMatrix[[i, All, 2]], valL[[i]]], {i, Length[valL]}];

headerRow = Prepend[valNoise, "L $\\backslash$ noise"];
accTag = ToString[Round[targetAccuracy*100]];


(* --- 6. Output ----------------------------------------------------------- *)

(* Readable version, to check the numbers directly in the notebook. *)
ShowTable[matrix_, title_] :=
 Print[
  Style[title <> "  (target accuracy: " <> accTag <> "%)", Bold, 14],
  "\n",
  Grid[Prepend[matrix, headerRow /. "L $\\backslash$ noise" -> "L \\ noise"],
   Frame -> All,
   Background -> {{LightGray, None}, {LightGray, None}},
   Alignment -> Center,
   Spacings -> {1.2, 1}]];

ShowTable[matrixClasses, "Total classes required"];
ShowTable[matrixIntervals, "Specific class intervals [Lmin, Lmax]"];

(* LaTeX version, ready to paste into the paper. *)
Print["\n--- TOTAL CLASSES (LaTeX) ---"];
Print[GenerateLaTeXCode[matrixClasses, headerRow,
   "tab:classes_" <> accTag, "Total classes required"]];

Print["\n--- SPECIFIC INTERVALS (LaTeX) ---"];
Print[GenerateLaTeXCode[matrixIntervals, headerRow,
   "tab:intervals_" <> accTag,
   "Specific class intervals $[L_{\\min}, L_{\\max}]$"]];

Print["\n--- Done. ---"];
