(* ::Package:: *)

(* ==========================================================================
   Example signal and outcome distributions
   --------------------------------------------------------------------------
   Produces two figures for one signal length at a time:
     - an example of the physical data stream: clean pulse train, corrupted
       signal, and the half-maximum threshold used for the filtration;
     - the outcome distributions of the topological filter at three noise
       levels, as in Figs. 6-10 of the paper.

   The three histograms share the same vertical scale, so that their heights
   can be compared directly on the printed page, and each panel carries its
   noise level in the title.

   Requirements: Mathematica 12 or later. No external packages.
   Runtime: a few minutes.

   Note on randomness: no random seed is fixed, so each run draws a fresh set
   of noise realisations.

   ========================================================================== *)

ClearAll["Global`*"];


(* --- 1. Parameters ------------------------------------------------------ *)
targetL = 20;            (* signal length shown in both figures             *)
numTrials = 1000;        (* Monte Carlo trials per noise level              *)

sigmaVal = 0.8;          (* standard deviation of each Gaussian pulse       *)
stepRes = 0.1;           (* sampling resolution, in spatial units           *)
dist = 5;                (* distance d between two consecutive peaks        *)
minWidthSamples = 8;     (* W_min = 0.8 units = 8 samples at step 0.1       *)
amplitudeThreshold = 0.5;(* half-maximum level (HML)                        *)

valNoise = {0.51, 0.55, 0.60};   (* noise levels of the three panels        *)
exampleNoise = 0.60;             (* noise level of the example signal       *)


(* --- 2. Signal generation and topological filter ------------------------ *)

GenerateCleanSignal[L_, maxSpace_] :=
  Quiet[
   Module[{xArray, peakPositions},
    xArray = Range[0, maxSpace, stepRes];
    peakPositions = Range[dist, L*dist, dist];
    Total[Table[Exp[-((xArray - pos)^2)/(2*sigmaVal^2)], {pos, peakPositions}]]],
   General::munfl];

GenerateNoise[noise_, maxSpace_, numPts_] :=
  ArrayResample[RandomReal[{-noise, noise}, maxSpace + 1], {numPts}];

CountValidPulses[noisySignal_] :=
  Count[Split[UnitStep[noisySignal - amplitudeThreshold]],
   l : {1 ..} /; Length[l] >= minWidthSamples];

(* One full Monte Carlo run: returns the list of estimated counts. *)
RunTrials[L_, noise_] :=
 Module[{maxSpace, numPts, clean},
  maxSpace = L*dist + 5;
  numPts = Round[maxSpace/stepRes] + 1;
  clean = GenerateCleanSignal[L, maxSpace];
  Table[CountValidPulses[clean + GenerateNoise[noise, maxSpace, numPts]],
   {numTrials}]];


(* --- 3. Example of the physical data stream ----------------------------- *)
Print["--- Example signal, L = ", targetL, ", s = ", exampleNoise, " ---"];

examplePlot =
 Module[{maxSpace, numPts, xGrid, clean, noisy, count},
  maxSpace = targetL*dist + 5;
  numPts = Round[maxSpace/stepRes] + 1;
  xGrid = Range[0, maxSpace, stepRes];
  clean = GenerateCleanSignal[targetL, maxSpace];
  noisy = clean + GenerateNoise[exampleNoise, maxSpace, numPts];
  count = CountValidPulses[noisy];
  ListLinePlot[
   {Transpose[{xGrid, clean}], Transpose[{xGrid, noisy}]},
   PlotStyle -> {{Black, AbsoluteThickness[2.5]}, {Red, AbsoluteThickness[1.2]}},
   PlotRange -> {{-1, maxSpace + 1}, {-1.2, 2.2}},
   GridLines -> {{}, {amplitudeThreshold}},
   GridLinesStyle -> Directive[Dashed, Black, AbsoluteThickness[1.5]],
   Frame -> True, FrameStyle -> Thick,
   FrameLabel -> {"Space (units)", "Amplitude"},
   PlotLabel ->
    Style[StringForm["L = ``, s = `` (counted: ``)",
      targetL, exampleNoise, count], Bold, 15],
   ImageSize -> 850, AspectRatio -> 1/4]];

Print[examplePlot];
Print["Black: clean pulse train. Red: corrupted signal. \
Dashed line: threshold at 0.5."];


(* --- 4. Outcome distributions at three noise levels --------------------- *)
Print["\n--- Outcome distributions, L = ", targetL, " ---"];

allData = Table[
   PrintTemporary["Computing s = ", noise, " ..."];
   RunTrials[targetL, noise],
   {noise, valNoise}];

(* Common vertical scale: the tallest bar of the three panels sets the height
   of all of them, so the plots can be compared directly. *)
yMax = 1.05*Max[Map[Max[Counts[#]] &, allData]];

xMin = targetL - 6;
xMax = targetL + 6;
ticksX = {targetL - 4, targetL - 2, targetL, targetL + 2, targetL + 4};

histograms = Table[
   Histogram[allData[[i]], {1}, "Count",
    ChartStyle -> EdgeForm[Thin],
    PlotRange -> {{xMin, xMax}, {0, yMax}},
    GridLines -> {{targetL}, {}},
    GridLinesStyle -> Directive[Red, Thick, Dashed],
    Frame -> True,
    PlotLabel ->
     Style[StringForm["L = ``, s = ``", targetL, valNoise[[i]]], 22, Bold],
    ImageSize -> 340,
    AspectRatio -> 0.85,
    (* The padding is the same on every panel, including the ones without
       tick labels, so that all the frames have identical size and sit on
       the same baseline. *)
    ImagePadding -> {{75, 15}, {55, 45}},
    FrameTicks -> {{If[i == 1, Automatic, None], None}, {ticksX, None}},
    FrameLabel -> {None, If[i == 1, Style["Count", 18, Bold], None]},
    LabelStyle -> {FontSize -> 18, Bold, Black}],
   {i, Length[valNoise]}];

Print[Row[Riffle[histograms, Spacer[10]]]];
Print["Red dashed line: true number of pulses. \
All panels share the same vertical scale."];

Print["\n--- Done. ---"];
