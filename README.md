# Topological data analysis for optical communications

Code and data accompanying the paper *"<paper title>"* by E. M. Del Regno,
E. Hyry and M. Ornigotti (Tampere University).

The paper studies a counting algorithm based on zero-dimensional persistent
homology for highly noisy optical channels. Information is encoded in the
total number of pulses received rather than in their individual positions,
which makes the transmission more resilient to noise. This repository
contains the Mathematica code that produces every figure and table of the
paper, together with the complete set of results, including the tables that
were left out of the article for reasons of space.

## Contents

### Code

All scripts are written for Mathematica 12 or later and require no external
packages. Each one is self-contained: parameters are collected at the top of
the file and can be changed without touching the rest of the code. Every
script comes with a PDF version of the evaluated notebook, so that the code
and its output can be read without running Mathematica.

| File | What it does | Paper |
|---|---|---|
| `tolerance_classes.wl` | Monte Carlo determination of the optimal tolerance classes over a grid of signal lengths and noise levels. Prints the results both as readable tables and as ready-to-paste LaTeX code. | Section 5 |
| `distribution_plots.wl` | Example of a corrupted pulse train, and outcome distributions of the topological filter at three noise levels. | Sections 4 and 6 |
| `pnm_vs_ook.wl` | Performance comparison between Pulse Number Modulation and On-Off Keying, including the single-bit reliability test. | Section 9 |
| `predictive_scaling.wl` | Calibration of the predictive model for the mean drift and the standard deviation, validation on unseen signal lengths, and extrapolation to very long pulse trains. | Section 10 |

### Results

| File | Contents |
|---|---|
| `figures.pdf` | The outcome distributions produced by `distribution_plots.wl`, collected for every signal length considered. |
| `tables.pdf` | The tolerance tables produced by `tolerance_classes.wl`, for all four target accuracies. |

## How to run

Open a script in Mathematica and evaluate the whole notebook. The parameters
that are meant to be changed are grouped in the first section of each file:

- `tolerance_classes.wl` computes one target accuracy per run. Set
  `targetAccuracy` to `0.95`, `0.90`, `0.85` or `0.80`; captions and labels
  of the generated LaTeX follow automatically.
- `distribution_plots.wl` produces the figures for one signal length per run.
  Set `targetL` to the desired value.
- `pnm_vs_ook.wl` runs at `noiseLevel = 0.6` by default. Set it to `0.7` for
  the extreme-noise comparison.

Runtimes range from a few minutes to considerably longer for the full
tolerance grid. To try a script quickly, shorten the lists of signal lengths
and noise levels at the top of the file.

## Reproducibility

All results come from Monte Carlo simulations with 1000 trials per
configuration. No random seed is fixed, so re-running a script does not
return the exact numbers printed in the paper, but values consistent with
them within the statistical fluctuation of the sampling.

## Citation

If you use this code, please cite the paper:

```
<bibtex entry>
```
