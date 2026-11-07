# Longitudinal Effects of Mantra-Based Meditation on Functional Connectivity

**Society for Neuroscience (SfN) 2026**

This repository contains MATLAB analysis code and supplementary materials for an EEG functional connectivity study investigating longitudinal changes associated with mantra-based meditation.

## Study Overview

The analysis included 31 participants:
- HK group: 16 participants
- SA group: 15 participants

Functional connectivity was quantified using phase-locking value (PLV) derived from RID-Rihaczek time-frequency representations.

The analyses examined a 12-electrode frontoparietal subnetwork across pre- and post-training resting-state and meditation recordings.

## Analysis Code

- [Alpha-band RID-PLV computation](code/RID_PLV_alpha.m) — 8–10 Hz and 10–12 Hz
- [Gamma-band RID-PLV computation](code/RID_PLV_gamma.m) — 30–50 Hz

The scripts use corrected RID frequency indexing and produce participant-level connectivity matrices.

## Dataset

The EEG recordings are available through the L-FAME Dataset Benchmark:

[Access the L-FAME dataset on Hugging Face](https://huggingface.co/datasets/L-FAME-Dataset-Benchmark/L-FAME)

The analyses presented at SfN used a subset of 31 participants. The original recordings are not duplicated in this repository.

See [dataset documentation](data/README.md) for additional information.

## Supplementary Figures

Additional figures will be available in the `figures/` directory.

## Requirements

- MATLAB
- EEGLAB
- `rid_rihaczek4.m`
- MATLAB Parallel Computing Toolbox

The scripts were developed for an HPC environment. Input paths and dependencies must be configured for the user's environment.
