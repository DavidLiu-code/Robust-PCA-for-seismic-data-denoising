# Robust PCA Seismic Denoising

This repository contains two single-file MATLAB workflows for reproducing and extending the f-x-y domain Robust Principal Component Analysis (RPCA) / Principal Component Pursuit (PCP) idea used by Cheng, Chen, and Sacchi (2015). The goal is to attenuate erratic noise in seismic data, especially sparse, high-amplitude, non-Gaussian outliers that contaminate only part of the spatial domain.

## Reference Paper

The main reference for this repository is:

> Cheng, J., Chen, K., & Sacchi, M. D. (2015). **Application of Robust Principal Component Analysis (RPCA) to suppress erratic noise in seismic records**. *SEG Technical Program Expanded Abstracts 2015*, 4646-4651.

- DOI: <https://doi.org/10.1190/segam2015-5869427.1>
- PDF: <https://saig.physics.ualberta.ca/lib/exe/fetch.php?media=publications%3Aconfpro%3Aseg2015_cheng_2.pdf>

## Motivation

Many conventional seismic random-noise attenuation methods work best when the noise is approximately Gaussian or can be smoothly suppressed in a transform domain. Field data, however, often contain erratic noise: dead or bad traces, impulsive interference, acquisition-system artifacts, or locally strong contamination. This type of noise is usually sparse in space but can have very large amplitude, so least-squares methods and fixed-rank f-x-y PCA/eigenimage filters can be strongly biased by these outliers.

The key idea in Cheng, Chen, and Sacchi (2015) is that, for each frequency slice, coherent seismic events can often be represented by a low-rank matrix, while erratic noise is better represented as a sparse matrix. Each f-x-y slice is therefore modeled as:

$$
\mathbf{D}(f) = \mathbf{L}(f) + \mathbf{S}(f)
$$

where:

- `D(f)` is the input f-x-y slice at frequency `f`;
- `L(f)` is the low-rank component, interpreted as coherent seismic signal;
- `S(f)` is the sparse component, interpreted as erratic noise or strong outlier contamination.

## Algorithm Overview

### f-x-y Frequency Slices

The input data are stored as `[nt, nx, ny]`. The workflows first apply an FFT along the time axis, transforming the data from the t-x-y domain to the f-x-y domain. For each positive frequency `f`, the code extracts a spatial matrix:

$$
\mathbf{D}_f = \mathbf{D}(f, :, :)
$$

Two-dimensional input `[nt, nx]` is internally reshaped to `[nt, nx, 1]`. Three-dimensional input directly forms `nx` by `ny` frequency slices.

### PCA / Eigenimage Baseline

Conventional f-x-y eigenimage filtering applies SVD to each frequency slice:

$$
\mathbf{D}_f = \mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^{*}
$$

and keeps only the first `r` singular values:

$$
\mathbf{L}_f =
\mathbf{U}_{(:,1:r)}
\boldsymbol{\Sigma}_{(1:r,1:r)}
\mathbf{V}_{(:,1:r)}^{*}
$$

This is effective for relatively stationary random noise. However, when a few traces or local spatial samples contain strong outliers, the dominant singular vectors can be distorted. The repository therefore keeps PCA/eigenimage filtering as a baseline for comparison with RPCA.

### RPCA / PCP Decomposition

RPCA is implemented through the Principal Component Pursuit formulation:

$$
\begin{aligned}
\min_{\mathbf{L},\mathbf{S}}\quad
& \lVert \mathbf{L} \rVert_* + \lambda \lVert \mathbf{S} \rVert_1 \\
\text{subject to}\quad
& \mathbf{D} = \mathbf{L} + \mathbf{S}.
\end{aligned}
$$

where:

- $\lVert \mathbf{L} \rVert_*$ is the nuclear norm, the sum of singular values, which promotes a low-rank `L`;
- $\lVert \mathbf{S} \rVert_1$ is the element-wise L1 norm, which promotes a sparse `S`;
- `lambda` controls how much energy is assigned to the sparse noise component.

The MATLAB code solves this problem with an inexact augmented Lagrange multiplier (IALM) style iteration. Each iteration mainly applies two thresholding operations:

- singular-value soft-thresholding to estimate the low-rank matrix `L`;
- element-wise soft-thresholding to estimate the sparse matrix `S`.

After estimating `L_f` for each positive frequency, the workflows restore the negative-frequency side using conjugate symmetry and then apply an inverse FFT back to the time domain. The main outputs are:

- `rpca_denoised`: the low-rank reconstruction;
- `rpca_sparse_noise`: the sparse noise separated by RPCA;
- `removed_noise = noisy - rpca_denoised`: the total removed component, which is often useful for QC.

## What This Repository Adds

The reference paper focuses on separating low-rank coherent signal and sparse erratic noise in the f-x-y domain. This repository follows that core idea and adds several practical details that make the workflows more stable for synthetic tests and easier to apply to real data.

### 1. Smooth Spectral Taper

The workflows use a raised-cosine spectral taper instead of hard frequency cutoffs:

```matlab
freq_taper.flow_stop  = 0;
freq_taper.flow_pass  = 2;
freq_taper.fhigh_pass = 45;
freq_taper.fhigh_stop = 75;
```

Frequencies between `flow_pass` and `fhigh_pass` are fully preserved. The low- and high-frequency edges are tapered smoothly, which helps reduce Gibbs/ringing artifacts caused by abrupt truncation.

### 2. Time-Axis Zero Padding

Before the FFT, the workflows optionally zero-pad the time axis using `pad_factor`. This gives smoother frequency sampling and reduces circular leakage and edge ringing after inverse FFT.

### 3. Rank Projection After RPCA

Pure PCP already separates `D = L + S`, but real data may still contain weak incoherent residual noise in `L`. The engineering workflow can optionally apply rank projection after RPCA:

```matlab
opts.use_rank_projection = true;
opts.rank_mode = 'adaptive';
```

This is not required by the original PCP formulation. It is a practical enhancement: RPCA first removes strong sparse outliers, and the subsequent low-rank projection suppresses weaker residual random energy.

### 4. Sliding-Window Overlap-Add

`rpca_engineering_workflow.m` supports sliding-window processing:

```matlab
USE_SLIDING_WINDOW = true;
win_length_s = 0.50;
win_overlap  = 0.50;
```

Each time window is processed independently in the f-x-y domain, and the final result is reconstructed with sine-window overlap-add weights. This helps handle nonstationary events, changing spectra, and time-localized noise behavior in longer records.

### 5. Adaptive Lambda and Rank

The engineering workflow supports frequency-dependent `lambda(f)` and `rank(f)`:

```matlab
opts.lambda_adaptive = true;
opts.lambda_freq_boost = 0.40;
opts.rank_mode = 'adaptive';
opts.rank_min = 2;
opts.rank_max = 7;
```

The motivation is that low frequencies are often stronger and more coherent, while high frequencies are more easily contaminated by noise. Allowing the sparse penalty and rank cap to vary with frequency avoids forcing one rigid parameter set on the full bandwidth.

### 6. Synthetic Parameter Sweep

When `INPUT_MODE = 'synthetic'`, the script has access to the clean reference and can compute:

$$
Q =
10 \log_{10}
\frac{
\lVert \mathbf{d}_{\mathrm{clean}} \rVert_F^2
}{
\lVert \mathbf{d}_{\mathrm{estimate}} - \mathbf{d}_{\mathrm{clean}} \rVert_F^2
}
$$

It can also scan `lambda` and rank parameters automatically. For real data, where no clean reference is available, Q metrics are skipped and QC relies on wiggle plots, image displays, removed-noise panels, rank statistics, iteration counts, and residual diagnostics.

## Repository Structure

```text
.
├── src/
│   ├── rpca_tapered_reproduction.m
│   └── rpca_engineering_workflow.m
├── .gitignore
└── README.md
```

## MATLAB Scripts

| Script | Purpose |
| --- | --- |
| `src/rpca_tapered_reproduction.m` | A compact tapered f-x-y RPCA/PCP reproduction workflow. It includes synthetic data generation, smooth spectral tapering, a PCA baseline, RPCA, and a synthetic parameter sweep. This is the best entry point for understanding the main algorithm. |
| `src/rpca_engineering_workflow.m` | A more practical workflow with sliding-window overlap-add, adaptive frequency-dependent lambda, adaptive rank projection, MAT input, and optional SEG-Y I/O hooks. This is the preferred entry point for testing real data. |

## Requirements

- MATLAB, preferably R2018b or newer.
- No extra toolbox is required for the default synthetic examples.
- SEG-Y input/output depends on `segyread` / `segywrite` support in your MATLAB environment.

## Quick Start

In MATLAB, change to the repository root and run the engineering workflow:

```matlab
cd('path/to/Robust PCA')
run('src/rpca_engineering_workflow.m')
```

For the lighter reproduction workflow:

```matlab
run('src/rpca_tapered_reproduction.m')
```

The default configuration generates synthetic 3D seismic data, adds Gaussian noise and erratic noise, and compares:

- smooth spectral taper reconstruction;
- f-x-y PCA/eigenimage filtering;
- enhanced f-x-y RPCA/PCP denoising;
- RPCA sparse noise and removed-noise QC.

## Using Your Own MAT Data

In the engineering workflow, modify:

```matlab
INPUT_MODE = 'mat';
USER_DATA_MAT = 'your_data.mat';
DATA_VAR_NAME = 'data';
dt = 0.002;
```

The array referenced by `DATA_VAR_NAME` should be:

- `[nt, nx]`: a 2D gather;
- `[nt, nx, ny]`: a 3D volume or a set of y slices.

Two-dimensional input is automatically reshaped to `[nt, nx, 1]`.

## Key Parameters

- `freq_taper`: smooth spectral taper used to avoid Gibbs/ringing artifacts from hard frequency cutoffs.
- `pad_factor`: zero-padding factor along the time axis before FFT.
- `USE_SLIDING_WINDOW`, `win_length_s`, `win_overlap`: sliding-window processing controls.
- `opts.lambda_scale`: base sparse penalty for RPCA/PCP.
- `opts.lambda_adaptive`: enables frequency-dependent lambda.
- `opts.use_rank_projection`, `opts.rank_mode`: low-rank projection strategy after RPCA.
- `DO_SYNTHETIC_SWEEP`: parameter sweep switch for synthetic data; recommended to disable for real data.

## Outputs

By default, the workflows save MATLAB result files:

- `rpca_tapered_reproduction_results.mat`
- `rpca_engineering_workflow_results.mat`

These result files, SEG-Y files, and large MAT data products are ignored by Git by default. If you want to publish example data, prepare a small anonymized sample and explicitly allow it in `.gitignore`.

## Notes

- This repository is a research reproduction and engineering testbed, not a figure-by-figure or number-by-number replication of the original paper.
- Synthetic experiments can compute `Q_full` and `Q_tapered_ref`; real data without a clean reference skip Q metrics.
- `rpca_engineering_workflow.m` is better suited for real-data workflows, while `rpca_tapered_reproduction.m` is better for algorithm comparison and parameter understanding.
- The scripts intentionally keep local functions in single MATLAB files so they are easy to run and move. If the project grows, the functions can be split into a `functions/` directory or a MATLAB package.
