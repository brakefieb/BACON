# BACON: Bayesian Clustering of *n*-gons

## Overview

<img src="images/MFM.png" alt="Hover Title" title="Hover Title" width="700"/>

## Mixture of Finite Mixtures

```r
bacon_mfm(NumericMatrix A, NumericMatrix L, double w_A, double w_L, bool doub_dirich, bool ddirch_A, bool ddirch_L, int Kmax_0, bool est_sr, bool est_s, bool est_r, double alpha_s, double beta_s, int iter)
```

- `A`: `m` x `n` matrix of relative interior angles
- `L`: `m` x `n` matrix of relative side lengths
- `w_A`: Weight for `A`
- `w_L`: Weight for `L`
- `doub_dirich`: Boolean flag indicating whether to use the double truncated Dirichlet distribution for both `A` and `L`
- `ddirich_A`: Boolean flag indicating whether to use the (single) truncated Dirichlet distribution for `A`
- `ddirich_L`: Boolean flag indicating whether to use the (single) truncated Dirichlet distribution for `L`
- `K_max_0`: Maximum number of clusters
- `est_sr`: Boolean flag for jointly updating `s` and `r`
- `est_s, est_r`: Boolean flag for separately updating `s` and `r`
- `alpha_s, beta_s`: Hyperparameters for prior on `s`
- `iter`: Number of MCMC iterations
