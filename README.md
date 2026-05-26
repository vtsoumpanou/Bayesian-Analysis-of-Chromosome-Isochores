# Bayesian Analysis of Chromosome Isochores

A Bayesian statistical analysis of GC content in chromosomes to detect isochore boundaries using model selection and Gibbs sampling.

## About

This project analyzes CG base counts across five chromosome segments (chr1 through chr5), each consisting of 100 consecutive windows of 5,000 DNA bases. The objective is to determine whether each chromosome is best described by:

- Model 1 (M1): A single isochore with constant CG proportion across all windows
- Model 2 (M2): Two isochores separated by an unknown changepoint t

A Bayesian framework is adopted throughout, using conjugate Beta priors (Jeffreys prior) for binomial likelihoods. The analysis includes:

1. Model selection using Bayes factors and posterior model probabilities for all five chromosomes
2. Gibbs sampling for full posterior inference on the changepoint model for chromosome 3

## Files

| File | Description |
|------|-------------|
| `Isochore_analysis.R` | Main script: Model specification, Bayes factors, model selection for all 5 chromosomes |
| `gibbs_isochores.R` | Gibbs sampler implementation for chr3 (changepoint estimation) |
| `Report_Bayes1_Tsoumpanou.pdf` | Full report: Model selection results for all chromosomes |
| `Report_Bayes2_Tsoumpanou.pdf` | Full report: Gibbs sampling analysis for chromosome 3 |

## Key Findings

### Model Selection Results (All Chromosomes)

| Chromosome | Posterior (M1) | Posterior (M2) | Log Bayes Factor |
|------------|----------------|----------------|------------------|
| chr1 | ~0.000 | ~1.000 | Strongly positive |
| chr2 | ~0.000 | ~1.000 | Strongly positive |
| chr3 | ~0.000 | ~1.000 | Strongly positive |
| chr4 | ~0.000 | ~1.000 | Strongly positive |
| chr5 | ~0.000 | ~1.000 | Strongly positive |

Conclusion: Decisive Bayesian evidence supporting two isochores for all five chromosomes.

### Gibbs Sampling Results (Chromosome 3)

- Changepoint location (t): Posterior concentrated around window 13
- Interpretation: Lower-GC region (windows 1-13) followed by higher-GC region (windows 14-100)
- p1 - p2: Posterior distribution entirely negative, excludes zero -> strong evidence against single isochore

## Dataset Description

Each chromosome dataset contains 100 observations. For each window of 5,000 DNA bases, the recorded value is the number of C or G bases (CG count).

| Chromosome | Sample size | Distribution |
|------------|-------------|--------------|
| chr1 | 100 | See data in script |
| chr2 | 100 | See data in script |
| chr3 | 100 | See data in script |
| chr4 | 100 | See data in script |
| chr5 | 100 | See data in script |

The data can be modeled as:

$$X_i \sim \text{Binomial}(n = 5000, \theta_i)$$

where:
- $X_i$ = number of CG bases in window i
- n = 5000 = number of bases per window
- $\theta_i$ = underlying CG proportion

## Statistical Methods

### Model 1: Single Isochore

- All windows share a single parameter theta
- Prior: $\theta \sim \text{Beta}(0.5, 0.5)$ (Jeffreys prior)
- Marginal likelihood computed analytically

### Model 2: Two Isochores with Changepoint

- Changepoint $t \in \{1, \ldots, 99\}$
- $\theta_1$ for windows $1, \ldots, t$
- $\theta_2$ for windows $t+1, \ldots, 100$
- Prior: Uniform over $t$
- Evidence obtained by summing over all $t$

### Gibbs Sampler for Chromosome 3

- Parameters: $p_1, p_2, t$
- Priors:
  - $p_1, p_2 \sim \text{Beta}(0.5, 0.5)$ (Jeffreys prior)
  - $t \sim \text{Discrete Uniform}(1, 99)$
- Sampling scheme:
  1\. Sample $p_1$ from Beta full conditional
  2\. Sample $p_2$ from Beta full conditional
  3\. Sample $t$ using normalized posterior probabilities (log-scale for stability)

## Requirements

### Install R and RStudio

1. Download R from https://cran.r-project.org/
2. Download RStudio Desktop (free) from https://posit.co/download/rstudio-desktop/
3. Install R first, then RStudio

### Required Packages

Base R only (no additional packages required for these scripts). The following base functions are used:

- `lbeta()` — Log Beta function for numerically stable marginal likelihood
- `cumsum()` — Cumulative sums for changepoint calculations
- `rbeta()` — Beta distribution sampling for p1 and p2
- `sample()` — Discrete uniform sampling for changepoint t

## Model Specifications

### Joint Posterior (up to constant)

$$
p(p_1, p_2, t \mid \mathbf{x}) \propto \text{Beta}(p_1 \mid a_1 + S_t,  b_1 + n t - S_t) \times \text{Beta}(p_2 \mid a_2 + (S - S_t),  b_2 + n (l - t) - (S - S_t)) \times \mathbb{I}(t \in \{1, \dots, l-1\})
$$

### Conditional Posteriors

$$
\begin{aligned}
p_1 \mid \mathbf{x}, t &\sim \text{Beta}(a_1 + S_t,  b_1 + n t - S_t) \\
p_2 \mid \mathbf{x}, t &\sim \text{Beta}(a_2 + (S - S_t),  b_2 + n (l - t) - (S - S_t)) \\
P(t \mid p_1, p_2, \mathbf{x}) &\propto \left(\frac{p_1}{1-p_1}\right)^{S_t} \times \left(\frac{p_2}{1-p_2}\right)^{S_{t}^{\prime}}
\end{aligned}
$$

Where:
- $S_t = \sum_{i=1}^{t} x_i$ (cumulative sum up to changepoint t)
- $S_{t}^{\prime} = \sum_{i=t+1}^{l} x_i$ (sum after changepoint t)
- $S = S_t + S_{t}^{\prime} = \sum_{i=1}^{l} x_i$ (total sum)
- $l = 100$ (number of windows)
- $n = 5000$ (bases per window)
- $\mathbb{I}(\cdot)$ is the indicator function

Note: Numerical stability is achieved using log-scale calculations in the implementation.

## Numerical Stability

The implementation uses:

- Log-sum-exp trick for marginal likelihood calculations in Model_2
- Log-scale probability calculations for sampling t in the Gibbs sampler
- lbeta() function for stable Beta function computations

## Limitations

- Model assumes at most one changepoint (does not detect multiple isochore boundaries)
- Gibbs sampler implemented for chr3 only (can be extended to other chromosomes)
- Jeffreys prior used for objectivity but may not reflect biological knowledge

## Future Work

- Extend to multiple changepoint detection
- Compare with other model selection criteria (DIC, TIC)
- Visualize posterior predictive distributions

## License

This project is for educational purposes as part of coursework at the National and Kapodistrian University of Athens.


## References

- Bayesian Inference Notes, Loukia Meligkotsidou

## Contact

Vasiliki Tsoumpanou
Email: vasiliki.tsoumpanou@gmail.com  
GitHub: [vtsoumpanou](https://github.com/vtsoumpanou)
