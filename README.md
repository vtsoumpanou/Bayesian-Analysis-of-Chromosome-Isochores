# Bayesian Analysis of Chromosome Isochores



Course: Bayesian Inference  

Institution: National and Kapodistrian University of Athens  

Department: Department of Mathematics  

Student: Tsoumpanou Vasiliki  

Fall Semester 2025-2026



A Bayesian statistical analysis of GC content in chromosomes to detect isochore boundaries using model selection and Gibbs sampling.

## 

## About

This project analyzes CG base counts across five chromosome segments (chr1 through chr5), each consisting of 100 consecutive windows of 5,000 DNA bases. The objective is to determine whether each chromosome is best described by:



\- \*\*Model 1 (M1):\*\* A single isochore with constant CG proportion across all windows

\- \*\*Model 2 (M2):\*\* Two isochores separated by an unknown changepoint t



A Bayesian framework is adopted throughout, using conjugate Beta priors (Jeffreys prior) for binomial likelihoods. The analysis includes:



1\. Model selection using Bayes factors and posterior model probabilities for all five chromosomes

2\. Gibbs sampling for full posterior inference on the changepoint model for chromosome 3

## 

## Files



| File | Description |

|------|-------------|

| `Isochore\_analysis.R` | Main script: Model specification, Bayes factors, model selection for all 5 chromosomes |

| `gibbs\_isochores.R` | Gibbs sampler implementation for chr3 (changepoint estimation) |

| `Report\_Bayes1\_Tsoumpanou.pdf` | Full report: Model selection results for all chromosomes |

| `Report\_Bayes2\_Tsoumpanou.pdf` | Full report: Gibbs sampling analysis for chromosome 3 |



## Key Findings





Model Selection Results (All Chromosomes)



| Chromosome | Posterior (M1) | Posterior (M2) | Log Bayes Factor |

|------------|----------------|----------------|------------------|

| chr1 | \~0.000 | \~1.000 | Strongly positive |

| chr2 | \~0.000 | \~1.000 | Strongly positive |

| chr3 | \~0.000 | \~1.000 | Strongly positive |

| chr4 | \~0.000 | \~1.000 | Strongly positive |

| chr5 | \~0.000 | \~1.000 | Strongly positive |



Conclusion: Decisive Bayesian evidence supporting two isochores for all five chromosomes.



\### Gibbs Sampling Results (Chromosome 3)



\- Changepoint location (t):\*\* Posterior concentrated around window 13

\- Interpretation: Lower-GC region (windows 1-13) followed by higher-GC region (windows 14-100)

\- p1 - p2: Posterior distribution entirely negative, excludes zero -> strong evidence against single isochore



## 

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



X\_i \~ Binomial(n = 5000, theta\_i)



where:

\- X\_i = number of CG bases in window i

\- n = 5000 = number of bases per window

\- theta\_i = underlying CG proportion



## Statistical Methods



Model 1: Single Isochore



\- All windows share a single parameter theta

\- Prior: theta \~ Beta(0.5, 0.5) (Jeffreys prior)

\- Marginal likelihood computed analytically



Model 2: Two Isochores with Changepoint



\- Changepoint t in {1, ..., 99}

\- theta\_1 for windows 1..t

\- theta\_2 for windows t+1..100

\- Prior: Uniform over t

\- Evidence obtained by summing over all t



Gibbs Sampler for Chromosome 3



\- Parameters: p1, p2, t

\- Priors:

&#x20; - p1, p2 \~ Beta(0.5, 0.5) (Jeffreys prior)

&#x20; - t \~ Discrete Uniform(1, 99)

\- Sampling scheme:

&#x20; 1. Sample p1 from Beta full conditional

&#x20; 2. Sample p2 from Beta full conditional

&#x20; 3. Sample t using normalized posterior probabilities (log-scale for stability)



## Requirements



\### Install R and RStudio



1\. Download R from https://cran.r-project.org/

2\. Download RStudio Desktop (free) from https://posit.co/download/rstudio-desktop/

3\. Install R first, then RStudio



\### Required Packages



Base R only (no additional packages required for these scripts). The following base functions are used:



\- lbeta() for log Beta functions

\- cumsum() for cumulative sums

\- rbeta() for Beta sampling

\- sample() for discrete sampling



## Usage



\### Clone the repository



git clone https://github.com/vtsoumpanou/bayesian-isochore-analysis.git

cd bayesian-isochore-analysis



\### Run Model Selection (All Chromosomes)



source("Isochore\_analysis.R")



\# Run model selection for all chromosomes

data <- list(chr1, chr2, chr3, chr4, chr5)

lapply(data, Model\_Selection)



\# Compute log Bayes factors

lapply(data, log\_BF)



\# Visualize CG content patterns

lapply(data, plots)

lapply(data, data\_hist)



\### Run Gibbs Sampler (Chromosome 3 Only)



source("gibbs\_isochores.R")



\# Run Gibbs sampler

draws <- gibbs(data = chr3, ndraw = 1000)



\# Extract parameters

p1 <- draws\[, 1]

p2 <- draws\[, 2]

t <- draws\[, 3]



\# Trace plots

par(mfrow = c(1, 3))

plot(p1, type = "l", main = "Draws p1")

plot(p2, type = "l", main = "Draws p2")

plot(t, type = "l", main = "Draws t")



\# Histograms

par(mfrow = c(2, 2))

hist(t, probab = TRUE, main = "Posterior of Changepoint t")

hist(p1, probab = TRUE, main = "Posterior of p1")

hist(p2, probab = TRUE, main = "Posterior of p2")

hist(p1 - p2, probab = TRUE, main = "Posterior of p1 - p2")



\## Model Specifications



\### Joint Posterior (up to constant)



p(p1, p2, t | x) proportional to Beta(p1 | a1 + S\_t, b1 + n\*t - S\_t) x Beta(p2 | a2 + (S - S\_t), b2 + n\*(l-t) - (S - S\_t)) x I(t in {1,...,l-1})



\### Conditional Posteriors



\- p1 | x, t \~ Beta(a1 + S\_t, b1 + n\*t - S\_t)

\- p2 | x, t \~ Beta(a2 + (S - S\_t), b2 + n\*(l-t) - (S - S\_t))

\- P(t | p1, p2, x) proportional to (p1/(1-p1))^S\_t x (p2/(1-p2))^(S'\_t)



Note: Numerical stability is achieved using log-scale calculations in the implementation.



## Project Structure



bayesian-isochore-analysis/

├── Isochore\_analysis.R

├── gibbs\_isochores.R

├── Report\_Bayes1\_Tsoumpanou.pdf

├── Report\_Bayes2\_Tsoumpanou.pdf

├── README.md

└── .gitignore





## Numerical Stability



The implementation uses:



\- Log-sum-exp trick for marginal likelihood calculations in Model\_2

\- Log-scale probability calculations for sampling t in the Gibbs sampler

\- lbeta() function for stable Beta function computations



## Limitations



\- Model assumes at most one changepoint (does not detect multiple isochore boundaries)

\- Gibbs sampler implemented for chr3 only (can be extended to other chromosomes)

\- Jeffreys prior used for objectivity but may not reflect biological knowledge



## Future work



\- Extend to multiple changepoint detection

\- Compare with other model selection criteria (DIC, TIC)

\- Visualize posterior predictive distributions



## License



This project is for educational purposes as part of coursework at the National and Kapodistrian University of Athens.





## References



\- Bayesian Inference Notes, Loukia Meligkotsidou

## 

## Contact



Vasiliki Tsoumpanou
Email: vasiliki.tsoumpanou@gmail.com  
GitHub: [your-username](https://github.com/your-username)



