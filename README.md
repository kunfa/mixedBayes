<!-- README.md is generated from README.Rmd. Please edit that file -->

# mixedBayes

> Bayesian Longitudinal Regularized Quantile Mixed Model
<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/mixedBayes)](https://cran.r-project.org/package=mixedBayes)
[![Codecov test
coverage](https://codecov.io/gh/kunfa/mixedBayes/branch/master/graph/badge.svg)](https://app.codecov.io/gh/kunfa/mixedBayes?branch=master)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/mixedBayes)](https://www.r-pkg.org:443/pkg/mixedBayes)
[![R-CMD-check](https://github.com/kunfa/mixedBayes/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kunfa/mixedBayes/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

In longitudinal studies, the same subjects are measured repeatedly over time, leading to correlations among the repeated measurements. Properly accounting for the intra-cluster correlations in the presence of data heterogeneity and long tailed distributions of the disease phenotype is challenging, especially in the context of high dimensional regressions. In this package, we developed a Bayesian quantile mixed effects model with spike- and -slab priors dissect important gene - environment interactions under longitudinal genomics studies. An efficient Gibbs sampler has been developed to facilitate fast computation. The Markov chain Monte Carlo algorithms of the proposed and alternative methods are efficiently implemented in 'C++'. The development of this software package and the associated statistical methods have been partially supported by an Innovative Research Award from Johnson Cancer Research Center, Kansas State University.

## How to install

  - To install from github, run these two lines of code in R

<!-- end list -->

    install.packages("devtools")
    devtools::install_github("kunfa/mixedBayes")

  - Released versions of mixedBayes are available on CRAN
    [(link)](https://cran.r-project.org/package=mixedBayes), and can be
    installed within R via

<!-- end list -->

    install.packages("mixedBayes")

## Examples

#### Example.1 (default method: robust group selection under random intercept and slope model)

    library(mixedBayes)
    data(data)
    
    fit = mixedBayes(y,e,X,g,w,k,structure=c("group"))
    fit$coefficient
    b = selection(fit,sparse=TRUE)
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)
#### Example.2 (alternative: robust individual selection under random intercept and slope model)

    fit = mixedBayes(y,e,X,g,w,k,structure=c("individual"))
    fit$coefficient
    b = selection(fit,sparse=TRUE)
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)

#### Example.3 (alternative: non-robust group selection)

    fit = mixedBayes(y,e,X,g,w,k,robust=FALSE, structure=c("group"))
    fit$coefficient
    b = selection(fit,sparse=TRUE)
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)
#### Example.4 (alternative: robust group selection under random intercept model)
    fit = mixedBayes(y,e,X,g,w,k,slope=FALSE, structure=c("group"))
    fit$coefficient    
    b = selection(fit,sparse=TRUE)
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)
## News

### mixedBayes 0.1.3 \[2024-09-25\]

- Updated with a more efficient Gibbs sampler in the longitudinal setting.
- Updated the inputs for the main mixedBayes function.
- Updated the documentation.
### mixedBayes 0.1.4 \[2025-01-07\]
- Updated README and DESCRIPTION file
  
## Methods

This package provides implementation for methods proposed in

  - Fan, K., Jiang, Y., Ma, S., Wang, W. and Wu, C. (2024+). Robust Sparse Bayesian Regression for Longitudinal Gene-Environment Interactions.(Under Review)
