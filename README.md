<!-- README.md is generated from README.Rmd. Please edit that file -->

# mixedBayes

> Bayesian Longitudinal Regularized Quantile Mixed Model


In longitudinal studies, the same subjects are measured repeatedly over time, leading to correlations among the repeated measurements. Properly accounting for the intra-cluster correlations in the presence of data heterogeneity and long tailed distributions of the disease phenotype is challenging, especially in the context of high dimensional regressions. Here, we aim at developing novel Bayesian regularized quantile mixed effect models to tackle these challenges. We have proposed a Bayesian variable selection in the mixed effect models for longitudinal genomics studies. To dissect important gene - environment interactions, our model can simultaneously identify important main and interaction effects on the individual and group level, which have been facilitated by imposing the spike- and -slab priors through Laplacian shrinkage in the Bayesian quantile hierarchical models. The within - subject dependence among data can be accommodated by incorporating the random effects. An efficient Gibbs sampler has been developed to facilitate fast computation. The Markov chain Monte Carlo algorithms of the proposed and alternative methods are efficiently implemented in 'C++'. The development of this software package and the associated statistical methods have been partially supported by an Innovative Research Award from Johnson Cancer Research Center, Kansas State University.

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
    
    fit = mixedBayes(y,e,C,g,w,k,structure=c("group"))
    fit$coefficient
   
    
    ## Compute TP and FP
    b = selection(fit,sparse=TRUE)
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)

#### Example.2 (alternative: robust individual selection under random intercept and slope model)

    fit = mixedBayes(y,e,C,g,w,k,structure=c("individual"))
    b = selection(fit,sparse=TRUE)
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)

#### Example.3 (alternative: non-robust group selection)

    fit = mixedBayes(y,e,C,g,w,k,robust=FALSE, structure=c("group"))
    b = selection(fit,sparse=TRUE)
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)
#### Example.4 (alternative: robust group selection under random intercept model)

    fit = mixedBayes(y,e,C,g,w,k,slope=FALSE, structure=c("group"))
    b = selection(fit,sparse=TRUE)
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)
## Methods

This package provides implementation for methods proposed in

  - Fan, K., Jiang, Y., Ma, S., Wang, W. and Wu, C. (2023+). Robust Bayesian variable selection of Gene-environment interactions in the longitudinal study.(submitted)
