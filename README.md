<!-- README.md is generated from README.Rmd. Please edit that file -->

# mixedBayes

> Bayesian Longitudinal Regularized Quantile Mixed Model
<!-- badges: start -->

[![CRAN](https://www.r-pkg.org/badges/version/mixedBayes)](https://cran.r-project.org/package=mixedBayes)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/mixedBayes)](https://www.r-pkg.org:443/pkg/mixedBayes)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/last-month/mixedBayes)](https://www.r-pkg.org:443/pkg/mixedBayes)

<!-- badges: end -->

With high-dimensional omics features, repeated measure ANOVA leads to longitudinal gene-environment interaction studies that have intra-cluster correlations, outlying observations and structured sparsity arising from the ANOVA design. In this package, we have developed robust sparse Bayesian mixed effect models tailored for the above studies ([Fan et al. (2025)](https://doi.org/10.1093/jrsssc/qlaf027)). An efficient Gibbs sampler has been developed to facilitate fast computation. The Markov chain Monte Carlo algorithms of the proposed and alternative methods are efficiently implemented in 'C++'. The development of this software package and the associated statistical methods have been partially supported by an Innovative Research Award from Johnson Cancer Research Center, Kansas State University.

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
    
## Load the package

    library(mixedBayes)
    
## Data

    The example data set `data` simulated under random intercept-and-slope model included in the package can be loaded by

    data(data)

## Examples

#### Example.1 (default method: robust sparse bi-level selection under random intercept -and- slope model)
    
    fit = mixedBayes(y,e,X,g,w,k,structure=c("bi-level"))
    
    # Estimated coefficients(posterior median)
    fit$coefficient
    
    # Identification
    b = selection(fit,sparse=TRUE)
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)
    
    # Prediction
    prediction=predict_mixedBayes(fit,y,X,e,g,w,k,slope=TRUE,loss = "L1")
    prediction
    
#### Example.2 (alternative: robust sparse individual level selection under random intercept -and- slope model)

    fit = mixedBayes(y,e,X,g,w,k,structure=c("individual"))
    
    # Estimated coefficients(posterior median)
    fit$coefficient
    
    # Identification
    b = selection(fit,sparse=TRUE)
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)
    
    # Prediction
    prediction=predict_mixedBayes(fit,y,X,e,g,w,k,slope=TRUE,loss = "L1")
    prediction

#### Example.3 (alternative: non-robust sparse bi-level selection under random intercept -and- slope model)

    fit = mixedBayes(y,e,X,g,w,k,robust=FALSE, structure=c("bi-level"))
    
    # Estimated coefficients(posterior median)
    fit$coefficient
    
    # Identification
    b = selection(fit,sparse=TRUE)
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)
    
    # Prediction
    prediction=predict_mixedBayes(fit,y,X,e,g,w,k,slope=TRUE,loss = "L2")
    prediction
    
#### Example.4 (alternative: robust sparse bi-level selection under random intercept model)

    fit = mixedBayes(y,e,X,g,w,k,slope=FALSE, structure=c("bi-level"))
    
    # Estimated coefficients(posterior median)
    fit$coefficient
    
    # Identification
    b = selection(fit,sparse=TRUE)
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)
    
    # Prediction
    prediction=predict_mixedBayes(fit,y,X,e,g,w,k,slope=FALSE,loss = "L1")
    prediction

  
## Methods

This package provides implementation for methods proposed in

  - Fan, K., Jiang, Y., Ma, S., Wang, W. and Wu, C. (2025). Robust Sparse Bayesian Regression for Longitudinal Gene-Environment Interactions. [Journal of the Royal Statistical Society Series C: Applied Statistics, 74(5), 1372–1394.](https://doi.org/10.1093/jrsssc/qlaf027)
