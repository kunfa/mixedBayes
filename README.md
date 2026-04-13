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
    
## Data generation under random intercept-and-slope model with t(2) error

    library(MASS)
    library(mixedBayes)
    Data <- function (n,p,k,q,quant){
    sigma2=1
    # Generate genetic factors
    sig = matrix(0,p,p)
    for (i in 1:p)
    {
    for(j in 1:p)
    {
      sig[i,j] = 0.5^abs(i-j)
    }
    }
  
    g = mvrnorm(n,rep(0,p),sig)
    g = as.matrix(g)
    g = scale(g)
  
    # True main effects for genetic variables 
    beta_true = rep(0,p)
    beta_true[c(2,4,7,12)] = runif(4,0.4,0.8)
  
    # Generate environmental factors (a group of 3 dummy variables)
    group <- sample(0:3, size = n, replace = TRUE)
  
    e1 <- as.numeric(group == 1)
    e2 <- as.numeric(group == 2)
    e3 <- as.numeric(group == 3)
  
    e <- cbind(e1, e2, e3)
    e <- scale(e)
  
    # Fixed effects for environmental main effects
    alpha1= runif(q,0.4,0.8)
  
    # Fixed effects for time (polynomial terms)
    alpha2= runif(3,0.4,0.8)
  
    # Interaction terms
    w = c()
  
    for (i in 1:ncol(g))
    {
    
    w = cbind(w,g[,i]*e)
    
    }
  
    # True interaction effects
    eta_true <- rep(0, p * q)
    nz_block <- list(4:6, 10:12, 16:18, 19:21, 34:36, 43:45)
    for (idx in nz_block) {
     eta_true[idx] <- runif(q, 0.4, 0.8)
    }
  
    # Combine all coefficients
    betas_true <- c(beta_true, eta_true)
  
    #  Random effects design 
    k_1 <- 1:k
    k_1 <- k_1 - mean(k_1)
    z <- cbind(1, k_1)
  
    # fixed-effect design for time
    xi <- cbind(1, 1:k, (1:k)^2)
    X <- do.call(rbind, replicate(n, xi, simplify = FALSE))
  
    # Generate response 
    y <- matrix(0, n, k)
    for (i in 1:n) {
    ata <- runif(2,0,1)
    ei <- matrix(rep(e[i, ], each = k), nrow = k)
    gi <- matrix(rep(g[i, ], each = k), nrow = k)
    wi <- matrix(rep(w[i, ], each = k), nrow = k)
    error <- rt(k, 2)
    error <- error -quantile(error, probs = quant)
    y[i, ] <- ei %*% alpha1 +xi %*% alpha2 +gi %*% beta_true +wi %*% eta_true +z %*% ata + error
    }
  
    # Convert to long format
    y = reformat(k,y,type="r");g = reformat(k,g,type="d");e = reformat(k,e,type="d");w = reformat(k,w,type="d")
    dat <- list(y = y, g = g, e  = e, w  = w ,X  = X, coef = betas_true)
    return (dat)
    }

#### Example.1 (proposed method: robust sparse bi-level selection under random intercept -and- slope model)
    
    fit = mixedBayes(y,e,X,g,k,structure="bilevel")
    
    # Estimated coefficients(posterior median)
    fit$coefficient
    
    # Identification
    beta_est = selection(fit,sparse=TRUE)
    coeff = data$coef
    index = which(coeff!=0) # true active predictors
    pos = which(beta_est != 0) # selected active predictors
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    
    # Prediction
    prediction=predict_mixedBayes(fit,y,X,e,g,k,slope=TRUE,loss = "L1")
    prediction
    
#### Example.2 (alternative: robust sparse individual level selection under random intercept -and- slope model)

    fit = mixedBayes(y,e,X,g,k,structure="individual")
    
    # Estimated coefficients(posterior median)
    fit$coefficient
    
    # Identification
    b = selection(fit,sparse=TRUE)
    coeff = data$coef
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)
    
    # Prediction
    prediction=predict_mixedBayes(fit,y,X,e,g,k,slope=TRUE,loss = "L1")
    prediction

#### Example.3 (alternative: non-robust sparse bi-level selection under random intercept -and- slope model)

    fit = mixedBayes(y,e,X,g,k,robust=FALSE, quant = NULL,structure = "bilevel")
    
    # Estimated coefficients(posterior median)
    fit$coefficient
    
    # Identification
    b = selection(fit,sparse=TRUE)
    coeff = data$coef
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)
    
    # Prediction
    prediction=predict_mixedBayes(fit,y,X,e,g,k,slope=TRUE,loss = "L2")
    prediction
    
#### Example.4 (alternative: robust bi-level selection under random intercept -and- slope model)

    fit = mixedBayes(y,e,X,g,k,robust=TRUE,sparse = FALSE,structure = "bilevel")
    
    # Estimated coefficients(posterior median)
    fit$coefficient
    
    # Identification
    b = selection(fit,sparse=FALSE)
    coeff = data$coef
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)
    
    # Prediction
    prediction=predict_mixedBayes(fit,y,X,e,g,k,slope=TRUE,loss = "L1")
    prediction
    
## Data generation under random intercept model with t(2) error

    library(MASS)
    library(mixedBayes)
    Data <- function (n,p,k,q,quant){
    sigma2=1
    # Generate genetic factors
    sig = matrix(0,p,p)
    for (i in 1:p)
    {
    for(j in 1:p)
    {
      sig[i,j] = 0.5^abs(i-j)
    }
    }
  
    g = mvrnorm(n,rep(0,p),sig)
    g = as.matrix(g)
    g = scale(g)
  
    # True main effects for genetic variables 
    beta_true = rep(0,p)
    beta_true[c(2,4,7,12)] = runif(4,0.4,0.8)
  
    # Generate environmental factors (a group of 3 dummy variables)
    group <- sample(0:3, size = n, replace = TRUE)
  
    e1 <- as.numeric(group == 1)
    e2 <- as.numeric(group == 2)
    e3 <- as.numeric(group == 3)
  
    e <- cbind(e1, e2, e3)
    e <- scale(e)
  
    # Fixed effects for environmental main effects
    alpha1= runif(q,0.4,0.8)
  
    # Fixed effects for time (polynomial terms)
    alpha2= runif(3,0.4,0.8)
  
    # Interaction terms
    w = c()
  
    for (i in 1:ncol(g))
    {
    
    w = cbind(w,g[,i]*e)
    
    }
  
    # True interaction effects
    eta_true <- rep(0, p * q)
    nz_block <- list(4:6, 10:12, 16:18, 19:21, 34:36, 43:45)
    for (idx in nz_block) {
     eta_true[idx] <- runif(q, 0.4, 0.8)
    }
  
    # Combine all coefficients
    betas_true <- c(beta_true, eta_true)
  
    #  Random effects design 
   
    z <- rep(1, k)
  
    # fixed-effect design for time
    xi <- cbind(1, 1:k, (1:k)^2)
    X <- do.call(rbind, replicate(n, xi, simplify = FALSE))
  
    # Generate response 
    y <- matrix(0, n, k)
    for (i in 1:n) {
    ata <- runif(1,0,1)
    ei <- matrix(rep(e[i, ], each = k), nrow = k)
    gi <- matrix(rep(g[i, ], each = k), nrow = k)
    wi <- matrix(rep(w[i, ], each = k), nrow = k)
    error <- rt(k, 2)
    error <- error -quantile(error, probs = quant)
    y[i, ] <- ei %*% alpha1 +xi %*% alpha2 +gi %*% beta_true +wi %*% eta_true +z * ata + error
    }
  
    # Convert to long format
    y = reformat(k,y,type="r");g = reformat(k,g,type="d");e = reformat(k,e,type="d");w = reformat(k,w,type="d")
    dat <- list(y = y, g = g, e  = e, w  = w ,X  = X, coef = betas_true)
    return (dat)
    }

    
#### Example.1 (proposed method: robust sparse bi-level selection under random intercept model)

    fit = mixedBayes(y,e,X,g,k,slope=FALSE, structure="bilevel")
    
    # Estimated coefficients(posterior median)
    fit$coefficient
    
    # Identification
    b = selection(fit,sparse=TRUE)
    coeff = data$coef
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)
    
    # Prediction
    prediction=predict_mixedBayes(fit,y,X,e,g,k,slope=FALSE,loss = "L1")
    prediction
    
#### Example.2 (alternative: robust sparse individual level selection under random intercept model)

    fit = mixedBayes(y,e,X,g,k,slope=FALSE, structure="individual")
    
    # Estimated coefficients(posterior median)
    fit$coefficient
    
    # Identification
    b = selection(fit,sparse=TRUE)
    coeff = data$coef
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)
    
    # Prediction
    prediction=predict_mixedBayes(fit,y,X,e,g,k,slope=FALSE,loss = "L1")
    prediction

#### Example.3 (alternative: non-robust sparse bi-level selection under random intercept model)

    fit = mixedBayes(y,e,X,g,k,slope=FALSE,robust=FALSE, quant = NULL,structure = "bilevel")
    
    # Estimated coefficients(posterior median)
    fit$coefficient
    
    # Identification
    b = selection(fit,sparse=TRUE)
    coeff = data$coef
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)
    
    # Prediction
    prediction=predict_mixedBayes(fit,y,X,e,g,k,slope=FALSE,loss = "L2")
    prediction
    
#### Example.4 (alternative: robust bi-level selection under random intercept model)

    fit = mixedBayes(y,e,X,g,k,slope=FALSE,robust=TRUE,sparse = FALSE,structure = "bilevel")
    
    # Estimated coefficients(posterior median)
    fit$coefficient
    
    # Identification
    b = selection(fit,sparse=FALSE)
    coeff = data$coef
    index = which(coeff!=0)
    pos = which(b != 0)
    tp = length(intersect(index, pos))
    fp = length(pos) - tp
    list(tp=tp, fp=fp)
    
    # Prediction
    prediction=predict_mixedBayes(fit,y,X,e,g,k,slope=FALSE,loss = "L1")
    prediction
  

## Methods

This package provides implementation for methods proposed in

  - Fan, K., Jiang, Y., Ma, S., Wang, W. and Wu, C. (2025). Robust Sparse Bayesian Regression for Longitudinal Gene-Environment Interactions. [Journal of the Royal Statistical Society Series C: Applied Statistics, 74(5), 1372–1394.](https://doi.org/10.1093/jrsssc/qlaf027)
