#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BL (arma::vec y, arma::mat e, arma::mat g, arma:: mat w, unsigned int q, unsigned int k, int maxSteps, arma::vec hatAlpha, arma::vec hatBeta, arma::vec hatEta,arma::mat hatAta, arma::mat z, arma::vec hatInvSigM0, arma::vec hatInvTauSq0, arma::vec hatInvTauSqStar,double hatLambdaSq0, double hatLambdaSqStar, double hatSigmaSq, double a0, double b0, double aStar, double bStar, double hatPhiSq, double alpha, double gamma, double alpha1, double gamma1, int progress)
{
  unsigned int n = g.n_rows, m = g.n_cols, c = z.n_cols,p = w.n_cols, n1 = n/k;
  arma::mat gsAlpha(maxSteps, q),
  gsBeta(maxSteps, m),
  gsEta(maxSteps, p),
  gsInvTauSq0(maxSteps, m),
  gsAta(maxSteps,n1*c),
  gsInvTauSqStar(maxSteps, p);

  arma::vec gsLambda0(maxSteps),
  gsLambdaStar(maxSteps),
  gsSigmaSq(maxSteps),
  gsPhiSq(maxSteps);

  arma::mat tBmBm = e.t()*e, tB0B0 = g.t()*g, tBrBr = w.t()*w;
  arma::vec tB0B0Diag = tB0B0.diag(),tBrBrDiag = tBrBr.diag();
  arma::mat invSigM0 = arma::diagmat(hatInvSigM0);

  arma::mat varM;
  arma::vec res, meanM, muInvTauSq0, muInvTauSqStar; // mu_m, mu_alpha,
  double temp0, tempS, varR0,varRs, B0jtRes,BrjtRes, meanR0, meanRs, lInvTauSq0, lInvTauSqStar;


  for (int t = 0; t < maxSteps; t++) {
    // m|y, r0, r.star

   varM = arma::inv(tBmBm/hatSigmaSq + invSigM0);
    res = y - (g * hatBeta + w * hatEta);
    for(unsigned int i=0;i<n1;i++){
      res.subvec((i*k), (i*k+k-1)) -= z*hatAta.col(i);
    }
    meanM = varM * (e.t() * res/hatSigmaSq);
    hatAlpha = mvrnormCpp(meanM, varM);
    res -= e * hatAlpha;
    gsAlpha.row(t) = hatAlpha.t();

    // ata|

    for(unsigned int i=0;i<n1;i++){
      res.subvec((i*k), (i*k+k-1)) += z * hatAta.col(i);
      arma::mat tzz =  z.t()*z;
      arma::mat A = tzz/hatSigmaSq;
      arma::vec B = z.t()*res.subvec((i*k), (i*k+k-1))/hatSigmaSq;
      arma:: mat T(c,c);
      T = T.eye();
      arma::mat invhatPhiSq = 1/hatPhiSq*T;
      arma::mat varAta = arma::inv(A+invhatPhiSq);
      arma::vec meanAta = varAta*B;
      hatAta.col(i) = mvrnormCpp(meanAta, varAta);
      res.subvec((i*k), (i*k+k-1)) -= z * hatAta.col(i);
    }

    gsAta.row(t) = arma::vectorise(hatAta).t();

    for(unsigned int j=0; j<m; j++){
      temp0 = 1/(tB0B0Diag(j) + hatInvTauSq0(j));
      varR0 = hatSigmaSq * temp0;
      res += g.col(j) * hatBeta(j);
      B0jtRes = arma::as_scalar(g.col(j).t() * res);
      meanR0 = temp0 * B0jtRes;
      hatBeta(j) = R::rnorm(meanR0, sqrt(varR0));
      res -= g.col(j) * hatBeta(j);
    }
    gsBeta.row(t) = hatBeta.t();
    for(unsigned int j=0; j<p; j++){
      tempS = 1/(tBrBrDiag(j) + hatInvTauSqStar(j));
      varRs = hatSigmaSq * tempS;
      res += w.col(j) * hatEta(j);
      BrjtRes = arma::as_scalar(w.col(j).t() * res);
      meanRs = tempS * BrjtRes;
      hatEta(j) = R::rnorm(meanRs, sqrt(varRs));
      res -= w.col(j) * hatEta(j);
    }
    gsEta.row(t) = hatEta.t();


    // sigma.sq|
    double shapeSig = alpha + (n+m+p)/2;

    double rateSig = gamma + 0.5*(arma::accu(arma::square(res)) +
                                  arma::accu(square(hatBeta) % hatInvTauSq0) +

                                  arma::accu(square(hatEta) % hatInvTauSqStar));
    hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
    gsSigmaSq(t) = hatSigmaSq;

    // invTAUsq.0|lambda, r0
    lInvTauSq0 = hatLambdaSq0;
    muInvTauSq0 = sqrt(hatLambdaSq0 * hatSigmaSq / square(hatBeta));
    for(unsigned int j = 0; j < m; j++){
      hatInvTauSq0(j) = rinvgaussian(muInvTauSq0(j), lInvTauSq0);
    }
    gsInvTauSq0.row(t) = hatInvTauSq0.t();


    // invTAUsq.star|lambda.star, r.star
    lInvTauSqStar =  hatLambdaSqStar;
    muInvTauSqStar = sqrt(hatLambdaSqStar * hatSigmaSq / square(hatEta));
    for(unsigned int j = 0; j<p; j++){
      hatInvTauSqStar(j) = rinvgaussian(muInvTauSqStar(j), lInvTauSqStar);
    }
    gsInvTauSqStar.row(t) = hatInvTauSqStar.t();

    // lambda0|invTAUsq.0
    double shape = a0 + m;
    double rate = b0 + arma::accu(1/hatInvTauSq0)/2;
    hatLambdaSq0 = R::rgamma(shape, 1/rate);
    gsLambda0(t) = hatLambdaSq0;


    // lambda.star|invTAUsq.star
    double shapeS = aStar +p;
    double rateS = bStar + arma::accu(1/hatInvTauSqStar)/2;
    hatLambdaSqStar = R::rgamma(shapeS, 1/rateS);
    gsLambdaStar(t) = hatLambdaSqStar;

    //phi;
    double shapePhi, ratePhi;
    shapePhi = alpha1 + n1*c/2;
    double diff;
    diff=0;
    for(unsigned int i=0;i<n1;i++){
      diff= diff+0.5*(arma::accu(square(hatAta.col(i))));
    }

    ratePhi = gamma1 + diff;
    hatPhiSq = 1/R::rgamma(shapePhi, 1/ratePhi);
    gsPhiSq(t) = hatPhiSq;


    if(progress != 0 && t % progress == 0){
      Rcpp::checkUserInterrupt();
      Rcpp::Rcout << "Iteration: " << t << std::endl;
      Rcpp::Rcout << "  mse    : " << arma::accu(arma::square(res))/n << std::endl;
      Rcpp::Rcout << "  sigmaSq: " << hatSigmaSq << std::endl;
    }
  }

  return Rcpp::List::create(Rcpp::Named("GS.alpha") = gsAlpha,

                            Rcpp::Named("GS.beta") = gsBeta,
                            Rcpp::Named("GS.ata") = gsAta,

                            Rcpp::Named("GS.eta") = gsEta,
                            Rcpp::Named("GS.invTAUsq.0") = gsInvTauSq0,

                            Rcpp::Named("GS.invTAUsq.star") = gsInvTauSqStar,
                            Rcpp::Named("GS.lambda.sq.0") = gsLambda0,

                            Rcpp::Named("GS.lambda.sq.star") = gsLambdaStar,
                            Rcpp::Named("GS.sigma.sq") = gsSigmaSq);
}


