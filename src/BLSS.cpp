#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include<vector>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BLSS(arma::vec y, arma::mat e, arma::mat g, arma:: mat w, arma:: mat z, unsigned int q,unsigned int k, int maxSteps, arma::vec hatAlpha, arma::vec hatBeta, arma::mat hatAta, arma::vec hatEta, arma::vec hatInvSigM0, arma::vec hatInvTauSq0, arma::vec hatInvTauSqStar, double hatPi0, double hatPiStar,double hatLambdaSq0, double hatLambdaSqStar, double hatSigmaSq, double Phi1Sq, double Phi2Sq, double a0, double b0, double aStar, double bStar, double alpha, double gamma, double alpha1,double gamma1, double mu0, double muStar, double nu0, double nuStar, int progress)
{
  unsigned int n = g.n_rows, m = g.n_cols, c = z.n_cols,p = w.n_cols, n1 = n/k;
  arma::mat gsAlpha(maxSteps, q),
  gsBeta(maxSteps, m),
  gsEta(maxSteps, p),
  gsInvTauSq0(maxSteps, m),
  gsAta(maxSteps,n1*c),
  gsInvTauSqStar(maxSteps, p),
  gsL0(maxSteps, m),
  gsLS(maxSteps, p);

  arma::vec gsLambda0(maxSteps),
  gsLambdaStar(maxSteps),
  gsSigmaSq(maxSteps),
  gsPi0(maxSteps),
  gsPhi1Sq(maxSteps),
  gsPhi2Sq(maxSteps),
  gsPiStar(maxSteps)

    ;

  arma::mat tBmBm = e.t()*e, tB0B0 = g.t()*g, tBrBr = w.t()*w;
  arma::vec tB0B0Diag = tB0B0.diag(),tBrBrDiag = tBrBr.diag();
  arma::mat invSigM0 = arma::diagmat(hatInvSigM0);

  arma::mat varM;
  arma::vec res, meanM, muInvTauSq0, muInvTauSqStar; // mu_m, mu_alpha,
  double temp0, tempS, varR0,varRs, B0jtRes,BrjtRes, meanR0, meanRs,l0, lS, u,u1, lInvTauSq0, lInvTauSqStar;


  for (int t = 0; t < maxSteps; t++) {
    // m|y, r0, r.star
    arma::mat Zblock(n, c * n1, arma::fill::zeros);
    for (unsigned int i = 0; i < n1; i++) {
      Zblock.submat(i * k, c * i, i * k + k - 1, c * i + c - 1) = z;
    }

    varM = arma::inv_sympd(tBmBm/hatSigmaSq + invSigM0);
    res = y - (g * hatBeta + w * hatEta+ Zblock * arma::vectorise(hatAta));
    meanM = varM * (e.t() * res/hatSigmaSq);
    hatAlpha = mvrnormCpp(meanM, varM);
    res -= e * hatAlpha;
    gsAlpha.row(t) = hatAlpha.t();

    // ata|
    arma::vec z0 = z.col(0);  // k x 1
    arma::vec z1 = z.col(1);  // k x 1

    for (unsigned int i = 0; i < n1; i++) {


      arma::vec rblock = res.subvec(i*k, i*k + k - 1);

      // ----- ata0 -----
      double old_ata0 = hatAta(0, i);
      rblock += z0 * old_ata0;

      double t00 = arma::dot(z0, z0) / hatSigmaSq;
      double b0  = arma::dot(z0, rblock) / hatSigmaSq;
      double var0  = 1 / (t00 + 1 / Phi1Sq);
      double mean0 = var0 * b0;

      hatAta(0, i) = R::rnorm(mean0, std::sqrt(var0));
      rblock -= z0 * hatAta(0, i);

      // ----- ata1 -----
      double old_ata1 = hatAta(1, i);
      rblock += z1 * old_ata1;

      double t11 = arma::dot(z1, z1) / hatSigmaSq;
      double b1  = arma::dot(z1, rblock) / hatSigmaSq;
      double var1  = 1 / (t11 + 1 / Phi2Sq);
      double mean1 = var1 * b1;

      hatAta(1, i) = R::rnorm(mean1, std::sqrt(var1));
      rblock -= z1 * hatAta(1, i);

      // write updated block back
      res.subvec(i*k, i*k + k - 1) = rblock;
    }

    gsAta.row(t) = arma::vectorise(hatAta).t();

    for(unsigned int j=0; j<m; j++){
      temp0 = 1/(tB0B0Diag(j) + hatInvTauSq0(j));
      varR0 = hatSigmaSq * temp0;
      res += g.col(j) * hatBeta(j);
      B0jtRes = arma::as_scalar(g.col(j).t() * res);
      meanR0 = temp0 * B0jtRes;
      l0 = hatPi0/(hatPi0 + (1-hatPi0)*sqrt(hatInvTauSq0(j)*temp0)*exp(0.5/hatSigmaSq*temp0*pow(B0jtRes,2)));
      gsL0(t,j) = l0;
      u = R::runif(0, 1);
      if(u<l0){
        hatBeta(j) = 0;
      }else{
        hatBeta(j) = R::rnorm(meanR0, sqrt(varR0));

      }
      res -= g.col(j) * hatBeta(j);
    }
    gsBeta.row(t) = hatBeta.t();

    for(unsigned int j=0; j<p; j++){
      tempS = 1/(tBrBrDiag(j) + hatInvTauSqStar(j));
      varRs = hatSigmaSq * tempS;
      res += w.col(j) * hatEta(j);
      BrjtRes = arma::as_scalar(w.col(j).t() * res);
      meanRs = tempS * BrjtRes;
      lS = hatPiStar/(hatPiStar + (1-hatPiStar)*sqrt(hatInvTauSqStar(j)*tempS)*exp(0.5/hatSigmaSq*tempS*pow(BrjtRes,2)));
      gsLS(t,j) = lS;
      u1 = R::runif(0, 1);
      if(u1<lS){
        hatEta(j) = 0;
      }else{
        hatEta(j) = R::rnorm(meanRs, sqrt(varRs));

      }
      res -= w.col(j) * hatEta(j);
    }
    gsEta.row(t) = hatEta.t();


    // invTAUsq.0|lambda, r0
    lInvTauSq0 = hatLambdaSq0;
    muInvTauSq0 = sqrt(hatLambdaSq0 * hatSigmaSq / square(hatBeta));
    for(unsigned int j = 0; j < m; j++){
      if(hatBeta(j) == 0){
        hatInvTauSq0(j) = 1/R::rgamma(1, 2/lInvTauSq0);
      }else{
        hatInvTauSq0(j) = rinvgaussian(muInvTauSq0(j), lInvTauSq0);
      }
    }
    gsInvTauSq0.row(t) = hatInvTauSq0.t();


    // invTAUsq.star|lambda.star, r.star
    lInvTauSqStar = hatLambdaSqStar;

    muInvTauSqStar = arma::sqrt(hatLambdaSqStar * hatSigmaSq / square(hatEta));
    for(unsigned int j = 0; j<p; j++){
      if(hatEta(j) == 0){
        hatInvTauSqStar(j) = 1/R::rgamma(1, 2/lInvTauSqStar);
      }else{
        hatInvTauSqStar(j) = rinvgaussian(muInvTauSqStar(j), lInvTauSqStar);
      }
    }
    gsInvTauSqStar.row(t) = hatInvTauSqStar.t();


    // lambda0|invTAUsq.0
    double shape = a0 + m;
    double rate = b0 + arma::accu(1/hatInvTauSq0)/2;
    hatLambdaSq0 = R::rgamma(shape, 1/rate);
    gsLambda0(t) = hatLambdaSq0;


    // lambda.star|invTAUsq.star
    double shapeS = aStar + p;
    double rateS = bStar + arma::accu(1/hatInvTauSqStar)/2;
    hatLambdaSqStar = R::rgamma(shapeS, 1/rateS);
    gsLambdaStar(t) = hatLambdaSqStar;

    // sigma.sq|
    double shapeSig = alpha + n/2 + arma::accu(hatEta != 0)/2 + arma::accu(hatBeta != 0)/2 ;

    double rateSig = gamma + 0.5*(arma::accu(arma::square(res)) +
                                  arma::accu(square(hatBeta) % hatInvTauSq0) +

                                  arma::accu(square(hatEta) % hatInvTauSqStar));
    hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
    gsSigmaSq(t) = hatSigmaSq;

    //phi1sq, phi2sq;
    double diff1 = 0.5 * arma::accu( arma::square(hatAta.row(0)) );
    double diff2 = 0.5 * arma::accu( arma::square(hatAta.row(1)) );

    double shape1 = alpha1 + n1 / 2;
    double rate1  = gamma1 + diff1;
    Phi1Sq = 1 / R::rgamma(shape1, 1 / rate1);

    double shape11 = alpha1 + n1 / 2;
    double rate11  = gamma1 + diff2;
    Phi2Sq = 1 / R::rgamma(shape11, 1 / rate11);

    gsPhi1Sq(t) = Phi1Sq;
    gsPhi2Sq(t) = Phi2Sq;


    // pi.0|
    double shape1_0 = mu0 + arma::accu(hatBeta == 0);
    double shape2_0 = nu0 + arma::accu(hatBeta != 0);
    hatPi0 = R::rbeta(shape1_0, shape2_0);
    gsPi0(t) = hatPi0;


    // pi.star|
    double shape1_s = muStar + arma::accu(hatEta == 0);
    double shape2_s = nuStar + arma::accu(hatEta != 0);
    hatPiStar = R::rbeta(shape1_s, shape2_s);
    gsPiStar(t) = hatPiStar;

    if(progress != 0 && t % progress == 0){
      Rcpp::checkUserInterrupt();
      Rcpp::Rcout << "Iteration: " << t << std::endl;
      Rcpp::Rcout << "  mse    : " << arma::accu(arma::square(res))/n << std::endl;
      Rcpp::Rcout << "  sigmaSq: " << hatSigmaSq << std::endl;
      // Rcpp::Rcout << "  pi0    : " << hatPi0 << std::endl;
      // Rcpp::Rcout << "  piStar : " << hatPiStar << std::endl;
      // Rcpp::Rcout << "  piZeta : " << hatPiZeta << std::endl;
    }
  }
  return Rcpp::List::create(Rcpp::Named("GS.alpha") = gsAlpha,

                            Rcpp::Named("GS.beta") = gsBeta,
                            Rcpp::Named("GS.ata") = gsAta,

                            Rcpp::Named("GS.eta") = gsEta

  );

}
