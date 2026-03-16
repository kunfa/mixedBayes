#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BGL_1 (arma::vec y, arma::mat e, arma::mat g, arma:: mat w, unsigned int q,unsigned int o,unsigned int k, int maxSteps, arma::vec hatM, arma::vec hatR0, arma::vec hatRStar,arma::mat hatAta, arma::mat z, arma::vec hatInvSigM0, arma::vec hatInvTauSq0, arma::vec hatInvTauSqStar,double hatLambdaSq0, double hatLambdaSqStar, double hatSigmaSq, double a0, double b0, double aStar, double bStar, double Phi1Sq, double alpha, double gamma, double alpha1, double gamma1, int progress)
{
  unsigned int L = q-o, n = g.n_rows,  s = g.n_cols, c = z.n_cols, n1 = n/k;
  arma::mat gsM(maxSteps, q),
  gsR0(maxSteps, s),
  gsRStar(maxSteps, s*L),
  gsInvTauSq0(maxSteps, s),
  gsAta(maxSteps,n1*c),
  gsInvTauSqStar(maxSteps, s);

  arma::vec gsLambda0(maxSteps),
  gsLambdaStar(maxSteps),
  gsSigmaSq(maxSteps),
  gsPhi1Sq(maxSteps);

  arma::mat tBmBm = e.t()*e, tB0B0 = g.t()*g;
  arma::vec tB0B0Diag = tB0B0.diag();

  arma::mat invSigM0 = arma::diagmat(hatInvSigM0);

  arma::mat Xr, varM, varRs, tempS, matRStar;
  arma::vec res, BrjtRes, meanM,  meanAlpha, meanRs, tRsRs, repInvTau, muInvTauSq0, muInvTauSqStar; // mu_m, mu_alpha,
  double temp0, meanR0, varR0, B0jtRes, lInvTauSq0, lInvTauSqStar;

    std::vector<arma::mat> tBrBr(s);
  for(unsigned int j=0; j<s; j++){
    Xr = w.cols((j*L), (j*L+L-1));
    tBrBr[j] = Xr.t()*Xr;
  }

  for (int t = 0; t < maxSteps; t++) {
    // m|y, r0, r.star
     varM = arma::inv_sympd(tBmBm/hatSigmaSq + invSigM0);
    arma::mat Zblock(n, c * n1, arma::fill::zeros);
    for (unsigned int i = 0; i < n1; i++) {
      Zblock.submat(i * k, c * i, i * k + k - 1, c * i + c - 1) = z;
    }
    res = y - (g * hatR0 + w * hatRStar+ Zblock * arma::vectorise(hatAta));

    meanM = varM * (e.t() * res/hatSigmaSq);
    hatM = mvrnormCpp(meanM, varM);
    res -= e * hatM;
    gsM.row(t) = hatM.t();

    // ata|

    arma::vec z0 = z.col(0);  // k x 1

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


      // write updated block back
      res.subvec(i*k, i*k + k - 1) = rblock;
    }

    gsAta.row(t) = arma::vectorise(hatAta).t();

    for(unsigned int j=0; j<s; j++){
      temp0 = 1/(tB0B0Diag(j) + hatInvTauSq0(j));
      varR0 = hatSigmaSq * temp0;
      res += g.col(j) * hatR0(j);
      B0jtRes = arma::as_scalar(g.col(j).t() * res);
      meanR0 = temp0 * B0jtRes;
      hatR0(j) = R::rnorm(meanR0, sqrt(varR0));
      res -= g.col(j) * hatR0(j);


      tempS = tBrBr[j];
      tempS.diag() += hatInvTauSqStar(j);
      tempS = arma::inv_sympd(tempS);
      varRs = hatSigmaSq * tempS;
      res += w.cols((j*L), (j*L+L-1)) * hatRStar.subvec((j*L), (j*L+L-1));
      BrjtRes = w.cols((j*L), (j*L+L-1)).t() * res;
      meanRs = tempS * BrjtRes;
      hatRStar.subvec((j*L), (j*L+L-1)) = mvrnormCpp(meanRs, varRs);
      res -= w.cols((j*L), (j*L+L-1)) * hatRStar.subvec((j*L), (j*L+L-1));
    }
    gsR0.row(t) = hatR0.t();
    gsRStar.row(t) = hatRStar.t();


    // sigma.sq|
    double shapeSig = alpha + (n+s+s*L)/2;
    repInvTau = arma::vectorise(arma::repelem(hatInvTauSqStar.t(), L, 1), 0);
    double rateSig = gamma + 0.5*(arma::accu(arma::square(res)) +
                                  arma::accu(square(hatR0) % hatInvTauSq0) +

                                  arma::accu(square(hatRStar) % repInvTau));
    hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
    gsSigmaSq(t) = hatSigmaSq;

    // invTAUsq.0|lambda, r0
    lInvTauSq0 = hatLambdaSq0;
    muInvTauSq0 = sqrt(hatLambdaSq0 * hatSigmaSq / square(hatR0));
    for(unsigned int j = 0; j < s; j++){
      hatInvTauSq0(j) = rinvgaussian(muInvTauSq0(j), lInvTauSq0);
    }
    gsInvTauSq0.row(t) = hatInvTauSq0.t();


    // invTAUsq.star|lambda.star, r.star
    lInvTauSqStar = L * hatLambdaSqStar;
    matRStar = arma::reshape(hatRStar, L, s);
    tRsRs = sum(square(matRStar), 0).t();
    muInvTauSqStar = sqrt(L * hatLambdaSqStar * hatSigmaSq / tRsRs);
    for(unsigned int j = 0; j<s; j++){
      hatInvTauSqStar(j) = rinvgaussian(muInvTauSqStar(j), lInvTauSqStar);
    }
    gsInvTauSqStar.row(t) = hatInvTauSqStar.t();

    // lambda0|invTAUsq.0
    double shape = a0 + s;
    double rate = b0 + arma::accu(1/hatInvTauSq0)/2;
    hatLambdaSq0 = R::rgamma(shape, 1/rate);
    gsLambda0(t) = hatLambdaSq0;


    // lambda.star|invTAUsq.star
    double shapeS = aStar + s*(L+1)/2;
    double rateS = bStar + L*arma::accu(1/hatInvTauSqStar)/2;
    hatLambdaSqStar = R::rgamma(shapeS, 1/rateS);
    gsLambdaStar(t) = hatLambdaSqStar;

    //phi1sq;
    double diff1 = 0.5 * arma::accu( arma::square(hatAta.row(0)) );

    double shape1 = alpha1 + n1 / 2;
    double rate1  = gamma1 + diff1;
    Phi1Sq = 1 / R::rgamma(shape1, 1 / rate1);

    gsPhi1Sq(t) = Phi1Sq;

    if(progress != 0 && t % progress == 0){
      Rcpp::checkUserInterrupt();
      Rcpp::Rcout << "Iteration: " << t << std::endl;
      Rcpp::Rcout << "  mse    : " << arma::accu(arma::square(res))/n << std::endl;
      Rcpp::Rcout << "  sigmaSq: " << hatSigmaSq << std::endl;
    }
  }

  return Rcpp::List::create(Rcpp::Named("GS.alpha") = gsM,

                            Rcpp::Named("GS.beta") = gsR0,

                            Rcpp::Named("GS.eta") = gsRStar,
                            Rcpp::Named("GS.ata") = gsAta)
                            ;
}


