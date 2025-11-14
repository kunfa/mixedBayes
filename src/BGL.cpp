#include<RcppArmadillo.h>
#include<Rmath.h>
#include<stdio.h>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;


// [[Rcpp::export()]]
Rcpp::List BGL (arma::vec y, arma::mat e, arma::mat g, arma:: mat w, unsigned int q,unsigned int o,unsigned int k, int maxSteps, arma::vec hatM, arma::vec hatR0, arma::vec hatRStar,arma::mat hatAta, arma::mat z, arma::vec hatInvSigM0, arma::vec hatInvTauSq0, arma::vec hatInvTauSqStar,double hatLambdaSq0, double hatLambdaSqStar, double hatSigmaSq, double a0, double b0, double aStar, double bStar, double hatPhiSq, double alpha, double gamma, double alpha1, double gamma1, int progress)
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
  gsPhiSq(maxSteps);

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
     varM = arma::inv(tBmBm/hatSigmaSq + invSigM0);
    res = y - (g * hatR0 + w * hatRStar);
    for(unsigned int i=0;i<n1;i++){
      res.subvec((i*k), (i*k+k-1)) -= z*hatAta.col(i);
    }
    meanM = varM * (e.t() * res/hatSigmaSq);
    hatM = mvrnormCpp(meanM, varM);
    res -= e * hatM;
    gsM.row(t) = hatM.t();

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
      tempS = arma::inv(tempS);
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

  return Rcpp::List::create(Rcpp::Named("GS.alpha") = gsM,

                            Rcpp::Named("GS.beta") = gsR0,

                            Rcpp::Named("GS.eta") = gsRStar,
                            Rcpp::Named("GS.ata") = gsAta,
                            Rcpp::Named("GS.invTAUsq.0") = gsInvTauSq0,

                            Rcpp::Named("GS.invTAUsq.star") = gsInvTauSqStar,
                            Rcpp::Named("GS.lambda.sq.0") = gsLambda0,

                            Rcpp::Named("GS.lambda.sq.star") = gsLambdaStar,
                            Rcpp::Named("GS.sigma.sq") = gsSigmaSq);
}


