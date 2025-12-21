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
Rcpp::List BGLSS_1(arma::vec y, arma::mat e, arma::mat g, arma:: mat w, arma:: mat z, unsigned int q, unsigned int o,unsigned int k, int maxSteps, arma::vec hatM, arma::vec hatR0, arma::mat hatAta, arma::vec hatRStar, arma::vec hatInvSigM0, arma::vec hatInvTauSq0, arma::vec hatInvTauSqStar, double hatPi0, double hatPiStar,double hatLambdaSq0, double hatLambdaSqStar, double hatSigmaSq, double hatPhiSq, double a0, double b0, double aStar, double bStar, double alpha, double gamma, double alpha1,double gamma1, double mu0, double muStar, double nu0, double nuStar, int progress)
{
  unsigned int L = q-o, n = g.n_rows,  s = g.n_cols, c = z.n_cols, n1 = n/k;
  arma::mat gsM(maxSteps, q),
  gsR0(maxSteps, s),
  gsRStar(maxSteps, s*L),
  gsRstRs(maxSteps, s),
  gsInvTauSq0(maxSteps, s),
  gsAta(maxSteps,n1*c),
  gsInvTauSqStar(maxSteps, s),
  gsL0(maxSteps, s),
  gsLS(maxSteps, s);

  arma::vec gsLambda0(maxSteps),
  gsLambdaStar(maxSteps),
  gsSigmaSq(maxSteps),
  gsPi0(maxSteps),
  gsPhiSq(maxSteps),
  gsPiStar(maxSteps)
  ;

  arma::mat tBmBm = e.t()*e, tB0B0 = g.t()*g;
  arma::vec tB0B0Diag = tB0B0.diag();
arma::mat invSigM0 = arma::diagmat(hatInvSigM0);

  arma::mat Xr, varM, varRs, tempS, matRStar;
  arma::vec res, meanM, BrjtRes, meanRs, tRsRs, repInvTauSq, residual, muInvTauSq0, muInvTauSqStar; // mu_m, mu_alpha,
  double temp0, varR0, B0jtRes, meanR0, l0, lS, u, lInvTauSq0, lInvTauSqStar;

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
      l0 = hatPi0/(hatPi0 + (1-hatPi0)*sqrt(hatInvTauSq0(j)*temp0)*exp(0.5/hatSigmaSq*temp0*pow(B0jtRes,2)));
      gsL0(k,j) = l0;
      u = R::runif(0, 1);
      if(u<l0){
        hatR0(j) = 0;
      }else{
        hatR0(j) = R::rnorm(meanR0, sqrt(varR0));

      }
      res -= g.col(j) * hatR0(j);


      tempS = tBrBr[j];
      tempS.diag() += hatInvTauSqStar(j);
      tempS = arma::inv(tempS);
      varRs = hatSigmaSq * tempS;
      res += w.cols((j*L), (j*L+L-1)) * hatRStar.subvec((j*L), (j*L+L-1));
      BrjtRes = w.cols((j*L), (j*L+L-1)).t() * res;
      meanRs = tempS * BrjtRes;
      lS = arma::as_scalar(hatPiStar/(hatPiStar+(1-hatPiStar)*pow(hatInvTauSqStar(j),(L/2))*sqrt(arma::det(tempS))*exp(0.5/hatSigmaSq*(BrjtRes.t()*tempS*BrjtRes))));
      gsLS(k, j) = lS;
      u = R::runif(0, 1);
      if(u<lS){
        hatRStar.subvec((j*L), (j*L+L-1)).zeros();
      }else{
        hatRStar.subvec((j*L), (j*L+L-1)) = mvrnormCpp(meanRs, varRs);

      }
      res -= w.cols((j*L), (j*L+L-1)) * hatRStar.subvec((j*L), (j*L+L-1));
    }
    gsR0.row(t) = hatR0.t();

    gsRStar.row(t) = hatRStar.t();

    // invTAUsq.0|lambda, r0
    lInvTauSq0 = hatLambdaSq0;
    muInvTauSq0 = sqrt(hatLambdaSq0 * hatSigmaSq / square(hatR0));
    for(unsigned int j = 0; j < s; j++){
      if(hatR0(j) == 0){
        hatInvTauSq0(j) = 1/R::rgamma(1, 2/lInvTauSq0);
      }else{
        hatInvTauSq0(j) = rinvgaussian(muInvTauSq0(j), lInvTauSq0);
      }
    }
    gsInvTauSq0.row(t) = hatInvTauSq0.t();


    // invTAUsq.star|lambda.star, r.star
    lInvTauSqStar = L * hatLambdaSqStar;
    matRStar = arma::reshape(hatRStar, L, s);
    tRsRs = arma::sum(arma::square(matRStar), 0).t();
    muInvTauSqStar = arma::sqrt(L * hatLambdaSqStar * hatSigmaSq / tRsRs);
    for(unsigned int j = 0; j<s; j++){
      if(tRsRs(j) == 0){
        hatInvTauSqStar(j) = 1/R::rgamma((L+1)/2, 2/lInvTauSqStar);
      }else{
        hatInvTauSqStar(j) = rinvgaussian(muInvTauSqStar(j), lInvTauSqStar);
      }
    }
    gsInvTauSqStar.row(t) = hatInvTauSqStar.t();
    gsRstRs.row(t) = tRsRs.t();

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

    // sigma.sq|
    double shapeSig = alpha + n/2 + L*arma::accu(tRsRs != 0)/2 + arma::accu(hatR0 != 0)/2 ;
    repInvTauSq = arma::vectorise(arma::repelem(hatInvTauSqStar.t(), L, 1), 0);
    double rateSig = gamma + 0.5*(arma::accu(arma::square(res)) +
                                  arma::accu(square(hatR0) % hatInvTauSq0) +

                                  arma::accu(square(hatRStar) % repInvTauSq));
    hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
    gsSigmaSq(t) = hatSigmaSq;

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

    // pi.0|
    double shape1_0 = mu0 + arma::accu(hatR0 == 0);
    double shape2_0 = nu0 + arma::accu(hatR0 != 0);
    hatPi0 = R::rbeta(shape1_0, shape2_0);
    gsPi0(t) = hatPi0;


    // pi.star|
    double shape1_s = muStar + arma::accu(tRsRs == 0);
    double shape2_s = nuStar + arma::accu(tRsRs != 0);
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
  return Rcpp::List::create(Rcpp::Named("GS.alpha") = gsM,

                            Rcpp::Named("GS.beta") = gsR0,
                            Rcpp::Named("GS.ata") = gsAta,

                            Rcpp::Named("GS.eta") = gsRStar,
                            Rcpp::Named("GS.tRsRs") = gsRstRs,
                            Rcpp::Named("GS.invTAUsq.0") = gsInvTauSq0,

                            Rcpp::Named("GS.invTAUsq.star") = gsInvTauSqStar,
                            Rcpp::Named("GS.pi0") = gsPi0,

                            Rcpp::Named("GS.pi.star") = gsPiStar,
                            Rcpp::Named("GS.lambda.sq.0") = gsLambda0,

                            Rcpp::Named("GS.lambda.sq.star") = gsLambdaStar,
                            Rcpp::Named("GS.sigma.sq") = gsSigmaSq,
                            Rcpp::Named("GS.l0") = gsL0,

                            Rcpp::Named("GS.lS") = gsLS
                             );

}
