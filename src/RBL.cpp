#include<RcppArmadillo.h>
#include<stdio.h>
#include<vector>
#include"BVCUtilities.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;
//using namespace R;

// [[Rcpp::export()]]

Rcpp::List RBL(arma::vec y, arma::mat e, arma::mat g, arma:: mat w,int maxSteps, int q, int k,arma::vec hatBeta, arma:: vec hatEta, arma::vec hatAlpha, arma::mat hatAta, arma::mat z,double hatTau, arma::vec hatV, arma::vec hatSg1,arma::vec hatSg2,arma::mat invSigAlpha0, double hatEtaSq1, double hatEtaSq2, double xi1, double xi2, double r1,double r2,double hatPhiSq,double a, double b,double alpha1, double gamma1,int progress)
{
  unsigned int n = g.n_rows, m = g.n_cols,p = w.n_cols,c = z.n_cols, n1 = n/k;
  arma::mat gsAlpha(maxSteps, q),
  gsBeta(maxSteps,m),
  gsAta(maxSteps,n1*c),
  gsEta(maxSteps,p),
  gsV(maxSteps, n),
  gsSg1(maxSteps, m),
  gsSg2(maxSteps, p)
    ;

  arma::vec gsEtaSq1(maxSteps),
  gsEtaSq2(maxSteps),
  gsPhiSq(maxSteps),
  gsTau(maxSteps)
    ;
  arma::mat temp1;

  double meanb,meane;
  double varb,vare;
  double XgXgoV1,RXgoV1,XgXgoV2,RXgoV2;

  arma::vec muV, muS1,  muS2, REoV(q), meanAlpha, res,meanAta,RZoV(c);
  double lambV, xi1Sq = std::pow(xi1, 2), xi2Sq = std::pow(xi2, 2);
  arma::mat varAlpha, tEEoV(q,q),varAta, tZZoV(c,c);

  for (int t = 0; t < maxSteps; t++) {


    // alpha|
    res = y - g*hatBeta-w*hatEta-xi1*hatV;
    for(unsigned int i=0;i<n1;i++){
      res.subvec((i*k), (i*k+k-1)) -= z*hatAta.col(i);
    }
    tEEoV = (e.each_col()/hatV).t() * e;
    REoV = arma::sum(e.each_col()% (res/hatV), 0).t();
    varAlpha = arma::inv_sympd(tEEoV*hatTau/xi2Sq+invSigAlpha0);
    meanAlpha = varAlpha* REoV * hatTau / xi2Sq;
    hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
    res -= e * hatAlpha;
    gsAlpha.row(t) = hatAlpha.t();

    // ata|

    for(unsigned int i=0;i<n1;i++){
      res.subvec((i*k), (i*k+k-1)) += z * hatAta.col(i);
      tZZoV = (z.each_col()/hatV.subvec((i*k), (i*k+k-1))).t() * z;
      RZoV = arma::sum(z.each_col()% (res.subvec((i*k), (i*k+k-1))/hatV.subvec((i*k), (i*k+k-1))), 0).t();

      temp1 = tZZoV*hatTau/xi2Sq;
      temp1.diag()+=1/hatPhiSq;
      varAta = arma::inv(temp1);
      meanAta = varAta* RZoV * hatTau / xi2Sq;
      hatAta.col(i) = mvrnormCpp(meanAta, varAta);
      res.subvec((i*k), (i*k+k-1)) -= z * hatAta.col(i);
    }

    gsAta.row(t) = arma::vectorise(hatAta).t();

    //v|
    res += xi1*hatV;
    lambV = hatTau*xi1Sq/xi2Sq + 2*hatTau;
    muV = arma::sqrt((xi1Sq+2*xi2Sq) / arma::square(res));
    for(unsigned int i=0;i<n;i++){
      bool flag = true;
      while(flag){
        hatV(i) = 1/rinvGauss(muV(i), lambV);
        if(hatV(i)<=0 || std::isinf(hatV(i)) || std::isnan(hatV(i))){
          if(progress != 0){
            Rcpp::Rcout << "hatV(i) <= 0 or nan or inf" << std::endl;
            Rcpp::checkUserInterrupt();
          }
        }else{
          flag = false;
        }
      }
    }
    res -= xi1*hatV;
    gsV.row(t) = hatV.t();

    //s1|
    muS1 = std::sqrt(hatEtaSq1)/ arma::abs(hatBeta);
    for(unsigned int j = 0; j<m; j++){
      bool flag = true;
      while(flag){
        hatSg1(j) = 1/rinvGauss(muS1(j), hatEtaSq1);
        if(hatSg1(j)<=0 || std::isinf(hatSg1(j)) || std::isnan(hatSg1(j))){
          if(progress != 0) Rcpp::Rcout << "hatSg1(j): " << hatSg1(j) << std::endl;
          Rcpp::checkUserInterrupt();
        }else{
          flag = false;
        }
      }
    }
    gsSg1.row(t) = hatSg1.t();


    //s2|
    muS2 = std::sqrt(hatEtaSq2)/ arma::abs(hatEta);
    for(unsigned int j = 0; j<p; j++){
      bool flag = true;
      while(flag){
        hatSg2(j) = 1/rinvGauss(muS2(j), hatEtaSq2);
        if(hatSg2(j)<=0 || std::isinf(hatSg2(j)) || std::isnan(hatSg2(j))){
          if(progress != 0) Rcpp::Rcout << "hatSg2(j): " << hatSg2(j) << std::endl;
          Rcpp::checkUserInterrupt();
        }else{
          flag = false;
        }
      }
    }
    gsSg2.row(t) = hatSg2.t();


    // Beta|

    for(unsigned int j=0; j<m; j++){
      res += g.col(j) * hatBeta(j);
      XgXgoV1 = arma::as_scalar((g.col(j)/hatV).t() * g.col(j));
      varb = 1/(XgXgoV1*hatTau/xi2Sq + 1/hatSg1(j));

      RXgoV1 = arma::sum(g.col(j) % (res/hatV));
      meanb = varb * RXgoV1 * hatTau / xi2Sq;
      hatBeta(j) = R::rnorm(meanb, sqrt(varb));
      res -= g.col(j) * hatBeta(j);
    }

    gsBeta.row(t) = hatBeta.t();


    // eta|

    for(unsigned int j=0; j<p; j++){
      res += w.col(j) * hatEta(j);
      XgXgoV2 = arma::as_scalar((w.col(j)/hatV).t() * w.col(j));
      vare = 1/(XgXgoV2*hatTau/xi2Sq + 1/hatSg2(j));

      RXgoV2 = arma::sum(w.col(j) % (res/hatV));
      meane = vare * RXgoV2 * hatTau / xi2Sq;
      hatEta(j) = R::rnorm(meane, sqrt(vare));
      res -= w.col(j) * hatEta(j);
    }

    gsEta.row(t) = hatEta.t();

    //etasq1;

    double shape2 = m+1;
    double rate2 = arma::accu(hatSg1)/2 + r1;
    hatEtaSq1 = R::rgamma(shape2, 1/rate2);
    gsEtaSq1(t) = hatEtaSq1;

    //etasq2;

    double shape21 = p+1;
    double rate21 = arma::accu(hatSg2)/2 + r2;
    hatEtaSq2 = R::rgamma(shape21, 1/rate21);
    gsEtaSq2(t) = hatEtaSq2;


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

    //tau|

    double shape = a + 3*n/2;
    double ResSqoV;
    ResSqoV = arma::accu(arma::square(res)/hatV);
    double rate = b + arma::accu(hatV)+ResSqoV/(2*xi2Sq);
    hatTau = R::rgamma(shape, 1/rate);
    gsTau(t) = hatTau;


  }

  return Rcpp::List::create(
    Rcpp::Named("GS.alpha") = gsAlpha,
    Rcpp::Named("GS.beta") = gsBeta,
    Rcpp::Named("GS.eta") = gsEta,
    Rcpp::Named("GS.v") = gsV,
    Rcpp::Named("GS.s1") = gsSg1,
    Rcpp::Named("GS.s2") = gsSg2,
    Rcpp::Named("GS.eta21.sq") = gsEtaSq1,
    Rcpp::Named("GS.eta22.sq") = gsEtaSq2,
    Rcpp::Named("GS.tau") = gsTau

  );

}
