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

Rcpp::List RBGLSS_1(arma::vec y, arma::mat e, arma::mat g, arma:: mat w, int maxSteps, unsigned int q, unsigned int o,unsigned int k, arma::vec hatBeta, arma:: mat hatEta, arma::vec hatAlpha, arma:: mat hatAta, arma:: mat z, double hatTau, arma::vec hatV, arma::vec hatSg1,arma::vec hatSg2,arma::mat invSigAlpha0,double hatPi1,double hatPi2, double hatEtaSq1, double hatEtaSq2,double xi1, double xi2, double r1,double r2,double Phi1Sq, double a, double b, double alpha1, double gamma1,double sh1, double sh0, int progress)
{
  unsigned int n = g.n_rows, L = q-o, m = g.n_cols, p = w.n_cols, n1 = n/k;
  arma::mat gsAlpha(maxSteps, q),
  gsBeta(maxSteps,m),
  gseta(maxSteps,p),
  gsV(maxSteps, n),
  gsAta(maxSteps,n1),
  gsSg1(maxSteps, m),
  gsSg2(maxSteps, m)
    ;

  arma::vec gsEtaSq1(maxSteps),
  gsEtaSq2(maxSteps),
  gsTau(maxSteps),
  gsPi1(maxSteps),
  gsPhi1Sq(maxSteps),
  gsPi2(maxSteps)
    ;
  arma::mat temp,temp1;

  double meanb;
  double varb;
  double XgXgoV1,RXgoV1;
  arma::vec meane;
  arma::mat varcove;
  arma::vec muV, muS1,tBgBg, REoV(q), meanAlpha, res;
  double muS2;
  double lambV, xi1Sq = std::pow(xi1, 2), xi2Sq = std::pow(xi2, 2),lj_1, lg_2,u1,u2;
  arma::mat XgXgoV2(L,L), varAlpha, tEEoV(q,q);
  arma::vec RXgoV2(L);

  for (int t = 0; t < maxSteps; t++) {

    // alpha|
    arma::mat Zblock(n, n1, arma::fill::zeros);
    for (unsigned int i = 0; i < n1; ++i) {
      Zblock.submat(i*k, i, i*k + k - 1, i).fill(1.0);
    }
    res = y - g*hatBeta-w*arma::vectorise(hatEta)-xi1*hatV - Zblock*arma::vectorise(hatAta);
    tEEoV = (e.each_col()/hatV).t() * e;
    REoV = arma::sum(e.each_col()% (res/hatV), 0).t();
    varAlpha = arma::inv_sympd(tEEoV*hatTau/xi2Sq+invSigAlpha0);
    meanAlpha = varAlpha* REoV * hatTau / xi2Sq;
    hatAlpha = mvrnormCpp(meanAlpha, varAlpha);
    res -= e * hatAlpha;
    gsAlpha.row(t) = hatAlpha.t();

    // ata|
    res+= Zblock*arma::vectorise(hatAta);

    for (unsigned int i = 0; i < n1; i++) {


      arma::vec rblock = res.subvec(i*k, i*k+k-1);
      arma::vec vblock = hatV.subvec(i*k, i*k+k-1);

      // ===== ata0 =====
      double tZZ0 = arma::as_scalar(
        (z.col(0) / vblock).t() * z.col(0)
      );
      double RZ0 = arma::sum(z.col(0) % (rblock / vblock));

      double var0  = 1.0 / (tZZ0 * hatTau / xi2Sq + 1.0 / Phi1Sq);
      double mean0 = var0 * (RZ0 * hatTau / xi2Sq);

      hatAta(0, i) = R::rnorm(mean0, std::sqrt(var0));

    }
    res-= Zblock*arma::vectorise(hatAta);

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
      if(hatBeta(j) == 0){
        hatSg1(j) = R::rexp(2/hatEtaSq1);
      }else{
        bool flag = true;
        while(flag){
          hatSg1(j) = 1/rinvGauss(muS1(j), hatEtaSq1);
          if(hatSg1(j)<=0 || std::isinf(hatSg1(j)) || std::isnan(hatSg1(j))){
            if(progress != 0){
              Rcpp::Rcout << "hatSg1(j)： " << hatSg1(j) << std::endl;
              Rcpp::checkUserInterrupt();
            }
          }else{
            flag = false;
          }
        }
      }

    }
    gsSg1.row(t) = hatSg1.t();


    //s2|
    tBgBg = arma::sum(arma::square(hatEta), 0).t();
    for(unsigned int j = 0; j<m; j++){
      if(tBgBg(j) == 0){
        hatSg2(j) = R::rgamma((L+1)/2, 2/hatEtaSq2);
      }else{
        muS2 = std::sqrt(hatEtaSq2/tBgBg(j));
        bool flag = true;
        while(flag){
          hatSg2(j) = 1/rinvGauss(muS2, hatEtaSq2);
          if(hatSg2(j)<=0 || std::isinf(hatSg2(j)) || std::isnan(hatSg2(j))){
            if(progress != 0){
              Rcpp::Rcout << "hatSg2(j): " << hatSg2(j) << " muS2: " << muS2 << " hatEtaSq2: " << hatEtaSq2 << std::endl;
              Rcpp::checkUserInterrupt();
            }
          }else{
            flag = false;
          }
        }
      }
    }
    gsSg2.row(t) = hatSg2.t();

    // Beta|

    for(unsigned int j=0;j<m;j++){
      res += g.col(j) * hatBeta(j);
      XgXgoV1 = arma::as_scalar((g.col(j)/ hatV).t() * g.col(j));
      varb = 1/(XgXgoV1 *hatTau/xi2Sq + 1/hatSg1(j));
      RXgoV1 = arma::sum(g.col(j) % res/ hatV)*hatTau / xi2Sq;
      meanb = varb* RXgoV1;
      double lj_temp_1 = std::sqrt(hatSg1(j))*std::exp(-0.5*varb*pow(RXgoV1,2))/std::sqrt(varb);
      lj_1 = hatPi1/(hatPi1+(1-hatPi1)*lj_temp_1);
      u1 = R::runif(0, 1);
      if(u1<lj_1){
        hatBeta(j) = R::rnorm(meanb, sqrt(varb));

      }else{
        hatBeta(j) = 0;
      }
      res -= g.col(j) * hatBeta(j);
    }

    gsBeta.row(t) = hatBeta.t();


    // eta|

    for(unsigned int j=0;j<m;j++){
      res += w.cols(j*L, j*L+L-1) * hatEta.col(j);
      XgXgoV2 = (w.cols((j*L),(j*L+L-1)).each_col()/hatV).t()*w.cols((j*L),(j*L+L-1));
      temp = XgXgoV2*hatTau/xi2Sq;
      temp.diag() += 1/hatSg2(j);
      varcove = arma::inv(temp);
      RXgoV2 = arma::sum(w.cols((j*L),(j*L+L-1)).each_col()%(res/hatV),0).t();
      RXgoV2 *= hatTau/xi2Sq;
      meane = varcove*RXgoV2;
      double lg_temp_2 = arma::as_scalar(arma::exp(-0.5*(RXgoV2.t()*varcove*RXgoV2)))*std::sqrt(arma::det(temp))*std::pow(hatSg2(j), L/2);
      lg_2 = hatPi2/(hatPi2+(1-hatPi2)*lg_temp_2);
      u2 = R::runif(0, 1);
      if(u2<lg_2){
        hatEta.col(j) = mvrnormCpp(meane, varcove);

      }else{
        hatEta.col(j).zeros();
      }
      res -= w.cols(j*L, j*L+L-1) * hatEta.col(j);
    }

    gseta.row(t) = arma::vectorise(hatEta).t();

    //etasq1;

    double shape2 = m+1;
    double rate2 = arma::accu(hatSg1)/2 + r1;
    hatEtaSq1 = R::rgamma(shape2, 1/rate2);
    gsEtaSq1(t) = hatEtaSq1;

    //etasq2;

    double shape21 = (m+m*L)/2+1;
    double rate21 = arma::accu(hatSg2)/2 + r2;
    hatEtaSq2 = R::rgamma(shape21, 1/rate21);
    gsEtaSq2(t) = hatEtaSq2;

    //phi1sq;
    double diff1 = 0.5 * arma::accu( arma::square(hatAta.row(0)) );

    double shape1 = alpha1 + n1 / 2.0;
    double rate1  = gamma1 + diff1;
    Phi1Sq = 1.0 / R::rgamma(shape1, 1.0 / rate1);

    gsPhi1Sq(t) = Phi1Sq;


    //tau|

    double shape = a + 3*n/2;
    double ResSqoV;
    ResSqoV = arma::accu(arma::square(res)/hatV);
    double rate = b + arma::accu(hatV)+ResSqoV/(2*xi2Sq);
    hatTau = R::rgamma(shape, 1/rate);
    gsTau(t) = hatTau;

    //pi1|
    double shapep1 = sh1 + arma::accu(hatBeta != 0);
    double shapep2 = sh0 + arma::accu(hatBeta == 0);
    hatPi1 = R::rbeta(shapep1, shapep2);
    gsPi1(t) = hatPi1;

    //pi2|
    double shape12 = sh1 + arma::accu(tBgBg != 0);
    double shape22 = sh0 + arma::accu(tBgBg == 0);
    hatPi2 = R::rbeta(shape12, shape22);
    gsPi2(t) = hatPi2;

  }

  return Rcpp::List::create(
    Rcpp::Named("GS.alpha") = gsAlpha,
    Rcpp::Named("GS.beta") = gsBeta,
    Rcpp::Named("GS.eta") = gseta,
    Rcpp::Named("GS.ata") = gsAta
  );

}
