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

Rcpp::List RBLSS_1(arma::mat y, arma:: mat e, arma::mat C, arma::mat g, arma:: mat w, arma:: vec z,int maxSteps, int n, int k,arma::vec hatBeta, arma:: vec hatEta, arma::vec hatAlpha, double hatTau, arma::vec hatV, arma::vec hatSg1,arma::vec hatSg2,arma::vec hatAta, arma::vec invSigAlpha0,double hatPi1,double hatPi2, double hatEtaSq1, double hatEtaSq2,double xi1, double xi2, double r1,double r2,double hatPhiSq,double a, double b, double alpha1,double gamma1, double sh1, double sh0, int progress)
{
  unsigned int q = e.n_cols-1,m = g.n_cols,p = w.n_cols;
  arma::mat gsAlpha(maxSteps, q+3),
  gsBeta(maxSteps,m),
  gseta(maxSteps,p),
  gsAta(maxSteps,n),
  gsV(maxSteps, n*k),
  gsSg1(maxSteps, m),
  gsSg2(maxSteps, p)
    ;

  arma::vec gsEtaSq1(maxSteps),
  gsEtaSq2(maxSteps),
  gsTau(maxSteps),
  gsPhiSq(maxSteps),
  gsPi1(maxSteps),
  gsPi2(maxSteps)
    ;


  arma::mat mat0,mat1,mat2,mat3,ei,gi,wi;


  double meanAlpha;
  double varAlpha;

  double meanb;
  double varb;

  double meane,varcove;

  arma::vec muV, muS1, muS2;
  double lambV, xi1Sq = std::pow(xi1, 2), xi2Sq = std::pow(xi2, 2),lj_1, lj_2,u1,u2;
  double RZoV;
  double tZZoV;
  for (int t = 0; t < maxSteps; t++) {


    mat0 = arma::repelem(e,k,1);
    mat1 = arma::repelem(g,k,1);
    mat2 = arma::repelem(w,k,1);


    // alpha|
    for(unsigned int j=0;j<q+3;j++){
      arma::vec res1, res11;
      double A0;
      A0 =0;
      double B0;
      B0=0;
      for(int i=0;i<n;i++){
        ei = mat0.rows((i*k),(i*k+k-1));
        ei.insert_cols(2, C);
        gi = mat1.rows((i*k),(i*k+k-1));
        wi = mat2.rows((i*k),(i*k+k-1));
        double tWWoV;
        tWWoV = arma::as_scalar((ei.col(j)/ hatV.subvec((i*k),(i*k+k-1))).t() * ei.col(j));
        A0 = A0+tWWoV;
        res1 = y.row(i).t()-ei*hatAlpha-gi*hatBeta-wi*hatEta-z*hatAta(i)-xi1*hatV.subvec((i*k),(i*k+k-1));
        res11 = res1+ei.col(j)*hatAlpha(j);
        double RWoV;
        RWoV = arma::sum(ei.col(j) % (res11/ hatV.subvec((i*k),(i*k+k-1))));
        B0 = B0+ RWoV;

      }

      varAlpha = 1/(A0*hatTau/xi2Sq + invSigAlpha0(j));
      meanAlpha = varAlpha*B0*hatTau / xi2Sq;

      hatAlpha(j) = R::rnorm(meanAlpha,sqrt(varAlpha));

    }

    gsAlpha.row(t) = hatAlpha.t();

    // ata|
    arma::vec res;
    for(int i=0;i<n;i++){
      ei = mat0.rows((i*k),(i*k+k-1));
      ei.insert_cols(2, C);
      gi = mat1.rows((i*k),(i*k+k-1));
      wi = mat2.rows((i*k),(i*k+k-1));

      tZZoV = arma::as_scalar((z/hatV.subvec((i*k),(i*k+k-1))).t() * z);
      res = y.row(i).t()-ei*hatAlpha-gi*hatBeta-wi*hatEta-xi1*hatV.subvec((i*k),(i*k+k-1));
      RZoV = arma::sum(z% (res/hatV.subvec((i*k),(i*k+k-1))));
      double varAta;
      varAta= 1/(tZZoV*hatTau/xi2Sq+1/hatPhiSq);
      double meanAta;
      meanAta= varAta* RZoV * hatTau / xi2Sq;
      hatAta(i) = R::rnorm(meanAta, sqrt(varAta));
    }


    gsAta.row(t) = hatAta.t();

    //v|
    arma::vec resv;
    for(int i=0;i<n;i++){
      ei = mat0.rows((i*k),(i*k+k-1));
      ei.insert_cols(2, C);
      gi = mat1.rows((i*k),(i*k+k-1));
      wi = mat2.rows((i*k),(i*k+k-1));
      resv = y.row(i).t()-ei*hatAlpha-gi*hatBeta-wi*hatEta-z*hatAta(i);
      lambV = hatTau*xi1Sq/xi2Sq + 2*hatTau;
      muV = arma::sqrt((xi1Sq+2*xi2Sq) / arma::square(resv));
      arma::vec v(k);
      for(int k0=0;k0<k;k0++){
        bool flag = true;
        while(flag){
          v(k0) = 1/rinvGauss(muV(k0), lambV);
          if(v(k0)<=0 || std::isinf(v(k0)) || std::isnan(v(k0))){
            if(progress != 0) Rcpp::Rcout << "v(k0) <= 0 or nan or inf" << std::endl;
            Rcpp::checkUserInterrupt();
          }else{
            flag = false;
          }
        }
      }
      hatV.subvec((i*k),(i*k+k-1)) = v;
    }


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
    muS2 = std::sqrt(hatEtaSq2)/ arma::abs(hatEta);
    for(unsigned int j = 0; j<p; j++){
      if(hatEta(j) == 0){
        hatSg2(j) = R::rexp(2/hatEtaSq2);
      }else{
        bool flag = true;
        while(flag){
          hatSg2(j) = 1/rinvGauss(muS2(j), hatEtaSq2);
          if(hatSg2(j)<=0 || std::isinf(hatSg2(j)) || std::isnan(hatSg2(j))){
            if(progress != 0){
              Rcpp::Rcout << "hatSg2(j)： " << hatSg2(j) << std::endl;
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
      arma::vec res2, res22;
      double A1;
      A1=0;
      double B1;
      B1=0;

      for(int i=0;i<n;i++){
        ei = mat0.rows((i*k),(i*k+k-1));
        ei.insert_cols(2, C);
        gi = mat1.rows((i*k),(i*k+k-1));
        wi = mat2.rows((i*k),(i*k+k-1));
        double XgXgoV1;
        XgXgoV1 = arma::as_scalar((gi.col(j)/ hatV.subvec((i*k),(i*k+k-1))).t() * gi.col(j));
        A1 = A1+XgXgoV1;
        res2 = y.row(i).t()-ei*hatAlpha-gi*hatBeta-wi*hatEta-z*hatAta(i)-xi1*hatV.subvec((i*k),(i*k+k-1));
        res22 = res2+gi.col(j)*hatBeta(j);
        double RXgoV1;
        RXgoV1 = arma::sum(gi.col(j) % (res22/ hatV.subvec((i*k),(i*k+k-1))));
        B1 = B1+ RXgoV1;

      }

      varb = 1/(A1*hatTau/xi2Sq + 1/hatSg1(j));
      meanb = varb*B1*hatTau / xi2Sq;
      double lj_temp_1 = std::sqrt(hatSg1(j))*std::exp(-0.5*varb*pow(B1*hatTau/xi2Sq,2))/std::sqrt(varb);
      lj_1 = hatPi1/(hatPi1+(1-hatPi1)*lj_temp_1);
      u1 = R::runif(0, 1);
      if(u1<lj_1){
        hatBeta(j) = R::rnorm(meanb, sqrt(varb));
      }else{
        hatBeta(j) = 0;
      }

    }

    gsBeta.row(t) = hatBeta.t();


    // eta|

    for(unsigned int j=0;j<p;j++){
      arma::vec res3, res33;
      double A2;
      A2=0;
      double B2;
      B2=0;

      for(int i=0;i<n;i++){
        ei = mat0.rows((i*k),(i*k+k-1));
        ei.insert_cols(2, C);
        gi = mat1.rows((i*k),(i*k+k-1));
        wi = mat2.rows((i*k),(i*k+k-1));
        double XgXgoV2;
        XgXgoV2 = arma::as_scalar((wi.col(j)/ hatV.subvec((i*k),(i*k+k-1))).t() * wi.col(j));
        A2 = A2+XgXgoV2;
        res3 = y.row(i).t()-ei*hatAlpha-gi*hatBeta-wi*hatEta-z*hatAta(i)-xi1*hatV.subvec((i*k),(i*k+k-1));
        res33 = res3+wi.col(j)*hatEta(j);
        double RXgoV2;
        RXgoV2 = arma::sum(wi.col(j) % (res33/ hatV.subvec((i*k),(i*k+k-1))));
        B2 = B2+ RXgoV2;

      }

      varcove = 1/(A2*hatTau/xi2Sq + 1/hatSg2(j));
      meane = varcove*B2*hatTau / xi2Sq;
      double lj_temp_2 = std::sqrt(hatSg2(j))*std::exp(-0.5*varcove*pow(B2*hatTau/xi2Sq,2))/std::sqrt(varcove);
      lj_2 = hatPi2/(hatPi2+(1-hatPi2)*lj_temp_2);
      u2 = R::runif(0, 1);
      if(u2<lj_2){
        hatEta(j) = R::rnorm(meane, sqrt(varcove));
      }else{
        hatEta(j) = 0;
      }

    }
    gseta.row(t) = hatEta.t();



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


    // phi.sq|
    double shapePhi, ratePhi;
    shapePhi = alpha1 + n/2;

    ratePhi = gamma1 + 0.5*(arma::accu(square(hatAta)));

    hatPhiSq = 1/R::rgamma(shapePhi, 1/ratePhi);
    gsPhiSq(t) = hatPhiSq;

    //tau|

    double shape = a + 3*n*k/2;
    double rest;
    arma::vec restt;
    rest = 0;
    double f;
    f=0;
    for(int i=0;i<n;i++){
      ei = mat0.rows((i*k),(i*k+k-1));
      ei.insert_cols(2, C);
      gi = mat1.rows((i*k),(i*k+k-1));
      wi = mat2.rows((i*k),(i*k+k-1));
      restt = y.row(i).t()-ei*hatAlpha-gi*hatBeta-wi*hatEta-z*hatAta(i)-xi1*hatV.subvec((i*k),(i*k+k-1));
      double ResSqoV;
      ResSqoV = arma::accu(arma::square(restt)/hatV.subvec((i*k),(i*k+k-1)));
      rest = rest+ResSqoV/(2*xi2Sq);
      f = f+ arma::accu(hatV.subvec((i*k),(i*k+k-1)));
    }
    double rate = b + f + rest/(2*xi2Sq);
    hatTau = R::rgamma(shape, 1/rate);
    gsTau(t) = hatTau;

    //pi1|
    double shapep1 = sh1 + arma::accu(hatBeta != 0);
    double shapep2 = sh0 + arma::accu(hatBeta == 0);
    hatPi1 = R::rbeta(shapep1, shapep2);
    gsPi1(t) = hatPi1;

    //pi2|
    double shape12 = sh1 + arma::accu(hatEta != 0);
    double shape22 = sh0 + arma::accu(hatEta == 0);
    hatPi2 = R::rbeta(shape12, shape22);
    gsPi2(t) = hatPi2;

  }

  return Rcpp::List::create(
    Rcpp::Named("GS.alpha") = gsAlpha,
    Rcpp::Named("GS.beta") = gsBeta,
    Rcpp::Named("GS.eta") = gseta,
    Rcpp::Named("GS.ata") = gsAta,
    Rcpp::Named("GS.v") = gsV,
    Rcpp::Named("GS.s1") = gsSg1,
    Rcpp::Named("GS.s2") = gsSg2,
    Rcpp::Named("GS.eta21.sq") = gsEtaSq1,
    Rcpp::Named("GS.eta22.sq") = gsEtaSq2,
    Rcpp::Named("GS.phi.sq") = gsPhiSq,
    Rcpp::Named("GS.tau") = gsTau,
    Rcpp::Named("GS.pi1") = gsPi1,
    Rcpp::Named("GS.pi2") = gsPi2

  );

}
