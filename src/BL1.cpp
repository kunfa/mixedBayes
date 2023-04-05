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

Rcpp::List BL_1(arma::mat y, arma:: mat e, arma:: mat C, arma::mat g, arma:: mat w, arma:: vec z,int maxSteps, int n, int k,arma::vec hatBeta, arma:: vec hatEta, arma::vec hatAlpha, arma::vec hatAta,arma::vec hatInvTauSq1, arma::vec hatInvTauSq2, arma::vec invSigAlpha0, double hatLambdaSqStar1, double hatLambdaSqStar2, double hatSigmaSq, double hatPhiSq,double aStar, double bStar, double alpha, double gamma, double alpha1,double gamma1)
{
  unsigned int q = e.n_cols-1,m = g.n_cols,p = w.n_cols;
  arma::mat gsAlpha(maxSteps, q+3),
  gsBeta(maxSteps,m),
  gseta(maxSteps,p),
  gsAta(maxSteps,n),
  gsInvTauSq1(maxSteps,m),
  gsInvTauSq2(maxSteps, p);


  arma::vec
    gsLambdaStar1(maxSteps),
    gsLambdaStar2(maxSteps),
    gsPhiSq(maxSteps),
    gsSigmaSq(maxSteps);


  arma::mat mat0,mat1,mat2,ei,gi,wi;


  double meanAlpha;
  double varAlpha;

  double meanb;
  double varb;



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

        A0 = A0+arma::as_scalar(ei.col(j).t()*ei.col(j))/hatSigmaSq;
        res1 = y.row(i).t()-ei*hatAlpha-gi*hatBeta-wi*hatEta-z*hatAta(i);
        res11 = res1+ei.col(j)*hatAlpha(j);
        B0 = B0+arma::as_scalar(ei.col(j).t()*res11)/hatSigmaSq;

      }


      meanAlpha = arma::as_scalar(B0/(invSigAlpha0(j)+A0));


      varAlpha = arma::as_scalar(1/(invSigAlpha0(j)+A0));

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
      double tzz;
      tzz = arma::accu(arma::square(z));
      double A;
      A = tzz/hatSigmaSq;
      res = y.row(i).t()-ei*hatAlpha-gi*hatBeta-wi*hatEta;
      double B;
      B = arma::as_scalar(z.t()*res)/hatSigmaSq;
      double varAta;
      varAta = 1/(A+1/hatPhiSq);
      double meanAta;
      meanAta = arma::as_scalar(B/(A+1/hatPhiSq));
      hatAta(i) = R::rnorm(meanAta, std::sqrt(varAta));
    }


    gsAta.row(t) = hatAta.t();


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
        arma::vec tggDiag = sum(square(gi),0).t();
        A1 = A1+tggDiag(j);
        res2 = y.row(i).t()-ei*hatAlpha-gi*hatBeta-wi*hatEta-z*hatAta(i);
        res22 = res2+gi.col(j)*hatBeta(j);
        B1 = B1+ arma::as_scalar(gi.col(j).t()*res22);

      }
      meanb = arma::as_scalar(B1/(hatInvTauSq1(j)+A1));

      varb = arma::as_scalar(hatSigmaSq/(hatInvTauSq1(j)+A1));


      hatBeta(j)=R::rnorm(meanb,sqrt(varb));

    }

    gsBeta.row(t) = hatBeta.t();


    // eta|

    double meane,varcove;

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
        arma::vec twwDiag = sum(square(wi),0).t();
        A2 = A2+twwDiag(j);
        res3 = y.row(i).t()-ei*hatAlpha-gi*hatBeta-wi*hatEta-z*hatAta(i);
        res33 = res3+wi.col(j)*hatEta(j);
        B2 = B2+ arma::as_scalar(wi.col(j).t()*res33);

      }


      meane = arma::as_scalar(B2/(hatInvTauSq2(j)+A2));
      varcove = arma::as_scalar(hatSigmaSq/(hatInvTauSq2(j)+A2));

      hatEta(j)=R::rnorm(meane,sqrt(varcove));


    }

    gseta.row(t) = hatEta.t();


    // sigma.sq|
    double shapeSig, rateSig;
    shapeSig = alpha + n*k/2 + m/2+ p/2;
    double ress;
    ress=0;
    for(int i=0;i<n;i++){
      ei = mat0.rows((i*k),(i*k+k-1));
      ei.insert_cols(2, C);
      gi = mat1.rows((i*k),(i*k+k-1));
      wi = mat2.rows((i*k),(i*k+k-1));
      arma::vec res4;
      res4 = y.row(i).t()-ei*hatAlpha-gi*hatBeta-wi*hatEta-z*hatAta(i);
      ress = ress+arma::accu(arma::square(res4));
    }
    rateSig = gamma + 0.5*(ress +arma::accu(square(hatBeta)%hatInvTauSq1)+ arma::accu(square(hatEta) % hatInvTauSq2));

    hatSigmaSq = 1/R::rgamma(shapeSig, 1/rateSig);
    gsSigmaSq(t) = hatSigmaSq;

    // phi.sq|
    double shapePhi, ratePhi;
    shapePhi = alpha1 + n/2;

    ratePhi = gamma1 + 0.5*(arma::accu(square(hatAta)));

    hatPhiSq = 1/R::rgamma(shapePhi, 1/ratePhi);
    gsPhiSq(t) = hatPhiSq;

    // invTAUsq.star1|

    double lInvTauSq1;
    lInvTauSq1   = hatLambdaSqStar1;
    arma::vec muInvTauSq1;
    muInvTauSq1 = arma::sqrt(hatLambdaSqStar1 * hatSigmaSq / arma::square(hatBeta));
    for(unsigned int j = 0; j<m; j++){

      hatInvTauSq1(j) = rinvgaussian(muInvTauSq1(j), lInvTauSq1);
    }



    gsInvTauSq1.row(t) = hatInvTauSq1.t();

    // invTAUsq.star2|

    double lInvTauSq2;
    lInvTauSq2 = hatLambdaSqStar2;
    arma::vec muInvTauSq2;
    muInvTauSq2= arma::sqrt(hatLambdaSqStar2 * hatSigmaSq / arma::square(hatEta));
    for(unsigned int j = 0; j<p; j++){

      hatInvTauSq2(j) = rinvgaussian(muInvTauSq2(j), lInvTauSq2);
    }



    gsInvTauSq2.row(t) = hatInvTauSq2.t();

    // lambda.star1|
    double shapeS1;
    shapeS1= aStar + m;
    double rateS1;
    rateS1= bStar + arma::accu(1/hatInvTauSq1)/2;
    hatLambdaSqStar1 = R::rgamma(shapeS1, 1/rateS1);
    gsLambdaStar1(t) = hatLambdaSqStar1;

    // lambda.star2|
    double shapeS2;
    shapeS2= aStar + p;
    double rateS2;
    rateS2= bStar + arma::accu(1/hatInvTauSq2)/2;
    hatLambdaSqStar2 = R::rgamma(shapeS2, 1/rateS2);
    gsLambdaStar2(t) = hatLambdaSqStar2;


  }

  return Rcpp::List::create(
    Rcpp::Named("GS.alpha") = gsAlpha,
    Rcpp::Named("GS.beta") = gsBeta,
    Rcpp::Named("GS.eta") = gseta,
    Rcpp::Named("GS.ata") = gsAta,
    Rcpp::Named("GS.invTAUsq1") = gsInvTauSq1,
    Rcpp::Named("GS.invTAUsq2") = gsInvTauSq2,
    Rcpp::Named("GS.lambda.sq1") = gsLambdaStar1,
    Rcpp::Named("GS.lambda.sq2") = gsLambdaStar2,
    Rcpp::Named("GS.phi.sq") = gsPhiSq,
    Rcpp::Named("GS.sigma.sq") = gsSigmaSq

  );

}
