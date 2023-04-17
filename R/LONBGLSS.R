LONBGLSS <- function(y,e,C,g,w,z,k,max.steps,sparse, structure){
  n = nrow(g)
  m = ncol(g)
  p = ncol(w)
  q = ncol(e)
  c = ncol(z)
  o = ncol(C)
  hatAta=rep(1,n*c)
  hatBeta = rep(1,m)
  hatEta = rep(1,p)
  hatAlpha = rep(1,q+o)
  hatInvTauSq1=rep(1,m)
  hatInvTauSq21=rep(1,p)
  hatInvTauSq22=rep(1,m)
  hatPiEta=1/2
  hatPiBeta=1/2
  invSigAlpha0 = rep(10^-3,q+o)
  hatLambdaSqStar1=1
  hatLambdaSqStar2=1
  hatSigmaSq=1
  hatPhiSq=1
  aStar=1
  bStar=1
  alpha=1
  gamma=1
  alpha1=1
  gamma1=1
  mu0=1
  nu0=1
  
  if(sparse){
    fit=switch (structure,
                "group" = BGLSS(y,e,C,g,w,z,max.steps,n,k,hatBeta,hatEta,hatAlpha,hatAta,hatInvTauSq1,hatInvTauSq22,hatPiEta,hatPiBeta,invSigAlpha0,hatLambdaSqStar1
                                ,hatLambdaSqStar2,hatSigmaSq,hatPhiSq,aStar,bStar,alpha,gamma,alpha1,gamma1,mu0,nu0),
                "individual" = BLSS(y,e,C,g,w,z,max.steps,n,k,hatBeta,hatEta,hatAlpha,hatAta,hatInvTauSq1,hatInvTauSq21,hatPiEta,hatPiBeta,invSigAlpha0,hatLambdaSqStar1
                                    ,hatLambdaSqStar2,hatSigmaSq,hatPhiSq,aStar,bStar,alpha,gamma,alpha1,gamma1,mu0,nu0)
    )
  }else{
    fit=switch (structure,
                "group" = BGL(y,e,C,g,w,z,max.steps,n,k,hatBeta,hatEta,hatAlpha,hatAta,hatInvTauSq1,hatInvTauSq22,invSigAlpha0,hatLambdaSqStar1
                              ,hatLambdaSqStar2,hatSigmaSq,hatPhiSq,aStar,bStar,alpha,gamma,alpha1,gamma1),
                "individual" = BL(y,e,C,g,w,z,max.steps,n,k,hatBeta,hatEta,hatAlpha,hatAta,hatInvTauSq1,hatInvTauSq21,invSigAlpha0,hatLambdaSqStar1
                                  ,hatLambdaSqStar2,hatSigmaSq,hatPhiSq,aStar,bStar,alpha,gamma,alpha1,gamma1)
    )
  }
  out = list( GS.beta = fit$GS.beta,
              GS.eta = fit$GS.eta)
  out
}