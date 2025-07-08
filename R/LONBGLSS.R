LONBGLSS <- function(y,e,X,g,w,z,k,max.steps,sparse, structure){
  n = nrow(g)
  m = ncol(g)
  p = ncol(w)
  c = ncol(z)
  E = cbind(e,X)
  q = ncol(E)
  o = ncol(X)
  n1 = n/k
  hatAta=matrix(c(rep(1,n1*c)),nrow=c)
  hatBeta = rep(1,m)
  hatEta = rep(1,p)
  hatAlpha = rep(1,q)
  hatInvTauSq1=rep(1,m)
  hatInvTauSq21=rep(1,p)
  hatInvTauSq22=rep(1,m)
  hatPiEta=1/2
  hatPiBeta=1/2
  invSigAlpha0 = rep(10^-3,q)
  hatLambdaSqStar1=1
  hatLambdaSqStar2=1
  hatSigmaSq=1
  hatPhiSq=1
  a0=aStar=1
  b0=bStar=1
  alpha=1
  gamma=1
  alpha1=1
  gamma1=1
  mu0=mu1=1
  nu0=nu1=1
  debugging=FALSE

  progress = ifelse(debugging, 10^(floor(log10(max.steps))-1), 0)

  if(sparse){
    fit=switch (structure,
                "bi-level" = BGLSS(y,E,g,w,z,q,o,k,max.steps,hatAlpha,hatBeta,hatAta,hatEta,invSigAlpha0,hatInvTauSq1,hatInvTauSq22,hatPiBeta,hatPiEta,hatLambdaSqStar1
                                ,hatLambdaSqStar2,hatSigmaSq,hatPhiSq,a0,b0,aStar,bStar,alpha,gamma,alpha1,gamma1,mu0,mu1,nu0,nu1,progress),
                "individual" = BLSS(y,E,g,w,z,q,k,max.steps,hatAlpha,hatBeta,hatAta,hatEta,invSigAlpha0,hatInvTauSq1,hatInvTauSq21,hatPiBeta,hatPiEta,hatLambdaSqStar1
                                    ,hatLambdaSqStar2,hatSigmaSq,hatPhiSq,a0,b0,aStar,bStar,alpha,gamma,alpha1,gamma1,mu0,mu1,nu0,nu1,progress)
    )
  }else{
    fit=switch (structure,
                "bi-level" = BGL(y,E,g,w,q,o,k,max.steps,hatAlpha,hatBeta,hatEta,hatAta,z,invSigAlpha0,hatInvTauSq1,hatInvTauSq22,hatLambdaSqStar1
                              ,hatLambdaSqStar2,hatSigmaSq,a0,b0,aStar,bStar,hatPhiSq,alpha,gamma,alpha1,gamma1,progress),
                "individual" = BL(y,E,g,w,q,k,max.steps,hatAlpha,hatBeta,hatEta,hatAta,z,invSigAlpha0,hatInvTauSq1,hatInvTauSq21,hatLambdaSqStar1
                                  ,hatLambdaSqStar2,hatSigmaSq,a0,b0,aStar,bStar,hatPhiSq,alpha,gamma,alpha1,gamma1,progress)
    )
  }
  out = list(GS.gamma1 = fit$GS.alpha[,1:(q-o)], GS.gamma0 = fit$GS.alpha[,-(1:(q-o))],GS.gamma2 = fit$GS.beta,
              GS.gamma3 = fit$GS.eta,GS.alpha = fit$GS.ata)
  out
}
