LONRBGLSS <- function(y,e,X,g,w,z,k,quant,max.steps,sparse, structure){

  n = nrow(g)
  m = ncol(g)
  p = ncol(w)
  E = cbind(e,X)
  o = ncol(X)
  q = ncol(E)
  c = ncol(z)
  n1 = n/k
  hatTau=1
  hatV = rep(1,n)
  hatSg1 = rep(1,m)
  hatSg21= rep(1,p)
  hatSg22 = rep(1,m)
  xi1=(1-2*quant)/(quant*(1-quant))
  xi2 = sqrt(2/(quant*(1-quant)))
  hatEtaSq1=1
  hatEtaSq2=1
  r1=1
  r2=1
  a=1
  b=1
  hatAta=matrix(c(rep(1,n1*c)),nrow=c)
  hatBeta = rep(1,m)
  hatEta1 = rep(1,p)
  hatEta2 = matrix(c(rep(1,p)),nrow=q-o)
  hatAlpha = rep(1,q)
  invSigAlpha0 = diag(rep(10^-3,q))
  alpha1=1
  gamma1=1
  hatPhiSq=1
  progress=0
  hatPi1=1/2
  hatPi2=1/2
  sh1=1
  sh0=1
  debugging=FALSE

  progress = ifelse(debugging, 10^(floor(log10(max.steps))-1), 0)
  if(sparse){
    fit=switch (structure,
                "bi-level" = RBGLSS(y,E,g,w,max.steps,q,o,k,hatBeta,hatEta2,hatAlpha,hatAta,z,hatTau,hatV,hatSg1,hatSg22,invSigAlpha0,
                                 hatPi1,hatPi2,hatEtaSq1,hatEtaSq2,xi1,xi2,r1,r2,hatPhiSq,a,b,alpha1,gamma1,sh1,sh0,progress),
                "individual" = RBLSS(y,E,g,w,max.steps,q,k,hatBeta,hatEta1,hatAlpha,hatAta,z,hatTau,hatV,hatSg1,hatSg21,invSigAlpha0,
                                     hatPi1,hatPi2,hatEtaSq1,hatEtaSq2,xi1,xi2,r1,r2,hatPhiSq,a,b,alpha1,gamma1,sh1,sh0,progress)
    )
  }else{
    fit=switch (structure,
                "bi-level" =  RBGL(y,E,g,w,max.steps,q,o,k,hatBeta,hatEta2,hatAlpha,hatAta,z,hatTau,hatV,hatSg1,hatSg22,invSigAlpha0,
                                hatEtaSq1,hatEtaSq2,xi1,xi2,r1,r2,hatPhiSq,a,b,alpha1,gamma1,progress),
                "individual" = RBL(y,E,g,w,max.steps,q,k,hatBeta,hatEta1,hatAlpha,hatAta,z,hatTau,hatV,hatSg1,hatSg21,invSigAlpha0,
                                   hatEtaSq1,hatEtaSq2,xi1,xi2,r1,r2,hatPhiSq,a,b,alpha1,gamma1,progress)
    )
  }

  out = list(GS.gamma1 = fit$GS.alpha[,1:(q-o)], GS.gamma0 = fit$GS.alpha[,-(1:(q-o))],GS.gamma2 = fit$GS.beta,
             GS.gamma3 = fit$GS.eta,GS.alpha = fit$GS.ata)
  out

}
