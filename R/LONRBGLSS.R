LONRBGLSS <- function(y,e,C,g,w,z,k,quant,max.steps,sparse, structure){
  
  n = nrow(g)
  m = ncol(g)
  p = ncol(w)
  q = ncol(e)
  o = ncol(C)
  c = ncol(z)
  hatTau=1
  hatV = matrix(c(rep(1,n*k)),nrow=k)
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
  hatAta=matrix(c(rep(1,n*c)),nrow=c)
  hatBeta = rep(1,m)
  hatEta1 = rep(1,p)
  hatEta2 = matrix(c(rep(1,p)),nrow=q)
  hatAlpha = rep(1,q+o)
  invSigAlpha0 = diag(rep(10^-3,q+o))
  alpha1=1
  gamma1=1
  hatPhiSq=1
  progress=0
  hatPi1=1/2
  hatPi2=1/2
  sh1=1
  sh0=1
  
  if(sparse){
    fit=switch (structure,
                "group" = RBGLSS(y,e,C,g,w,z,max.steps,n,k,hatBeta,hatEta2,hatAlpha,hatTau,hatV,hatSg1,hatSg22,hatAta,invSigAlpha0,
                                 hatPi1,hatPi2,hatEtaSq1,hatEtaSq2,xi1,xi2,r1,r2,hatPhiSq,a,b,alpha1,gamma1,sh1,sh0,progress),
                "individual" = RBLSS(y,e,C,g,w,z,max.steps,n,k,hatBeta,hatEta1,hatAlpha,hatTau,hatV,hatSg1,hatSg21,hatAta,invSigAlpha0,
                                     hatPi1,hatPi2,hatEtaSq1,hatEtaSq2,xi1,xi2,r1,r2,hatPhiSq,a,b,alpha1,gamma1,sh1,sh0,progress)
    )
  }else{
    fit=switch (structure,
                "group" =  RBGL(y,e,C,g,w,z,max.steps,n,k,hatBeta,hatEta2,hatAlpha,hatTau,hatV,hatSg1,hatSg22,hatAta,invSigAlpha0,
                                hatEtaSq1,hatEtaSq2,xi1,xi2,r1,r2,hatPhiSq,a,b,alpha1,gamma1,progress),
                "individual" = RBL(y,e,C,g,w,z,max.steps,n,k,hatBeta,hatEta1,hatAlpha,hatTau,hatV,hatSg1,hatSg21,hatAta,invSigAlpha0,
                                   hatEtaSq1,hatEtaSq2,xi1,xi2,r1,r2,hatPhiSq,a,b,alpha1,gamma1,progress)
    )
  }
  
  out = list( GS.beta = fit$GS.beta,
              GS.eta = fit$GS.eta)
  out
  
}