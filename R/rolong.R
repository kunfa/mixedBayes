rolong <- function(y,e,g,w,k, iterations=10000, burn.in=NULL, slope=TRUE, robust=TRUE, quant=0.5, sparse=TRUE, structure=c("group","individual"))
{
  structure = match.arg(structure)
  this.call = match.call()

  if(iterations<1) stop("iterations must be a positive integer.")
  if(is.null(burn.in)){
    BI = floor(iterations/2)
  }else if(burn.in>=1){
    BI = as.integer(burn.in)
  }else{
    stop("burn.in must be a positive integer.")
  }
  if(iterations<=BI) stop("iterations must be larger than burn.in.")
  if(slope){
    z = cbind(rep(1,k),c(1:k))
  if(robust){
    out = LONRBGLSS(y,e,g,w,z,k,quant,iterations,sparse, structure)
  }else{
    out = LONBGLSS(y,e,g,w,z,k,iterations,sparse, structure)
  }
  }else{
    z = rep(1,k)
    if(robust){
      out = LONRBGLSS_1(y,e,g,w,z,k,quant,iterations,sparse, structure)
    }else{
      out = LONBGLSS_1(y,e,g,w,z,k,iterations,sparse, structure)
    }
  }

  coeff.main = apply(out$GS.beta[-(1:BI),,drop=FALSE], 2, stats::median);
  coeff.interaction = apply(out$GS.eta[-(1:BI),,drop=FALSE], 2, stats::median);


  coefficient = list(main=coeff.main, inter=coeff.interaction)

  fit = list(call = this.call, posterior = out, coefficient=coefficient, burn.in = BI, iterations=iterations)

  class(fit)=c("rolong", class(out))
  fit
}
