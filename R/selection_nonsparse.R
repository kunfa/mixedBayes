selection_nonsparse = function(obj,burn.in=obj$burn.in){
  sg1 = obj$posterior$GS.beta
  sg2 = obj$posterior$GS.eta
  q_t1=c()
  for(j in 1:ncol(sg1)){
    t1=as.matrix(sg1[,j])
    t1 = t1[-c(1:burn.in),]
    q_t1 = as.matrix(stats::quantile(t1,c(0.025,0.975)))
  }
  main = apply(q_t1, 2, fun)
 q_t2=c()
  for(j in 1:ncol(sg2)){
    t2=as.matrix(sg2[,j])
    t2 = t2[-c(1:burn.in),]
    q_t2 = as.matrix(stats::quantile(t2,c(0.025,0.975)))
  }
 inter = apply(q_t2, 2, fun)
 inde = c(main,inter)
 inde
}
