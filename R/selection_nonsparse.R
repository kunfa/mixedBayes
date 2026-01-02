selection_nonsparse = function(obj){
  sg1 = obj$posterior$GS.gamma2
  sg2 = obj$posterior$GS.gamma3
  q_1=c()
  for(j in 1:ncol(sg1)){
    t1=as.matrix(sg1[,j])
    q_t1 = as.matrix(stats::quantile(t1,c(0.025,0.975)))
    q_1 = cbind(q_1,q_t1)

  }
  main = apply(q_1, 2, fun)
  q_2=c()
  for(j in 1:ncol(sg2)){
    t2=as.matrix(sg2[,j])
    q_t2 = as.matrix(stats::quantile(t2,c(0.025,0.975)))
    q_2 = cbind(q_2,q_t2)
  }
 inter = apply(q_2, 2, fun)
 index = c(main,inter)
 index
}
