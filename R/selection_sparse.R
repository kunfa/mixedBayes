selection_sparse = function(obj,burn.in=obj$burn.in){
  sg1 = obj$posterior$GS.beta
  sg1[which(sg1!=0)]=1
  sg2 = obj$posterior$GS.eta
  sg2[which(sg2!=0)]=1
  main=c()
  for(j in 1:ncol(sg1)){
    t1=as.matrix(sg1[,j])
    t1 = t1[-c(1:burn.in),]
    q_t1 = mpm(t1)
    main = c(main,q_t1)
  }
  inter=c()
  for(j in 1:ncol(sg2)){
    t2=as.matrix(sg2[,j])
    t2 = t2[-c(1:burn.in),]
    q_t2 = mpm(t2)
    inter = c(inter,q_t2)
  }
  inde = c(main,inter)
  inde
}
