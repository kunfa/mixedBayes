selection_sparse = function(obj){
  sg1 = obj$posterior$GS.gamma2
  sg1[which(sg1!=0)]=1
  sg2 = obj$posterior$GS.gamma3
  sg2[which(sg2!=0)]=1
  main=c()
  for(j in 1:ncol(sg1)){
    t1=as.matrix(sg1[,j])
    q_t1 = mpm(t1)
    main = c(main,q_t1)
  }
  inter=c()
  for(j in 1:ncol(sg2)){
    t2=as.matrix(sg2[,j])
    q_t2 = mpm(t2)
    inter = c(inter,q_t2)
  }
  index = c(main,inter)
  index
}
