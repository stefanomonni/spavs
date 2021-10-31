###########################################################################################################################
#random intialization
#creates k groups and assign randomly the y and x samples  a random number of X from 0 to pmax  either p or floor(p/N) is p>N


conf_initialize=function(k, p, q, N){

  yalloc=c(1:k, sample(1:k, q-k, replace=TRUE))
  yalloc= sample(yalloc)
  
  pmax=ifelse(p<N, p, floor(p/N))
  
  conf=list()    
  for(i in 1:k)  conf[[i]]=list(x= sample(1:p, sample(0:pmax,1),replace=FALSE), y =which(yalloc==i) )
  
  return(conf)
}




#creates k groups and assign randomly the y and x samples  a random number of X from 1 to p
#conf_initialize=function(k, p, q){

#  yalloc=c(1:k, sample(1:k, q-k, replace=TRUE))
#  yalloc= sample(yalloc)
    
#  conf=list()    
#  for(i in 1:k)  conf[[i]]=list(x= sample(1:p, sample(1:p,1),replace=FALSE), y =which(yalloc==i) )
  
#  return(conf)
#}




