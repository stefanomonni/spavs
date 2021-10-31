
starting_conf=function(initial.configuration,  initial.configuration.params, X, Y, params, tchains){

 ic=match.arg(initial.configuration, c("random", "predetermined", "data"))

 if(params$q<3 && ic=="data"){
   cat(" the number of responses is 2, the data determined configuration is not implemented ")
   cat(" choose the  random option or load a configuration ")
   stop()     
 }
 
 switch(ic,                       
        random=starting_conf_random(initial.configuration.params, X, Y, params, tchains),
        predetermined=starting_conf_upload(initial.configuration.params, X, Y, params, tchains),
        data=starting_config_data(X, Y, params, tchains)  )

} 
 
 
  
starting_conf_random=function(initial.configuration.params, X, Y, params, tchains){
       
 conf=list()
 struct=list()
 logpost=rep(NA, tchains)
 logpostvector=list()

 p=params$p
 q=params$q
 N=params$N
 
 icp=initial.configuration.params 
 if(!is.numeric(icp)  ) {
     cat("initial.configuration.params  must be: \n")
     cat("a numeric vector or a number specifying the number of  components of the  random starting  configuration(s) \n")
     cat("when the  option initial.configuration=random")
     stop()
   } 
 
 
  icp=round(icp)
  if(any(icp<1) || any(icp>q))  stop("number of components must be between 1 and the number of responses ")
   
  if(!(length(icp) ==1 || length(icp)==tchains)){
     stop("initial.configurations.params must be a natural number, or a vector of natural numbers of length tchains ")
  }
       
  if(length(icp)==1){   
     
      for (tc in 1:tchains){
          tmp=conf_initialize(k=icp, p=p, q=q, N=N)
          tmpstr=conf_structure(tmp)
          tmplp= log_posterior(X=X, Y=Y, conf=tmp, conf_struct=tmpstr,  params=params)
            
          conf[[tc]]= tmp
          struct[[tc]]=tmpstr
          logpostvector[[tc]]= tmplp$logpcomp
          logpost[tc]=  tmplp$logpconf          
       }
  }
  else{
        for (tc in 1:tchains){
          tmp=conf_initialize(k=icp[tc], p=p, q=q, N=N)
          tmpstr=conf_structure(tmp)
          conf[[tc]]= tmp
          struct[[tc]]=tmpstr
          tmplp= log_posterior(X=X, Y=Y, conf=tmp, conf_struct=tmpstr,  params=params)
          logpostvector[[tc]]= tmplp$logpcomp
          logpost[tc]=  tmplp$logpconf         
        }
  }
   
   
  return(list(conf=conf, struct=struct, logpostvector=logpostvector, logpost=logpost))            
}
 
 
 
 
starting_conf_upload=function(initial.configuration.params, X, Y, params, tchains){

 if(!is.character(initial.configuration.params) ) {
     cat("initial.configuration.params must be  the name (with path) of a file from which to read the initial configuration(s): \n")
     cat("if initial.configuration=predetermined  is chosen\n")     
     stop()
 } 

      
 conf=list()
 struct=list()
 logpost=rep(NA, tchains)
 logpostvector=list()

 p=params$p
 q=params$q
 
 tmp=confdata_file_to_list(initial.configuration.params)  
             

 if(length(tmp)!=tchains){ 
         stop('number of loaded configurations is not equal to the number of tempering chains')
 }
        
 for(tc in 1:tchains){
            conf[[tc]]= tmp[[tc]]
            struct[[tc]]=conf_structure(tmp[[tc]])
            tmplp= log_posterior(X=X, Y=Y, conf=conf[[tc]], conf_struct= struct[[tc]], params=params)
            logpostvector[[tc]]= tmplp$logpcomp
            logpost[tc]=  tmplp$logpconf        
 }
         
 cat("Warning: The software does not check  your initial configuration(s) is in the right form. Use at your own risk \n")


 return(list(conf=conf, struct=struct, logpostvector=logpostvector, logpost=logpost))
}
 
 

starting_config_data=function(X, Y, params, tchains){
       
 q=params$q  
 tY=t(Y) 
 
 if(q==3){ 
    ktry=2 
 } else if(q==4){
          ktry=2:3
 }else if (q %in% 5:13) {
        ktry=2:4
 }else{ 
       kmin=2; kmax=ceiling(q/4)
       kspace = ceiling((kmax-kmin)/10)
       ktry = seq(kmin, kmax, by=kspace)
 }
 
 kwss <- sapply(ktry, function(k)  
       kmeans(tY, k, nstart=20)$tot.withinss)
 
 relD <- sapply(2:length(kwss), function(i) 
    round(-(kwss[i]-kwss[i-1])/kwss[i], 2))
 
 ksel = ktry[max(which(relD/sum(relD)>0.1))]
  
  # fit cluster 
 kclu <- kmeans(tY, ksel, nstart=20)$cluster
 Yalloc <- sapply(1:ksel, function(k) which(kclu==k))
  
  # association with X
 mean.clu <- sapply(1:ksel, function(k) 
    apply(as.matrix(Y[,Yalloc[[k]]]), 1, mean))
 
 Xassoc <- sapply(1:ksel, function(k) {
 unip <- apply(X, 2, function(j) summary(lm(mean.clu[,k] ~ j))$coef[2,4])
 unip.fdr <- p.adjust(unip, method="fdr")
            which(unip.fdr < 0.05)})
  
   
 tmp=list()   
 for(i in 1:length(Xassoc)){  
           tmp[[i]]=list(x=as.numeric(Xassoc[[i]]), y=as.numeric(Yalloc[[i]]) )  
 }  

 tmpstr=conf_structure(tmp)  
 
 conf=list()
 struct=list()
 logpost=rep(NA, tchains)
 logpostvector=list() 
   
 for (tc in 1:tchains){
          conf[[tc]]= tmp
          struct[[tc]]=tmpstr
          tmplp= log_posterior(X=X, Y=Y, conf=tmp, conf_struct=tmpstr,  params=params)
          logpostvector[[tc]]= tmplp$logpcomp
          logpost[tc]=  tmplp$logpconf         
 }
    
 return(list(conf=conf, struct=struct, logpostvector=logpostvector, logpost=logpost))

}
 
