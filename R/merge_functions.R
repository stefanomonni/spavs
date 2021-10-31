
##########################################################################################################
##############################  Type  1   merge 
MergeType1=function(invtemp, conf, struct, logpostvector, logpost,  allowed_comp, mselect, X, Y, params, pm){
    
   rho=params$rho
   p=params$p
     
   newconf=conf
   newlogpostvector=logpostvector

  
   comp=allowed_comp[sample.int(length(allowed_comp),1)]
   
   xold= conf[[comp]]$x
   xoldcomplement=setdiff(1:p, xold) 
 
   newconf[[comp]]$x=c(xoldcomplement[sample.int(length(xoldcomplement), 1)], xold)             
   m=struct$xcard[comp]
   n=struct$ycard[comp] 
   
   logpostcomp= logpostvector[comp]         
   newlogpostcomp=  log_marginal(newconf[[comp]], X, Y, params) +log(rho)*n*(m+1)     
   newlogpostvector[comp]=  newlogpostcomp 
   newstruct=  conf_structure(newconf) 
   
  
   #number of ++ conf is the same or one more if the select comp was (0,n)     same as length(which(newstruct$xcard>0))
   pospos=length(which(struct$xcard!=0)) +as.numeric(length(xold)==0)     
          
   lambdaratio=   log(mselect)+ log(p-m)  -log(m+1)  -log(pospos)
   lambdaratio=lambdaratio +log(pm)-log(length(allowed_moves(newstruct,p)$moveset))   #addition for the moves randomly picked
        
   lpaccpt=  newlogpostcomp- logpostcomp                     
   lu=log(runif(1))           
   new=0  
   if(lu< lambdaratio+ invtemp*lpaccpt){        
      new=1          
      conf=newconf  
      struct=newstruct  
      logpostvector=newlogpostvector
      logpost=sum(newlogpostvector)        
    }
    
  return(list(conf=conf, struct=struct, logpostvector=logpostvector, logpost=logpost, new=new))
       
}


##########################################################################################################
##########################################################################################################
############################ MOVE TYPE 2 MERGE 
MergeType2=function(invtemp,conf, struct, logpostvector, logpost, K, X, Y, params, pm){   
   
    rho=params$rho
    p=params$p
       
    index=sample(1:K, 2, replace=FALSE) 
        
    comp1=conf[[index[1]]]
    comp2=conf[[index[2]]]
    logpostcomp1=logpostvector[index[1]]
    logpostcomp2=logpostvector[index[2]]

    newy=c(comp1$y, comp2$y)
    newx=unique(c(comp1$x, comp2$x))
    n1=struct$ycard[index[1]]
    n2=struct$ycard[index[2]]   
    m1=struct$xcard[index[1]]
    m2=struct$xcard[index[2]]
    
    m=length(newx)
    cnumber=m1+m2-m
             
    newconf=list()
    newlogpostvector=rep(NA, K-2)
     
    #move the components that have not changed and attach the new one as last 
    tmp=1
    for(i in  setdiff(1:K, index)){    
      tmplist=list()
      tmplist$x= conf[[i]]$x
      tmplist$y=conf[[i]]$y
       
      newconf[[tmp]]=tmplist
      newlogpostvector[tmp]= logpostvector[i]
      tmp=tmp+1    
   }
  
   tmplist=list()
   tmplist$x= newx      
   tmplist$y=newy   
   newconf[[tmp]]=tmplist
   
   newlogpostcomp=  log_marginal(tmplist, X, Y, params) + log(rho)*(n1+n2)*m    
   newlogpostvector[tmp]=newlogpostcomp
   
   lpaccpt= newlogpostcomp-logpostcomp1-logpostcomp2
      
   n_new_ngeq2=length(which(struct$ycard>1))+1 -as.integer(n1>1) -as.integer(n2 >1)
         
   newstruct=  conf_structure(newconf)   
    
   qsplit=qsplit(nvar_tosplit=m, commonvars=cnumber, psplit=params$psplit)
   lambdaratio=lchoose(K, 2)+log(2) -log(n1+n2-1) -lchoose(n1+n2, n1) +log(qsplit)-log(n_new_ngeq2)    
   lambdaratio=lambdaratio +log(pm)-log(length(allowed_moves(newstruct,p)$moveset))   #addition for the moves randomly picked
     
   lu=log(runif(1))      
    
   new=0
   if(lu< lambdaratio+invtemp*lpaccpt){    
       new=1
       conf=newconf  
       struct=newstruct   
       logpostvector=newlogpostvector
       logpost=sum(newlogpostvector) 
    }
    
  return(list(conf=conf, struct=struct, logpostvector=logpostvector, logpost=logpost, new=new))  
  
}

