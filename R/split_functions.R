##########################################################################################################
## Type 1 split move
SplitType1=function(invtemp, conf, struct, logpostvector, logpost,  allowed_comp, mselect, X, Y, params, pm){   
    
   rho=params$rho
   p=params$p
        
   newconf=conf
   newlogpostvector=logpostvector  
 
   comp=allowed_comp[sample.int(length(allowed_comp),1)]
   xold= conf[[comp]]$x   
    
   m=struct$xcard[comp]
   n=struct$ycard[comp]
       
   index= sample(1:m,1)
   xnew=xold[setdiff(1:m, index)]            
   newconf[[comp]]$x=xnew
    
   logpostcomp= logpostvector[comp]
   newlogpostcomp=  log_marginal(newconf[[comp]], X, Y, params) +log(rho)*n*(m-1)  
 
   newlogpostvector[comp]=  newlogpostcomp 
     
   newstruct=  conf_structure(newconf)  
   
   #conf with length  p those that were p before minus eventually 1 if one split was p-length
   #k-nmber_p   are the number of condfiguration with fewer than p exes so we can add in the reverse
   #struct$card-number_p!=length(which(newstruct$xcard<p))) 
   number_p=length(which(struct$xcard==p)) -as.numeric(m==p)   
         
   lambdaratio=   log(mselect)+ log(m) -log(p-m+1) -log(struct$card-number_p)      
   lambdaratio=lambdaratio +log(pm)-log(length(allowed_moves(newstruct,p)$moveset))   #addition for the moves randomly picked
          
   new=0
   lpaccpt=  newlogpostcomp- logpostcomp 
   lu=log(runif(1))      
    
   if(lu< lambdaratio+invtemp*lpaccpt){     
       new=1
       conf=newconf  
       struct=newstruct  
       logpostvector=newlogpostvector
       logpost=sum(newlogpostvector)
   }
     
   return(list(conf=conf, struct=struct, logpostvector=logpostvector, logpost=logpost, new=new) )  
           
}


##########################################################################################################
##########################################################################################################
SplitType2=function(invtemp, conf, struct, logpostvector, logpost,  allowed_comp, mselect, X, Y, params, pm){ 

   rho=params$rho
   p=params$p
   newconf=conf
   newlogpostvector=logpostvector
   
   comp=sample2(allowed_comp, 1)   

   xold=conf[[comp]]$x
   yold=conf[[comp]]$y
    
   m=struct$xcard[comp]
   n=struct$ycard[comp]
  
    
   n1=sample(1:(n-1), 1)   
   ynew=sample2(yold, n1, replace=FALSE)
   ynewcomp=setdiff(yold, ynew)
   
   objsplit=split_vars(xold, params$psplit)   
   
   m1=length(objsplit$x1)
   m2=length(objsplit$x2)
    
   newconf[[comp]]$x=objsplit$x1
   newconf[[comp]]$y=ynew      
   newcomponent=list()
   newcomponent$x=  objsplit$x2  
   newcomponent$y=  ynewcomp     
   newconf[[struct$card+1]]=newcomponent    
   newstruct=  conf_structure(newconf) 
   
   qsplit=objsplit$qsplit 
   lambdaratio=log(mselect) + log(n-1)+lchoose(n, n1)
   lambdaratio=lambdaratio-lchoose(struct$card+1,  2)-log(2) -log(qsplit)
   lambdaratio=lambdaratio +log(pm)-log(length(allowed_moves(newstruct,p)$moveset))   #addition for the moves randomly picked
     
   logpostcomp= logpostvector[comp]
    
   newlogpostcomp1=  log_marginal(newconf[[comp]], X, Y, params) + log(rho)*n1*m1 
   newlogpostcomp2= log_marginal(newcomponent, X, Y, params) +log(rho)*(n-n1)*m2 
        
   newlogpostvector[comp]=  newlogpostcomp1
   newlogpostvector=c(newlogpostvector, newlogpostcomp2)
   
   lpaccpt= newlogpostcomp1+newlogpostcomp2-logpostcomp
      
   lu=log(runif(1))      
   new=0
    
   if(lu< lambdaratio+invtemp*lpaccpt){    
       new=1   
       conf=newconf  
       struct=newstruct   
       logpostvector=newlogpostvector
       logpost=sum(newlogpostvector)  
   }
  
   return(list(conf=conf, struct=struct,  logpostvector=logpostvector, logpost=logpost, new=new))   
}   


##########################################################################################################
##########################################################################################################
qsplit=function(nvar_tosplit, commonvars, psplit){   
 
     if(nvar_tosplit==0) return(1)                 #no variables to split so one possibilit for splitting empty groups
    
     if(psplit==0.5){
        pchosen=1/(nvar_tosplit+1)
     }   
     else{    
       fl=floor(nvar_tosplit/2)   
       if(commonvars <=fl){
         pchosen= psplit/(fl+1)
       }else{            
       pchosen= (1-psplit)/(nvar_tosplit-fl)  
      }        
    }  
  prob=pchosen*0.5^{nvar_tosplit-commonvars +1} /choose(nvar_tosplit,commonvars)    
  return(prob)       
}


##########################################################################################################
##########################################################################################################
split_vars=function(xvar, psplit) { 
  
  m=length(xvar)  
     
  if(m>0){  
     
     if(psplit==0.5){
         cr=sample(0:m,1)
         pchosen=1/(m+1)
      }else{
        fl=floor(m/2)
        tmp=rbinom(1, size=1, prob=psplit)            
        if(tmp==1){
          cr=sample(0:fl,1)
          pchosen=psplit/(fl+1)
        }else{
            cr=sample((fl+1):m,1)
            pchosen=(1-psplit)/(m-fl)
        }
     } 
     
    commonvars=sample2(x=xvar, size=cr) 
    uniqvars=setdiff(xvar, commonvars) 
    
    indicator= rbinom(length(uniqvars),1,prob=1/2)   
    whichgroup=sample(0:1,1)     
    xnew1=c(commonvars, uniqvars[which(indicator==whichgroup)])
    xnew2=c(commonvars, uniqvars[which(indicator!=whichgroup)])       
    qsplit= pchosen*(0.5)^{m-cr+1} /choose(m, cr)
              
  }else{  #m==0  
      xnew1=numeric(0)
      xnew2=numeric(0)
      qsplit=1
 }      
        
 return(list(x1= xnew1,  x2= xnew2, qsplit=qsplit) )                    
    
}  


