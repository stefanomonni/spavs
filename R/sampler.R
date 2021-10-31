##########################################################################################################
##########################################################################################################

sampler=function(n_iter_t, invtemp,  configuration, structure, logpostvector, logpost, proposed, accepted,  X, Y, params, print, thin,  outputfile){
    
  conf=configuration
  struct=structure
  p=params$p 
     
  for(iter in 1:n_iter_t){   
     tmp=allowed_moves(struct,p)
     whichmove= sample(tmp$moveset,1)
 
     if(whichmove=="m1"){           
          obj=MergeType1(invtemp, conf, struct, logpostvector, logpost, allowed_comp=tmp$comps$m1, mselect=tmp$number$m1, X, Y, params,    
                         pm=length(tmp$moveset))
          proposed$m1=proposed$m1+1
          accepted$m1=accepted$m1+obj$new         
     }else if(whichmove=="s1"){  
          obj=SplitType1(invtemp, conf, struct, logpostvector, logpost, allowed_comp=tmp$comps$s1, mselect=tmp$number$s1, X, Y, params, 
                        pm=length(tmp$moveset) ) 
          proposed$s1=proposed$s1+1
          accepted$s1=accepted$s1+obj$new  
     }else if(whichmove=="s2"){
          obj=SplitType2(invtemp, conf, struct, logpostvector, logpost, allowed_comp=tmp$comps$s2, mselect=tmp$number$s2, X, Y, params,  
                        pm=length(tmp$moveset)) 
          proposed$s2=proposed$s2+1
          accepted$s2=accepted$s2+obj$new           
     }else{
          obj=MergeType2(invtemp, conf, struct, logpostvector, logpost, K=tmp$number$m2, X, Y, params, pm=length(tmp$moveset))   
          proposed$m2=proposed$m2+1
          accepted$m2=accepted$m2+obj$new    
    }

    conf=obj$conf
    
    struct=obj$struct
    logpostvector=obj$logpostvector
    logpost=obj$logpost
    
    if(print){
      if((iter%%thin)==0){
   
       filename=paste(outputfile, "_conf.txt", sep="") 
       confvec=conf_to_vector(conf)              
       write(file=filename, x=confvec, ncolumns=length(confvec),  append=TRUE)
       
       filename=paste(outputfile, "_logpost_size.txt", sep="") 
       write(file=filename, x=c(logpost, struct$card), ncolumns=2,  append=TRUE)             
       }
    }    
 }

 output=list(conf=conf, struct=struct, logpostvector=logpostvector, logpost=logpost, proposed=proposed, accepted=accepted)
  
 return(output)
} 

##########################################################################################################

