    
 parallel_sampler=function(tc, n_iter_t, inversetemp,  conf, struct, logpostvector, logpost, proposed, accepted, X, Y, params, printf, thin,   outputfile){
      
     filename= paste(outputfile, paste("_t", tc, sep="" ), sep="")
     res=sampler(n_iter_t, invtemp =inversetemp[tc], configuration=conf[[tc]], structure=struct[[tc]], logpostvector=logpostvector[[tc]],    logpost  =logpost[tc], proposed=proposed[[tc]], accepted=accepted[[tc]], X, Y,  params,  print=printf[tc], thin,  
     outputfile=filename)
    
     return(res)
        
  }
    
