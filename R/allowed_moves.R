allowed_moves=function(struct, p){ 
     allowed=list()
     number=list()
     
     m1=0; s1=0; m2=0; s2=0;
     
     
     allowed$m1=which(struct$xcard<p)
     number$m1=length(allowed$m1)     
     if(number$m1>0) m1=1
      
      
     allowed$s1=which(struct$xcard>0)
     number$s1=length(allowed$s1)
     if(number$s1>0) s1=1
   
     
     allowed$s2=which(struct$ycard>1)
     number$s2=length(allowed$s2)
     if(number$s2>0) s2=1
             
     number$m2= struct$card 
     allowed$m2=1:number$m2
     if(number$m2>1)  m2=1
     
     moveset=c(rep("s1", s1), rep("m1", m1), rep("s2", s2), rep("m2", m2))
     
     return(list(comps=allowed, number=number, moveset=moveset) ) 
    
}



