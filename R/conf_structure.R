
##########################################################################################################
#returns the number of components and for each component the  number of Y and X not which ones
conf_structure=function(conf){

 k=length(conf) 
 xl=rep(0, k)
 yl=rep(0,k)
 
 for(i in 1:k){
   xl[i]=length(conf[[i]]$x)
   yl[i]=length(conf[[i]]$y)
 } 
  
 return(list(card= length(conf), xcard=xl,  ycard=yl ) )

}


#faster?
#conf_structure2=function(conf){

# tmp=mapply(conf, FUN=function(comp){return(c(length(comp$x), length(comp$y)))}) 
# return(list(card=length(conf), xcard=tmp[1,], ycard=tmp[2, ]))

 
#} 
