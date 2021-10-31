
################# WRITE and READ FUNCTIONS and POSTERIOR ANALYSIS

write_matrix_rates=function(accepted, proposed, tempered_accepted,  tempered_prop){

 k=length(accepted)
 mat= cbind(proposed[[1]], accepted[[1]])  
 nm=c(1,1)
 vec=matrix(nrow=1, ncol=2*k, 0)
 i=2 

 while( i<k+1){
   mat=cbind(mat, cbind(proposed[[i]], accepted[[i]]))
   vec[1, (2*i-3):(2*i-2)]= c(tempered_prop[i-1],  tempered_accepted[i-1])   
   nm=c(nm,i,i)
   i=i+1
 } 
 
 mat=rbind(mat, vec)
 colnames(mat)=nm
 rownames(mat)[5]=c("tswap")

 return(mat)
}


#write the configuration  as a numeric vectors
# a comma (,) separates x and y of a component and  a semicolon (;)  one comp from another
# 5 , 3 ; 7 , 2 ; , 1 ;  would be a 3-component conf with C1= [5][3] C2=[5][2]  C3 [][1] 
#so length(which(is.na(vec)) gives quickly the number of components, components are only (m, n) with n>0 and m>=0
conf_to_vector=function(conf){ 
    vec= append(append(append(conf[[1]]$x, "," ),conf[[1]]$y), ";")  
    i=2  
    while(i < length(conf)+1){  
     tmp= append(append(append(conf[[i]]$x, "," ),conf[[i]]$y), ";")
     vec=append(vec, tmp )
     i=i+1
     }
    return(vec) 
}





#this is to reload the  output samples from a file back into R as lists
#need to check how config (0, n) work since now this file set $x=numeric(0) need to check whether logposterior has

#' Sampled configurations as R lists.
#'
#' \code{conf.file.to.list} is used to load into R the sampled configurations saved by the function \code{\link{spavs}} in a file.

#' @param filename the file  where the sampled configurations are saved. 
#                   The file must be (in the same format as) the output file "..._conf.txt" of the function 'spavs'.
#' @param  burn.in the number of  sampled configurations  that are omitted (i.e., the first \code{burn.in} configurations are not loaded).
#'         If configurations are written out  at every iteration, this is the standard burn-in.
#' @return   The output is a list whose elements are lists. The i-th list, list[[[i]]], represents a configuration, whose elements are 
#'           its components:  list[[[i]]][[[j]]]$x is the vector of the indices of the predictors and  list[[[i]]][[[j]]]$y that of the responses,
#'           in the j-th component of the i-th configuration.

#' @export
#' @examples   conf.file.to.list("run1_t1_conf.txt")
       
conf.file.to.list=function(filename, burn.in=0){

  connection = file(filename, "r") 
  data=readLines(connection, n = -1L) 
  l= length(data)
  
  stp=as.integer(burn.in) +1
  if(stp<1) {cat("all saved configurations will be used ")
             stp=1
            }  
  if(stp>l) stop(" the parameter burn.in may not be larger than the total number of sampled configurations") 
    
  conf=list() 
  index=0
  for(i in stp:l){
    index=index+1
    conf[[index]]=list()
    tmp=unlist(strsplit(data[[i]], ", | ;"))
       
    for(j in 1:(length(tmp)/2)){            
      conf[[index]][[j]]=list( x =as.numeric(unlist(strsplit(trimws(tmp[2*j-1]), split=" +"))),
                               y =as.numeric(unlist(strsplit(trimws(tmp[2*j]), split=" +")))
                             )
    }  
  }   
 close(connection)
 
 return(conf)
}              


#'  Posterior Pair-wise Probabilities
#'
#' \code{pairwise.probs}  returns  estimates of   pairwise posterior probabilities 
#'  of any two variables belonging to the same component.
#'  Specifically, it outputs: the matrix PXX, whose \eqn{i,j} element gives the relative incidence of the covariates \eqn{X_i, X_j} occurring in the same component, across components and configurations;
#'  the matrix PYY, whose \eqn{i,j} element  estimates the probability that the responses \eqn{Y_i, Y_j} are in the same component;  
#'  the matrix PXY, whose \eqn{i,j} element  estimates the probability that the predictor \eqn{X_i} and response \eqn{Y_j} are in the same component.               

 
#' @param filename the file  where the sampled configurations are saved. 
#'                  The file must be (in the same format as) the file with suffix "_conf.txt" output by  \code{\link{spavs}}.
#' @param p,q  the number of regressors and responses, respectively.
#' @param  burn.in  the number of  sampled configurations  that are omitted in the computation of the pairwise probabilities
#'                  (i.e., the first \code{burn.in} configurations are not used to compute the probabilities).
#'                 If configurations are written out  at every iteration, this is the standard burn-in.
#' @return    The output is a list with elements: 
#'
#'     $samples  number of configurations on which the inference is based
#'               (number of configurations saved in the input file minus burn.in) 
#'
#'     $total_comp total number of components on which the inference is based
#'                
#'     $pxx the matrix PXX 
#'
#'     $pyy the matrix PYY 
#'
#'     $pxy the matrix PXY
#'
#' @export
#' @examples   pairwise.probs("run1_t1_conf.txt", p=12, q=7)

pairwise.probs= function(filename, p, q, burn.in=0){
  
  dt=conf.file.to.list(filename, burn.in)
  l=length(dt)    
  YY=matrix(ncol=q, nrow=q, 0)
  XY=matrix(nrow=p, ncol=q, 0)
  XX=matrix(nrow=p, ncol=p, 0)
  comp_number=0 
  for(i in 1:l){
     lst=dt[[i]]
     llst=length(lst)
     comp_number=comp_number+llst
     for(j in 1:llst){  
      x=lst[[j]]$x
      y=lst[[j]]$y 
      XX[x,x]=XX[x,x]+1
      XY[x,y]=XY[x,y]+1
      YY[y,y]=YY[y,y]+1
     }   
 }  
 
 res=list()
 res$samples=l
 res$total_comp= comp_number
 res$pxx=XX/comp_number
 res$pxy=XY/l
 res$pyy=YY/l
 return(res) 
}
