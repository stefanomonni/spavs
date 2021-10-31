##################################################### checks input data and parameters ################################
check_input=function(Y, X,  h0, h, nu, sigma_0.square, rho,  alpha_0, beta_0, psplit, tempering, inv.temp, tchains, tswaps, mws, thin, cores){

 if(is.null(Y)){stop("Y is missing") }
 if(is.null(X)){stop("X is missing") }   
 X=as.matrix(X)
 Y=as.matrix(Y) 
 if(mode(X)!="numeric") stop("X matrix is not numeric")
 if(mode(Y)!="numeric") stop("Y matrix is not numeric")
 p=ncol(X)
 q=ncol(Y)    
 N=nrow(X)
 if(N!=nrow(Y)){stop("Y and X must have same number of rows")}
 
  
 if(any(c(h0, h, nu,  sigma_0.square, rho) <=0))  stop(" 'h0', 'h', 'nu', 'sigma0.square', 'rho', must be positive ") 
 if(rho>1) stop("'rho' must be in (0, 1)")   
 if(is.null(alpha_0)) alpha_0=rep(0, q)
 if(is.null(beta_0))  beta_0=rep(0, p) 
 if(length(alpha_0)!=q) stop(" 'alpha_0' must be a vector of length equal to the number of responses ")
 if(length(beta_0)!=p) stop("'beta_0'  must be a vector of length equal to the number of predictors ")
 if(psplit>=1 || psplit<=0)  stop("'p.split' must be in the open interval (0,1) ")
 
 tchains=as.integer(tchains)
 tswaps=as.integer(tswaps)
 mws=as.integer(mws)
 thin=as.integer(thin)
 cores <- as.integer(cores)
   
 if(tempering){
   if(tchains<2) stop("When tempering=TRUE, 'tchains' must be >=2")
 }  
 
 if(tchains<1)   stop("Number of chains must be greater than or equal to one ") 
 if(!is.numeric(inv.temp)) stop("'inv.temp' must be numeric, a number or a vector of length tchains ")
 
 if(tempering){ 
    if(length(inv.temp)==1){
         if(inv.temp<=0 || inv.temp>1/tchains) stop("If 'inv.temp' is a  non-negative number, it must be between [0, 1/tchains] ")   
        
         inversetemp= 1 -(0:(tchains-1) )*inv.temp           
    }
    else if(length(inv.temp)==tchains){
         if(any(inv.temp<=0) || any(inv.temp>1))  stop("If 'inv.temp' is a vector, its entries must be in (0,1] ") 
         if(!any(inv.temp==1))  stop("If 'inv.temp' is  a vector, one of its entries must be equal to 1 ")      
         if(length(unique(inv.temp))!=tchains) stop("If 'inv.temp' is a vector, its entries must be different from each other ") 
         
          inversetemp=  sort(inv.temp, decreasing=TRUE)
    }else{   
         stop("'inv.temp' must be numeric, a number or a vector of length tchains when 'tempering=TRUE' ")
   }  
   
   if(tswaps==0 && tchains >1) stop("When tempering=TRUE,  'tswaps' should be >=1")
         
 } else{    
      inversetemp= rep(1, tchains)
 }
 
 
 
 if(tswaps<0)    stop(" 'tswaps' must be a non negative integer ")
 if(mws<0)       stop(" 'mws' must be a non negative integer ")
 if(thin>mws)    cat("Sampled configurations will not be printed to a file since thin>mws ")
 
 if(cores<1)  stop(" 'cores' in 'parallel.params' must be >= 1 ")
 if(cores>tchains) cores=tchains

 params=list(h0=h0, h=h, nu=nu, sigma02=sigma_0.square, rho=rho, q=q, p=p, N=N, alpha0=alpha_0, beta0=beta_0, psplit=psplit)   
 temp.params=list(inversetemp=inversetemp, tchains= tchains, tswaps=tswaps, mws=mws)
                                      
 return(list(Y=Y, X=X, params=params, temp.params=temp.params, thin=thin, cores=cores))

}
