##########################################################################################################
#the log-marginal for a component.  does not include log(rho)*n*m
log_marginal=function(comp, X, Y, params){
   h0=params$h0                    
   h=params$h                      
   nu=params$nu                    
   sigma_0.square=params$sigma02   
   N=params$N
   nk=length(comp$y)
   mk=length(comp$x)
   
   nnh=N*nk*0.5
   dum=1/(N+1/h0) 


   alpha0=params$alpha0[comp$y]
   YC=  as.matrix(Y[ ,comp$y])
   ycalpha=apply(YC,2,sum) +alpha0/h0
   Omega=  sum(alpha0^2)/h0+  sum(YC^2)-dum*sum(ycalpha^2) 
    
   loglik=lgamma(nnh+ sigma_0.square)-lgamma(sigma_0.square)-nnh*log(2*pi)-mk*.5*log(h)-nk*.5*log(N*h0+1)+sigma_0.square*log(nu)  
    
   if(mk>0){   
      beta0=params$beta0[comp$x]                           
      XC=  as.matrix(X[ ,comp$x])   
      tXC= t(XC) 
      A= diag(mk)/h +nk*tXC%*%(diag(N) -dum*matrix(1, ncol=N, nrow=N))%*%XC      
      V=-dum*sum(ycalpha)*apply(XC,2, sum)+apply(tXC%*%YC,1, sum) +beta0/h        
      Omega= Omega +sum(beta0^2)/h  -drop(V%*%solve(A)%*%V)    
      loglik= loglik -.5*determinant(A, log=TRUE)$modulus[1]      
   }
             
  
   loglik=loglik-(sigma_0.square+nnh)*log(nu+Omega*.5)
     
   return(loglik) 
}

##########################################################################################################
#write posterior of the configurations (up to constant) and each component contribution
log_posterior=function(X, Y, conf, conf_struct,  params){
  rho=params$rho  
  k=conf_struct$card  
  mvec=conf_struct$xcard
  nvec=conf_struct$ycard
 
  lpcomp=rep(NA, k)
  for(i in 1:k)  lpcomp[i]= log_marginal(comp=conf[[i]], X=X, Y=Y, params=params) +mvec[i]*nvec[i]*log(rho)   
     
  return(list(logpconf=sum(lpcomp), logpcomp=lpcomp))

}


#old one using mvtnorm #not used anymore
#log_marginal2=function(comp, X, Y, params){
#   h0=params$h0                    
#   h=params$h                      
#   nu=params$nu                    
#   sigma_0.square=params$sigma02   
#   N=params$N
#   alpha0=params$alpha0
#   beta0=params$beta0
                        
#   nk=length(comp$y)
#   mk=length(comp$x)
#   XC=  as.matrix(X[ ,comp$x])
#   YC=  Y[ ,comp$y]
 
#   theta0=c(alpha0[comp$y],beta0[comp$x])  

   #to avoid problem when nk==1 and mk==0, specify ncol 
#   HO=diag(c(rep(h0, nk), rep(h, mk)), ncol=nk+mk)
    
#   oneN=matrix(ncol=1, nrow=N, 1)
#   onenk=matrix(ncol=1, nrow=nk,1)  
#   W1=kronecker(oneN, diag(nk))
#   W2=kronecker(XC, onenk )  
#   W=cbind(W1, W2)    
#   Sigma= W%*%HO%*%t(W)+diag(nk*N)
#   Sigma=Sigma*(nu/sigma_0.square)
   
#   delta= W%*%theta0
   #place the matrix Y as a vector one: one row after the other 
   # (Y_11, Y_12, ..., Y_1n Y_21,..)  
#   y=as.vector(t(YC))  
   
#   loglik=mvtnorm::dmvt(x=y, sigma=Sigma,   df=2*sigma_0.square, delta=delta,   log = TRUE, type = "shifted", checkSymmetry = FALSE)   
#   return(loglik)
#}



