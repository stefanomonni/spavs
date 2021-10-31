##########################################################################################################################################
#' Posterior sampling of configurations of sets of responses and associated regressors.
#'
#' \code{spavs}  is the main function of the package. It runs a Markov chain or multiple chains corresponding either to different temperatures when parallel tempering is used or to different initializations when tempering is not used.  At each iteration,  a move is  proposed that is chosen with equal probability among those allowed. A move can be a merge or a split move of one of two kinds (type 1 and of type 2). The details of the algorithm are described in the paper by \href{https://doi.org/10.1214/09-BA416}{S. Monni and M. G. Tadesse (2009)}.  As far as possible, the input statistical parameters of the \code{spavs} function follow the notation in that paper.
#' This function does not return any R objects. Instead the MCMC output (i.e., the sampled configurations)  and other information are written into output files as detailed below.
#' If needed for analysis, the sampled configurations can then be loaded in R from these output files as lists, employing the
#' function \code{\link{conf.file.to.list}}.  Estimates of posterior probabilities can be computed with the function \code{\link{pairwise.probs}}.

#' @param   Y   the  matrix  of outcomes: an \eqn{n*q} matrix, with \eqn{q}  the number of reponses and n the number of samples (the transpose of the matrix Y in the paper).
#'
#' @param   X   the  matrix of predictors: an \eqn{n*p} matrix,  with \eqn{p} the number of predictors (the transpose of the matrix X in the paper).
          
#' @param   h0,h   hyperparameters that control the variance of the normal priors of the  regression coefficients.
#'                 \code{h0} controls the variance of the intercepts  and \code{h} the variance of the regression coefficients of the regressors X.   

#' @param   sigma0.square,nu  the shape and scale hyperparameters for the inverse-Gamma prior of the component residuals variances \eqn{sigma_k^2}. 
#'                    The prior mode is \eqn{nu/(sigma0.square+1)}.         
#' @param   rho  a real number \eqn{0< rho < 1}. It controls the prior probability:  \eqn{P(C)=rho^{m*n}}, for an (m, n) configuration C.     

#' @param   alpha0   a length-q numeric vector. If  NULL (default) the entries are set to 0.
#'                  It is the mean of the normal prior of  the intercepts of the mean function of  \eqn{(Y_1, \ldots, Y_q)}.
#'
#' @param   beta0     a length-p numeric vector. If NULL (default) the entries are set to 0.
#'                It is the mean of the normal prior of the vector of the regression coefficients for predictors \eqn{(X_1, \ldots, X_p)}.
#'
#' @param   p.split    controls the number of common predictors in the  components resulting from  a split move of type 2.
#'                     \code{p.split} must be in the interval \eqn{(0, 1)}. If the component to be split has \eqn{m} regressors,
#'                     the number  of common regressors in the two resulting components is \eqn{k}, with probability,  when \eqn{m>0},
#'                     \deqn{P(k)=p.split/(1+\floor(m/2)),} for \eqn{k=0,1, \ldots, \floor(m/2)}
#'                     \deqn{P(k)=(1-p.split)/(m-\floor(m/2)),} for \eqn{k=\floor(m/2), \ldots, m},
#'                     and, when \eqn{m=0},  \eqn{k=0} with probability 1.
#' @param  tempering   TRUE/FALSE specifies whether or not tempering is to be employed.      
#' @param   inv.temp   a numeric vector of either length 1 or length equal to \code{tchains}.
#'                     It is considered only when \code{tempering=TRUE}.
#'                     If it has length 1, it must be  a positive number in the interval  \eqn{(0, 1/tchains]}, 
#'                     in which case it represents the difference between inverse temperatures of the chains, so                      
#'                     that the tempering schedule is  \eqn{1/T_i= 1- (i-1)*inv.temp},  for   \eqn{i=1, \ldots, tchains}.
#'                     If it has length equal to \code{tchains}, it represents  the vector of the inverse temperatures in the chains: 
#'                     one entry must then be equal to one, and the others must be unique values in \eqn{(0, 1)}.
#'                     The program will order the entries  in decreasing order.

#' @param   tchains    an integer greater than 0  indicating the number of  chains. 
#'                     If \code{tempering=TRUE} the chains are run at different temperatures,  otherwise no tempering swaps 
#'                     are attempted.
#' @param   tswaps     total number of temperature swaps (moves that exchange configurations at different temperatures),
#'                     when \code{tempering=TRUE} in which case it must be \eqn{>=1}. The parameter can be used even
#'                     when  \code{tempering=FALSE} and no temperature swaps occur. The total number of iterations in each chain 
#'                     is given by \eqn{(tswaps+1)*mws}.
#' @param   mws        number of moves in each chain within attempts at temperature-transfers.
#'                     If tempering is not implemented, one can either set \eqn{tswaps=0}, so that \code{mws}  is the number of iterations,
#'                     or one can use \eqn{tswaps >0}, in which case the total number of iterations is given by \eqn{(tswaps+1)*mws}.
#'                      This latter option is preferable when one wants to monitor the acceptance rates, which, with the option        \code{print.monitoring=TRUE}, will appear on the screen.

#' @param  parallel.params   a list with three entries,  \code{cores}, \code{preschedule}, \code{affinity.list}.
#'               They are the parameters for the parallelization,  with different chains run in parallel.
#'               The code is parallelized with the function \code{parallel::mclapply} and these three values correspond to \code{mc.cores}, 
#'               \code{mc.preschedule},  and \code{affinity.list} in that function. See \code{\link[parallel]{mclapply}}.
#'               In particular, \code{cores} is the number of cores to use. Ideally it should be  equal to \code{tchains},
#'               and there is no gain in selecting it larger than \code{tchains}.
#'              If \code{preschedule} and \code{affinity.list} are not specified, the default values of \code{mclapply} are used.
#'              All other parameters in \code{mclapply}  are the default ones. 
#' @param  init.conf  one of the strings:  "random", predetermined", "data".
#'         If \code{init.conf="random"}, the initial configuration is randomly chosen, with the initial number of groups specified in 
#'                 \code{init.conf.params}.
#'             If \code{init.conf="predetermined"}, the initial configuration is loaded from a file, specified in 
#'                   \code{init.conf.params}.
#'                The file must have one configuration per line  and the number of lines must be equal to the number of tempering chains.
#'                Each initial configuration, one per tempering chain, must appear in the same form as the software output: e.g.,
#'                       5 , 3 ; 7 , 2 ; , 1 ; see below the description of the output files.
#'                There is no check that the file is correct other than possible errors thrown by the function used to read
#'                the file, which is the function  \code{\link{conf.file.to.list}}. 
#'             If \code{init.conf="data"}, the starting configuration (the same for all  chains) 
#'                is determined using k-means on the Y data.
#'                This option is available only when the number of responses \eqn{q \ge 3}.

#' @param init.conf.params
#'             If  \code{init.conf="random"}, this parameter must be  an integer or a vector of integers of length equal to \code{tchains} (number of tempering chains). 
#'             It denotes the number of components of the initial configuration in each chain.
#'             If more than one chain is run,  and the vector has length 1, all initial configurations have the same number of components equal to this integer, but they are chosen randomly. 
#'             If \code{init.conf="predetermined"}, it must be the name of the file that contains the initial configuration(s) from which to start.
#'             If \code{init.conf="data"}, this parameter is ignored.                     

#' @param    thin   an integer indicating the thinning used for the sampler
#'                   (namely, the configuration of the chain at \eqn{T=1}  are output every \code{thin} iterations).
#'           If one wants to have any sampled configurations saved in the file, choose \code{thin} such that \code{thin} \eqn{\le}  \code{mws}.

#' @param  outputfile   a string, the prefix to the names of the 3 files where the sampled configurations are written. 
#'        
#' @param   print.monitoring  TRUE/FALSE.  If TRUE,  print on screen the acceptance rates after each temperature swap.   
#'       
#' @return none.     The routine writes the samples and other information in  files and does not save R objects.
#'                   Use \code{\link{conf.file.to.list}} to load the saved configurations into R in the form of lists.
#'        There are three types of output files, all with the same prefix specified in the option  outputfile. They are:
#'
#'        -"outputfile_tk_conf.txt" will contain the sampled configurations
#'         of the k-th chain.  If 'tempering=TRUE', only one file for T=1 is
#'         printed.
#'         Each configuration is in one line.
#'         Components are separated by a semicolon. A comma separates
#'         the indices of the regressors and responses of a component.
#'         E.g., the line:  5 , 3 ; 5 7 , 2 ; , 1 ; represents:
#'         C1= X_5, Y_3;  C2=X_5, X_7, Y_2;   C3=Y_1.
#'
#'        -"outputfile_logpost_size.txt",  a two-column file
#'          with each row having logpost and size of the sampled configuration 
#'          saved in the corresponding row of the "outputfile_t_k_conf.txt" file.
#'
#'        -"outputfile_ar.txt" tabulates the acceptance rates for the moves in 
#'          all chains and, with tempering, for the swaps between thermal chains.
#'          Each row  has the number of attempted and accepted moves of a specific
#'          type, for all different chains in order of increasing temperature 
#'          (the first column being thus the T_1=1 chain).
#'          The last row gives proposed and accepted numbers of swaps between 
#'          contiguous tempering chains, in order of increasing temperature: 
#'          the first two columns refer to the swap between the T_1=1 and T_2
#'          chains, the second two to the swap between the T_2 and T_3 chains... 
#'
#'@references 
#'S. Monni, M. G. Tadesse, A Stochastic Partitioning Method to Associate High-dimensional Responses and Covariates, Bayesian Analysis (2009) 4, pp. 413-436. 
#' \href{https://doi.org/10.1214/09-BA416}{https://doi.org/10.1214/09-BA416}

#'
#' @export



#' @examples
#'# four chains at temperatures 1, 0.9, 0.85, 0.81, 
#'# with a parallelization using 4 cores  
#'spavs(Y=Y, X=X, h0=2, h=2, nu=1, sigma0.square=3, rho=0.3, tchains=4, tempering=TRUE, inv.temp=c(1,0.9, 0.85, 0.81), init.conf="random", init.conf.params=c(3,3,1,4), parallel.params=list(cores=4))
#'
#' 

#' @examples
#'# four chains at the same temperature T=1 (no tempering) 
#'# with a parallelization using 4 cores.
#'# this is equivalent to 4 different runs
#'spavs(Y=Y, X=X, h0=2, h=2, nu=1, sigma0.square=3, rho=0.3, tchains=4, tempering=FALSE, init.conf="random", init.conf.params=c(3,3,1,4), parallel.params=list(cores=4))
#'
#'
#' @examples
#'# one  chain 
#'spavs(Y=Y, X=X, h0=2.2, h=2, nu=2, sigma0.square=1.3, rho=0.01, p.split=0.8, tempering=FALSE, tchain=1, init.conf="data", outputfile="run1")
#'
  
             
spavs=function(Y, X, h0, h, nu, sigma0.square, rho,  alpha0=NULL, beta0=NULL, p.split=1/2, tempering=FALSE,  inv.temp=1, tchains=1, tswaps=0, mws=30000, parallel.params=list(cores=1, preschedule=TRUE, affinity.list=NULL), init.conf=c("random", "predetermined", "data"), init.conf.params, thin=200, outputfile="output", print.monitoring=FALSE){

  tmp=check_input(Y=Y, X=X,  h0=h0, h=h, nu=nu, sigma_0.square=sigma0.square, rho=rho,  alpha_0=alpha0, beta_0=beta0, psplit=p.split, tempering=tempering, inv.temp=inv.temp, tchains=tchains, tswaps=tswaps, mws=mws,  thin=thin, cores=parallel.params$cores)
 
  mc.ps=ifelse(is.null(parallel.params$preschedule), TRUE, parallel.params$preschedule)
  mc.aff.list=parallel.params$affinity.list   


  X=tmp$X
  Y=tmp$Y 
  params=tmp$params 
  temp.params=tmp$temp.params  
  tchains=temp.params$tchains
  inversetemp=temp.params$inversetemp                                                          
  tswaps=temp.params$tswaps                                    
  mws=temp.params$mws    
  
  thin=tmp$thin  
  cores=tmp$cores
  outputfile=as.character(outputfile)


  tmp=starting_conf(initial.configuration=init.conf,  initial.configuration.params= init.conf.params,  X=X, Y=Y, params=params, tchains=tchains)
  logpost= tmp$logpost             
  logpostvector=tmp$logpostvector
  conf=tmp$conf
  struct= tmp$struct
  

  proposed=list()
  accepted=list()       
  for(tc in 1:tchains){     
    proposed[[tc]]=list(m1=0, m2=0, s1=0, s2=0)  
    accepted[[tc]]= list(m1=0, m2=0, s1=0, s2=0)   
   } 
   
   tempered_prop=rep(0, tchains-1)
   tempered_accepted=rep(0, tchains-1)    
      
   #we sample only from the chain at Temp =1  if tempering=TRUE
   printf=rep(TRUE, tchains)
  
   if(tempering)  printf[2:tchains]= FALSE
    
 
   
 ######################################################################################################################################### 
 ####sampler starts 
 for(iter  in 1:(tswaps+1)) {   
            
   tmp= parallel::mclapply(1:tchains, parallel_sampler, n_iter_t= mws,  inversetemp,  conf,  struct, logpostvector, logpost, proposed, accepted, X, Y, params, printf, thin,  outputfile,  mc.cores = cores, mc.preschedule=mc.ps, affinity.list=mc.aff.list)   
   
     
  for(tc in 1:tchains){ 
      logpostvector[[tc]]= tmp[[tc]]$logpostvector     
      logpost[tc]=  tmp[[tc]]$logpost  
      proposed[[tc]]= tmp[[tc]]$proposed
      accepted[[tc]]=tmp[[tc]]$accepted
      conf[[tc]]=tmp[[tc]]$conf
      struct[[tc]]=tmp[[tc]]$struct  
    }
    
    rm(tmp)
           
    #attempts the (tchains -1) swaps if there is tempering         
     tc=ifelse(tempering, tchains, 0)
     
     while(tc >1){          
       tempered_prop[tc-1]= tempered_prop[tc-1]+1
       deltainvtemp= inversetemp[tc-1] -inversetemp[tc] 
      
       lprobaccept= deltainvtemp*(logpost[tc]-logpost[tc-1])
       
       lu=log(runif(1))          
       if(lu < lprobaccept){  
         tempered_accepted[tc-1]= tempered_accepted[tc-1]+1
        
         tmp=conf[[tc]]
         conf[[tc]]= conf[[tc-1]]
         conf[[tc-1]]=tmp
        
         tmp=struct[[tc]]
         struct[[tc]]=struct[[tc-1]]
         struct[[tc-1]]=tmp
        
         tmp=logpostvector[[tc]]
         logpostvector[[tc]] =logpostvector[[tc-1]]
         logpostvector[[tc-1]]=tmp
        
         tmp=logpost[tc]
         logpost[tc]=logpost[tc-1]
         logpost[tc-1]=tmp
       }               
    tc=tc-1   
   }       #end swaps
   
   if(print.monitoring) print(write_matrix_rates(accepted, proposed, tempered_accepted,  tempered_prop))
 }  
 
   filename=paste(outputfile, "_ar.txt", sep="")
   matrates=  write_matrix_rates(accepted, proposed, tempered_accepted,  tempered_prop)
   write.table(file=filename, x=matrates, quote=FALSE, col.names=TRUE) 
   
} 


