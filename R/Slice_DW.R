




Dir_data_inf_DW_slice<-function(data,MaxIter=100 ){
     start_Time =Sys.time()
      ###############Data averages############################################
        theta=data;
        n=nrow(theta);
        theta_bar=apply(theta,2, mean);
        nu=length(theta_bar)/2-1
      ########################################################################
      ########################################################################
        mu_n=theta_bar/norm_vec(theta_bar)
        R_n=n*norm_vec(theta_bar)
      ##### Initial value for slice sampler ###################################
      #########################################################################
        w=5;v=1;mu=mu_n; kappa=KAPPA_INITIAL(data)
        x=1
        #print(kappa)
      ##### storing different components#######################################
        x_all = log_v_all = w_all = kappa_all = mu_all = NULL
      #########################################################################
       #browser()
      for(McLen in 1:MaxIter){
          #### Sampling x #######################################
                  #print(paste0("nnn=",(n-1)*log(w)))
                  #if( (n-1)*log(w) < -100 ){x=0}
                  #if((n-1)*log(w) > -100) {  x <- runif(n = 1, min = 0,max = w^(n-1))}
                  x_pow_1_by_n_1<-sample_x_pow_1_by_n_1(w, n)



          ######%%%%%%%%%%%%%%%%#################################
          ##### sampling log_v  #################################
                  log_v_upper <- (R_n*kappa*(1+ inner_prod(mu,mu_n)   ))
                  log_v      <- log_v_upper-rexp(n=1, rate=1)
          #######v=runif(n = 1, min = 0,max = v_upper )
          #######################################################################
          ######### sampling mu ######################################### hh#####
                  mu<-sample_mu(mu_n, kappa, R_n, log_v)
          #######################################################################
                  Exp_rate <- (alt_besselI(kappa, nu = nu)-1) ## nu=d/2-1; here d=length(mu)
                  #print(Exp_rate)
                  E        <- rexp(n=1, rate = Exp_rate);
                  M        <- E+w
                  #print(M)
          ########################################################################
                # w<-runif(n = 1, min = x^(1/(n-1)),   max = M)
                  #print(x)
                  w<-rTruncatedExp(n=1, rate = 1, lowerLim = x_pow_1_by_n_1, UpperLim = M)
          #########################################################################
          #########################################################################

                   seq_k           <- 1:100;
                   seq_lambda_k    <- generate_Lambda_Seq(seq_k, nu)
                       #u_k_upper_bound <- ifelse(seq_lambda_k>0, exp(-w*seq_lambda_k*kappa^(2*seq_k)), 1)
                       u_k_upper_bound=exp(-w*seq_lambda_k*kappa^(2*seq_k))
                   u_k_seq         <- apply(matrix(u_k_upper_bound,ncol=1), 1, function(yy){runif(1, 0,as.numeric(yy))})
                   #print(which(u_k_seq==min(u_k_seq)))
                   #Sys.sleep(.1)
                   kappa_upper     <- min(    (-  log(u_k_seq)/(w*seq_lambda_k) )^(1/(2*seq_k))) ### N in the notation of the paper
                   kappa_lower     <- max(log_v/( R_n+R_n* inner_prod(mu,mu_n)  ), 0)
                          if(kappa_lower> kappa_upper){kappa_upper=kappa_lower+.00001}
                  #kappa<-rexptr(n = 1,lambda = R_n, range = c( kappa_lower, kappa_upper))
                  kappa=rTruncatedExp(n = 1,rate=R_n,  lowerLim =kappa_lower, UpperLim =kappa_upper  )



                 #print(kappa)
           ##################################################################################################
           ############ Storing #############################################################################
               x_all     = c(x_all,x);
               log_v_all = c(log_v_all, log_v);
               w_all     = c(w_all, w);
               kappa_all = c(kappa_all, kappa)
               mu_all    = cbind(mu_all, mu)

      }
        Run_Time= Sys.time()-start_Time
 MC=list(mu=mu_all, kappa=kappa_all, x=x_all, log_v=log_v_all, w=w_all )
 Vonf_MC=list(MC=MC, data=data, Run_Time=Run_Time, method="Damien Walker Slice Sampler")
 return(Vonf_MC)
}



generate_Lambda_Seq<-function(k_seq, nu ){
  val1=apply(matrix(k_seq, nrow=length(k_seq)),1,  function(k){val=ifelse(k<=100, 0.5^(2*k+nu)/(gamma(k+1+nu)*gamma(k+1)), 0);return(val)} )
  return(val1)
}



 sample_x_pow_1_by_n_1<-function(w,n){
   #browser()
   #x <- runif(n = 1, min = 0,max = w^(n-1))
   log_x=(n-1)*log(w)-rexp(n=1, rate=1)
   x_pow_n_1=exp(log_x/(n-1))
   return(x_pow_n_1)
 }



 sample_mu<-function(mu_n, kappa, R_n, log_v){
    #browser()
   # two domensional case
   pert<-acos( log_v/(R_n*kappa)-1)
   mu_rad=atan(mu_n[2]/ mu_n[1])
   theta_lower<-mu_rad-pert; theta_upper<-mu_rad+pert
   mu_sample<-runif(n = 1, min = theta_lower, max = theta_upper)
   mu_coordinate=c(cos(mu_sample), sin(mu_sample))
   return(mu_coordinate)
 }


 log_runif<-function(n=1,Upper ){
   #Lower=0; Upper>0
   if(Upper<=0){return(NULL)}
   log_val=log(Upper)-rexp(n=1,rate=1)
   return(log_val)
 }



 ##################

 #
 # rTruncatedExp<-function(n=1, rate=1,lowerLim=0, UpperLim=Inf, eps=1e-100 ){
 #   a=lowerLim;b=UpperLim
 #   #(   exp(-rate*a)-exp(-rate*x)  )/(   exp(-rate*a)-exp(-rate*b)    )=u
 #   #-log(exp(-rate*a)-log(u)*Norm_const)/rate=(x)
 #   #exp(-rate*a)-log(u)*Norm_const=exp(-rate*x)
 #   #browser()
 #   u=runif(n)
 #   Norm_const= exp(-rate*a)-exp(-rate*b)
 #   if(is.na(Norm_const)){browser()}
 #   if(Norm_const<eps){return(discreteTruckExp(rate,lowerLim,UpperLim, prec_len=100))}
 #   x=-log(exp(-rate*a)-(u)*Norm_const)/rate
 #   return(x)
 # }
 #
 # discreteTruckExp<-function(rate,lowerLim,UpperLim, prec_len=100){
 #   seq_rng=seq(lowerLim, UpperLim, length.out = prec_len)
 #   log_d_val= -seq_rng*rate
 #   d_val=exp(log_d_val-max(log_d_val))
 #   sampled_bin=sample(x = seq_rng, size = 1, prob =d_val )
 #
 #   val=runif(1,sampled_bin, sampled_bin+ (UpperLim-lowerLim)/prec_len   )
 #   if(is.na(val)){browser()}
 #   return(val)
 # }
 #



 ##############

 #' @title  Random generator for a Truncated Exponential distribution.
 #
 #' @description Simulate random number from a truncated Exponential distribution.
 #
 #' @details
 #' It provide a way to simulate from a truncated Exponential
 #' distribution with given pameter \eqn{\lambda} and the range \eqn{range}.
 #' This will be used during the posterior sampling in th Gibbs sampler.
 #'
 #' @param n       int, optional \cr
 #'                number of simulations.
 #'
 #' @param lambda  double, optional \cr
 #'                parameter of the distribution.
 #'
 #' @param range   array_like, optional \cr
 #'                domain of the distribution, where we truncate our
 #'                Exponential. \eqn{range(0)} is the min of the range
 #'                and \eqn{range(1)} is the max of the range.
 #'
 #'
 #' @return \code{rexptr} returns the simulated value of the
 #' distribution:
 #'
 #' \item{u}{  double \cr
 #'            it is the simulated value of the truncated Exponential
 #'            distribution. It will be a value in \eqn{(range(0),
 #'            range(1))}.
 #'             }
 #'
 #'
 #' @references
 #' \itemize{
 #'
 #'  \item   [1] Y. Wang, A. Canale, D. Dunson.
 #'          "Scalable Geometric Density Estimation" (2016).\cr
 #'          The implementation of rgammatr is inspired to the Matlab
 #'          implementation of rexptrunc by Ye Wang.
 #' }
 #'
 #' @author L. Rimella, \email{lorenzo.rimella@hotmail.it}
 #'
 #' @examples
 #' #Simulate a truncated Exponential with parameters 0.5 in the range
 #' #5,Inf.
 #' #Set the range:
 #' range<- c(1,Inf)
 #'
 #' #Simulate the truncated Gamma
 #' set.seed(123)
 #' vars1<-rexptr(1000,0.5,range)
 #'
 #' #Look at the histogram
 #' hist(vars1,freq=FALSE,ylim=c(0,2),xlim = c(0,5))
 #' lines(density(vars1))
 #'
 #' #Compare with a non truncated Exponential
 #' set.seed(123)
 #' vars2<-rexp(1000,0.5)
 #'
 #'
 #' #Compare the two results
 #' lines(density(vars2),col='red')
 #'
 #' #Observation: simulate without range is equivalent to simulate from
 #' #rexp(1000,0.5)
 #'
 #' @import stats
 #'
 #' @export

 rexptr <- function(n=1, lambda=1, range=NULL)
 {
    #*************************************************************************
    #***        Author:      L. Rimella <lorenzo.rimella@hotmail.it>       ***
    #***                                                                   ***
    #***        Supervisors: M. Beccuti                                    ***
    #***                     A. Canale                                     ***
    #***                                                                   ***
    #*************************************************************************

    if(is.null(range))
    {
       warning("No range provided in the truncated Exponential. A simple Exponential random variable is simulated.")
       return(rexp(n, rate= lambda))
    }

    a1= range[1]
    a2= range[2]
    smallvalue = 1e-8

    cdf1 = pexp(a1, rate= lambda)
    cdf2 = pexp(a2, rate= lambda)

    if(cdf2-cdf1< smallvalue)
    {
       u = a1
    }

    # Otherwise we simulate from a truncated Exp according to the
    # transformation:
    # Ftrunc(x)= \frac{FExp(x)-FExp(a1)}{FExp(a2)-FExp(a1)}

    else
    {
       u = qexp(cdf1+ runif(n)*(cdf2- cdf1), rate= lambda);
    }

    return(u)

 }



 ###### app ######
 #
 # theta=rvmf(n=1000, mu=c(1,0,0), k = 5)
 # n=nrow(theta);theta_bar=apply(theta,2, mean);nu=length(theta_bar)/2-1
 #
 # plot(f, xlim=c(0, 1))
 #
 # f<-function(kappa){
 #
 #   A=log(besselI(x=n*kappa*norm_vec(theta_bar), nu=nu, expon.scaled = TRUE))+n*kappa*norm_vec(theta_bar) + nu*n*log(kappa)
 #   B=n*log(besselI(x=kappa, nu=nu, expon.scaled = TRUE)) + n*kappa  + nu*log( n*kappa*norm_vec(theta_bar))
 #   return((A-B-300))
 # }
 #
 #
 # g<-function(kappa){
 #
 #   A=n*kappa*norm_vec(theta_bar)* inner_prod(theta_bar, mu)
 #   B=n*log(besselI(x=kappa, nu=nu, expon.scaled = TRUE))+n*kappa
 #   return((A-B))
 # }
