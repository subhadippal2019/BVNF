
#' Fits a Data Augmentation algorithm for Vonmises Distribution in general dimension
#'   Generates MCMC samples from Posterior of mode and concentration parameter.
#'
#' @param Y:  n\times d  matrix containing n directional data of dimension d.
#' norm of each row of the matrix Y is 1.
#' @param kappa_start: A positive number, starting value for the MCMC algorithm.
#' If set NUll the Default procedure to get the initial value will be used.
#' @param MCSamplerSize: Number of MCMC samples required to be generated.
#' @return MCMC samples from for  Posterior of mode and concentration parameter of Vonmises distribution.
#' @examples
#' library(Rfast)
#' library(benchmarkme)
#' library(Bessel)
#' library(gsl)
#' data=rvmf(n =10, mu=c(1,0,0,0,0),k = 10)
#' MC_OBject=VONF_pD_MPG_INVCDF_DA_POSTERIOR(Y=data, MCSamplerSize=50)
#' @export

VONF_pD_MPG_INVCDF_DA_POSTERIOR<-function(Y,kappa_start=NULL,MCSamplerSize=5000, K=100,eps_accuracy=.00000001){
  #K is the number of terms for  calculating CDF
  Start_Time=Sys.time()
  if(is.null(kappa_start)){(kappa_start=KAPPA_INITIAL(Y))}
  if(is.character(kappa_start)){kappa_start=KAPPA_INITIAL(Y)}
  beta_prior_kappa=0;alpha_prior_kappa=1;
  #datasummary
  Y_bar=apply(Y,2,'mean');Y_SUM= apply(Y,2,'sum');n=dim(Y)[1]
  nu=length(Y_bar)/2-1;
  #initialization
  T=replicate(n,-1); kappa=kappa_start;McSample_kappa=replicate(MCSamplerSize,-1);McSample_mu=replicate(MCSamplerSize,0*Y_bar)
  #browser()
  #Sampling loop
  j_nu_0=bessel_zero_Jnu(nu,1:K);
  J_nuPlus1=BesselJ(j_nu_0,(nu+1));

  for(iter in 1:MCSamplerSize){
    # for(i in 1:n){
    #
    #   T[i]= PG_randomSample(a = kappa,nu=nu, K = K, j_nu_0 = j_nu_0, J_nuPlus1 =J_nuPlus1   );#PG_mod_randomSample: Modified Polya Gamma INVersion of CDF
    # }

    T=replicate(n,PG_randomSample(a = kappa,nu=nu, K = K, j_nu_0 = j_nu_0, J_nuPlus1 =J_nuPlus1 , eps_accuracy = eps_accuracy  ))
    #### mu: mean direction sampler #####################
    Post_mean_dir=Y_bar/norm_vec(Y_bar);  Post_concentration=n*kappa*norm_vec(Y_bar);
    mu=rvmf(n=1, mu=Post_mean_dir, k=Post_concentration)
    #### kappa: concentration sampler #####################
    #beta_Kappa_post=#calculate_beta_kappa_post(W,Z,mu,Y_SUM)
    ExtGamma_a= sum(T);ExtGamma_b= as.double((mu)%*%(Y_SUM)-beta_prior_kappa ) ;

    kappa=rExtendedGamma(1 ,a=ExtGamma_a,b =ExtGamma_b )
    #print(kappa)

    #################################################
    #### storing MC sample###########################
    #################################################
    McSample_kappa[iter]=kappa;McSample_mu[,iter]=mu;
  } # End iter loop
  Run_Time=Sys.time()-Start_Time
  #Sys_info=get_sys_details()


######################
  function_name=as.character(match.call()[1])
  function_def0<-  capture.output(print(get( function_name)));
  function_def0[1]=paste0( function_name,"<-",function_def0[1])
  function_def=paste( function_def0 , collapse = "\n")
  #function_def1=dput(get(function_name))
  #########

  lst=list(McSample_kappa=McSample_kappa,
           McSample_mu=McSample_mu,
           Y=Y,
           Method="PG Augmentation InversionOfCDS: function:VONF_pD_MPG_INVCDF_DA_POSTERIOR",
           Run_Time=Run_Time,
          # System_info=Sys_info,
           call=match.call(),
           function_def=function_def
           )

  return(lst)
  #return(list(McSample_kappa=McSample_kappa,McSample_mu=McSample_mu))
}# End funtion







#' Fits a Data Augmentation algorithm for Vonmises Distribution in general dimension
#'   Generates MCMC samples from Posterior of mode and concentration parameter.
#'
#' @param Y:  n\times 3  matrix containing n directional data of dimension 3.
#' norm of each row of the matrix Y is 1.
#' @param kappa_start: A positive number, starting value for the MCMC algorithm.
#' If set NUll the Default procedure to get the initial value will be used.
#' @param MCSamplerSize: Number of MCMC samples required to be generated.
#' @return MCMC samples from for  Posterior of mode and concentration parameter of Vonmises distribution.
#' @examples
#' library(Rfast)
#' library(benchmarkme)
#' library(Bessel)
#' library(gsl)
#' data=rvmf(n =10, mu=c(1,0,0),k = 10)
#' MC_OBject=VONF_3D_MPG_ES_DA_POSTERIOR(Y=data, MCSamplerSize=50)
#' @export


VONF_3D_MPG_DA_POSTERIOR<-function(Y,kappa_start=NULL,MCSamplerSize=5000, K=100){
  #K is the number of terms for  calculating CDF
  Start_Time=Sys.time()
  if(is.null(kappa_start)){(kappa_start=KAPPA_INITIAL(Y))}
  if(is.character(kappa_start)){kappa_start=KAPPA_INITIAL(Y)}
  beta_prior_kappa=0;alpha_prior_kappa=1;
  #datasummary
  Y_bar=apply(Y,2,'mean');Y_SUM= apply(Y,2,'sum');n=dim(Y)[1]
  nu=length(Y_bar)/2-1;
  if(nu!=0.5){ print("This is not 3 dimensional data");return(NULL)}
  #initialization
  T=replicate(n,-1); kappa=kappa_start;McSample_kappa=replicate(MCSamplerSize,-1);McSample_mu=replicate(MCSamplerSize,0*Y_bar)

  #Sampling loop
  #j_nu_0=bessel_zero_Jnu(nu,1:K);
  #J_nuPlus1=besselJ(j_nu_0,(nu+1));

  for(iter in 1:MCSamplerSize){
    # for(i in 1:n){
    #
    #   T[i]= PG_randomSample(a = kappa,nu=nu, K = K, j_nu_0 = j_nu_0, J_nuPlus1 =J_nuPlus1   );#PG_mod_randomSample: Modified Polya Gamma INVersion of CDF
    # }

    #T=replicate(n,PG_randomSample(a = kappa,nu=nu, K = K, j_nu_0 = j_nu_0, J_nuPlus1 =J_nuPlus1 , eps_accuracy = .00001  ))
    #T=replicate(n, MPG_sample_nu_Half(alpha = kappa ))
    T=rMPG_nu_Half(n = n, alpha = kappa)
    #### mu: mean direction sampler #####################
    Post_mean_dir=Y_bar/norm_vec(Y_bar);  Post_concentration=n*kappa*norm_vec(Y_bar);
    mu=rvmf(n=1, mu=Post_mean_dir, k=Post_concentration)
    #### kappa: concentration sampler #####################
    #beta_Kappa_post=#calculate_beta_kappa_post(W,Z,mu,Y_SUM)
    ExtGamma_a= sum(T);ExtGamma_b= as.double((mu)%*%(Y_SUM)-beta_prior_kappa ) ;

    kappa=rExtendedGamma(1 ,a=ExtGamma_a,b =ExtGamma_b )
    #print(kappa)

    #################################################
    #### storing MC sample###########################
    #################################################
    McSample_kappa[iter]=kappa;McSample_mu[,iter]=mu;
  } # End iter loop
  Run_Time=Sys.time()-Start_Time
 # Sys_info=get_sys_details()

  #########
  #browser()
  function_name=as.character(match.call()[1])
  function_def0<-  capture.output(print(get( function_name)));
  function_def0[1]=paste0( function_name,"<-",function_def0[1])
  function_def=paste( function_def0 , collapse = "\n")
  #function_def1=dput(get(function_name))
  #########
  lst=list(
    McSample_kappa=McSample_kappa,
    McSample_mu=McSample_mu,
    Y=Y,
    Method="PG Augmentation 3D_MPG_Sampler: function:VONF_3D_MPG_ES_DA_POSTERIOR",
    Run_Time=Run_Time,
    #System_info=Sys_info,
    function_def=function_def)

  return(lst)
  #return(list(McSample_kappa=McSample_kappa,McSample_mu=McSample_mu))
}# End funtion


