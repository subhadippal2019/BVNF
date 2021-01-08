
#############################################################################################################
##############################################################################################################
############## Geometric Augmentation ########################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################


#' Fits a Data Augmentation algorithm for Vonmises Distribution in 3 dimensions.
#'   Generates MCMC samples from Posterior of mode and concentration parameter.
#'
#' @param Y:  n \times 3 matrix containing n directional data of dimension d.
#' norm of each row of the matrix Y is 1.
#' @param kappa_start: A positive number, starting value for the MCMC algorithm.
#' If set NUll the Default procedure to get the initial value will be used.
#' @param MCSamplerSize: Number of MCMC samples required to be generated.
#' @return MCMC samples from for  Posterior of mode and concentration parameter of Vonmises distribution.
#' @examples
#' library(Rfast)
#' library(benchmarkme)
#' data = rvmf(n =10, mu=c(1,0,0),k = 10)
#' MC_OBject = VONF_pD_MPG_INVCDF_DA_POSTERIOR( Y=data, MCSamplerSize=50 )
#' @export
VONF_3D_GEO_DA_POSTERIOR<-function(Y,MCSamplerSize=100, kappa_start=NULL){
  # Initially this function was called VONF_GEO_DA_Mu_kappa_post_sampler
  #browser()
  start_time=Sys.time()
  if(is.null(kappa_start)){(kappa_start=KAPPA_INITIAL(Y))}
  if(is.character(kappa_start)){kappa_start=KAPPA_INITIAL(Y)}
  #kappa_start=20
  beta_prior_kappa=0;alpha_prior_kappa=1;
  #datasummary
  Y_bar=apply(Y,2,'mean');Y_SUM= apply(Y,2,'sum');n=dim(Y)[1]

  #initialization
  Z=replicate(n,-1); kappa=kappa_start;McSample_kappa=replicate(MCSamplerSize,-1);McSample_mu=replicate(MCSamplerSize,0*Y_bar)

  #Sampling loop
  for(iter in 1:MCSamplerSize){
    for(i in 1:n){    Z[i]=rgeom(n=1, prob=(1-exp(-2*kappa)))  }
    #browser()
    #Z=replicate(n = n, expr = rgeom(n=1, prob=(1-exp(-2*kappa))))
    #### mu: mean direction sampler #####################
    Post_mean_dir=Y_bar/norm_vec(Y_bar);  Post_concentration=n*kappa*norm_vec(Y_bar);
    mu=rvmf(n=1, mu=Post_mean_dir, k=Post_concentration)
    #### kappa: concentration sampler #####################
    #beta_Kappa_post= calculate_beta_kappa_post_GEO(Z,mu,Y_SUM)
    beta_Kappa_post=sum(2*Z+1)-(mu)%*%(Y_SUM)+beta_prior_kappa
    alpha_kappa_post=n+alpha_prior_kappa
    kappa=rgamma(n = 1,shape = alpha_kappa_post,rate = beta_Kappa_post)


    #################################################
    #### storing MC sample###########################
    #################################################
    McSample_kappa[iter]=kappa;McSample_mu[,iter]=mu;
  } # End iter loop

            Run_Time= Sys.time()-start_time;
            Sys_info=get_sys_details()
            #########
            #browser()
            function_name=as.character(match.call()[1])
            function_def0<-  capture.output(print(get( function_name)))
            function_def0[1]=paste0( function_name,"<-",function_def0[1])
            function_def=paste( function_def0 , collapse = "\n")
            #function_def1=dput(get(function_name))

  #########
  lst=list(McSample_kappa=McSample_kappa,
           McSample_mu=McSample_mu,
           Y=Y,
           Method="Geometric RV Augmentation: function:VONF_3D_GEO_DA_POSTERIOR",
           Run_Time=Run_Time,
           System_info=Sys_info,
           call=match.call(),
           function_def=(function_def))
  return(lst)
}# End funtion










#########################################################################################
#########################################################################################
#########################################################################################
#Geometric GIG Mixture Augmentation ( Nmodified Polya Gamma) ############################
#########################################################################################
#########################################################################################
#########################################################################################

#' @export
VONF_3D_GEOGIG_DA_POSTERIOR<-function(Y,kappa_start=NULL,MCSamplerSize=5000 ){
  #kappa_start=20
  start_time=Sys.time()
  if(is.null(kappa_start)){(kappa_start=KAPPA_INITIAL(Y))}
  if(is.character(kappa_start)){kappa_start=KAPPA_INITIAL(Y)}
  beta_prior_kappa=0;alpha_prior_kappa=1;
  #datasummary
  Y_bar=apply(Y,2,'mean');Y_SUM= apply(Y,2,'sum');n=dim(Y)[1]

  #initialization
  Z=W=replicate(n,-1); kappa=kappa_start;McSample_kappa=replicate(MCSamplerSize,-1);McSample_mu=replicate(MCSamplerSize,0*Y_bar)

  #Sampling loop
  #browser()
  for(iter in 1:MCSamplerSize){
    for(i in 1:n){
         Z[i]=rgeom(n=1, prob=(1-exp(-2*kappa)));
         param_chi= kappa*(2*Z[i]+1)^2/2;
        W[i]=rgig(n=1,lambda = -.5,psi = 2*kappa, chi =param_chi);
    }


    #### mu: mean direction sampler #####################
    Post_mean_dir=Y_bar/Norm(Y_bar);  Post_concentration=n*kappa*norm_vec(Y_bar);
    mu=rvmf(n=1, mu=Post_mean_dir, k=Post_concentration)
    #### kappa: concentration sampler #####################
    #beta_Kappa_post=#calculate_beta_kappa_post(W,Z,mu,Y_SUM)
    beta_Kappa_post= sum(W)+.25*sum((2*Z+1)^2/W)-(mu)%*%(Y_SUM) +beta_prior_kappa
    alpha_kappa_post=3*n/2+alpha_prior_kappa
    kappa=rgamma(n = 1,shape = alpha_kappa_post,rate = beta_Kappa_post)


    #################################################
    #### storing MC sample###########################
    #################################################
    McSample_kappa[iter]=kappa;McSample_mu[,iter]=mu;
  } # End iter loop
  Run_Time= Sys.time()-start_time;
  Sys_info=get_sys_details()
  #########
  #browser()
  function_name=as.character(match.call()[1])
  function_def0<-  capture.output(print(get( function_name)))
  function_def0[1]=paste0( function_name,"<-",function_def0[1])
  function_def=paste( function_def0 , collapse = "\n")
  #function_def1=dput(get(function_name))

  lst=list(McSample_kappa=McSample_kappa,
           McSample_mu=McSample_mu,
           Y=Y,
           Method="GIG_PG_RV Augmentation: function:VONF_3D_GEOGIG_DA_POSTERIOR",
           Run_Time=Run_Time,
           System_info=Sys_info,
           call=match.call(),
           function_def=function_def
           )
  return(lst)
  #return(list(McSample_kappa=McSample_kappa,McSample_mu=McSample_mu))
}# End funtion


