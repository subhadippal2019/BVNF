
VONF_pD_MPG_INVCDF_DA_POSTERIOR<-function(Y,kappa_start=NULL,MCSamplerSize=5000, K=100){
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
  
  #Sampling loop
  j_nu_0=bessel_zero_Jnu(nu,1:K); 
  J_nuPlus1=besselJ(j_nu_0,(nu+1));
  
  for(iter in 1:MCSamplerSize){
    # for(i in 1:n){
    #   
    #   T[i]= PG_randomSample(a = kappa,nu=nu, K = K, j_nu_0 = j_nu_0, J_nuPlus1 =J_nuPlus1   );#PG_mod_randomSample: Modified Polya Gamma INVersion of CDF
    # }
     
    T=replicate(n,PG_randomSample(a = kappa,nu=nu, K = K, j_nu_0 = j_nu_0, J_nuPlus1 =J_nuPlus1 , eps_accuracy = .00001  ))
    #### mu: mean direction sampler #####################
    Post_mean_dir=Y_bar/Norm(Y_bar);  Post_concentration=n*kappa*Norm(Y_bar);
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
  lst=list(McSample_kappa=McSample_kappa,McSample_mu=McSample_mu, Y=Y, Method="PG Augmentation InversionOfCDS: function:VONF_pD_MPG_INVCDF_DA_POSTERIOR",Run_Time=Run_Time)
  
  return(lst)
  #return(list(McSample_kappa=McSample_kappa,McSample_mu=McSample_mu))
}# End funtion




VONF_3D_MPG_ES_DA_POSTERIOR<-function(Y,kappa_start=NULL,MCSamplerSize=5000, K=100){
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
    T=replicate(n, MPG_sample_nu_Half(alpha = kappa ))
    #### mu: mean direction sampler #####################
    Post_mean_dir=Y_bar/Norm(Y_bar);  Post_concentration=n*kappa*Norm(Y_bar);
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
  lst=list(McSample_kappa=McSample_kappa,McSample_mu=McSample_mu, Y=Y, Method="PG Augmentation InversionOfCDS: function:VONF_pD_MPG_INVCDF_DA_POSTERIOR",Run_Time=Run_Time)
  
  return(lst)
  #return(list(McSample_kappa=McSample_kappa,McSample_mu=McSample_mu))
}# End funtion


