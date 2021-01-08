

###########################################################################
###########################################################################
############################# Accesory functions ##########################
###########################################################################





###########################################################################
###########################################################################
############################# Technical Functions ##########################
###########################################################################

#library(gsl)

  AllTerms_CDF_log<-function(t,a,nu,j_nu_0,J_nuPlus1){
          #j_nu_0=bessel_zero_Jnu(nu,1:K)
           #J_nuPlus1=besselJ(j_nu_0,(nu+1));
          #log(besselI(a,nu,TRUE))+a actually returns log(besselI(a,nu,FALSE))
          #Note for Future: to make this function efficint. In stead of loop, compute the privious function as a vector operation.
          #Note for Future: One can incorporate the j_nu_0 and J_nuPlus1 here as required, insead of taking it as output.
      Generic_log_constantTerm=log(besselI(a,nu,TRUE))+a     -nu*log(a) + log(2);
      termK=(nu+1)*log( j_nu_0)-log(  abs(J_nuPlus1)  )- log( a^2+(j_nu_0)^2)  -( a^2+(j_nu_0)^2)*t
      termK=termK+Generic_log_constantTerm
           #MAx_term_log=max(termK)
     return(termK)
  }

  ##############################################
  ##############################################

  log_pos<-function(x){
            #browser()
    if(x>0){log_pos=log(x)}
    if(x<=0){log_pos=log(-x);
           #log_pos=-Inf;
           # print('negative rgument to log. see function log_pos')
    }
   return(log_pos)
  }

  ##############################################
  ##############################################


  log_survival_nu<-function(t,a,nu,j_nu_0,J_nuPlus1, iflog=TRUE){
           #log_survival=log(1-cdf)
    K=length(j_nu_0);
     Terms=AllTerms_CDF_log(t,a,nu,j_nu_0,J_nuPlus1) #Terms=AllTerms_log(t,a,nu,K)
    Max=max(Terms);  Terms_adj=Terms-Max
            #################
            # while exponenting we want to take out the maximum and add it back to log.
            ################
     sign_of_Terms= sign(J_nuPlus1);

     cumSum=cumsum(exp(Terms_adj)*sign_of_Terms)
     cumavg1=(cumSum[2:K]+cumSum[1:(K-1)])/2
     cumavg2=(cumavg1[2:(K-1)]+cumavg1[1:(K-2)])/2
            #################
            # while exponenting we want to rak out the maximum and add it back to log.
            #################
            # browser()
      value=log_pos(cumavg2[K-2])+Max

     if(!iflog){ value=exp(value)}
     return(value)
  }


##############################################
##############################################
##############################################
##############################################

##############################################
##############################################
##############################################
##############################################





#To Show the convergence for the new method for computing term by term log CDF
#library(ggplot2)
t_conv_log_CDF_plot<-function(t,a,nu,K,iflog=TRUE){
  #shows the term profiles of the CDF
  j_nu_0=bessel_zero_Jnu(nu,1:K); J_nuPlus1=besselJ(j_nu_0,(nu+1));
  Terms=AllTerms_CDF_log(t, a, nu,j_nu_0,J_nuPlus1);Max=max(Terms); Terms_adj=Terms-Max
  browser()
   sign_of_Terms= sign(J_nuPlus1);

  cumSum=cumsum(exp(Terms_adj)*sign_of_Terms)
  cumavg1=(cumSum[2:K]+cumSum[1:(K-1)])/2
  cumavg2=(cumavg1[2:(K-1)]+cumavg1[1:(K-2)])/2
  sh=ggplot(data=NULL, mapping=(aes(x=1:(K-2), y=cumSum[1:(K-2)]))) + geom_point(col='red')+geom_line()
  sh=sh+geom_point( aes(x=1:(K-2), y=cumavg1[1:(K-2)]),col='blue')
  sh=sh+geom_point( aes(x=1:(K-2), y=cumavg2),col='purple')

  return(sh)

}



###########################################################
############################################################
###########################################################
##################### density ##############################
########################################
########################################
########################################
########################################








All_DensityTerms_nu_log<-function(t,a,nu,j_nu_0,J_nuPlus1){
  #j_nu_0=bessel_zero_Jnu(nu,1:K)
  #J_nuPlus1=besselJ(j_nu_0,(nu+1));
  #Note for Future: to make this function efficint. In stead of loop, compute the privious function as a vector operation.
  #Note for Future: One can incorporate the j_nu_0 and J_nuPlus1 here as required, insead of taking it as output.
  Generic_log_constantTerm=log(besselI(a,nu,TRUE))+a     -nu*log(a) + log(2);

  termK=(nu+1)*log( j_nu_0)-log(  abs(J_nuPlus1)  )  -( a^2+(j_nu_0)^2)*t
  termK=termK+Generic_log_constantTerm;
  #MAx_term_log=max(termK)
  return(termK)
}



########################################################################
############################################################################

dMPG<-function(t,a,nu,iflog=TRUE, NumofTerms=200){
  K=NumofTerms
  t=t*(t>.0001)+.0001*(t<=.0001)
  j_nu_0=bessel_zero_Jnu(nu,1:K);
  J_nuPlus1=besselJ(j_nu_0,(nu+1));
  dMPG_single<-function(t1, a, nu, iflog=iflog){
    return(log_f_nu(t=t1,a=a,nu=nu,j_nu_0=j_nu_0,J_nuPlus1=J_nuPlus1,iflog=iflog))
  }
  val=apply(matrix(t, ncol=1), MARGIN = 1, FUN = function(xvec){dMPG_single(xvec,a=a, nu=nu, iflog=iflog )})

  return(val)
}

log_f_nu<-function(t,a,nu,j_nu_0,J_nuPlus1,iflog=TRUE){
  K=length(j_nu_0);
  Terms=All_DensityTerms_nu_log(t,a,nu,j_nu_0,J_nuPlus1)
  Max=max(Terms);  Terms_adj=Terms-Max

  sign_of_Terms= sign(J_nuPlus1);

  cumSum=cumsum(exp(Terms_adj)*sign_of_Terms)
  cumavg1=(cumSum[2:K]+cumSum[1:(K-1)])/2
  cumavg2=(cumavg1[2:(K-1)]+cumavg1[1:(K-2)])/2
  value=log_pos(cumavg2[K-2])+Max
  if(!iflog){
    value=exp(value)
    #print("Inside log_f_nu")
  }
  return(value)
}

########################################
########################################
########################################
########################################



########################################
########################################
#################NEWTON RAPHSON########
########################################
t_start<-function(u,a,nu,j_nu_0_1,J_nuPlus1_1){
  #browser()
  #log(besselI(a,nu,TRUE))+a actually returns log(besselI(a,nu,FALSE))
  ##Note for Future:  To make this function efficient take the constants log(besselI(a,nu,TRUE))+a and add at th end with all terms
  Right_part1= log(besselI(x=a,nu=nu,expon.scaled=TRUE))+a     -nu*log(a) + log(2)+(nu+1)*log( j_nu_0_1)-log(  abs(J_nuPlus1_1)  )- log( a^2+(j_nu_0_1)^2)

  t_start=(Right_part1-log(1-u))/( a^2+(j_nu_0_1)^2)

  return(t_start)
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

PG_randomSample<-function(a,nu,K=200,j_nu_0=NULL, J_nuPlus1=NULL){
  #Generates a Sample from the appropriate distribution.
   #K=200;
  Total_NR_steps=30;
  u=runif(1)
  val_iter=1;t=.01
  ################################################
      if(is.null(j_nu_0)){
          j_nu_0=bessel_zero_Jnu(nu,1:K);
          J_nuPlus1=besselJ(j_nu_0,(nu+1));
      }
      if(is.null(J_nuPlus1)){
        J_nuPlus1=besselJ(j_nu_0,(nu+1));
      }
  ###############################################
  ###############################################
  ############ starting point ###################

  t=t_start(u,a,nu,j_nu_0[1],J_nuPlus1[1])
  if(t==-Inf){t=.01}
  t=round(t,10)
  #print(paste("start=",t))
  val_iter=t;NR_steps=1;STOP_FLAG=1
 # browser()
  while((NR_steps <= Total_NR_steps)*(STOP_FLAG)){
      log_survival=log_survival_nu(t,a,nu,j_nu_0,J_nuPlus1)
      t=t+( log_survival-log(1-u) )* exp( log_survival-log_f_nu(t,a,nu,j_nu_0,J_nuPlus1)   )
     # t=round(t,10)
      val_iter[NR_steps+1]=t
      #print(t)
       if( abs(t-val_iter[NR_steps])<.0001){STOP_FLAG=0}
        NR_steps=NR_steps+1;
  }
  return(t)
}


######################################################################################
######################################################################################
######################################################################################




