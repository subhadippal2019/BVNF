

#' rMPG_nu_Half(alpha) generates random sample from the Modified Ploya Gamma distribution.
#' @param n: Number of samples to be generated.
#' @param alpha: The parameter of the distribution, a positive number.
#' @return The sample from the MPG(nu=.5, \code{alpha}) dostribution.
#' @examples
#' library(ghyp)
#' rMPG_nu_Half(10,1)
#' rMPG_nu_Half(200,100)
#' rMPG_nu_Half(5,10000)
#' @export
rMPG_nu_Half<-function(n=1, alpha){


        rMPG_nu_Half_single<-function(alpha){
              # This MPG sampler is only available for dealing with three dimensional directional data.
              Accept=-1
              #############################################################################
              #############################################################################
              #browser()
              while(Accept!=1){
                #####
                X<-mix_Exp_Gig_sampler(alpha)
                density_f=ch_a_n(n = 0,t = X,alpha = alpha) #=S_0
                U=runif(1)*density_f
                S_n=density_f;n=0;flag=1;#Accept=-1
                ##############################################################
                ##############################################################
                while(   (flag==1)*(n<100 ) ){
                  n=n+1
                  if( !n%%2){ # if n even
                      S_n = S_n + ch_a_n(n = n,t = X,alpha = alpha)
                      if(U>S_n){flag=0;  Accept=0}
                  }
                  ############################
                  if( n%%2){ # if n odd
                    S_n = S_n - ch_a_n(n = n,t = X,alpha = alpha)
                    if(U<S_n){flag=0;  Accept=1}
                  }
                }
                ##############################################################
                ##############################################################
                if(Accept==1){return(X)}
              }#End While
              #############################################################################
              #############################################################################
        }
  #dummy_var=matrix(rep(1, n), ncol=1)
  #val=apply(dummy_var, MARGIN = 1, FUN =function(x){return(rMPG_nu_Half_single(alpha))} )
  val=replicate(n = n , rMPG_nu_Half_single(alpha))
  return(val)
} # End of functrion


###############################################################################################################
###############################################################################################################





sampling_gig_segment_selector<-function(alpha){
  if(alpha>100){return(1)}
  t_star=0.1005082
      #nu=1.5
  #browser()
      exp_rate=pi^2+alpha^2
      log_p_e=log(2*pi^2)- exp_rate*t_star  -log(exp_rate)
      log_p_g= 1.5*log(2*alpha)-.5*log(pi) + log_Bessel_K(alpha, nu=1.5)+
               log(pgig(t_star,lambda=-1.5, chi=.5, psi=2*alpha^2))
            M=max(log_p_e, log_p_g)
            p_g=exp( log_p_g-M ); p_e=exp( log_p_e-M );


            p_gig= p_g/(p_e+p_g)
       #print(c(log_p_e,log_p_g,p_gig))
      val=rbinom(n=1, size=1,prob = p_gig )
      #browser()
      return(val)
}





ch_a_n<-function(n=0,t, alpha){
  val=log_ch_a_n(n=n,t = t, alpha = alpha)
  return(exp(val))
}

ch_a_n_OLD<-function(n=0,t, alpha){ # this one also works efficiently after introducing ch_t =exp(log_ch_t)
  t_star=0.1005082; exp_rate=pi^2+alpha^2
  if(t<=t_star){
    log_ch_t<-  (log_sinh(alpha)-log(alpha))-log(2) -2.5*log(t)-.5*log(pi)-.25/t-t*alpha^2
    ch_t=exp(log_ch_t)
    #ch_t=(sinh(alpha)/alpha)*  ( 1/(2*t^2*sqrt(pi*t)) )  * exp(-.25/t -t*alpha^2)
    #browser()
        if(!n%%2){
          a_n=(n+1)^2 * exp(-(n+1)^2/(4*t) +  1/(4*t))
        }
        if(n%%2){
          a_n=2*t*exp(    -n^2/(4*t)+1/(4*t)    )
        }
    val_ch_a_n=a_n*ch_t
  }

  if(t>t_star){
    ch_t=( sinh(alpha)/alpha) * (   2*pi^2*exp(-exp_rate*t)    )
    a_n=((n+1)^2) * exp(   -(   (n+1)^2-1   )*pi^2*t  )
    val_ch_a_n=ch_t*a_n
  }

  return(val_ch_a_n)

}





log_ch_a_n<-function(n=0,t, alpha){
  t_star=0.1005082; exp_rate=pi^2+alpha^2
  if(t<=t_star){
    log_ch_t<-  (log_sinh(alpha)-log(alpha))-log(2) -2.5*log(t)-.5*log(pi)-.25/t-t*alpha^2
    #ch_t=exp(log_ch_t)
    #ch_t=(sinh(alpha)/alpha)*  ( 1/(2*t^2*sqrt(pi*t)) )  * exp(-.25/t -t*alpha^2)
    #browser()
    if(!n%%2){
      log_a_n=2*log(n+1) +(-(n+1)^2/(4*t) +  1/(4*t))
      #a_n=(n+1)^2 * exp(-(n+1)^2/(4*t) +  1/(4*t))
    }
    if(n%%2){
      log_a_n=log(2*t)+(    -n^2/(4*t)+1/(4*t)    )
      #a_n=2*t*exp(    -n^2/(4*t)+1/(4*t)    )
    }
    log_val_ch_a_n=log_a_n+log_ch_t
  }

  if(t>t_star){
    log_ch_t=(log_sinh(alpha)-log(alpha)) + log(   2*pi^2) +  ( -exp_rate*t)
    #ch_t=( sinh(alpha)/alpha) * (   2*pi^2*exp(-exp_rate*t)    )
    log_a_n=(2*log(n+1)) + (   -(   (n+1)^2-1   )*pi^2*t  )
    #a_n=((n+1)^2) * exp(   -(   (n+1)^2-1   )*pi^2*t  )
    log_val_ch_a_n=log_ch_t+log_a_n
  }

  return(log_val_ch_a_n)

}


##########################################################################################
##########################################################################################

#' @export
sample_truncated_Exp<-function(n, alpha){
  t_star=0.1005082; exp_rate=pi^2+alpha^2
  #val=t_star+rexp(n=n, rate =pi^2+alpha^2 )
  val=rTruncatedExp(n = n, rate = pi^2+alpha^2, lowerLim =t_star, upperLim = Inf )
  return(val)
}

##########################################################################################
##########################################################################################
    sample_truncated_gig<-function( alpha){
      t_star=0.1005082;  flag=0
      mode_dist=(-2.5+sqrt(2.5^2+ alpha^2))/(2*alpha^2)

    ##############################################################
    ##############################################################
          #if(ub>=.4){
          if(mode_dist<t_star){
                iter=1
                while((flag==0)*(iter<10)){
                            iter=iter+1
                            X=rgig(n = 1,lambda = -1.5, chi = 1/2, psi = 2*alpha^2)
                            if(X<=t_star){flag=1; val=X} #inf if
                }#End while
          }#End if
    ##############################################################
    ##############################################################
      #browser()
            if(flag==0){ # includes the case ub<.3
                    ub=pgig(t_star, lambda = -1.5,chi = 1/2, psi = 2*alpha^2)
                    U = runif(1)*ub
                    val= qgig(p = U,lambda = -1.5, chi = 1/2, psi = 2*alpha^2  )
            }#endif

      return(val)
    }

##########################################################################################
##########################################################################################


sample_truncated_gig_0<-function( alpha){
  # this is an alternative function to sample_truncated_gig; not as efficient as sample_truncated_gig
  t_star=0.1005082;
  #browser()

        ub=pgig(t_star, lambda = -1.5,chi = 1/2, psi = 2*alpha^2)
        U = runif(1)*ub
        val= qgig(p = U,lambda = -1.5, chi = 1/2, psi = 2*alpha^2  )
  return(val)
}

##########################################################################################
##########################################################################################


  mix_Exp_Gig_sampler<-function(alpha){
        gig_selector=sampling_gig_segment_selector(alpha)
            if(gig_selector){    val=sample_truncated_gig( alpha)      }
            if(!gig_selector){   val=sample_truncated_Exp(n=1, alpha)  }
        return(val)
  }

###############################################################################################################
###############################################################################################################

#' Generates a sample from the truncated Exponential random variable
#'  rTruncatedExp(n=1, rate=1,lowerLim=0, upperLim=Inf, eps=1e-100 )
#'
#' @param n number of samples to be generated.
#' @param rate rate of the Exponential Distribution.
#' @param lowerLim Lower limit of the support. A number greater than 0.
#' @param upperLim UpperLimit of the Support. A number greater then \code{lowerLim}
#' @param eps=1e-100 Incase a discrete sampler is to be deployed, \code{eps} pertains to the bin size.
#' @return Samples from the truncated Exponential distribution
#' @examples
#' rTruncatedExp(n=1, rate=2)
#' rTruncatedExp(n=1, rate=200, lowerLim=20, upperLim=21)
#' @export
 rTruncatedExp<-function(n=1, rate=1,lowerLim=0, upperLim=Inf, eps=1e-100 ){
         a=lowerLim;b=upperLim
         #(   exp(-rate*a)-exp(-rate*x)  )/(   exp(-rate*a)-exp(-rate*b)    )=u
         #-log(exp(-rate*a)-log(u)*Norm_const)/rate=(x)
         #exp(-rate*a)-log(u)*Norm_const=exp(-rate*x)
         #browser()
         u=runif(n)
         Norm_const= exp(-rate*a)-exp(-rate*b)
         if(is.na(Norm_const)){browser()}
         if(Norm_const<eps){return(discreteTruckExp(rate,a,b, prec_len=100))}
         x=-log(exp(-rate*a)-(u)*Norm_const)/rate
   return(x)
 }



 discreteTruckExp<-function(rate,lowerLim,UpperLim, prec_len=100){
         seq_rng=seq(lowerLim, UpperLim, length.out = prec_len)
         log_d_val= -seq_rng*rate
         d_val=exp(log_d_val-max(log_d_val))
         sampled_bin=sample(x = seq_rng, size = 1, prob =d_val )

         val=runif(1,sampled_bin, sampled_bin+ (UpperLim-lowerLim)/prec_len   )
         if(is.na(val)){browser()}
   return(val)
 }
 ###############################################################################################################
 ###############################################################################################################



 AcceptanceRate_sampler_nu_half_1<-function(x){
   AcceptanceRate_sampler_nu_half<-function(alpha){
     t_star=0.1005082
     #nu=1.5
     #browser()
     Const=log_sinh(alpha)-log(alpha)
     exp_rate=pi^2+alpha^2
     log_p_e=log(2*pi^2)- exp_rate*t_star  -log(exp_rate)
     log_p_g= 1.5*log(2*alpha)-.5*log(pi) + log_Bessel_K(alpha, nu=1.5)+
     log(pgig(t_star,lambda=-1.5, chi=.5, psi=2*alpha^2))
     p_g=exp( log_p_g +Const); p_e=exp( log_p_e+Const );
     #print(c(p_g, p_e, p_e+p_g))
     #p_gig= p_g/(p_e+p_g)
     #print(c(log_p_e,log_p_g,p_gig))
     #val=rbinom(n=1, size=1,prob = p_gig )
     #browser()
     val=p_g+p_e
     return(val)
   }
   return(apply(matrix(x, ncol=1), MARGIN = 1, AcceptanceRate_sampler_nu_half))
 }



#' Plots the density and the proposal density
#' @examples
#' plot_GIG_MPG_Proposal_h_densities(3)
#' plot_GIG_MPG_Proposal_h_densities(5)
#' @export
 plot_GIG_MPG_Proposal_h_densities<-function(alpha=5 , UpperLim=NULL, LowerLim=NULL){
   nu=.5



   f_mpg<-function(x){
     if(length(x)>1){
       return(apply(matrix(x, ncol=1), MARGIN = 1, function(y){log_f_nu(t=y, a=alpha, nu=nu, j_nu_0 = j_0, J_nuPlus1 = J_1, iflog = FALSE)}))
     }
     return(log_f_nu(t=x, a=alpha, nu=nu, j_nu_0 = j_0, J_nuPlus1 = J_1, iflog = FALSE))
   }


   f_gig<-function(x){
     val=dgig(x =x, lambda =-(nu+1), chi = .5, psi = 2*alpha^2 )
     return(val)
   }

   h_fun_single<-function(t){
     t_star=0.1005082
     Const=log_sinh(alpha)-log(alpha)
     if( t<0.1005082){val= Const-1/(4*t)-alpha^2*t -log(2*t^2*sqrt(pi*t))}
     if( t>=0.1005082){ val= Const+log(2*pi^2)-(pi^2+alpha^2)*t   }
     return(val)
   }

   log_h_fun<-function(t){
     return(apply(matrix(t, ncol=1), MARGIN = 1, FUN = h_fun_single) )
   }

   #browser()

   #nu=.5
   #alpha=5;
   j_0=bessel_zero_Jnu(nu = nu, s = 1:1000)
   J_1=besselJ(j_0,nu=1+nu);
   if(is.null(UpperLim)){
     UpperLim=.05+.75/alpha
   }
   if(is.null(LowerLim)){LowerLim=.006}
   #UpperLim=.25*(nu==0)+ .15*(nu==.5)+ .1*(nu>.5)*(nu<4)+ .06*(nu>4)
   library(ggplot2)
   nu_name=nu
   if(nu==.5){nu_name="half"}
   #browser()
   print(alpha)



   #
   #
   #
   # #nu=1.5
   # #browser()
   # Const=log_sinh(alpha)-log(alpha)
   # exp_rate=pi^2+alpha^2
   # log_p_e=log(2*pi^2)- exp_rate*t_star  -log(exp_rate)
   # log_p_g= 1.5*log(2*alpha)-.5*log(pi) + log_Bessel_K(alpha, nu=1.5)+
   #   log(pgig(t_star,lambda=-1.5, chi=.5, psi=2*alpha^2))
   # p_g=exp( log_p_g +Const); p_e=exp( log_p_e+Const );
   #
   #
   #




   p=ggplot(data.frame(x = rnorm(100)), aes(x)) +
     stat_function(fun = function(x){f_mpg(x = x )}, colour = "white", size=.3)+
     stat_function(fun = function(x){exp(log_h_fun(t = x))}, colour = "black", size=.3)+
     stat_function(fun = function(x){f_gig(x = x)}, colour = "red", size=.3)+
     xlim(LowerLim, UpperLim)+xlab("Support of the Distributions")+ylab("Density")

   p=p+theme(
     panel.background = element_rect(fill = "gray",
                                     colour = "gray",
                                     size = 0.05, linetype = "solid"),
     panel.grid.major = element_blank(),
     #element_line(size = 0.001, linetype = 'solid',
     #colour = "#f2f2f2"),
     panel.grid.minor = element_blank()
     #element_line(size = 0.001, linetype = 'solid',
     #colour = "gray")
   )
   p
   #Fileloc="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/Codes/VonFCodes1.5/densityComparison plots/"
   #FileName=paste0(Fileloc,"plot_GIG_MPG_densities_nu_", nu_name, "alpha", alpha , ".pdf")
   #ggsave(file=FileName, height = 3,width=4)
   #dev.off()
   return(p)
 }


