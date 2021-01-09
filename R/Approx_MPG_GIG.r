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




log_f_mpg<-function(x, alpha,nu,  j_0, J_1){
  if(length(x)>1){
    return(apply(matrix(x, ncol=1), MARGIN = 1, function(y){log_f_nu(t=y, a=alpha, nu=nu, j_nu_0 = j_0, J_nuPlus1 = J_1, iflog = TRUE)}))
  }
  return(log_f_nu(t=x, a=alpha, nu=nu,j_nu_0 = j_0, J_nuPlus1 = J_1, iflog = TRUE))

}


log_f_gig<-function(x, alpha, nu){
  val=dgig(x =x, lambda =-(nu+1), chi = .5, psi = 2*alpha^2 , log = TRUE)
  return(val)
}




f_diff<-function(x, alpha, j_0, J_1 ){
  #ff1=f_mpg(x,alpha )
  #val=(ff1-f_gig(x,alpha))
  return(exp(log_f_mpg(x, alpha,j_0=j_0, J_1=J_1))-exp(log_f_gig(x=x, alpha=alpha)))
}





###################
#' @examples
#' x=seq(from = .01,to=.15,by = .0001)
#' alpha=seq(from=10, to = 30, by = 1)
#' pp=Plot_Approximation_MPG_GIG_Relative_diff(x=seq(from = .01,to=.15,by = .001),alpha= seq(from=10, to = 20, by = 1))
#' pp$p
#' pp$p1
#' @export
Plot_Approximation_MPG_GIG_Relative_diff<-function(x, alpha, nu=0){


#nu=0#### change nu for different results
j_0=bessel_zero_Jnu(nu = nu, s = 1:1000)
J_1=besselJ(j_0,nu=1+nu);


Alpha_all=alpha %x% replicate(length(x),1)
x_all=replicate(length(alpha),1)%x% x

data_mat=cbind(Alpha_all,x_all)
#diff_den1=apply(data_mat, MARGIN = 1, FUN = function(x_vec){ f_diff(x_vec[2], x_vec[1]) })
den1=apply(data_mat, MARGIN = 1, FUN = function(x_vec){ exp(log_f_mpg(x = x_vec[2],alpha =  x_vec[1],nu=nu, j_0 = j_0, J_1 = J_1)) })
den2=apply(data_mat, MARGIN = 1, FUN = function(x_vec){ exp(log_f_gig(x_vec[2], x_vec[1], nu) ) })

df=data.frame(x=x_all, alpha=Alpha_all, den1=den1, den2=den2, diff_den=den1-den2)

df$scaled_diff=df$diff_den*0
for(index in alpha){
  mx=max(df$den1[df$alpha==index])
  df$scaled_diff[df$alpha==index]=df$diff_den[df$alpha==index]/mx
}


p=ggplot(data=df, aes(x=x, y=scaled_diff, group=alpha)) +
  #geom_line(aes(color=alpha), size=.5)+scale_color_gradient(low="skyblue", high="#9999CC")+
  #geom_line(aes(color=alpha), size=.5)+scale_color_gradient(low="white", high="black")+
  geom_line(aes(color=alpha), size=.2)+scale_color_gradient(low="black", high="white")+
  ylim(-.015,.015)+
  xlab("Support of the Distributions")+
  ylab("Relative Difference Between the Densities ") +
  labs(color=expression(alpha))

p+theme(
  panel.background = element_rect(fill = "gray",
                                  colour = "gray",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.002, linetype = 'solid',
                                  colour = "white"),
  panel.grid.minor = element_line(size = 0.001, linetype = 'solid',
                                  colour = "lightgray")
)
#dev.off()
p


p1=ggplot(data=df, aes(x=x, y=diff_den, group=alpha)) +
  geom_line(aes(color=alpha))+scale_color_gradient(low="skyblue", high="#9999CC")+
  ylim(-1,1)+
  xlab("Support of the Distributions")+
  ylab("Difference Between the Densities") +
  labs(color=expression(alpha))

p1
return(list(p=p,p1=p1))

}











#' @export
plot_GIG_MPG_densities<-function(alpha=5, nu=0 , UpperLim=NULL, LowerLim=NULL){
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


  p=ggplot(data.frame(x = rnorm(100)), aes(x)) +
    stat_function(fun = function(x){f_mpg(x = x )}, colour = "white", size=.3)+
    stat_function(fun = function(x){f_gig(x = x)}, colour = "black", size=.2)+
    xlim(LowerLim, UpperLim)+xlab("Support of the Distributions")+ylab("Density")

  p=p+theme(
    panel.background = element_rect(fill = "gray",
                                    colour = "gray",
                                    size = 0.05, linetype = "solid"),
    panel.grid.major = element_line(size = 0.02, linetype = 'solid',
                                    colour = "#f2f2f2"),
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                    colour = "#f2f2f2")
  )
  p
  #Fileloc="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/Codes/VonFCodes1.5/densityComparison plots/"
  #FileName=paste0(Fileloc,"plot_GIG_MPG_densities_nu_", nu_name, "alpha", alpha , ".pdf")
  #ggsave(file=FileName, height = 3,width=4)
  #dev.off()
  return(p)
}

