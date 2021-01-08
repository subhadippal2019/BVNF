
#' @export
rExtendedGamma<-function(alpha, a, b){
  # alpha>0, a>0 and b = any real

  if(1-(a>0)*(alpha>0)){  print("The parameters a and alpha required to be positive."); return(NULL) }


  if(b<0){
    return(rExtendedGamma_negative_b(alpha, a, abs(b)))
  }
  if(b>0){
   return( rExtendedGamma_positive_b_alpha(alpha, a, b))
  }

  if(b==0){
    shape_par=alpha/2; rate_par=a;
    Exp_g<-  rgamma(1,shape=shape_par,rate=rate_par)
    g=sqrt(Exp_g)
    return(g)
  }

}

# alpha=1;a=400;b=+100
# samlpe_g= rExtendedGamma(alpha, a,b)
#
# for(iii in 1:100000){
#   samlpe_g[iii]=rExtendedGamma(alpha,a,b)
# }
#
# par(mfrow=c(1,2))
#
# plot(density(samlpe_g))
# plot(function(x){x^(alpha-1)*exp(-a*x^2+b*x)},xlim=c(0,.4))
#
