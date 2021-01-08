
#' calculates log(K(x, nu)) where K is the modified Bessel function of the Second kind.
#'
#' @param x A complex number.
#' @param nu degree of the Bessel function
#' @return The log(K(\code{x},\code{nu}))
#' @examples
#' log_Bessel_K(10,.5)
#' log_Bessel_K(10+3*1i,0)
#'  log_Bessel_K(10000,0)
#' @export
log_Bessel_K<-function(x, nu){ return( log(BesselK(x,nu =nu, expon.scaled = TRUE  ))-x ) }



#Norm<-function(x){  sqrt(sum(x^2)) }



KAPPA_INITIAL<-function(Y){
  Y_bar=apply(Y,2,'mean');n=dim(Y)[1]
  nu=length(Y_bar)/2-1;

  const_a= Norm(Y_bar);
  #kappa_lower_1= (nu+.5)*const_a/(1-const_a)
  kappa_lower= 2*(nu+.5)*const_a/(1-const_a^2) # Segura 2011 bound
  #kappa_upper=(nu+1)*const_a/(1-const_a)
  return(kappa_lower)
}







####################################################################
####################################################################
####################################################################
####################################################################

#' calculates log(sinh(x)) for small and large arguments.
#'
#' @param x A real number, a vector of real numbers.
#' @return The log(sinh(\code{x}))
#' @examples
#' log_sinh(10)
#' log_sinh(c(1000,1,10,10000))
#' @export
log_sinh<-function(x){
  log_sinh_single<-function(y){
    if(y<200){return(log(sinh(y)))}
    if(y>=200){return(y-log(2))}
  }
  val=apply(X = matrix(x, ncol=1),MARGIN = 1,FUN = log_sinh_single)
  return(val)
}





