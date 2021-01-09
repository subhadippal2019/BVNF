
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



norm_vec<-function(x){  sqrt(sum(x^2)) }


#' Calculates initial value of KAPPA_INITIAL(Y). Approximate posterior mode.
#'
#' @param Y:  n\times d  matrix containing n directional data of dimension d.
#' norm of each row of the matrix Y is 1.
#' @return Approximate posterior mode for the concentration parameter.
#' @examples
#' library(Rfast)
#' Y=rvmf(n =10, mu=c(1,0,0),k = 10)
#' KAPPA_INITIAL(Y)
#' @export
KAPPA_INITIAL<-function(Y){
  Y_bar=apply(Y,2,'mean');n=dim(Y)[1]
  nu=length(Y_bar)/2-1;

  const_a= norm_vec(Y_bar);
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




inner_prod<-function(x,y){(sum(x*y))}
#### Data Generation #####
alt_besselI<-function(x, nu){
  if(nu==0){
    val=gsl::bessel_I0(x)
  }
  if(nu>0){
    val=besselI(10,nu = nu)
  }
  return(val)
}




