
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
##################################
#alpha=10
