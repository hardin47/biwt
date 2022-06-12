#' Title
#'
#' @param p a number
#' @param a a number
#' @param c1 a number
#'
#' @return a number
#'
#' @examples
#' chi.int.p(1,2,3)
#' @export
chi.int.p <- function(p,a,c1){
return( exp(lgamma((p+a)/2)-lgamma(p/2))*2^{a/2}*stats::dchisq(c1^2,p+a)*2*c1 )
}
