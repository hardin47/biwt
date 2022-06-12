#' Internal function
#'
#' @param p a number
#' @param a a number
#' @param c1 a number
#'
#' @return a number
#'
#' @examples
#' chi.int(2,3,4)
#' @export
chi.int2 <- function(p,a,c1) {
return( exp(lgamma((p+a)/2)-lgamma(p/2))*2^{a/2}*(1-stats::pchisq(c1^2,p+a)) )
}
