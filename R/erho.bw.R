#' Internal function
#'
#' @param p a number
#' @param c1 a cutoff
#'
#' @return a number, The expected value of rho
#' @export
#'
#' @examples
#' erho.bw(2,3)
`erho.bw` <- function(p,c1) {
  # gives b0 = E(rho)
return(chi.int(p,2,c1)/2-chi.int(p,4,c1)/(2*c1^2)+
     chi.int(p,6,c1)/(6*c1^4)+c1^2*chi.int2(p,0,c1)/6)
}

