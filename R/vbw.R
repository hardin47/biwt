#' Internal function
#'
#' @param x a vector
#' @param c1 a cutoff
#'
#' @return a vector
#'
#' @examples
#' vbw(rnorm(10),3)
#' @export
vbw <- function(x,c1) {
  psibw(x,c1)*x
}

