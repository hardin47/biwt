#' Internal function
#'
#' @param x a vector
#' @param c1 a number
#'
#' @return a vector
#'
#' @examples
#' wtbw(rnorm(10),3)
#' @export
wtbw <- function(x, c1) {
  ivec <- (abs(x) > c1)
  return((1 - ivec) * (1 - (x / c1) ^ 2) ^ 2)
}
