#' Internal function
#'
#' @param x a vector
#' @param c1 a cutoff
#'
#' @return a vector
#'
#' @examples
#' psibw(rnorm(10),3)
#' @export
psibw <- function(x, c1) {
  ivec <- (abs(x) > c1)
  (1 - ivec) * (x * (1 - (x / c1) ^ 2) ^ 2)
}

