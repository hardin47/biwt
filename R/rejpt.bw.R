#' Internal function
#'
#' @param p a number
#' @param r the breakdown
#'
#' @return the asymptotic rejection point
#'
#' @examples
#' rejpt.bw(2,3)
#' @export
rejpt.bw <- function(p, r) {
  c1 <- 2 * p
  iter <- 1
  crit <- 100
  eps <- 1e-5
  while ((crit > eps) & (iter < 100)) {
    c1.old <- c1
    fc <- erho.bw(p, c1) - c1 ^ 2 * r / 6
    fcp <- erho.bw.p(p, c1) - c1 * r / 3
    c1 <- c1 - fc / fcp
    if (c1 < 0)
      c1 <- c1.old / 2
    crit <- abs(fc)
    iter <- iter + 1
  }
  return(c(c1, pchisq(c1 ^ 2, p), log10(1 - stats::pchisq(c1 ^ 2, p))))
}
