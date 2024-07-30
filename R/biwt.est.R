#' A function to compute the biweight mean vector and covariance matrix
#'
#' Compute a multivariate location and scale estimate based on Tukey's biweight weight function.
#'
#' @param x a \code{2 x n} matrix or data frame (\code{n} is the number of measurements)
#' @param r breakdown (\code{k/n} where \code{k} is the largest number of measurements that can be replaced with arbitrarily large values while keeping the estimates bounded). Default is \code{r = 0.2}.
#' @param med.init a (robust) initial estimate of the center and shape of the data.  The format is a list with components center and cov (as in the output of \code{covMcd()} from the rrcov package). Default is the minimum covariance determinant (MCD) on the data.
#'
#' @return A list with components
#' \item{biwt.mu}{the final estimate of center}
#' \item{biwt.sig}{the final estimate of shape}
#' @references Hardin, J., Mitani, A., Hicks, L., VanKoten, B.; A Robust Measure of Correlation Between Two Genes on a Microarray, \emph{BMC Bioinformatics, 8}:220; 2007.
#'
#' @importFrom MASS mvrnorm
#' @importFrom robustbase covMcd
#' @importFrom stats dchisq mad mahalanobis pchisq
#'
#' @examples
#' samp.data <- t(MASS::mvrnorm(30,mu=c(0,0),
#'                Sigma=matrix(c(1,.75,.75,1),ncol=2)))
#'
#' samp.bw <- biwt.est(samp.data)
#' samp.bw
#'
#' samp.bw.var1 <- samp.bw$biwt.sig[1,1]
#' samp.bw.var2 <- samp.bw$biwt.sig[2,2]
#' samp.bw.cov <- samp.bw$biwt.sig[1,2]
#'
#' samp.bw.cor <- samp.bw.cov / sqrt(samp.bw.var1 * samp.bw.var2)
#' samp.bw.cor
#'
#' # or:
#'
#' samp.bw.cor <- samp.bw$biwt.sig[1,2] / sqrt(samp.bw$biwt.sig[1,1]*samp.bw$biwt.sig[2,2])
#' samp.bw.cor
#'
#' ##############
#' # to speed up the calculations, use the median/mad for the initialization:
#' ##############
#'
#' samp.init <- list()
#' samp.init$cov <- diag(apply(samp.data, 1, stats::mad, na.rm=TRUE))
#' samp.init$center <- apply(samp.data, 1, median, na.rm=TRUE)
#' samp.init
#'
#' samp.bw <- biwt.est(samp.data, med.init = samp.init)
#' samp.bw.cor <- samp.bw$biwt.sig[1,2] / sqrt(samp.bw$biwt.sig[1,1]*samp.bw$biwt.sig[2,2])
#' samp.bw.cor
#' @export
biwt.est <- function (x, r = 0.2, med.init = robustbase::covMcd(x))
{
  p <- 2
  n <- dim(x)[2]
  c1 <- rejpt.bw(p = 2, r)[1]
  b0 <- erho.bw(p = 2, c1)[1]
  x <- t(x)
  d <- sqrt(mahalanobis(x, med.init$center, med.init$cov))
  k <- ksolve(d, p, c1, b0)
  if (is.na(k)) {
    med.init <- covMcd(x)
    d <- sqrt(mahalanobis(x, med.init$center, med.init$cov))
    k <- ksolve(d, p, c1, b0)
  }
  eps <- 1e-05
  crit <- 100
  iter <- 1
  while (crit > eps & iter < 100) {
    d <- d/k
    biwt.mu <- apply(wtbw(d, c1) * x, 2, sum, na.rm = TRUE)/
      sum(wtbw(d, c1), na.rm = TRUE)
    cent <- array(dim = c(n, p, p))
    for (i in 1:n) {
      cent[i, , ] <- (x[i, ] - biwt.mu) %*% t(x[i, ] - biwt.mu)
    }
    biwt.sig <- apply(cent * wtbw(d, c1), c(2, 3), sum, na.rm = TRUE)/
      sum(vbw(d, c1), na.rm = TRUE)
    d2 <- sqrt(mahalanobis(x, biwt.mu, biwt.sig))
    k <- ksolve(d2, p, c1, b0)
    crit <- max(abs(d - (d2/k)), na.rm = TRUE)
    d <- d2
    iter <- iter + 1
  }
  return(list(biwt.mu = biwt.mu, biwt.sig = biwt.sig))
}
