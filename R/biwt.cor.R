#' A function to compute the biweight mean vector and covariance matrix
#'
#' @param x a \code{g x n} matrix or data frame (\code{n} is the number of measurements, \code{g} is the number of observations (genes) )
#' @param r breakdown (\code{k/n} where \code{k} is the largest number of observations that can be replaced with arbitrarily large values while keeping the estimates bounded). Default is \code{r = 0.2}.
#' @param output a character string specifying the output format.  Options are "matrix" (default), "vector", or "distance".  See value below.
#' @param median a logical command to determine whether the initialization is done using the coordinate-wise median and MAD^2 (TRUE, default) or using the minimum covariance determinant (MCD) (FALSE).  Using the MCD is substantially slower. The MAD is the median of the absolute deviations from the median.  See R help file on \code{mad}.
#' @param full.init a logical command to determine whether the initialization is done for each pair separately (FALSE) or only one time at the beginning using the entire data matrix (TRUE, default).  Initializing for each pair separately is substantially slower.
#' @param absval a logical command to determine whether the distance should be measured as 1 minus the absolute value of the correlation (TRUE, default) or as 1 minus the correlation (FALSE).
#'
#' @importFrom MASS mvrnorm
#' @importFrom robustbase covMcd
#' @importFrom stats dchisq mad mahalanobis pchisq
#'
#' @return Specifying "vector" for the output argument returns a vector consisting of the lower triangle of the correlation matrix stored by columns in a vector, say \eqn{bwcor}. If \eqn{g} is the number of observations and \eqn{bwcor} is the correlation vector, then for \eqn{i < j <= g}, the biweight correlation between (rows) \eqn{i} and \eqn{j} is \eqn{bwcor[(j-1)*(j-2)/2 + i]}. The length of the vector is \eqn{g*(g-1)/2}, i.e., of order \eqn{g^2}.
#'
#' Specifying "matrix" for the output argument returns a matrix of the biweight correlations.
#'
#' Specifying "distance" for the output argument returns a matrix of the biweight distances (default is 1 minus absolute value of the biweight correlation).
#'
#' If there is too much missing data or if the initialization is not accurate, the function will compute the MCD for a given pair of observations before computing the biweight correlation (regardless of the initial settings given in the call to the function).
#'
#' The "vector" output option is given so that correlations can be stored as vectors which are less computationally intensive than matrices.
#'
#' Returns a list with components:
#' \item{corr}{a vector consisting of the lower triangle of the correlation matrix stored by columns in a vector, say bwcor. If \code{g} is the number of observations, i.e., then for \eqn{i < j \leq g}, the biweight correlation between (rows) \code{i} and \code{j} is bwcor[\eqn{g*(i-1) - i*(i-1)/2 + j-i}]. The dimension of the matrix is \eqn{g x g}.}
#' \item{corr.mat}{a matrix consisting of the lower triangle of the correlation matrix stored by columns in a vector, say bwcor. If \code{g} is the number of observations, i.e., then for \eqn{i < j \leq g}, the biweight correlation between (rows) \code{i} and \code{j} is bwcor[\eqn{g*(i-1) - i*(i-1)/2 + j-i}]. The length of the vector is \eqn{g*(g-1)/2}, i.e., of order \eqn{g^2}. }
#' \item{dist.mat}{a matrix consisting of the correlations converted to distances (either 1 - correlation or 1 - abs(correlation)).}
#' @examples
#'
#' # note that biwt.cor() takes data that is gxn where the
#' # goal is to find correlations or distances between each of the g items
#'
#' samp.data <- t(MASS::mvrnorm(30,mu=c(0,0,0),
#'                Sigma=matrix(c(1,.75,-.75,.75,1,-.75,-.75,-.75,1),ncol=3)))
#'
#' # To compute the 3 pairwise correlations from the sample data:
#'
#' samp.bw.cor <- biwt.cor(samp.data, output="vector")
#' samp.bw.cor
#'
#' # To compute the 3 pairwise correlations in matrix form:
#'
#' samp.bw.cor.mat <- biwt.cor(samp.data)
#' samp.bw.cor.mat
#'
#' # To compute the 3 pairwise distances in matrix form:
#'
#' samp.bw.dist.mat <- biwt.cor(samp.data, output="distance")
#' samp.bw.dist.mat
#'
#' # To convert the distances into an object of class `dist'
#'
#' as.dist(samp.bw.dist.mat)
#' @export
biwt.cor <- function (x,
                      r = 0.2,
                      output = "matrix",
                      median = TRUE,
                      full.init = TRUE,
                      absval = TRUE)
{
  if (full.init == TRUE) {
    rand.samp <- x[sample(1:nrow(x), 2), ]
    if (median != TRUE) {
      med.init <- robustbase::covMcd(t(rand.samp))
    }
    else {
      med.init <- list()
      med.init$cov <- diag(1, 2) * (apply(rand.samp, 1, mad, na.rm = TRUE)) ^
        2
      med.init$center <- c(1, 1) * apply(rand.samp, 1, median, na.rm = TRUE)
    }
  }
  corr <- c()
  g <- dim(x)[1]
  if (output == "matrix") {
    for (i in 1:g) {
      j <- 1
      while (j < i) {
        if (full.init != TRUE) {
          if (median != TRUE) {
            med.init <- robustbase::covMcd(cbind(x[i, ], x[j, ]))
          }
          else {
            med.init <- list()
            med.init$cov <- diag(1, 2) * (apply(cbind(x[i, ], x[j, ]), 2, mad, na.rm = TRUE)) ^
              2
            med.init$center <- apply(cbind(x[i, ], x[j, ]), 2, median, na.rm = TRUE)
          }
        }
        biwt <- biwt.est(rbind(x[i, ], x[j, ]), r, med.init)
        corr <- c(corr,
                  biwt$biwt.sig[1, 2] / sqrt(biwt$biwt.sig[1, 1] * biwt$biwt.sig[2, 2]))
        j <- j + 1
      }
    }
    corr.mat <- vect2diss(corr)
    diag(corr.mat) <- 1
    return(corr.mat)
  }
  if (output == "distance") {
    for (i in 1:g) {
      j <- 1
      while (j < i) {
        if (full.init != TRUE) {
          if (median != TRUE) {
            med.init <- robustbase::covMcd(cbind(x[i, ], x[j, ]))
          }
          else {
            med.init <- list()
            med.init$cov <- diag(1, 2) * (apply(cbind(x[i, ], x[j, ]), 2, mad, na.rm = TRUE)) ^
              2
            med.init$center <- apply(cbind(x[i, ], x[j, ]), 2, median, na.rm = TRUE)
          }
        }
        biwt <- biwt.est(rbind(x[i, ], x[j, ]), r, med.init)
        corr <- c(corr,
                  biwt$biwt.sig[1, 2] / sqrt(biwt$biwt.sig[1, 1] * biwt$biwt.sig[2, 2]))
        j <- j + 1
      }
    }
    if (absval == TRUE) {
      dist.mat <- vect2diss(1 - abs(corr))
    }
    else {
      dist.mat <- vect2diss(1 - corr)
    }
    diag(dist.mat) <- 0
    return(dist.mat)
  }
  if (output == "vector") {
    for (i in 1:g) {
      j <- 1
      while (j < i) {
        if (full.init != TRUE) {
          if (median != TRUE) {
            med.init <- robustbase::covMcd(cbind(x[i, ], x[j, ]))
          }
          else {
            med.init <- list()
            med.init$cov <- diag(1, 2) * (apply(cbind(x[i, ], x[j, ]), 2, mad, na.rm = TRUE)) ^
              2
            med.init$center <- apply(cbind(x[i, ], x[j, ]), 2, median, na.rm = TRUE)
          }
        }
        biwt <- biwt.est(rbind(x[i, ], x[j, ]), r, med.init)
        corr <- c(corr,
                  biwt$biwt.sig[1, 2] / sqrt(biwt$biwt.sig[1, 1] * biwt$biwt.sig[2, 2]))
        j <- j + 1
      }
    }
    return(corr)
  }
}

