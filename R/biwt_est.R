#' A function to compute the biweight mean vector and covariance matrix
#'
#' Compute a multivariate location and scale estimate based on Tukey's biweight weight function.
#'
#' @param x an \code{n x 2} matrix or data frame (\code{n} is the number of observations)
#' @param r breakdown (\code{k/n} where \code{k} is the largest number of observations that can be replaced with arbitrarily large values while keeping the estimates bounded)
#' @param med.init a (robust) initial estimate of the center and shape of the data.  form is a list with components center and cov.
#'
#' @return A robust measure of center and shape is computed using Tukey's biweight M-estimator.  The biweight estimates are essentially weighted means and covariances where the weights are calculated based on the distance of each measurement to the data center with respect to the shape of the data.  The estimates should be computed pair-by-pair because the weights should depend only on the pairwise relationship at hand and not the relationship between all the observations globally.
#'
#' Returns a list with components:
#' \item{biwt_mu}{the final estimate of location}
#' \item{biwt_sig}{the final estimate of scatter}
#' \item{biwt_NAid}{a logical of whether the given initialization was used (coded as 0) or whether a more precise initialization was used (namely, the pair by pair MCD, coded as 1).}
#' @references Hardin, J., Mitani, A., Hicks, L., VanKoten, B.; A Robust Measure of Correlation Between Two Genes on a Microarray, \emph{BMC Bioinformatics, 8}:220; 2007.
#'
#' @importFrom MASS mvrnorm
#' @importFrom robustbase covMcd
#' @importFrom stats dchisq mad mahalanobis pchisq
#'
#' @examples
#' samp_data <- MASS::mvrnorm(30,mu=c(0,0),Sigma=matrix(c(1,.75,.75,1),ncol=2))
#' r <- 0.2 # breakdown
#'
#' samp_mcd <- robustbase::covMcd(samp_data)
#' samp_bw <- biwt_est(samp_data, r, samp_mcd)
#'
#' samp_bw_var1 <- samp_bw$biwt_sig[1,1]
#' samp_bw_var2 <- samp_bw$biwt_sig[2,2]
#' samp_bw_cov <- samp_bw$biwt_sig[1,2]
#'
#' samp_bw_corr <- samp_bw_cov / sqrt(samp_bw_var1 * samp_bw_var2)
#'
#' # or:
#'
#' samp_bw_corr <- samp_bw$biwt_sig[1,2] / sqrt(samp_bw$biwt_sig[1,1]*samp_bw$biwt_sig[2,2])
#' samp_bw_corr
#'
#' ##############
#' # to speed up the calculations, use the median/mad for the initialization:
#' ##############
#'
#' samp_init <- list()
#' samp_init$cov <- diag(apply(samp_data, 2, stats::mad, na.rm=TRUE))
#' samp_init$center <- apply(samp_data, 2, median, na.rm=TRUE)
#'
#' samp_bw <- biwt_est(samp_data, r, samp_init)
#' samp_bw_corr <- samp_bw$biwt.sig[1,2] / sqrt(samp_bw$biwt.sig[1,1]*samp_bw$biwt.sig[2,2])
#' @export
biwt_est <- function(x, r, med.init){

if(!is.matrix(x)) x <- as.matrix(x)
p<-2
n <- dim(x)[1]
c1<-rejpt.bw(p=2,r)[1]
b0<-erho.bw(p=2,c1)[1]

NAid<-FALSE
d <- sqrt(stats::mahalanobis(x,med.init$center,med.init$cov))
k <- ksolve(d,p,c1,b0)

if(is.na(k)) {
NAid<-TRUE
med.init <- robustbase::covMcd(x)
d <- sqrt(stats::mahalanobis(x,med.init$center,med.init$cov))
k <- ksolve(d,p,c1,b0)}
# MCD is a more robust estimate of the center/shape
# than the median which is sometimes used


eps <- 1e-5
crit <- 100
iter <- 1
while (crit > eps & iter < 100) {
d <- d/k
biwt.mu <- apply(wtbw(d,c1)*x,2,sum,na.rm=TRUE) / sum(wtbw(d,c1),na.rm=TRUE)
cent <- array(dim=c(n,p,p))
for (i in 1:n){
cent[i,,] <- (x[i,] - biwt.mu)%*%t(x[i,]-biwt.mu)}
biwt.sig <- apply(cent*wtbw(d,c1),c(2,3),sum,na.rm=TRUE)/
sum(vbw(d,c1),na.rm=TRUE)


d2 <- sqrt(stats::mahalanobis(x,biwt.mu,biwt.sig))
k <- ksolve(d2,p,c1,b0)
crit <- max(abs(d-(d2/k)),na.rm=TRUE)
d <- d2
iter <-  iter+1}

return(list(biwt_mu=biwt.mu,biwt_sig=biwt.sig, biwt_NAid=as.numeric(NAid)))}

