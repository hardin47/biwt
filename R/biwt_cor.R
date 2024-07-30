#' A function to compute the biweight mean vector and covariance matrix
#'
#' @param x an \code{n x g} matrix or data frame (\code{n} is the number of measurements, \code{g} is the number of observations (genes) )
#' @param r breakdown (\code{k/n} where \code{k} is the largest number of observations that can be replaced with arbitrarily large values while keeping the estimates bounded)
#' @param median a logical command to determine whether the initialization is done using the coordinate-wise median and MAD (TRUE) or using the minimum covariance determinant (MCD)  (FALSE).  Using the MCD is substantially slower.
#' @param full_init a logical command to determine whether the initialization is done for each pair separately (FALSE) or only one time at the beginning using the entire data matrix (TRUE).  Initializing for each pair separately is substantially slower.
#'
#' @importFrom MASS mvrnorm
#' @importFrom robustbase covMcd
#' @importFrom stats dchisq mad mahalanobis pchisq
#'
#' @return Using \code{\link{biwt_est}} to estimate the robust covariance matrix, a robust measure of correlation is computed using Tukey's biweight M-estimator. The biweight correlation is essentially a weighted correlation where the weights are calculated based on the distance of each measurement to the data center with respect to the shape of the data.  The correlations are computed pair-by-pair because the weights should depend only on the pairwise relationship at hand and not the relationship between all the observations globally.  The biwt functions compute many pairwise correlations and create distance matrices for use in other algorithms (e.g., clustering).
#'
#' In order for the biweight estimates to converge, a reasonable initialization must be given. Typically, using TRUE for the median and full_init arguments will provide acceptable initializations. With particularly irregular data, the MCD should be used to give the initial estimate of center and shape. With data sets in which the observations are orders of magnitudes different, full_init=FALSE should be specified.
#'
#' Returns a list with components:
#' \item{biwt_corr}{a vector consisting of the lower triangle of the correlation matrix stored by columns in a vector, say bwcor. If \code{g} is the number of observations, i.e., then for \eqn{i < j \leq g}, the biweight correlation between (rows) \code{i} and \code{j} is bwcor[\eqn{g*(i-1) - i*(i-1)/2 + j-i}]. The length of the vector is \eqn{g*(g-1)/2}, i.e., of order \eqn{g^2}. }
#' \item{biwt_NAid}{a vector which is indexed in the same way as \code{biwt_corr}.  The entries represent whether the biweight correlation was possible to compute (will be NA if too much data is missing or if the initializations are not accurate).  0 if computed accurately, 1 if NA.}
#'
#' @examples
#'
#' samp_data <- MASS::mvrnorm(30,mu=c(0,0,0),Sigma=matrix(c(1,.75,-.75,.75,1,-.75,-.75,-.75,1),ncol=3))
#' r <- 0.2 # breakdown
#'
#' # To compute the 3 pairwise correlations from the sample data:
#'
#' samp_bw_cor <- biwt_cor(samp_data,r)
#' samp_bw_cor
#' @export
biwt_cor <- function(x, r, median=TRUE, full_init=TRUE){

if(!is.matrix(x)) x <- as.matrix(x)

if (full_init==TRUE){

	if(median!=TRUE){med.init <- robustbase::covMcd(x)}
	else	{med.init<-list()
		med.init$cov<-diag(1,2)*median(apply(x,1,stats::mad,na.rm=TRUE))
		med.init$center<-c(1,1)*median(apply(x,1,median,na.rm=TRUE))
		}
	}

	corr <- c()
	NAid <- c()
	g <- dim(x)[2]


for(i in 1:g){
	j <- 1
	while(j < i){

if (full_init !=TRUE){

	if (median!=TRUE) {med.init <- robustbase::covMcd(cbind(x[,i],x[,j]))}
	else		{med.init<-list()
			med.init$cov <- diag(1,2)*apply(cbind(x[,i],x[,j]),2,stats::mad,na.rm=TRUE)
			med.init$center <- apply(cbind(x[,i],x[,j]),2,median,na.rm=TRUE)}
	}

	biwt <- biwt_est(cbind(x[,i],x[,j]),r,med.init)
	corr <- c(corr,biwt$biwt_sig[1,2]/sqrt(biwt$biwt_sig[1,1]*biwt$biwt_sig[2,2]))
	NAid <- c(NAid,biwt$NAid)
	j<-j+1
	}

	}
return(list(biwt_corr=corr,biwt_NAid=NAid))}
