#' A function to compute the biweight mean vector and covariance matrix
#'
#' Compute a multivariate location and scale estimate based on Tukey's biweight weight function.
#'
#' @param x an \code{n x g} matrix or data frame (\code{n} is the number of measurements, \code{g} is the number of observations (genes) )
#' @param r breakdown (\code{k/n} where \code{k} is the largest number of observations that can be replaced with arbitrarily large values while keeping the estimates bounded)
#' @param median a logical command to determine whether the initialization is done using the coordinate-wise median and MAD (TRUE) or using the minimum covariance determinant (MCD)  (FALSE).  Using the MCD is substantially slower.
#' @param full_init a logical command to determine whether the initialization is done for each pair separately (FALSE) or only one time at the beginning using the entire data matrix (TRUE).  Initializing for each pair separately is substantially slower.
#' @param absval a logical command to determine whether the distance should be measured as 1 minus the absolute value of the correlation (TRUE) or simply 1 minus the correlation (FALSE)
#'
#' @importFrom MASS mvrnorm
#' @importFrom robustbase covMcd
#' @importFrom hopach dissmatrix
#' @importFrom stats dchisq mad mahalanobis pchisq
#'
#' @return Using \code{\link{biwt_est}} to estimate the robust covariance matrix, a robust measure of correlation is computed using Tukey's biweight M-estimator. The biweight correlation is essentially a weighted correlation where the weights are calculated based on the distance of each measurement to the data center with respect to the shape of the data.  The correlations are computed pair-by-pair because the weights should depend only on the pairwise relationship at hand and not the relationship between all the observations globally.  The biwt functions compute many pairwise correlations and create distance matrices for use in other algorithms (e.g., clustering).
#'
#' In order for the biweight estimates to converge, a reasonable initialization must be given. Typically, using TRUE for the median and full_init arguments will provide acceptable initializations. With particularly irregular data, the MCD should be used to give the initial estimate of center and shape. With data sets in which the observations are orders of magnitudes different, full_init=FALSE should be specified.
#'
#' Returns a list with components:
#' \item{biwt_dist_matrix}{a matrix of the biweight distances (default is 1 minus absolute value of the biweight correlation).}
#' \item{biwt_NAid_matrix}{a matrix representing whether the biweight correlation was possible to compute (will be NA if too much data is missing or if the initializations are not accurate).  0 if computed accurately, 1 if NA.}
#'
#' @references Hardin, J., Mitani, A., Hicks, L., VanKoten, B.; A Robust Measure of Correlation Between Two Genes on a Microarray, \emph{BMC Bioinformatics, 8}:220; 2007.
#'
#' @examples
#' samp_data <- MASS::mvrnorm(30,mu=c(0,0,0),Sigma=matrix(c(1,.75,-.75,.75,1,-.75,-.75,-.75,1),ncol=3))
#' r <- 0.2 # breakdown
#'
#' # To compute the 3 pairwise distances in matrix form:
#' samp_bw_dist_mat <- biwt_dist_matrix(samp_data, r)
#' samp_bw_dist_mat
#'
#' # To convert the distances into an element of class 'dist'
#' as.dist(samp_bw_dist_mat$biwt_dist_mat)
#' @export
biwt_dist_matrix <- function(x, r, median=TRUE, full_init=TRUE, absval=TRUE){

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

	if (median!=TRUE) {med.init<-robustbase::covMcd(cbind(x[,i],x[,j]))}
	else		{med.init<-list()
			med.init$cov <- diag(1,2)*apply(cbind(x[,i],x[,j]),2,stats::mad,na.rm=TRUE)
			med.init$center <- apply(cbind(x[,i],x[,j]),2,median,na.rm=TRUE)}
	}

	biwt <- biwt_est(cbind(x[,i],x[,j]), r, med.init)
	corr <- c(corr,biwt$biwt_sig[1,2]/sqrt(biwt$biwt_sig[1,1]*biwt$biwt_sig[2,2]))
	NAid <- c(NAid,biwt$biwt_NAid)
	j<-j+1
	}

	}

if(absval==TRUE){dist.mat <- hopach::dissmatrix(1 - abs(corr))}
else {dist.mat <- hopach::dissmatrix(1 - corr)}

diag(dist.mat) <- 0

NAid.mat <- hopach::dissmatrix(NAid)

return(list(biwt_dist_mat = dist.mat, biwt_NAid_mat=NAid.mat))}




