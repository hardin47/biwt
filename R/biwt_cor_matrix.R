#' A function to compute the biweight mean vector and covariance matrix
#'
#' Compute a multivariate location and scale estimate based on Tukey's biweight weight function.
#'
#' @param x an \code{n x g} matrix or data frame (\code{n} is the number of measurements, \code{g} is the number of observations (genes) )
#' @param r breakdown (\code{k/n} where \code{k} is the largest number of observations that can be replaced with arbitrarily large values while keeping the estimates bounded)
#' @param median a logical command to determine whether the initialization is done using the coordinate-wise median and MAD (TRUE) or using the minimum covariance determinant (MCD)  (FALSE).  Using the MCD is substantially slower.
#' @param full.init a logical command to determine whether the initialization is done for each pair separately (FALSE) or only one time at the beginning using the entire data matrix (TRUE).  Initializing for each pair separately is substantially slower.
#'
#' @return returns a list consisting of:
#' \item{biwt.corr.matrix}{a matrix of the biweight correlations.}
#' \item{biwt.NAid.matrix}{a matrix representing whether the biweight correlation was possible to compute (will be NA if too much data is missing or if the initializations are not accurate).  0 if computed accurately, 1 if NA.}
#' @references Hardin, J., Mitani, A., Hicks, L., VanKoten, B.; A Robust Measure of Correlation Between Two Genes on a Microarray, \emph{BMC Bioinformatics, 8}:220; 2007.
#'
#' @importFrom robustbase covMcd
#' @importFrom MASS mvrnorm
#' @importFrom hopach dissmatrix
#' @importFrom stats dchisq mad mahalanobis pchisq
#'
#' @examples
#'
#' samp_data <- MASS::mvrnorm(30,mu=c(0,0,0),Sigma=matrix(c(1,.75,-.75,.75,1,-.75,-.75,-.75,1),ncol=3))
#' r <- 0.2 # breakdown
#'
#' # To compute the 3 pairwise correlations in matrix form:
#'
#' samp_bw_cor_mat <- biwt_cor_matrix(samp_data,r)
#' samp_bw_cor_mat
#' @export
biwt_cor_matrix <- function(x, r, median=TRUE, full.init=TRUE){


if (full.init==TRUE){

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

if (full.init !=TRUE){

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

corr.mat <- hopach::dissmatrix(corr)
diag(corr.mat) <- 1

NAid.mat <- hopach::dissmatrix(NAid)

return(list(biwt_corr_mat = corr.mat, biwt_NAid_mat=NAid.mat))}




