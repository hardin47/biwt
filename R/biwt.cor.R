#' A function to compute the biweight mean vector and covariance matrix
#'
#' @param x an \code{n x g} matrix or data frame (\code{n} is the number of measurements, \code{g} is the number of observations (genes) )
#' @param r breakdown (\code{k/n} where \code{k} is the largest number of observations that can be replaced with arbitrarily large values while keeping the estimates bounded)
#' @param median a logical command to determine whether the initialization is done using the coordinate-wise median and MAD (TRUE) or using the minimum covariance determinant (MCD)  (FALSE).  Using the MCD is substantially slower.
#' @param full.init a logical command to determine whether the initialization is done for each pair separately (FALSE) or only one time at the beginning using the entire data matrix (TRUE).  Initializing for each pair separately is substantially slower.
#'
#' @importFrom MASS mvrnorm
#' @importFrom robustbase covMcd
#'
#' @return returns a list consisting of:
#' \item{biwt.corr}{a vector consisting of the lower triangle of the correlation matrix stored by columns in a vector, say bwcor. If \code{g} is the number of observations, i.e., then for \eqn{i < j \leq g}, the biweight correlation between (rows) \code{i} and \code{j} is bwcor[\code{g*(i-1) - i*(i-1)/2 + j-i}]. The length of the vector is \code{g*(g-1)/2}, i.e., of order \code{g^2}. }
#' \item{biwt.NAid}{a vector which is indexed in the same way as \code{biwt.corr}.  The entries represent whether the biweight correlation was possible to compute (will be NA if too much data is missing or if the initializations are not accurate).  0 if computed accurately, 1 if NA.}
#'
#' @examples
#'
#' samp.data <- MASS::mvrnorm(30,mu=c(0,0,0),Sigma=matrix(c(1,.75,-.75,.75,1,-.75,-.75,-.75,1),ncol=3))
#' r<-0.2 # breakdown
#'
#' # To compute the 3 pairwise correlations from the sample data:
#'
#' samp.bw.cor <- biwt.cor(samp.data,r)
#' samp.bw.cor
#' @export
biwt.cor <- function(x,r,median=T,full.init=T){

if (full.init==T){

	if(median!=T){med.init <- robustbase::covMcd(x)}
	else	{med.init<-list()
		med.init$cov<-diag(1,2)*median(apply(x,1,stats::mad,na.rm=T))
		med.init$center<-c(1,1)*median(apply(x,1,median,na.rm=T))
		}
	}

	corr <- c()
	NAid <- c()
	g <- dim(x)[2]


for(i in 1:g){
	j <- 1
	while(j < i){

if (full.init !=T){

	if (median!=T) {med.init <- robustbase::covMcd(cbind(x[,i],x[,j]))}
	else		{med.init<-list()
			med.init$cov <- diag(1,2)*apply(cbind(x[,i],x[,j]),2,stats::mad,na.rm=T)
			med.init$center <- apply(cbind(x[,i],x[,j]),2,median,na.rm=T)}
	}

	biwt <- biwt.est(cbind(x[,i],x[,j]),r,med.init)
	corr <- c(corr,biwt$biwt.sig[1,2]/sqrt(biwt$biwt.sig[1,1]*biwt$biwt.sig[2,2]))
	NAid <- c(NAid,biwt$NAid)
	j<-j+1
	}

	}
return(list(biwt.corr=corr,biwt.NAid=NAid))}