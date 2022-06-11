#' A function to compute the biweight mean vector and covariance matrix
#'
#' @param x an $n \times g$ matrix or data frame ($n$ is the number of measurements, $g$ is the number of observations (genes) )
#' @param r breakdown ($k/n$ where $k$ is the largest number of observations that can be replaced with arbitrarily large values while keeping the estimates bounded
#' @param median a logical command to determine whether the initialization is done using the coordinate-wise median and MAD (TRUE) or using the minimum covariance determinant (MCD)  (FALSE).  Using the MCD is substantially slower.
#' @param full.init  a logical command to determine whether the initialization is done for each pair separately (FALSE) or only one time at the beginning using the entire data matrix (TRUE).  Initializing for each pair separately is substantially slower.
#'
#' @return
#' @export
#'
#' @examples
`biwt.cor` <-
function(x,r,median=T,full.init=T){


if (full.init==T){

	if(median!=T){med.init <- covMcd(x)}
	else	{med.init<-list()
		med.init$cov<-diag(1,2)*median(apply(x,1,mad,na.rm=T))
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

	if (median!=T) {med.init<-covMcd(cbind(x[,i],x[,j]))}
	else		{med.init<-list()
			med.init$cov <- diag(1,2)*apply(cbind(x[,i],x[,j]),2,mad,na.rm=T)
			med.init$center <- apply(cbind(x[,i],x[,j]),2,median,na.rm=T)}
	}

	biwt <- biwt.est(cbind(x[,i],x[,j]),r,med.init)
	corr <- c(corr,biwt$biwt.sig[1,2]/sqrt(biwt$biwt.sig[1,1]*biwt$biwt.sig[2,2]))
	NAid <- c(NAid,biwt$NAid)
	j<-j+1
	}

	}
return(list(biwt.corr=corr,biwt.NAid=NAid))}
