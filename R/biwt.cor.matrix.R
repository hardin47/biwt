`biwt.cor.matrix` <-
function(x,r,median=T,full.init=T){

library(hopach)
library(robustbase)

if (full.init==T){

	if(median!=T){med.init <- robustbase::covMcd(x)}
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

	if (median!=T) {med.init<-robustbase::covMcd(cbind(x[,i],x[,j]))}
	else		{med.init<-list()
			med.init$cov <- diag(1,2)*apply(cbind(x[,i],x[,j]),2,mad,na.rm=T)
			med.init$center <- apply(cbind(x[,i],x[,j]),2,median,na.rm=T)}
	}

	biwt <- biwt.est(cbind(x[,i],x[,j]),r,med.init)
	corr <- c(corr,biwt$biwt.sig[1,2]/sqrt(biwt$biwt.sig[1,1]*biwt$biwt.sig[2,2]))
	NAid <- c(NAid,biwt$biwt.NAid)
	j<-j+1
	}

	}

corr.mat <- hopach::dissmatrix(corr)
diag(corr.mat) <- 1

NAid.mat <- hopach::dissmatrix(NAid)

return(list(biwt.corr.mat = corr.mat, biwt.NAid.mat=NAid.mat))}




