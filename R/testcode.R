
samp.data <- MASS::mvrnorm(30,mu=c(0,0,0),Sigma=matrix(c(1,.75,-.75,.75,1,-.75,-.75,-.75,1),ncol=3))

r<-0.2 # breakdown

biwt.cor(samp.data[,1:2],r=.2)$corr
biwt.cor(samp.data[,c(1,3)],r=.2)$corr
biwt.cor(samp.data[,c(2,3)],r=.2)$corr


biwt.cor(samp.data,r=.2)

biwt.cor.matrix(samp.data,r=.2)

biwt.dist.matrix(samp.data,r=.2,absval=T)
biwt.dist.matrix(samp.data,r=.2,absval=F)


