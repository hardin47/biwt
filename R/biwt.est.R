`biwt.est` <-
function(x,r,med.init){

p<-2
n <- dim(x)[1]
c1<-rejpt.bw(p=2,r)[1]
b0<-erho.bw(p=2,c1)[1]

NAid<-FALSE
d <- sqrt(mahalanobis(x,med.init$center,med.init$cov))
k <- ksolve(d,p,c1,b0)

if(is.na(k)) {  
NAid<-TRUE
med.init <- covMcd(x)
d <- sqrt(mahalanobis(x,med.init$center,med.init$cov))
k <- ksolve(d,p,c1,b0)}
# MCD is a more robust estimate of the center/shape
# than the median which is sometimes used


eps <- 1e-5
crit <- 100
iter <- 1
while (crit > eps & iter < 100) {
d <- d/k
biwt.mu <- apply(wtbw(d,c1)*x,2,sum,na.rm=T) / sum (wtbw(d,c1),na.rm=T)
cent <- array(dim=c(n,p,p))
for (i in 1:n){
cent[i,,] <- (x[i,] - biwt.mu)%*%t(x[i,]-biwt.mu)}
biwt.sig <- apply(cent*wtbw(d,c1),c(2,3),sum,na.rm=T)/
sum(vbw(d,c1),na.rm=T) 


d2 <- sqrt(mahalanobis(x,biwt.mu,biwt.sig))
k <- ksolve(d2,p,c1,b0)
crit <- max(abs(d-(d2/k)),na.rm=T)
d <- d2
iter <-  iter+1}

return(list(biwt.mu=biwt.mu,biwt.sig=biwt.sig, biwt.NAid=as.numeric(NAid)))}

