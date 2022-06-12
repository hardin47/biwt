#' Internal function
#'
#' @param d a vector
#' @param p a number
#' @param c1 a number
#' @param b0 a number
#'
#' @return a number
#'
#' @examples
#' ksolve(rnorm(20,.1,2),1,3,1)
#' @export
ksolve <- function(d,p,c1,b0){
    k <- 1
    iter <- 1
    crit <- 100
    eps <- 1e-10
    while ((crit > eps)&(iter<100)){
    k.old <- k
        fk <- mean(rhobw(d/k,c1),na.rm=T)-b0
        fkp <- -mean(psibw(d/k,c1)*d/k^2,na.rm=T)
if (fkp==0) {k<-NA
return(k)
stop("no values close enough")}
        k <- k - fk/fkp
        if (k < 0)  k <- k.old/2
        crit <- abs(k-k.old)
        iter <- iter+1    }
    return(k) }

