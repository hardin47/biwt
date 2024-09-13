#' Internal function
#'
#' @param v a vector
#'
#' @return a vector
#'
#' @examples
#' vect2diss_hop(rnorm(10))
#' @export
vect2diss_hop <- function (v){
 # This is exactly the dissmatrix function in the hopach package
 # I have re-typed it to avoid dependencies of the hopach package.
    if (!is.vector(v))
      stop("arg to dissmatrix() must be a vector")
    p <- (1 + sqrt(1 + 8 * length(v)))/2
    M <- matrix(0, nrow = p, ncol = p)
    count <- 1
    for (i in 1:(p - 1)) {
      M[i, (i + 1):p] <- v[count:(count + p - i - 1)]
      count <- count + p - i
    }
    return(M + t(M))
  }
