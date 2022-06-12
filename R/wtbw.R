#' Internal function
#'
#' @param x a number
#' @param c1 a cutoff
#'
#' @return a number
#' @export
#'
#' @examples
#' wtbw(2,3)
`wtbw` <- function(x,c1){
    ivec <- (abs(x)>c1)
    return((1-ivec)*(1-(x/c1)^2)^2)}

