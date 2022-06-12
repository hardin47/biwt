#' Internal function
#'
#' @param x a number
#' @param c1 a cutoff
#'
#' @return a number
#' @export
#'
#' @examples
#' psibw(2,3)
`psibw` <- function(x,c1){
ivec <- (abs(x)>c1)
    return((1-ivec)*(x*(1-(x/c1)^2)^2))}

