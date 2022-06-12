#' Internal function
#'
#' @param x a number
#' @param c1 a cutoff
#'
#' @return a number
#' @export
#'
#' @examples
#' vbw(2,3)
`vbw` <- function(x,c1){
  return(psibw(x,c1)*x)
}

