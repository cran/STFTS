#' Default Kernel Function
#'
#' Default kernel function used in calculating the estimated long-run covariance function.
#'
#' @param x The variable.
#' @return The kernel value.
#' @export
#' @examples
#' default_kernel(1)
#' default_kernel(0)

default_kernel<-function(x){
  result<-0
  if (x>=0 & x<0.1) result<-1
  else {if (x>=0.1 & x<1.1) result<-1.1-abs(x)}
  return(result)
}
