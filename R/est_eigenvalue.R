#' Eigenvalues of An Estimated Long-run Covariance Function
#'
#' Calculate the eigenvalues of an estimated long-run covariance function.
#'
#' @param este Estimated errors in the long-run covariance function, inputted in a matrix form with each row representing each observation of the functional error values on equidistant points of any prespecified interval.
#' @param K Kernel function in the estimation of the long-run covariance function.
#' @param h_power Power of sample size 'N' (valued in (0,1)) for the smoothing bandwidth in the kernel function.
#' @param estev Number of the largest eigenvalues chosen in the output.
#' @return A vector of first 'estev' largest eigenvalues in descending order.
#' @export
#' @examples
#' N<-100
#' EE<-matrix(rep(0,N*100),ncol=100)
#' set.seed(1)
#' for (i in 1:N) {
#' temp<-rnorm(100,0,1)
#' EE[i,1]<-temp[1]
#' for (j in 2:100) {
#' EE[i,j]<-EE[i,j-1]+temp[j]
#' }
#' }
#' est_eigenvalue(este=EE,K=default_kernel,h_power=1/2,estev=10)


est_eigenvalue<-function(este,K,h_power,estev){
  if (h_power<=0 | h_power>=1) stop("The value of 'h_power' exceeds the range (0,1).")

  estN<-nrow(este)
  estcol<-ncol(este)

  gamma<-rep(0,estN*estcol*estcol)
  dim(gamma)<-c(estN,estcol,estcol)
  for (i in 1:estN){
    for (j in i:estN){
      gamma[i,,]<-gamma[i,,]+(este[j,]%*%t(este[j-i+1,]))/estN
    }
  }

  cc<-matrix(rep(0,estcol*estcol),ncol=estcol)
  cc<-gamma[1,,]
  for (i in 2:estN){
    cc<-cc+K((i-1)/(estN^(h_power)))*(gamma[i,,]+t(gamma[i,,]))
  }

  eigens<-eigen(cc/estcol)$values[1:estev]

  return(eigens)

}
