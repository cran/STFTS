#' Test Statistic M_N in Functional Stationarity Test
#'
#' Calculate test statistic M_N for functional stationarity test, which was constructed in Horvath et al. (2014).
#'
#' @param X The functional time series being tested, inputted in a matrix form with each row representing each observation of the functional data values on equidistant points of any prespecified interval.
#' @return The value of test statistic M_N calculated for functional stationarity test.
#' @export
#' @references Horvath, L., Kokoszka, P., & Rice, G. (2014). Testing stationarity of functional time series. Journal of Econometrics, 179(1), 66-82.
#' @examples
#' N<-100
#' EX<-matrix(rep(0,N*100),ncol=100)
#' set.seed(1)
#' for (i in 1:N) {
#' temp<-rnorm(100,0,1)
#' EX[i,1]<-temp[1]
#' for (j in 2:100) {
#' EX[i,j]<-EX[i,j-1]+temp[j]
#' }
#' }
#' FST_teststat_MN(EX)

FST_teststat_MN<-function(X){
  n_col<-ncol(X)

  SN<-partial_sum(X)

  ZN<-SN-(1:n_col)%*%t(SN[n_col,])/n_col

  z2<-ZN-rep(1,n_col)%*%t(apply(ZN,2,sum)/n_col)

  MN<-mean(z2^2)
  return(MN)
}
