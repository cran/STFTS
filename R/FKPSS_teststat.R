#' Test Statistic in Functional KPSS Test
#'
#' Calculate test statistic R_N for functional KPSS test, which Was constructed in Kokoszka and Young (2016) and Chen and Pun (2019).
#'
#' @param X The functional time series being tested, inputted in a matrix form with each row representing each observation of the functional data values on equidistant points of any prespecified interval.
#' @return The value of test statistic R_N for functional KPSS test.
#' @export
#' @references Chen, Y., & Pun, C. S. (2019). A bootstrap-based KPSS test for functional time series. Journal of Multivariate Analysis, 174, 104535.
#' @references Kokoszka, P., & Young, G. (2016). KPSS test for functional time series. Statistics, 50(5), 957-973.
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
#' FST_teststat_TN(EX)

FKPSS_teststat<-function(X){

  N<-nrow(X)

  MX<-apply(X,2,sum)/N
  IX<-apply((1:N)*X,2,sum)
  ee<-X+(((1:N-(N+1)/2)*6/(N-1)-1))%*%t(MX)-((1:N-(N+1)/2)*12/(N*(N^2-1)))%*%t(IX)

  ZN<-partial_sum(ee)
  TN<-mean(ZN^2)

  return(TN)

}
