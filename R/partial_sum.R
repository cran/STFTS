#' Partial Sum
#'
#' Calculate the partial sum matrix of a functional time series.
#'
#' @param X The functional time series being calculated, inputted in a matrix form with each row representing each observation of the functional data values on equidistant points of any prespecified interval.
#' @return A square matrix of partial sums of the functional time series. The column (row) dimension of the output matrix is equal to the column dimension of the input matrix.
#' @export
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
#' partial_sum(EX)

partial_sum<-function(X){

  X_N<-nrow(X)
  X_ncol<-ncol(X)
  temp<-matrix(rep(0,X_ncol*X_ncol),ncol=X_ncol)

  if (floor(X_N/X_ncol)>0) {
    for (i in 1:floor(X_N/X_ncol)){
      temp[1,]<-temp[1,]+X[i,]/sqrt(X_N)
    }
  }
  if (X_ncol>1){
  for (i in 2:X_ncol){
    temp[i,]<-temp[i-1,]
    if (floor(X_N*(i-1)/X_ncol)<floor(X_N*i/X_ncol)) {
      for (j in (floor(X_N*(i-1)/X_ncol)+1):floor(X_N*i/X_ncol)){
        temp[i,]<-temp[i,]+X[j,]/sqrt(X_N)
      }
    }
  }
  }
  return(temp)

}
