#' Integrals of Squared Brownian Motion
#'
#' Generate a dataset of independent simulated values of \eqn{\int_0^1{W^2(t)}dt}, where \eqn{W} is a standard Brownian motion on [0,1].
#'
#' @param NUM Number of simulated values generated in the dataset.
#' @return A vector with length equals to NUM.
#' @importFrom e1071 rwiener
#' @export
#' @examples
#' dataset_bm(10)

dataset_bm<-function(NUM){
  bm<-rep(0,NUM)
  for (i in 1:NUM){
    W<-rwiener()
    bm[i]<-sum(W^2)/1000
  }
  return(bm)
}

