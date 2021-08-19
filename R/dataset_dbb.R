#' Integrals of Squared De-meaned Brownian Bridge
#'
#' Generate a dataset of independent simulated values of \eqn{\int_0^1{[B(t)-\int_0^1{B(x)dx}]^2}dt}, where \eqn{B} is a standard Brownian bridge on [0,1].
#'
#' @param NUM Number of simulated values generated in the dataset.
#' @return A vector with length equals to NUM.
#' @importFrom e1071 rwiener
#' @export
#' @examples
#' dataset_dbb(10)

dataset_dbb<-function(NUM){
  dbb<-rep(0,NUM)
  for (i in 1:NUM){
    W<-rwiener()
    V<-W-(1:1000)*W[1000]/1000
    V1<-V-rep(sum(V)/1000,1000)
    dbb[i]<-sum(V1^2)/1000
  }
  return(dbb)
}
