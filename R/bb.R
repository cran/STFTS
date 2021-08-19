#' Integrals of Squared Brownian Bridge
#'
#' Generate a dataset of independent simulated values of \eqn{\int_0^1{B^2(t)}dt}, where \eqn{B} is a standard Brownian bridge on [0,1].
#'
#' @param NUM Number of simulated values generated in the dataset.
#' @return A vector with length equals to NUM.
#' @importFrom e1071 rwiener
#' @export
#' @examples
#' dataset_bb(10)

dataset_bb<-function(NUM){
  bb<-rep(0,NUM)
  for (i in 1:NUM){
    W<-rwiener()
    V<-W-(1:1000)*W[1000]/1000
    bb[i]<-sum(V^2)/1000
  }
  return(bb)
}
