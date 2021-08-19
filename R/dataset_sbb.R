#' Integrals of Squared Second-level Brownian Bridge
#'
#' Generate a dataset of independent simulated values of \eqn{\int_0^1{[W(t)+(2x-3x^2)W(1)+(-6x+6x^2)\int_0^1{W(x)dx}]^2}dt}, where \eqn{W} is a standard Brownian motion on [0,1].
#'
#' @param NUM Number of simulated values generated in the dataset.
#' @return A vector with length equals to NUM.
#' @importFrom e1071 rwiener
#' @export
#' @examples
#' dataset_sbb(10)

dataset_sbb<-function(NUM){
  sbb<-rep(0,NUM)
  for (i in 1:NUM){
    W<-rwiener()
    V2<-W+(2*(1:1000)/1000-3*((1:1000)^2)/(1000^2))*W[1000]+(6*((1:1000)^2/(1000^2))-6*(1:1000)/1000)*sum(W)/1000
    sbb[i]<-sum(V2^2)/1000
  }
  return(sbb)
}

