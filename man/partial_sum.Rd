% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partial_sum.R
\name{partial_sum}
\alias{partial_sum}
\title{Partial Sum}
\usage{
partial_sum(X)
}
\arguments{
\item{X}{The functional time series being calculated, inputted in a matrix form with each row representing each observation of the functional data values on equidistant points of any prespecified interval.}
}
\value{
A square matrix of partial sums of the functional time series. The column (row) dimension of the output matrix is equal to the column dimension of the input matrix.
}
\description{
Calculate the partial sum matrix of a functional time series.
}
\examples{
N<-100
EX<-matrix(rep(0,N*100),ncol=100)
set.seed(1)
for (i in 1:N) {
temp<-rnorm(100,0,1)
EX[i,1]<-temp[1]
for (j in 2:100) {
EX[i,j]<-EX[i,j-1]+temp[j]
}
}
partial_sum(EX)
}
