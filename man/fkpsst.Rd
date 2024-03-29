% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fkpsst.R
\name{fkpsst}
\alias{fkpsst}
\title{Functional KPSS Test}
\usage{
fkpsst(
  X,
  ALPHA = 0.05,
  METHOD = "MC",
  K = default_kernel,
  h_power = 2/5,
  est_ev = nrow(X),
  MCNsim = 10000,
  sbb_set = NULL,
  M = 1000,
  b = ceiling((2 * nrow(X))^(1/3))
)
}
\arguments{
\item{X}{The functional time series being tested, inputted in a matrix form with each row representing each observation of the functional data values on equidistant points of any prespecified interval.}

\item{ALPHA}{Significance level of the test. The default value is 5\%.}

\item{METHOD}{Method to determine the critical value of the test statistic. The default value 'METHOD="MC"' represents the Monte Carlo method. 'METHOD="SB"' represents the simple bootstrap method and 'METHOD="MBB"' represents the moving block bootstrap method.}

\item{K}{Kernel function in the estimation of the long-run covariance function, which is only effective in the Monte Carlo method. The default function is 'default_kernel' function in this package.}

\item{h_power}{Power of sample size 'N' (valued in (0,1)) for the smoothing bandwidth, which is only effective in the Monte Carlo method. The default value is 2/5.}

\item{est_ev}{Number of the largest eigenvalues chosen to estimate the limiting distribution, which is only effective in the Monte Carlo method. The default value is the sample size 'N'.}

\item{MCNsim}{Number of Monte Carlo datasets generated in the Monte Carlo method, which is only effective in Monte Carlo method. The default value is 10000.}

\item{sbb_set}{A vector of independent simulated data generated from the function 'dataset_sbb', which is only effective and essential in Monte Carlo method.}

\item{M}{Number of bootstrap datasets generated in the bootstrap method, which is only effective in bootstrap methods. The default value is 1000.}

\item{b}{Block length used in the moving block bootstrap method, which is only effective in the moving block bootstrap method. The default value is ceiling((2N)^(1/3)), where 'N' is the sample size.}
}
\value{
The result of the test is presented with the value of test statistic and its p-value under the null hypothesis of trend stationarity.
}
\description{
KPSS test for functional time series with different methods on determining the critical values of the test statistic. The Monte Carlo method (for a large sample size) was constructed in Kokoszka and Young (2016) and the bootstrap-based methods (both the simple bootstrap method and moving block bootstrap method, for a small/moderate sample size) were constructed in Chen and Pun (2019).
}
\examples{
N<-100
EX<-matrix(rep(0,N*100),ncol=100)
EXI<-rep(0,100)
set.seed(1)
EXI[1]<-rnorm(1,0,1)
for (j in 2:100) {EXI[j]<-EXI[j-1]+rnorm(1,0,1)}
for (i in 1:N) {
temp<-rnorm(100,0,1)
EX[i,1]<-temp[1]
for (j in 2:100) {
EX[i,j]<-EX[i,j-1]+temp[j]
}
}
fkpsst(X=EX,METHOD="SB")
}
\references{
Chen, Y., & Pun, C. S. (2019). A bootstrap-based KPSS test for functional time series. Journal of Multivariate Analysis, 174, 104535.

Kokoszka, P., & Young, G. (2016). KPSS test for functional time series. Statistics, 50(5), 957-973.
}
