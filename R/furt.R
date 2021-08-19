#' Functional Unit Root Test
#'
#' Unit root test for functional time series with different methods on determining the critical values of the test statistic. The Monte Carlo method was constructed in Chen and Pun (2021), while the bootstrap-based methods have not been validated in the literature (although such an option is provided, please use them at your own risk).
#'
#' @param X The functional time series being tested, inputted in a matrix form with each row representing each observation of the functional data values on equidistant points of any prespecified interval.
#' @param ALPHA Significance level of the test. The default value is 5\%.
#' @param TYPE Type of hypothesis test. The default value 'TYPE="R"' represents the right-tailed test, which is used when the alternative hypothesis is trend stationarity model. 'TYPE="L"' represents the left-tailed test, which is used when the alternative hypothesis is simple stationarity model or AR(1) model.
#' @param METHOD Method to determine the critical value of the test statistic. The default value 'METHOD="MC"' represents the Monte Carlo method. 'METHOD="SB"' represents the simple bootstrap method and 'METHOD="MBB"' represents the moving block bootstrap method.
#' @param K Kernel function in the estimation of the long-run covariance function, which is only effective in the Monte Carlo method. The default function is 'default_kernel' function in this package.
#' @param h_power Power of sample size 'N' (valued in (0,1)) for the smoothing bandwidth, which is only effective in the Monte Carlo method. The default value is 2/5.
#' @param est_ev Number of the largest eigenvalues chosen to estimate the limiting distribution, which is only effective in the Monte Carlo method. The default value is the sample size 'N'.
#' @param MCNsim Number of Monte Carlo datasets generated in the Monte Carlo method, which is only effective in Monte Carlo method. The default value is 10000.
#' @param bm_set A vector of independent simulated data generated from the function 'dataset_bm', which is only effective and essential in Monte Carlo method.
#' @param M Number of bootstrap datasets generated in the bootstrap method, which is only effective in bootstrap methods. The default value is 1000.
#' @param b Block length used in the moving block bootstrap method, which is only effective in the moving block bootstrap method. The default value is ceiling((2N)^(1/3)), where 'N' is the sample size.
#' @return The result of the test is presented with the value of test statistic and its p-value under the null hypothesis of functional random walk.
#' @export
#' @references Chen, Y., & Pun, C. S. (2021). Functional Unit Root Test. Available at SSRN.
#' @examples
#' N<-100
#' EE<-matrix(rep(0,N*100),ncol=100)
#' EX<-matrix(rep(0,N*100),ncol=100)
#' set.seed(1)
#' for (i in 1:N) {
#' temp<-rnorm(100,0,1)
#' EE[i,1]<-temp[1]
#' for (j in 2:100) {
#' EE[i,j]<-EE[i,j-1]+temp[j]
#' }
#' }
#' EX[1,]<-EE[1,]
#' for (i in 2:N) {EX[i,]<-EX[i-1]+EE[i,]}
#' furt(X=EX,METHOD="SB")


furt<-function(X, ALPHA=0.05, TYPE="R", METHOD="MC",
               K=default_kernel, h_power=2/5, est_ev=nrow(X)-1, MCNsim=10000, bm_set=NULL,
               M=1000, b=ceiling((2*(nrow(X)-1))^(1/3))){
  if (h_power<=0 | h_power>=1) stop("The value of 'h_power' exceeds the range (0,1).")
  if (METHOD=="MC" & is.null(bm_set)) stop("You must generate a simluated dataset for Monte Carlo method from the function 'generate_bm'.")

  N<-nrow(X)
  n_col<-ncol(X)

  hat_eta<-X[2:N,]-X[1:(N-1),]

  ZN<-partial_sum(hat_eta)
  TN<-mean(ZN^2)

  if (METHOD=="MC"){

  lambda<-est_eigenvalue(hat_eta,K,h_power,est_ev)

  intV<-matrix(sample(bm_set,MCNsim*est_ev,replace=T),ncol=est_ev)

  monte=(intV%*%lambda)[,1]

  if (TYPE=="R") {

    furt_result<-list(statistic=c(R_N=TN),p.value=mean(TN<monte),sample.size=N,alternative="trend stationary",method="Monte-Carlo-based functional unit root test")

  } else {

    furt_result<-list(statistic=c(R_N=TN),p.value=mean(TN>monte),sample.size=N,alternative="local stationary & AR(1) stationary",method="Monte-Carlo-based functional unit root test")

  }

  } else {

    if (METHOD=="SB"){
      bs_teststat<-sapply(1:M,function(o){
        bs_ind<-sample(1:(N-1),N-1,replace=T)
        bs_hat_eta<-hat_eta[bs_ind,]

        bsZN<-partial_sum(bs_hat_eta)
        bsTN<-mean(bsZN^2)

        return(bsTN)

      }
      )
      if (TYPE=="R") {

        furt_result<-list(statistic=c(R_N=TN),p.value=mean(TN<bs_teststat),sample.size=N,alternative="trend stationary",method="Simple-bootstrap-based functional unit root test")

      } else {

        furt_result<-list(statistic=c(R_N=TN),p.value=mean(TN>bs_teststat),sample.size=N,alternative="local stationary & AR(1) stationary",method="Simple-bootstrap-based functional unit root test")

      }

    } else {
      bbs_teststat<-sapply(1:M,function(o){

        bs<-sample(1:(N-b),ceiling((N-1)/b),replace=T)
        bs_ind<-bs[1]:(bs[1]+b-1)
        for (i in 2:ceiling((N-1)/b)){
          bs_ind<-c(bs_ind,bs[i]:(bs[i]+b-1))
        }
        bs_hat_eta<-hat_eta[bs_ind[1:(N-1)],]

        bsZN<-partial_sum(bs_hat_eta)
        bsTN<-mean(bsZN^2)

        return(bsTN)

      }
      )
      if (TYPE=="R") {

        furt_result<-list(statistic=c(R_N=TN),p.value=mean(TN<bbs_teststat),sample.size=N,alternative="trend stationary",method="Moving-block-bootstrap-based functional unit root test")

      } else {

        furt_result<-list(statistic=c(R_N=TN),p.value=mean(TN>bbs_teststat),sample.size=N,alternative="local stationary & AR(1) stationary",method="Moving-block-bootstrap-based functional unit root test")

      }

    }


  }

  class(furt_result)<-"htest"
  return(furt_result)

}

