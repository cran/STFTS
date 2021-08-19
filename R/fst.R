#' Functional Stationarity Test
#'
#' Stationarity test for functional time series with different methods on determining the critical values of the test statistic. The Monte Carlo method was constructed in Horvath et al. (2014), while the bootstrap-based methods have not been validated in the literature (although such an option is provided, please use them at your own risk).
#'
#' @param X The functional time series being tested, inputted in a matrix form with each row representing each observation of the functional data values on equidistant points of any prespecified interval.
#' @param STAT Test statistic used in the stationarity test. The default value 'STAT=1' represents 'T_N' and 'STAT=2' represents 'M_N' in original paper.
#' @param ALPHA Significance level of the test. The default value is 5\%.
#' @param METHOD Method to determine the critical value of the test statistic. The default value 'METHOD="MC"' represents the Monte Carlo method. 'METHOD="SB"' represents the simple bootstrap method and 'METHOD="MBB"' represents the moving block bootstrap method.
#' @param K Kernel function in the estimation of the long-run covariance function, which is only effective in the Monte Carlo method. The default function is 'default_kernel' function in this package.
#' @param h_power Power of sample size 'N' (valued in (0,1)) for the smoothing bandwidth, which is only effective in the Monte Carlo method. The default value is 2/5.
#' @param est_ev Number of the largest eigenvalues chosen to estimate the limiting distribution, which is only effective in the Monte Carlo method. The default value is the sample size 'N'.
#' @param MCNsim Number of Monte Carlo datasets generated in the Monte Carlo method, which is only effective in Monte Carlo method. The default value is 10000.
#' @param bb_set A vector of independent simulated data generated from the function 'dataset_bb', which is only effective and essential in Monte Carlo method with 'STAT=1'.
#' @param dbb_set A vector of independent simulated data generated from the function 'dataset_dbb', which is only effective and essential in Monte Carlo method with 'STAT=2'.
#' @param M Number of bootstrap datasets generated in the bootstrap method, which is only effective in bootstrap methods. The default value is 1000.
#' @param b Block length used in the moving block bootstrap method, which is only effective in the moving block bootstrap method. The default value is ceiling((2N)^(1/3)), where 'N' is the sample size.
#' @return The result of the test is presented with the value of test statistic and its p-value under the null hypothesis of local stationarity.
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
#' fst(X=EX,METHOD="SB")


fst<-function(X, STAT=1, ALPHA=0.05, METHOD="MC",
              K=default_kernel, h_power=2/5, est_ev=nrow(X), MCNsim=10000, bb_set=NULL, dbb_set=NULL,
              M=1000, b=ceiling((2*nrow(X))^(1/3))){
  if (h_power<=0 | h_power>=1) stop("The value of 'h_power' exceeds the range (0,1).")
  if (METHOD=="MC" & STAT==1 & is.null(bb_set)) stop("You must generate a simluated dataset for Monte Carlo method from the function 'generate_bb'.")
  if (METHOD=="MC" & STAT==2 & is.null(dbb_set)) stop("You must generate a simluated dataset for Monte Carlo method from the function 'generate_dbb'.")

  N<-nrow(X)

  MeanX<-apply(X,2,sum)/N

  hat_eta<-X-rep(1,N)%*%t(MeanX)

  if (METHOD=="MC"){

  lambda<-est_eigenvalue(hat_eta,K,h_power,est_ev)

  if (STAT==1) {
    TN<-FST_teststat_TN(X)

    intV<-matrix(sample(bb_set,MCNsim*est_ev,replace=T),ncol=est_ev)
    monte=(intV%*%lambda)[,1]

    fst_result<-list(statistic=c(T_N=TN),p.value=mean(TN<monte),sample.size=N,method="Monte-Carlo-based functional stationarity test")

  }
  if (STAT==2) {
    MN<-FST_teststat_MN(X)

    intV<-matrix(sample(dbb_set,MCNsim*est_ev,replace=T),ncol=est_ev)
    monte=(intV%*%lambda)[,1]

    fst_result<-list(statistic=c(M_N=MN),p.value=mean(MN<monte),sample.size=N,method="Monte-Carlo-based functional stationarity test")

  }

  } else {

    if (METHOD=="SB"){
      if (STAT==1) {
      TN<-FST_teststat_TN(X)

      bs_teststat<-sapply(1:M,function(o){
        bs_ind<-sample(1:N,N,replace=T)
        bs_hat_eta<-hat_eta[bs_ind,]

        bsTN<-FST_teststat_TN(bs_hat_eta)
        return(bsTN)

      }
      )
      fst_result<-list(statistic=c(T_N=TN),p.value=mean(TN<bs_teststat),sample.size=N,method="Simple-bootstrap-based functional stationarity test")
      }

      if (STAT==2) {

        MN<-FST_teststat_MN(X)

        bs_teststat<-sapply(1:M,function(o){
          bs_ind<-sample(1:N,N,replace=T)
          bs_hat_eta<-hat_eta[bs_ind,]

          bsMN<-FST_teststat_MN(bs_hat_eta)
          return(bsMN)

        }
        )
        fst_result<-list(statistic=c(M_N=MN),p.value=mean(MN<bs_teststat),sample.size=N,method="Simple-bootstrap-based functional stationarity test")
      }

    } else {
      if (STAT==1){
      TN<-FST_teststat_TN(X)
      bbs_teststat<-sapply(1:M,function(o){

        bs<-sample(1:(N-b+1),ceiling(N/b),replace=T)
        bbs<-bs[1]:(bs[1]+b-1)
        for (i in 2:ceiling(N/b)){
          bbs<-c(bbs,bs[i]:(bs[i]+b-1))
        }
        bs_hat_eta<-hat_eta[bbs[1:N],]

        bsTN<-FST_teststat_TN(bs_hat_eta)
        return(bsTN)

      }
      )
      fst_result<-list(statistic=c(T_N=TN),p.value=mean(TN<bbs_teststat),sample.size=N,method="Moving-block-bootstrap-based functional stationarity test")
      }
      if (STAT==2){
        MN<-FST_teststat_MN(X)
        bbs_teststat<-sapply(1:M,function(o){

          bs<-sample(1:(N-b+1),ceiling(N/b),replace=T)
          bbs<-bs[1]:(bs[1]+b-1)
          for (i in 2:ceiling(N/b)){
            bbs<-c(bbs,bs[i]:(bs[i]+b-1))
          }
          bs_hat_eta<-hat_eta[bbs[1:N],]

          bsMN<-FST_teststat_MN(bs_hat_eta)
          return(bsMN)

        }
        )
        fst_result<-list(statistic=c(M_N=MN),p.value=mean(MN<bbs_teststat),sample.size=N,method="Moving-block-bootstrap-based functional stationarity test")
      }
    }

  }
  class(fst_result)<-"htest"
  return(fst_result)
}

