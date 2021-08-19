#' Functional KPSS Test
#'
#' KPSS test for functional time series with different methods on determining the critical values of the test statistic. The Monte Carlo method (for a large sample size) was constructed in Kokoszka and Young (2016) and the bootstrap-based methods (both the simple bootstrap method and moving block bootstrap method, for a small/moderate sample size) were constructed in Chen and Pun (2019).
#'
#' @param X The functional time series being tested, inputted in a matrix form with each row representing each observation of the functional data values on equidistant points of any prespecified interval.
#' @param ALPHA Significance level of the test. The default value is 5\%.
#' @param METHOD Method to determine the critical value of the test statistic. The default value 'METHOD="MC"' represents the Monte Carlo method. 'METHOD="SB"' represents the simple bootstrap method and 'METHOD="MBB"' represents the moving block bootstrap method.
#' @param K Kernel function in the estimation of the long-run covariance function, which is only effective in the Monte Carlo method. The default function is 'default_kernel' function in this package.
#' @param h_power Power of sample size 'N' (valued in (0,1)) for the smoothing bandwidth, which is only effective in the Monte Carlo method. The default value is 2/5.
#' @param est_ev Number of the largest eigenvalues chosen to estimate the limiting distribution, which is only effective in the Monte Carlo method. The default value is the sample size 'N'.
#' @param MCNsim Number of Monte Carlo datasets generated in the Monte Carlo method, which is only effective in Monte Carlo method. The default value is 10000.
#' @param sbb_set A vector of independent simulated data generated from the function 'dataset_sbb', which is only effective and essential in Monte Carlo method.
#' @param M Number of bootstrap datasets generated in the bootstrap method, which is only effective in bootstrap methods. The default value is 1000.
#' @param b Block length used in the moving block bootstrap method, which is only effective in the moving block bootstrap method. The default value is ceiling((2N)^(1/3)), where 'N' is the sample size.
#' @return The result of the test is presented with the value of test statistic and its p-value under the null hypothesis of trend stationarity.
#' @export
#' @references Chen, Y., & Pun, C. S. (2019). A bootstrap-based KPSS test for functional time series. Journal of Multivariate Analysis, 174, 104535.
#' @references Kokoszka, P., & Young, G. (2016). KPSS test for functional time series. Statistics, 50(5), 957-973.
#' @examples
#' N<-100
#' EX<-matrix(rep(0,N*100),ncol=100)
#' EXI<-rep(0,100)
#' set.seed(1)
#' EXI[1]<-rnorm(1,0,1)
#' for (j in 2:100) {EXI[j]<-EXI[j-1]+rnorm(1,0,1)}
#' for (i in 1:N) {
#' temp<-rnorm(100,0,1)
#' EX[i,1]<-temp[1]
#' for (j in 2:100) {
#' EX[i,j]<-EX[i,j-1]+temp[j]
#' }
#' }
#' fkpsst(X=EX,METHOD="SB")

fkpsst<-function(X, ALPHA=0.05, METHOD="MC",
                 K=default_kernel, h_power=2/5, est_ev=nrow(X), MCNsim=10000, sbb_set=NULL,
                 M=1000, b=ceiling((2*nrow(X))^(1/3))){
  if (h_power<=0 | h_power>=1) stop("The value of 'h_power' exceeds the range (0,1).")
  if (METHOD=="MC" & is.null(sbb_set)) stop("You must generate a simluated dataset for Monte Carlo method from the function 'generate_sbb'.")

  N<-nrow(X)

  MeanX<-apply(X,2,sum)/N
  StairsumX<-apply((1:N)*X,2,sum)

  hat_eta<-X+(((1:N-(N+1)/2)*6/(N-1)-1))%*%t(MeanX)-((1:N-(N+1)/2)*12/(N*(N^2-1)))%*%t(StairsumX)

  ZN<-partial_sum(hat_eta)
  TN<-mean(ZN^2)

  if (METHOD=="MC"){

    lambda<-est_eigenvalue(hat_eta,K,h_power,est_ev)

    intV<-matrix(sample(sbb_set,MCNsim*est_ev,replace=T),ncol=est_ev)
    monte=(intV%*%lambda)[,1]

    fkpsst_result<-list(statistic=c(R_N=TN),p.value=mean(TN<monte),sample.size=N,method="Monte-Carlo-based functional KPSS test")

  } else {

    hat_xi<-apply((1:N-(N+1)/2)*X,2,sum)*(12/(N*(N^2-1)))
    hat_mu<-apply(X,2,sum)/N-hat_xi*(N+1)/2

    if (METHOD=="SB"){
      bs_teststat<-sapply(1:M,function(o){
        bs_ind<-sample(1:N,N,replace=T)
        bs_hat_eta<-hat_eta[bs_ind,]

        bsX=rep(1,N)%*%t(hat_mu)+(1:N)%*%t(hat_xi)+bs_hat_eta

        bsTN<-FKPSS_teststat(bsX)
        return(bsTN)

      }
      )
      fkpsst_result<-list(statistic=c(R_N=TN),p.value=mean(TN<bs_teststat),sample.size=N,method="Simple-bootstrap-based functional KPSS test")

    } else {
      bbs_teststat<-sapply(1:M,function(o){

        bs<-sample(1:(N-b+1),ceiling(N/b),replace=T)
        bbs<-bs[1]:(bs[1]+b-1)
        for (i in 2:ceiling(N/b)){
          bbs<-c(bbs,bs[i]:(bs[i]+b-1))
        }
        bs_hat_eta<-hat_eta[bbs[1:N],]

        bsX=rep(1,N)%*%t(hat_mu)+(1:N)%*%t(hat_xi)+bs_hat_eta

        bsTN<-FKPSS_teststat(bsX)
        return(bsTN)

      }
      )
      fkpsst_result<-list(statistic=c(R_N=TN),p.value=mean(TN<bbs_teststat),sample.size=N,method="Moving-block-bootstrap-based functional KPSS test")
    }
  }

  class(fkpsst_result)<-"htest"
  return(fkpsst_result)

}
