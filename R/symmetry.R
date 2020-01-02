#' @title A symmetry test using R
#' @description A symmetry test using R by the nonparametric method.
#' @param d  the data vector
#' @param alpha  the significance for test
#' @return p the true significance
#' @return T the result of the test
#' @importFrom stats rnorm rchisq var dnorm pnorm
#' @examples
#' \dontrun{
#' d1<-rnorm(100)
#' d2<-rchisq(100,df=1)
#' Test1<-symmetry(d1,0.05)
#' Test1
#' Test2<-symmetry(d2,0.05)
#' Test2
#' }
#' @export
symmetry<-function(d,alpha){
  n<-length(d)
  v<-var(d)
  d<-d/sqrt(v)
  T=0
  for(i in 1:n){
    for(j in 1:n){
      T<-T+(dnorm(d[i]-d[j])-dnorm(d[i]+d[j]))/n/n
      }
  }
  c<-dnorm(0)/n
  sigma<-0
  for(i in 1:n){
    for(j in 1:n){
      sigma<-sigma+dnorm(d[i]-d[j])
    }
  }
  sigma<-sigma/(4*n)/sqrt(pi)
  Ttest<-n*(T-c)/sqrt(sigma)
  p<-(1-pnorm(Ttest))/2
  return(list(p=p,T= p>alpha))
}

