#' @title An independence test using R
#' @description An independence test using R by the nonparametric method.
#' @param d1  the first data vector
#' @param d2  the second data vector
#' @param alpha  the significance for test
#' @return p the true significance
#' @return T the result of the test
#' @importFrom stats rnorm dnorm pnorm var
#' @import Ball
#' @import GeneralizedHyperbolic
#' @import MASS
#' @import bootstrap 
#' @import microbenchmark 
#' @import rootSolve
#' @examples
#' \dontrun{
#' d1<-faithful$eruptions
#' d2<-faithful$waiting
#' Test<-independence(d1,d2,0.05)
#' Test$T
#' Test$p
#' }
#' @export
independence<-function(d1,d2,alpha){
  v1<-var(d1)
  v2<-var(d2)
  d1<-(d1-mean(d1))/sqrt(v1)
  d2<-(d2-mean(d2))/sqrt(v2)         
  n<-length(d1)
  T<-0
  f2<-f1<-numeric(n)           
  for(i in 1:n){
    for(j in 1:n){
      if(j!=i){
        f1[i]<-f1[i]+dnorm(d1[i]-d1[j])/(n-1)
      }
    }
  }
  for(i in 1:n){
    for(j in 1:n){
      if(j!=i){
        f2[i]<-f2[i]+dnorm(d2[i]-d2[j])/(n-1)
      }
    }
  }
  for(i in 1:n){                   
    for(j in 1:n){
      if(j!=i){
        T<-T+dnorm(d1[j]-d1[i])*dnorm(d2[j]-d2[i])/n/(n-1)
      }
    }
    T<-T+2*f1[i]*f2[i]/n
    for(j in 1:n){
      T<-T+f1[i]*f2[j]/n/n
    }
  }
  sigma<-0                     
  for(i in 1:n){
    for(j in 1:n){
      if(j!=i){
        sigma<-sigma+dnorm(d1[i]-d1[j])^2*dnorm(d2[i]-d2[j])^2
      }
    }
  }
  sigma<-sigma*2/n/n
  sigma<-sqrt(sigma)
  Ttest<-n*abs(T)/sigma                
  p<-(1-pnorm(Ttest))/2
  return(list(p=p,T= p>alpha))
}   
