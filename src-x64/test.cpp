#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//' @title MetropolisC
//' @description The Rcpp function in homework.
//' @param sigma the variance of the normal
//' @param x0 the initial value
//' @param N the number of the trying
//' @return the path and the number of the accept point
//' @useDynLib SC19025
//' @examples
//' \dontrun{
//' x<-Metropolisc(1,25,2000)
//' x
//' }
//' @export
// [[Rcpp::export]]
NumericVector MetropolisC(double sigma,double x0,int N) {
  NumericVector x(N+1);
  x[0] = x0;
  int k=0;
  double y=0;
  double u=0;
  
  for (int i=1;i<N;i++) {
    u = runif(1)[0];
    y = rnorm(1,x[i-1], sigma)[0];
    if(u <= exp(-abs(y)+abs(x[i-1]))){
      x[i] = y;
    }
    else {
      x[i] = x[i-1];
      k = k + 1;
    }
  }
  x[N]=k;
  return(x);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//





