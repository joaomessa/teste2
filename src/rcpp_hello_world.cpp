
#include <Rcpp.h>
#include "Cabecalho.h"
using namespace Rcpp;


// [[Rcpp::export]]
List rcpp_hello_world() {

    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;

    return z ;
}

// [[Rcpp::export]]
int signC(int x, int y){
  if (x+y>0) {
    return 1;
  } else if (x+y==0){
    return 0;
  } else {
    return -1;
  }
}

// [[Rcpp::export]]
double sumC(NumericVector x){
  int n = x.size();
  double total = 0;
  for(int i = 0; i < n; ++i){
    total += x[i];
  }
  return total;
}

// [[Rcpp::export]]
int add2(int j, int k){
  return teste::add(j,k);
}

