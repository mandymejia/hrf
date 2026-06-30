// Rcpp scaffold. Add new C++ kernels here as standalone .cpp files; run
// devtools::document() (which calls Rcpp::compileAttributes()) to regen
// src/RcppExports.cpp and R/RcppExports.R.

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void say_hello() {
  Rcout << "Hello from hrf C++ scaffold\n";
}
