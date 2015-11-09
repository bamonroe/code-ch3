#include <Rcpp.h>
using namespace Rcpp;

// pow takes any base, but only up to double exponent, this function fixes that

// [[Rcpp::export]]
NumericVector vpow(const NumericVector base, const NumericVector exp) {
  NumericVector out(base.size());
  std::transform(base.begin(), base.end(),
                 exp.begin(), out.begin(), ::pow);
  return out;
}
