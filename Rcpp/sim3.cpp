#include <Rcpp.h>
using namespace Rcpp;

// pow takes any base, but only up to double exponent, this function fixes that
NumericVector vpow(const NumericVector base, const NumericVector exp) {
  NumericVector out(base.size());
  std::transform(base.begin(), base.end(),
                 exp.begin(), out.begin(), ::pow);
  return out;
}

NumericVector crra(NumericVector x , NumericVector r) {
		return(	ifelse(r!= 1, vpow(x,(1-r)) / (1-r), log(x)));
}

NumericVector perH(NumericVector r, NumericVector mu,
					NumericVector A0, NumericVector A1, NumericVector B0, NumericVector B1,
					NumericVector pA0, NumericVector pA1, NumericVector pB0, NumericVector pB1,
					NumericVector max, NumericVector min, NumericVector c){

	int nn = A0.size();

	NumericVector ll(nn);

	NumericVector N0(nn);
	NumericVector N1(nn);
	NumericVector pA(nn);

	// Calculate the context
	NumericVector ctx = crra(max,r) - crra(min,r);

	// Calculate the utility of the lotteries
	NumericVector UA = (pA0 * crra(A0,r)) + (pA1 * crra(A1,r));
	NumericVector UB = (pB0 * crra(B0,r)) + (pB1 * crra(B1,r));

	// Re-base utility of B and add in context and fechner
	NumericVector UB1  = (UB/ctx/mu) - (UA/ctx/mu);

	// If we have no issues, this is the choice probability of A
	NumericVector PA = (1 / (1 + exp(UB1)));

	// Are we dealing with an insane number?
	// yes
	N0 = ifelse( UB > UA , 0 , 1 ); 
	// no, but are we making an insane number via exp?
	N1 = ifelse( UB1 > 709 , 0 , PA );

	// Check for the 2 issues, and return the probability of A
	pA = ifelse( is_na(UB1) , N0 , N1);

	// Making pB = 1-pA saves us the exponential calculations - it's faster
	NumericVector pB = 1 - pA;

	// Grab the choice probability for the chosen option
	ll = ifelse(c==0,pA,pB);

	return(ll);

}

// [[Rcpp::export]]
NumericVector simcpp3(NumericMatrix r, NumericMatrix mu, 
						NumericVector A0, NumericVector A1, NumericVector B0, NumericVector B1,
						NumericVector pA0, NumericVector pA1, NumericVector pB0, NumericVector pB1,
						NumericVector max, NumericVector min, NumericVector c){

	int ncol = r.ncol(), rnum = c.length();

	// Sim is a matrix of the probabilities associated with the chosen option
	// The rows will be averaged in R as the means routine is more accurate and is very fast
	NumericMatrix sim(rnum,ncol);

	for(int i=0;i<ncol;i++){
		sim(_,i) = perH( r(_,i), mu(_,i), A0, A1, B0, B1, pA0, pA1, pB0, pB1, max, min, c);
	}

	return(sim);

}
