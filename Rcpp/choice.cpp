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

// [[Rcpp::export]]
NumericMatrix getChoice(NumericVector r, NumericVector mu,
					NumericMatrix A,  NumericMatrix B,
					NumericMatrix pA, NumericMatrix pB,
					NumericVector max, NumericVector min, int tnum){

//IntegerMatrix getChoice(NumericVector r, NumericVector mu,

	int nn = r.length();

	NumericMatrix c(nn,tnum);

	NumericVector N0(nn);
	NumericVector N1(nn);
	NumericVector PA(nn);

	// Calculate the context
	NumericVector ctx = crra(max,r) - crra(min,r);

	// Calculate the utility of the lotteries
	NumericVector UA = (pA(_,0) * crra(A(_,0),r)) + (pA(_,1) * crra(A(_,1),r));
	NumericVector UB = (pB(_,0) * crra(B(_,0),r)) + (pB(_,1) * crra(B(_,1),r));

	// Re-base utility of B and add in context and fechner
	NumericVector UB1  = (UB/ctx/mu) - (UA/ctx/mu);

	// If we have no issues, this is the choice probability of A
	NumericVector pa = (1 / (1 + exp(UB1)));

	// Are we dealing with an insane number?
	// yes
	N0 = ifelse( UB > UA , 0 , 1 ); 
	// no, but are we making an insane number via exp?
	
	N1 = ifelse( UB1 > 709 , 0 , pa );

	// Check for the 2 issues, and return the probability of A
	PA = ifelse( is_na(UB1) , N0 , N1);

	//Making pB = 1-pA saves us the exponential calculations - it's faster
	// Don't actually need to calculate pB
	//	NumericVector pB = 1 - pA;
	

	NumericVector rand(nn);

	for(int i = 0; i < tnum ; i++){
        
		rand = runif(nn);

		// Grab the choice 
		c(_,i) = ifelse(PA>rand,0,1);
		//c(_,i) = PA;
	}

	return(c);

}


