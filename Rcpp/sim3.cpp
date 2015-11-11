# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// pow takes any base, but only up to double exponent, this function fixes that
NumericVector vpow(const NumericVector base, const NumericVector exp) {
  NumericVector out(base.size());
  std::transform(base.begin(), base.end(),
                 exp.begin(), out.begin(), ::pow);
  return out;
}

// pow takes any base, but only up to double exponent, this function fixes that
arma::vec avpow(const arma::vec base, const arma::vec exp) {
	arma::vec out(base.n_rows);
  std::transform(base.begin(), base.end(),
                 exp.begin(), out.begin(), ::pow);
  return out;
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
	NumericVector ctx = vpow(max,(1-r)/(1-r)) - vpow(min,(1-r)/(1-r));

	// Calculate the utility of the lotteries
	NumericVector UA = (pA0 * vpow(A0,(1-r))/(1-r)) + (pA1 * vpow(A1,(1-r))/(1-r));
	NumericVector UB = (pB0 * vpow(B0,(1-r))/(1-r)) + (pB1 * vpow(B1,(1-r))/(1-r));

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

// [[Rcpp::export]]
double MSL_EUT(NumericVector par, NumericMatrix h1, NumericMatrix h2,
					NumericVector A0, NumericVector A1, NumericVector B0, NumericVector B1,
					NumericVector pA0, NumericVector pA1, NumericVector pB0, NumericVector pB1,
					NumericVector max, NumericVector min, NumericVector c){
    

	double rm = par[0];
	double rs = exp(par[1]);
	double um = exp(par[2]);
	double us = exp(par[3]);

	double k = pow(um,2) / pow(us,2);
	double t = pow(us,2) / um;

	int rnum = h1.nrow();
	int h = h1.ncol();

	NumericMatrix r(rnum,h);
	NumericMatrix mu(rnum,h);

	NumericVector N0(rnum);
	NumericVector N1(rnum);
	NumericVector pA(rnum);

	NumericMatrix sim(rnum,h);

	for(int i = 0; i < h ; i++){

		r(_,i)  = qnorm(h1(_,i),rm,rs);
		mu(_,i) = qgamma(h2(_,i),k,t);

		// Calculate the context
		NumericVector ctx = vpow(max,(1-r(_,i)))/(1-r(_,i)) - vpow(min,(1-r(_,i)))/(1-r(_,i));

		// Calculate the utility of the lotteries
		NumericVector UA = (pA0 * vpow(A0,(1-r(_,i)))/(1-r(_,i))) + (pA1 * vpow(A1,(1-r(_,i)))/(1-r(_,i)));
		NumericVector UB = (pB0 * vpow(B0,(1-r(_,i)))/(1-r(_,i))) + (pB1 * vpow(B1,(1-r(_,i)))/(1-r(_,i)));

		// Re-base utility of B and add in context and fechner
		NumericVector UB1  = (UB/ctx/mu(_,i)) - (UA/ctx/mu(_,i));

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
		sim(_,i) = ifelse(c==0,pA,pB);

	}

	NumericVector sl(rnum);

	for(int i = 0; i < rnum ; i++){
		sl[i] = mean(sim(i,_));
	}

	sl = log(sl);
	return(-sum(sl));

}

NumericVector pw(NumericVector p , NumericVector a, NumericVector b){

     return( exp(-b * vpow(-log(p), a)) );
    
}

// [[Rcpp::export]]
double MSL_RDU(NumericVector par, NumericMatrix h1, NumericMatrix h2, NumericMatrix h3, NumericMatrix h4,
					NumericVector A0, NumericVector A1, NumericVector B0, NumericVector B1,
					NumericVector pA0, NumericVector pA1, NumericVector pB0, NumericVector pB1,
					NumericVector max, NumericVector min, NumericVector c){
    

	double rm = par[0];
	double rs = exp(par[1]);
	double um = exp(par[2]);
	double us = exp(par[3]);

	double am = exp(par[4]);
	double as = exp(par[5]);

	double bm = exp(par[6]);
	double bs = exp(par[7]);

	double uk = pow(um,2) / pow(us,2);
	double ut = pow(us,2) / um;

	double ak = pow(am,2) / pow(as,2);
	double at = pow(as,2) / am;

	double bk = pow(bm,2) / pow(bs,2);
	double bt = pow(bs,2) / bm;

	int rnum = h1.nrow();
	int h = h1.ncol();

	NumericMatrix r(rnum,h);
	NumericMatrix mu(rnum,h);
	NumericMatrix alpha(rnum,h);
	NumericMatrix beta(rnum,h);


	NumericVector wA1(rnum);
	NumericVector wA0(rnum);
	NumericVector wB1(rnum);
	NumericVector wB0(rnum);


	NumericVector N0(rnum);
	NumericVector N1(rnum);
	NumericVector pA(rnum);

	NumericMatrix sim(rnum,h);

	for(int i = 0; i < h ; i++){

		r(_,i)  = qnorm(h1(_,i),rm,rs);
		mu(_,i) = qgamma(h2(_,i),uk,ut);
		alpha(_,i) = qgamma(h3(_,i),ak,at);
		beta (_,i) = qgamma(h4(_,i),bk,bt);


		wA1 = pw(pA1,alpha(_,i),beta(_,i));
		wA0 = pw(pA0 + pA1,alpha(_,i),beta(_,i)) - pw(pA1,alpha(_,i),beta(_,i));

		wB1 = pw(pB1,alpha(_,i),beta(_,i));
		wB0 = pw(pB0 + pB1,alpha(_,i),beta(_,i)) - pw(pB1,alpha(_,i),beta(_,i));


		// Calculate the context
		NumericVector ctx = vpow(max,(1-r(_,i)))/(1-r(_,i)) - vpow(min,(1-r(_,i)))/(1-r(_,i));

		// Calculate the utility of the lotteries
		NumericVector UA = (wA0 * vpow(A0,(1-r(_,i)))/(1-r(_,i))) + (wA1 * vpow(A1,(1-r(_,i)))/(1-r(_,i)));
		NumericVector UB = (wB0 * vpow(B0,(1-r(_,i)))/(1-r(_,i))) + (wB1 * vpow(B1,(1-r(_,i)))/(1-r(_,i)));

		// Re-base utility of B and add in context and fechner
		NumericVector UB1  = (UB/ctx/mu(_,i)) - (UA/ctx/mu(_,i));

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
		sim(_,i) = ifelse(c==0,pA,pB);

	}

	NumericVector sl(rnum);

	for(int i = 0; i < rnum ; i++){
		sl[i] = mean(sim(i,_));
	}

	sl = log(sl);
	return(-sum(sl));

}




// [[Rcpp::export]]
double MSL_cpp3(arma::vec par, NumericMatrix h1, NumericMatrix h2,
					arma::vec prizes,
					arma::mat pA0, arma::mat pA1, arma::mat pB0, arma::mat pB1,
					arma::mat max, arma::mat min, arma::imat c){
    

	double rm = par[0];
	double rs = exp(par[1]);
	double um = exp(par[2]);
	double us = exp(par[3]);

	double k = pow(um,2) / pow(us,2);
	double t = pow(us,2) / um;

	int rnum = h1.nrow();
	int h = h1.ncol();

	NumericMatrix r(rnum,h);
	NumericMatrix mu(rnum,h);

	for(int i = 0; i < h ; i++){
		r(_,i)  = qnorm(h1(_,i),rm,rs);
		mu(_,i) = qgamma(h2(_,i),k,t);
	}

	// Change into arma
	arma::mat rr = as<arma::mat>(r);
	arma::mat mu2 = as<arma::mat>(mu);


	arma::mat UA0 = exp( log(prizes[0]) * rr ) / (1-rr); 
	arma::mat UA1 = exp( log(prizes[1]) * rr ) / (1-rr); 
	arma::mat UB0 = exp( log(prizes[2]) * rr ) / (1-rr); 
	arma::mat UB1 = exp( log(prizes[3]) * rr ) / (1-rr); 
			  max = exp(       log(max) % rr ) / (1-rr); 
			  min = exp(       log(min) % rr ) / (1-rr); 


	// Calculate the context
	arma::mat ctx = max - min;

	// Calculate the utility of the lotteries
	arma::mat UA = (pA0 % UA0) + (pA1 % UA1);
	arma::mat UB = (pB0 % UB0) + (pB1 % UB1);

	// Re-base utility of B and add in context and fechner
	arma::mat UBB  = (UB/ctx/mu2) - (UA/ctx/mu2);

	// If we have no issues, this is the choice probability of A
	arma::mat PA = (1 / (1 + exp(UBB)));

	// Inverse the choices
	arma::imat cA = 1 - c;

	// Element by element to get choice probs
	arma::mat sim = (cA % PA) + (c % (1-PA));

	// Get the logged means
	arma::vec sl = log( mean(sim,1) );

	double res = sum(sl);

	return(-res);

}

