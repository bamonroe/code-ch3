#include <math.h>       /* exp, pow */
#include <Rcpp.h>
#include <omp.h>
using namespace Rcpp;


double crra(double x,double r){
	return(pow(x,r) / r);
}

NumericVector perH(NumericVector x, NumericVector A0, NumericVector A1, NumericVector B0, NumericVector B1,
										NumericVector pA0, NumericVector pA1, NumericVector pB0, NumericVector pB1,
										NumericVector max, NumericVector min, NumericVector c){

	int n  = x.size();
	int nn = n/2;

	double ctx[nn];
	double UA[nn];
	double UB[nn];
	double UA1[nn];
	double UB1[nn];
	double pA[nn];
	double pB[nn];


	double r[nn];
	double mu[nn];

	NumericVector ll(nn);

#pragma omp parallel for 
	for(int i=0;i<nn;i++){
		r[i] = 1 - x[i];
		mu[i] = x[i+nn];

		ctx[i] = crra(max[i],r[i]) - crra(min[i],r[i]);

		UA[i] = (pA0[i] * crra(A0[i],r[i])) + (pA1[i] * crra(A1[i],r[i]));
		UB[i] = (pB0[i] * crra(B0[i],r[i])) + (pB1[i] * crra(B1[i],r[i]));

		UB1[i]  = (UB[i]/ctx[i]/mu[i]) - (UA[i]/ctx[i]/mu[i]);
		UA1[i]  = 0;

		if(  UB1[i] == NA_REAL ){
			if( UB[i] > UA[i]) {
				pA[i] = 0;
				pB[i] = 1;
			}
			else{
				pA[i] = 1;
				pB[i] = 0;
			}
		}
		else{
			if( UB1[i] > 709 ){
				pA[i] = 0;
				pB[i] = 1;
			}
			else{
				if( UB1[i] < -709 ) {
					pA[i] = 1;
					pB[i] = 0;
				}
				else{
					pA[i] = exp(UA1[i]) / (exp(UB1[i]) + 1);
					pB[i] = exp(UB1[i]) / (exp(UB1[i]) + 1);
				}
			}
		}
		if( c[i] == 0 ){
			ll[i] = pA[i];
		}
		else{
			ll[i] = pB[i];
		}
	}

	return(ll);

}

// [[Rcpp::export]]
NumericVector simpar(NumericMatrix DD, NumericVector A0, NumericVector A1, NumericVector B0, NumericVector B1,
										NumericVector pA0, NumericVector pA1, NumericVector pB0, NumericVector pB1,
										NumericVector max, NumericVector min, NumericVector c){

	int nrow = DD.nrow() , ncol = DD.ncol(), cnum = c.length();
	NumericMatrix sim(cnum,ncol);

	for(int i=0;i<ncol;i++){
		sim(_,i) = perH( DD(_,i), A0, A1, B0, B1, pA0, pA1, pB0, pB1, max, min, c);
	}

	return(sim);

}
