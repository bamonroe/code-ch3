#include <Rcpp.h>
using namespace Rcpp;


double crra(double x,double r){
	return(pow(x,r) / r);
}


NumericVector perH(NumericVector R, NumericVector MU, NumericVector A0, NumericVector A1, NumericVector B0, NumericVector B1,
										NumericVector pA0, NumericVector pA1, NumericVector pB0, NumericVector pB1,
										NumericVector max, NumericVector min, NumericVector c){

	int nn  = R.size();

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

	for(int i=0;i<nn;i++){
		r[i] = 1 - R[i];
		mu[i] = MU[i];

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
double MSL_cpp(NumericVector par, NumericMatrix h1, NumericMatrix h2,
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

	NumericMatrix sim(rnum,h);

	for(int i = 0; i < h ; i++){
		r(_,i)  = qnorm(h1(_,i),rm,rs);
		mu(_,i) = qgamma(h2(_,i),k,t);
		sim(_,i) = perH( r(_,i),mu(_,i), A0, A1, B0, B1, pA0, pA1, pB0, pB1, max, min, c);
	}

	NumericVector sl(rnum);

	for(int i = 0; i < rnum ; i++){
		sl[i] = mean(sim(i,_));
	}

	sl = log(sl);
	return(-sum(sl));


}

