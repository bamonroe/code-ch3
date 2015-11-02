#include <Rcpp.h>
using namespace Rcpp;

double crra(double x,double r){
	if(r==1)
	return(pow(x,(1-r)) / (1-r));
	else
	return(log(x));
}

NumericVector perH(NumericVector x, NumericVector A0, NumericVector A1, NumericVector B0, NumericVector B1,
										NumericVector pA0, NumericVector pA1, NumericVector pB0, NumericVector pB1,
										NumericVector max, NumericVector min){

	int nn  = A0.size();

	double ctx[nn];
	double UA[nn];
	double UB[nn];
	double UA1[nn];
	double UB1[nn];
	double pA[nn];
	double pB[nn];

	double Aerr[nn];
	double Berr[nn];

	double CEA[nn];
	double CEB[nn];
	double CEM[nn];

	double r = 1 - x[0];
	double mu = x[1];

	NumericVector res(nn*7);

	for(int i=0;i<nn;i++){

		ctx[i] = crra(max[i],r) - crra(min[i],r);

		UA[i] = (pA0[i] * crra(A0[i],r)) + (pA1[i] * crra(A1[i],r));
		UB[i] = (pB0[i] * crra(B0[i],r)) + (pB1[i] * crra(B1[i],r));

		Aerr[i] = 0;
		Berr[i] = 0;

		if( UB[i] > UA[i] )
			Aerr[i] = 1;
		else
			Berr[i] = 1;

		CEA[i] = pow((UA[i] * r), (1/r));
		CEB[i] = pow((UB[i] * r), (1/r));

		if( CEA[i] > CEB[i] )
			CEM[i] = CEA[i];
		else
			CEM[i] = CEB[i];

		UB1[i]  = (UB[i]/ctx[i]/r) - (UA[i]/ctx[i]/r);
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

		res[i]        = Aerr[i];
		res[i+nn]     = Berr[i];
		res[i+(2*nn)] = CEA[i];
		res[i+(3*nn)] = CEB[i];
		res[i+(4*nn)] = CEM[i];
		res[i+(5*nn)] = pA[i];
		res[i+(6*nn)] = pB[i];

	}

	return(res);

}

// [[Rcpp::export]]
NumericMatrix getRes(NumericMatrix sim, NumericVector A0, NumericVector A1, NumericVector B0, NumericVector B1,
										NumericVector pA0, NumericVector pA1, NumericVector pB0, NumericVector pB1,
										NumericVector max, NumericVector min){

	int nrow = sim.nrow() , ncol = sim.ncol(), cnum = A0.length();
	NumericMatrix RES((cnum*7),ncol);

	for(int i=0;i<ncol;i++){
		RES(_,i) = perH( sim(_,i), A0, A1, B0, B1, pA0, pA1, pB0, pB1, max, min);
	}

	return(RES);

}
