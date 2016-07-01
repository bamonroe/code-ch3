# include <math.h>
# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


NumericVector crra(NumericVector x , double r) {
	if(r!=1)
		return(pow(x,(1-r)) / (1-r));
	else
		return(log(x));
}


NumericVector perH(NumericVector x, NumericVector A0, NumericVector A1, NumericVector B0, NumericVector B1,
							NumericVector pA0, NumericVector pA1, NumericVector pB0, NumericVector pB1,
							NumericVector max, NumericVector min){

	int nn  = A0.size();

	NumericVector ctx(nn) ;
	NumericVector UA(nn)  ;
	NumericVector UB(nn)  ;
	NumericVector UA1(nn) ;
	NumericVector UB1(nn) ;

	NumericVector wA0(nn)  ;
	NumericVector wA1(nn)  ;
	NumericVector wB0(nn)  ;
	NumericVector wB1(nn)  ;

	NumericVector pA(nn)  ;
	NumericVector pB(nn)  ;

	NumericVector Aerr(nn);
	NumericVector Berr(nn);

	NumericVector PA(nn);

	NumericVector N0(nn);
	NumericVector N1(nn);

	NumericVector CEA(nn) ;
	NumericVector CEB(nn) ;
	NumericVector CEM(nn) ;

	ctx.fill(0);
	UA.fill(0);
	UB.fill(0);
	UA1.fill(0);
	UB1.fill(0);

	wA0.fill(0);
	wA1.fill(0);
	wB0.fill(0);
	wB1.fill(0);

	pA.fill(0);
	pB.fill(0);
	
	Aerr.fill(0);
	Berr.fill(0);
	
	PA.fill(0);
	
	N0.fill(0);
	N1.fill(0);
	
	CEA.fill(0);
	CEB.fill(0);
	CEM.fill(0);

	double r = x[0];
	double mu = x[1];

	if (x.size() == 2) {

		wA1 = pA1 ;
		wA0 = pA0 ;

		wB1 = pB1 ;
		wB0 = pB0 ;
	
	}
	else{

		double alpha = x[2];
		double beta = x[3];

		wA1 = exp( -1 * beta * pow( -1 * log(pA1), alpha) ) ;
		wA0 = exp( -1 * beta * pow( -1 * log(pA0 + pA1), alpha) )  - wA1 ;
		wB1 = exp( -1 * beta * pow( -1 * log(pB1), alpha) ) ;
		wB0 = exp( -1 * beta * pow( -1 * log(pB0 + pB1), alpha) )  - wB1 ;
	
	}

	NumericVector res(nn*7);
	res.fill(999);

	ctx = crra(max,r) - crra(min,r);

	UA = (wA0 * crra(A0,r)) + (wA1 * crra(A1,r));
	UB = (wB0 * crra(B0,r)) + (wB1 * crra(B1,r));

	Aerr = ifelse(UB > UA,1,0);
	Berr = ifelse(UA > UB,1,0);

	CEA = pow((UA * (1-r)), (1/(1-r)));
	CEB = pow((UB * (1-r)), (1/(1-r)));

	CEM = ifelse(CEA>CEB,CEA,CEB);

	UB1  = (UB/ctx/mu) - (UA/ctx/mu);

	PA = (1 / (1 + exp(UB1)));

	// Are we dealing with an insane number?
	// yes
	N0 = ifelse( UB > UA , 0 , 1 ); 
	// no, but are we making an insane number via exp?
	N1 = ifelse( UB1 > 709 , 0 , PA );

	pA = ifelse( is_na(UB1) , N0 , N1);

	pB = 1 - pA;

	for(int i = 0 ; i < nn ; i++){
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


// [[Rcpp::export]]
NumericMatrix DDcpp(NumericMatrix pat, NumericVector M, NumericMatrix Errors, 
				NumericMatrix Cert, NumericMatrix CEMax, NumericMatrix Probs){

	int snum = pat.ncol();
	int choices = pat.nrow();
	int cnum = Errors.ncol(); 
	double EE=0, PC=1, LPC=0, WC=0, WP=0;

	int choice[20];
	double Eset, E0=0, E1=0, E2=0, E3=0, E4=0, E5=0, E6=0, E7=0, E8=0, E9=0, E10=0;

	NumericMatrix Res( snum,16 );
	/* These are the 24 columns
	EE
	PC
	LPC
	WC
	WP
	E.0   
	E.1
	E.2
	E.3
	E.4
	E.5
	E.6
	E.7
	E.8
	E.9
	E.10 
	**Stopped here
	WP.1 
	WP.2 
	WP.3 
	WP.4 
	WP.5 
	WP.6 
	WP.7 
	WP.8 
	WP.9 
	WP.10
	*/

	// Loop through each pattern
	for(int i = 0 ; i < snum ; i++){

	EE=0, PC=1, LPC=0, WC=0, WP=0;
	E0=0, E1=0, E2=0, E3=0, E4=0, E5=0, E6=0, E7=0, E8=0, E9=0, E10=0;

		// Set up the choice vector
		for(int j = 0 ; j < choices ; j++){

			if(pat(j,i)==0){
				choice[j] = 1;
				EE  = EE + M[j];
				PC  = PC * M[50+j];
				LPC = LPC + std::log(M[50+j]);
				WC  = WC + ( M[20+j] - M[30+j] ) ; 
				WP  = WP + ( M[20+j] / M[40+j] ) ; 
			}
			else{
				choice[j] = 0;
				EE  = EE + M[j+10];
				PC  = PC * M[60+j];
				LPC = LPC + std::log(M[60+j]);
				WC  = WC + ( M[30+j] - M[20+j] ) ;
				WP  = WP + ( M[30+j] / M[40+j] ) ; 
			}

			choice[j+10] = pat(j,i);

		}

		//Set up the Eset
		for(int j = 0 ; j < cnum ; j++){
			Eset = 0 ;

			for(int k = 0 ; k < 20 ; k++){
				Eset = Eset + (choice[k] * Errors(k,j));
				//WPset = WPset + ( (choice[k] * Cert(k,j))/CEMax(k,j) )/10   ;
			}

			if(Eset ==0)
				E0 = E0 + 1;
			else if(Eset==1)
				E1 = E1 + 1;
			else if(Eset==2)
				E2 = E2 + 1;
			else if(Eset==3)
				E3 = E3 + 1;
			else if(Eset==4)
				E4 = E4 + 1;
			else if(Eset==5)
				E5 = E5 + 1;
			else if(Eset==6)
				E6 = E6 + 1;
			else if(Eset==7)
				E7 = E7 + 1;
			else if(Eset==8)
				E8 = E8 + 1;
			else if(Eset==9)
				E9 = E9 + 1;
			else if(Eset==10)
				E10 = E10 + 1;

		}

		Res(i,0) = EE;
		Res(i,1) = PC;
		Res(i,2) = LPC;
		Res(i,3) = WC;
		Res(i,4) = WP;
		Res(i,5) = E0;
		Res(i,6) = E1;
		Res(i,7) = E2;
		Res(i,8) = E3;
		Res(i,9) = E4;
		Res(i,10) = E5;
		Res(i,11) = E6;
		Res(i,12) = E7;
		Res(i,13) = E8;
		Res(i,14) = E9;
		Res(i,15) = E10;

	}

	return(Res);

}


