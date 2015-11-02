#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix DDcpp(NumericMatrix pat, NumericVector M, NumericMatrix Errors, 
				NumericMatrix Cert, NumericMatrix CEMax, NumericMatrix Probs){

	int snum = pat.ncol();
	int cnum = Errors.ncol(); 
	double EE=0,PC=1,WC=0,WP=0;

	int choice[20];
	double Eset,E0=0,E1=0,E2=0,E3=0,E4=0,E5=0,E6=0,E7=0,E8=0,E9=0,E10=0;

	NumericMatrix Res( snum,15 );
	/* These are the 24 columns
	EE
	PC
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

	EE=0,PC=1,WC=0,WP=0;
	E0=0,E1=0,E2=0,E3=0,E4=0,E5=0,E6=0,E7=0,E8=0,E9=0,E10=0;

		// Set up the chocie vector
		for(int j = 0 ; j < 10 ; j++){

			if(pat(j,i)==0){
				choice[j] = 1;
				EE = EE + M[j];
				PC = PC * M[50+j];
				WC = WC + ( M[20+j] - M[30+j] ) ; 
				WP = WP + ( M[20+j] / M[40+j] ) ; 
			}
			else{
				choice[j] = 0;
				EE = EE + M[j+10];
				PC = PC * M[60+j];
				WC = WC + ( M[30+j] - M[20+j] ) ;
				WP = WP + ( M[30+j] / M[40+j] ) ; 
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
		Res(i,2) = WC;
		Res(i,3) = WP;
		Res(i,4) = E0;
		Res(i,5) = E1;
		Res(i,6) = E2;
		Res(i,7) = E3;
		Res(i,8) = E4;
		Res(i,9) = E5;
		Res(i,10) = E6;
		Res(i,11) = E7;
		Res(i,12) = E8;
		Res(i,13) = E9;
		Res(i,14) = E10;

	}

	return(Res);

}



