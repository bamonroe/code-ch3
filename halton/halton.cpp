#include <math.h>  
#include <Rcpp.h>
using namespace Rcpp;

double hnum(int index,int base){

	double r = 0;
	double f = 1;
	int i = index;

	while(i >0){
		f = f / base;
		r = r + f * (i % base);
		i = floor(i / base);
	}
	return(r);
}

// [[Rcpp::export]]
NumericVector halton(int init, int H, int prime){

	int stop = H - init +1 ;
    
	NumericVector out(stop);

    int j;

    for(int i = 0 ; i < stop ; i++){

    	j = i + init;

		out[i] = hnum(j,prime);
        
	}
	return(out);

}


