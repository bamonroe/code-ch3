clear all

use "../sim/choice.dta" ,clear

* Define number of Halton replications (internal default is $H if not specified) -- see nrep($hrep) option below
global H "250"

* Define maximum number of iterations if some have convergence problems
global imax "30"


/* MODEL 6 ADDING FECHNER ERROR*/
mata: mata clear
mata:
function utilfunc(XMAT, BETA)
{

	// BETA = 1 x nrep matrix of coefficients
	// XMAT[1,1..4] = 1 x 4 matrix of probabilities for lottery 1
	// XMAT[2,1..4] = 1 x 4 matrix of probabilities for lottery 2
	// XMAT[1,5..8] = XMAT[2,5..8] = 1 x 4 matrix of prizes
	// Note that :* denotes element-by-element multiplication
	// Note that :^ denotes element-by-element exponentiation 

	R  = BETA[1,.]	
	MU = exp(BETA[2,.])	

	if(R[1,1]!=0) {

    	V1 = XMAT[1,1..2] :* (XMAT[1,3]:^R',XMAT[1,4]:^R') 
    	V2 = XMAT[2,1..2] :* (XMAT[2,3]:^R',XMAT[2,4]:^R') 
    	CTX = XMAT[1,5]:^R' :- XMAT[1,6]:^R'

    }
    else{

		L = R' :+ 1

    	V1 = XMAT[1,1..2] :* (ln(XMAT[1,3] :* L),ln(XMAT[1,4] :* L)) 
    	V2 = XMAT[2,1..2] :* (ln(XMAT[2,3] :* L),ln(XMAT[2,4] :* L)) 
    	CTX = ln(XMAT[1,5] :* L) :- ln(XMAT[1,6] :* L)
    }

/*
	printf("%f\n",V1[1,1])
	printf("%f\n",V1[2,1])
	printf("%f\n",V1[3,1])
	printf("%f\n",CTX[1,1])
	printf("%f\n",CTX[2,1])
	printf("%f\n",CTX[3,1])
	printf("\n")
	*/

	// Fechner error index, the extra beta
	F1 = colsum(V1') :/ CTX'
	F2 = colsum(V2') :/ CTX'

	/*

    	printf("%f\n",F1[1,16])
    	printf("%f\n",F1[1,15])
    	printf("%f\n",F1[1,14])
    	printf("%f\n",F1[1,13])
    	printf("%f\n",F1[1,12])
    	printf("%f\n",F1[1,11])
    	printf("%f\n",F1[1,10])
    	printf("%f\n",F1[1,9])
    	printf("%f\n",F1[1,8])
    	printf("%f\n",F1[1,7])
    	printf("%f\n",F1[1,6])
    	printf("%f\n",F1[1,5])
    	printf("%f\n",F1[1,4])
    	printf("%f\n",F1[1,3])
    	printf("%f\n",F1[1,2])
    	printf("%f\n",F1[1,1])
    	printf("\n")

    	*/

	F1 = F1 :/ MU
	F2 = F2 :/ MU
	PEU = F1 \ F2
	EV = exp(PEU) :/ colsum(exp(PEU))

	return(EV)	

}
end

/*
mata: mata clear
mata:
function utilfunc(XMAT, BETA)
{

	// BETA = 1 x nrep matrix of coefficients
	// XMAT[1,1..4] = 1 x 4 matrix of probabilities for lottery 1
	// XMAT[2,1..4] = 1 x 4 matrix of probabilities for lottery 2
	// XMAT[1,5..8] = XMAT[2,5..8] = 1 x 4 matrix of prizes
	// Note that :* denotes element-by-element multiplication
	// Note that :^ denotes element-by-element exponentiation 

	R  = BETA[1,.]	
	MU = exp(BETA[2,.])	

    	V1 = XMAT[1,1..2] :* (XMAT[1,3]:^R',XMAT[1,4]:^R') 
    	V2 = XMAT[2,1..2] :* (XMAT[2,3]:^R',XMAT[2,4]:^R') 

	// Fechner error index, the extra beta
	F1 = colsum(V1') :/ MU
	F2 = colsum(V2') :/ MU
	PEU = F1 \ F2
	EV = exp(PEU) :/ colsum(exp(PEU))

	return(EV)

}
end
*/

constraint 1 [sd1]_cons = 0
constraint 2 [sd2]_cons = 0

global choice "c"
global probs "p0 p1"
global prizes "pz0 pz1"
global context "max min"

* assume that both r and mu are non-random
//mixlognl choice prob0 prob1 prob2 prob3 prize0 prize1 prize2 prize3 , group(gid) id(id) kfix(0) krnd(2) search nrep($H) constraint(1 2)
//mixlognl $choice $probs $prizes $context, group(gid) id(id) kfix(0) krnd(2) search nrep($H) constraint(1 2)

* assume that both r is random and mu is non-random
//mixlognl $choice $probs $prizes $context, group(gid) id(id) kfix(0) krnd(2) search nrep($H) constraint(2)

* assume that both r and mu are random
mixlognl $choice $probs $prizes $context, group(gid) id(id) kfix(0) krnd(2) search nrep($H)

* assume that both r and mu are random, and allow for correlation
mixlognl $choice $probs $prizes $context, group(gid) id(id) kfix(0) krnd(2) search nrep($H) corr
mixlcov_nl
mixlcov_nl, sd

