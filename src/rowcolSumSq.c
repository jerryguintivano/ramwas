#include <R.h>
#include <Rinternals.h>

SEXP CrowSumsSq(SEXP input, SEXP outpt)
{
	input = PROTECT(coerceVector(input, REALSXP));
	outpt = PROTECT(coerceVector(outpt, REALSXP));
	int nrow = length(outpt);
	int ncol = length(input) / nrow;
	
	double *xinput = REAL(input);
	double *xoutpt = REAL(outpt);
	
	for(int j = 0; j < ncol; j++) 
		for(int i = 0; i < nrow; i++)
			xoutpt[i] += xinput[j*nrow+i] * xinput[j*nrow+i];
	UNPROTECT(2);
	return R_NilValue;
}

SEXP CcolSumsSq(SEXP input, SEXP outpt)
{
	input = PROTECT(coerceVector(input, REALSXP));
	outpt = PROTECT(coerceVector(outpt, REALSXP));
	int ncol = length(outpt);
	int nrow = length(input) / ncol;
	
	double *xinput = REAL(input);
	double *xoutpt = REAL(outpt);
	
	for(int j = 0; j < ncol; j++) {
		double sm = 0;
		for(int i = 0; i < nrow; i++)
			sm += xinput[j*nrow+i] * xinput[j*nrow+i];
		xoutpt[j] = sm;
	}
	UNPROTECT(2);
	return R_NilValue;
}
