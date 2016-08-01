#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

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

	
SEXP cover_frw_c(SEXP starts, SEXP cpgs, SEXP fragdist, SEXP ind1, SEXP ind2, SEXP out) {
	starts   = PROTECT(coerceVector(starts,    INTSXP));
	cpgs     = PROTECT(coerceVector(cpgs,      INTSXP));
	fragdist = PROTECT(coerceVector(fragdist, REALSXP));
	ind1     = PROTECT(coerceVector(ind1,      INTSXP));
	ind2     = PROTECT(coerceVector(ind2,      INTSXP));
	out      = PROTECT(coerceVector(out,      REALSXP));
	
	int *xind1 = INTEGER(ind1);
	int *xind2 = INTEGER(ind2);
	int nind = length(ind1);
	int *xcpgs = INTEGER(cpgs);
	double *xfragdist = REAL(fragdist);
	int *xstarts = INTEGER(starts);
	double *xout = REAL(out);
		
	for(int i = 0; i < nind; i++) {
		if(xind1[i] < xind2[i]) {
			int curcpg = xcpgs[i];
			double tempsum = 0;
			for(int j = xind1[i]; j < xind2[i]; j++) {
//				Rprintf("Adding read %d at position %d, distance to start %d, extra coverage %f\n", j, i, curcpg-xstarts[j], xfragdist[curcpg-xstarts[j]]);
				tempsum = tempsum + xfragdist[curcpg-xstarts[j]];
			}
			xout[i] = xout[i] + tempsum;
		}
	}
	
	UNPROTECT(6);
	return(R_NilValue);
}

SEXP cover_rev_c(SEXP starts, SEXP cpgs, SEXP fragdist, SEXP ind1, SEXP ind2, SEXP out) {
	starts   = PROTECT(coerceVector(starts,    INTSXP));
	cpgs     = PROTECT(coerceVector(cpgs,      INTSXP));
	fragdist = PROTECT(coerceVector(fragdist, REALSXP));
	ind1     = PROTECT(coerceVector(ind1,      INTSXP));
	ind2     = PROTECT(coerceVector(ind2,      INTSXP));
	out      = PROTECT(coerceVector(out,      REALSXP));
	
	int *xind1 = INTEGER(ind1);
	int *xind2 = INTEGER(ind2);
	int nind = length(ind1);
	int *xcpgs = INTEGER(cpgs);
	double *xfragdist = REAL(fragdist);
	int *xstarts = INTEGER(starts);
	double *xout = REAL(out);
	
	for(int i = 0; i < nind; i++) {
		if(xind1[i] < xind2[i]) {
			int curcpg = xcpgs[i];
			double tempsum = 0;
			for(int j = xind1[i]; j < xind2[i]; j++) {
//				Rprintf("Adding read %d at position %d, distance to start %d, extra coverage %f\n", j, i, xstarts[j]-curcpg, xfragdist[xstarts[j]-curcpg]);
				tempsum = tempsum + xfragdist[xstarts[j]-curcpg];
			}
			xout[i] = xout[i] + tempsum;
		}
	}
	
	UNPROTECT(6);
	return(R_NilValue);
}

static R_CallMethodDef callMethods[] = {
	{"CrowSumsSq", (DL_FUNC) &CrowSumsSq, 2},
	{"CcolSumsSq", (DL_FUNC) &CcolSumsSq, 2},
	{"cover_frw_c", (DL_FUNC) &cover_frw_c, 6},
	{"cover_rev_c", (DL_FUNC) &cover_rev_c, 6},
	{NULL, NULL, 0}
};

void R_init_ramwas(DllInfo *info)	{
	// Rprintf("Registering AAA\n");
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
