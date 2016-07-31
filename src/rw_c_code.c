#include <R.h>
#include <Rinternals.h>

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
