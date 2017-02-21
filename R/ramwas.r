# Replace "x" with "replacement" if x is NULL
.notnull = function(x, replacement){
    if(is.null(x)){ replacement }else{ x }
}

# addition of vectors of mismatching length
`%add%` = function(x, y){
    if(is.null(x))
        return(y);
    if(is.null(y))
        return(x);
    if( length(x) > length(y) ){
        length(y) = length(x);
        y[is.na(y)] = 0;
    } else {
        length(x) = length(y);
        x[is.na(x)] = 0;
    }
    return(x + y);
}

# chech if path is absolute
# (common alternatives are flawed)
isAbsolutePath = function( path ){
    if( path == "~" )
        return(TRUE);
    if( grepl("^~/", path) )
        return(TRUE);
    if( grepl("^.:(/|\\\\)", path) )
        return(TRUE);
    if( grepl("^(/|\\\\)", path) )
        return(TRUE);
    return(FALSE);
}

# Get full path to the "filename" assuming current directory is "path"
makefullpath = function(path, filename){
    if( is.null(path) )
        return(filename);
    if(isAbsolutePath(filename)){
        return(filename)
    } else {
        return( paste0(path, "/", filename) );
    }
}

# Delete file, no warning on NULL or missing file
.file.remove = function(x){
    if( !is.null(x) )
        if( file.exists(x) )
            file.remove(x);
}

# Find best N p-values, in unsorted vector
findBestNpvs = function(pv, n){
    if(n < 1)
        return(which(pv <= n));
    if(n > length(pv))
        stop('n > length(pv) in findBestNpvs() call');
    if(n == length(pv))
        return(seq_len(n));
    
    # Thresholds are chosen a reasonable for p-value input
    pvthr = c(10^((-100):0), .Machine$double.xmax);
    fi = findInterval(pv, pvthr);
    tab = cumsum(tabulate(fi));
    # Minimum threshold below which there are at least n p-values
    upperfi = which(tab >= n)[1];
    set1 = which(fi <= upperfi);
    cpgsetraw = set1[sort.list(pv[set1])[seq_len(n)]];
    cpgset = sort.int(cpgsetraw);
    return(cpgset);
}

# Standard BH p-value to q-value calculation
pvalue2qvalue = function(pv, n = length(pv)){
    ord = sort.list(pv);
    FDR = pv[ord] * n / seq_along(pv);
    FDR[length(FDR)] = min(FDR[length(FDR)], 1);
    FDR = rev(cummin(rev(FDR)));

    rez = double(length(pv));
    rez[ord] = FDR;
    return(rez)
}

# Caching environment
.ramwasEnv = new.env()

# Load an RDS file and cache it
# or load from cache
cachedRDSload = function(rdsfilename){
    if(is.null(rdsfilename))
        return(NULL);
    cachename = rdsfilename; #paste0(".ramwas.",rdsfilename);
    if( exists(x = cachename, envir = .ramwasEnv) ){
        # message("cachedRDSload: Using cache for: ", rdsfilename);
        return(base::get(x = cachename, envir = .ramwasEnv));
    } else {
        # message("cachedRDSload: Loading to cache: ", rdsfilename);
        data = readRDS(rdsfilename);
        base::assign(x = cachename, value = data, envir = .ramwasEnv);
        return(data);
    }
}

# Orthonormalize a set of covariates
orthonormalizeCovariates = function(cvrt, modelhasconstant = TRUE){
    if(any(sapply(lapply(cvrt, is.na), any)))
        stop("Missing values are not allowed in the covariates")
    if(modelhasconstant) {
        cvrtset = c(const = list(rep(1, nrow(cvrt))), cvrt);
    } else {
        cvrtset = cvrt;
    }
    if(is.list(cvrtset)) {
        factorset = which(sapply(cvrtset, class) %in% c("character","factor"));
        for( ind in factorset ){ # ind = 3
            fctr = factor(cvrtset[[ind]]);
            cvrtset[[ind]] = model.matrix(~fctr)[,-1];
            rm(fctr);
        }
        cvrtmat = matrix(unlist(cvrtset), nrow = nrow(cvrt));
    } else {
        cvrtmat = cvrtset;
    }
    cvrtqr = qr.Q(qr(cvrtmat));  ### tcrossprod(cvrtqr) - diag(nrow(cvrtqr))
    return(cvrtqr)
}

# find how samples in CpG score matrix
# match those in "covariates" parameter
# get the total number of CpGs along the way
.matchCovmatCovar = function( param ){
    cvsamples = param$covariates[[1]];

    fm = fm.open( paste0(param$dircoveragenorm, "/Coverage"), readonly = TRUE);
    fmsamples = rownames(fm);
    ncpgs = ncol(fm);
    close(fm);

    if(is.null(cvsamples))
        cvsamples = fmsamples;
    
    rowsubset = match(cvsamples, fmsamples, nomatch = 0L);
    if( any(rowsubset==0) )
        stop( paste("Unknown samples in covariate file:",
                    cvsamples[head(which(rowsubset==0))]) );

    if( length(cvsamples) == length(fmsamples) ){
        if( all(rowsubset == seq_along(rowsubset)) ){
            rowsubset = NULL;
        }
    }
    return(list(rowsubset = rowsubset, ncpgs = ncpgs, samplenames = fmsamples));
}

# Get covariates + PCs matrix for analysis
# orthonormalized unless normalize == FALSE
.getCovariates = function(param, rowsubset = NULL, normalize = TRUE, modelhasconstant){
    cvrtqr = param$covariates[ param$modelcovariates ];
    ### Reading PCs, add as coveriates
    if( param$modelPCs > 0 ){
        e = readRDS(paste0(param$dirpca,"/eigen.rds"));
        PCs = e$vectors[, seq_len(param$modelPCs), drop=FALSE];
        if(!is.null( rowsubset ))
            PCs = PCs[rowsubset,];
        mwascvrtqr = cbind(cvrtqr, PCs);
        rm(e);
    } else {
        mwascvrtqr = cvrtqr;
    }
    # stopifnot( all.equal( tcrossprod(mwascvrtqr), diag(nrow(mwascvrtqr))) );
    if(normalize){
        rez = t(orthonormalizeCovariates(mwascvrtqr, modelhasconstant));
    } else {
        if(modelhasconstant) {
            rez = t(cbind(rep(1, nrow(mwascvrtqr)),mwascvrtqr));
        } else {
            rez = t(mwascvrtqr); #;
        }
    }
    return(rez);
}

