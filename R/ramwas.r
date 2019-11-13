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
        stop("n > length(pv) in findBestNpvs() call");
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
    if( length(pv) == 0 )
        return(pv);
    if( is.unsorted(pv) ){
        ord = sort.list(pv, decreasing = TRUE);
        FDR = pv[ord] * n / (length(pv):1);
        FDR[1] = min(FDR[1], 1);
        FDR = cummin(FDR);
    
        rez = double(length(pv));
        rez[ord] = FDR;
        return(rez);
    } else {
        FDR = pv * n / seq_along(pv);
        FDR[length(FDR)] = min(FDR[length(FDR)], 1);
        FDR = rev(cummin(rev(FDR)));
        return(FDR);
    }
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
    # Prevent missing values
    if(any(vapply(lapply(cvrt, is.na), any, TRUE)))
        stop("Missing values are not allowed in the covariates");
    
    # Add a constant?
    if(modelhasconstant){
        cvrtset = c(const = list(rep(1, nrow(cvrt))), cvrt);
    } else {
        cvrtset = cvrt;
    }
    
    # Transform factor covariates into dummies, kill zero covariates
    if(is.list(cvrtset)){
        cvrtset = as.list(cvrtset);
        isfactorset = vapply(cvrtset, class, "") %in% c("character","factor");
        for( ind in seq_along(isfactorset) ){ # ind = 1
            if(isfactorset[ind]){
                fctr = factor(cvrtset[[ind]]);
                if(nlevels(fctr) >= 2) {
                    cvrtset[[ind]] = model.matrix(~fctr)[,-1];
                } else {
                    cvrtset[[ind]] = numeric(0);
                }
                rm(fctr);
            } else {
                # Kill pure zero covariates
                if(all(cvrtset[[ind]] == 0))
                    cvrtset[[ind]] = numeric(0);
            }
        }
        cvrtmat = matrix(unlist(cvrtset), nrow = nrow(cvrt));
    } else {
        cvrtmat = cvrtset;
    }
    
    # Orthonormalize the covariates
    cvrtqr = qr.Q(qr(cvrtmat));  ### tcrossprod(cvrtqr) - diag(nrow(cvrtqr))
    return(cvrtqr)
}

# Get covariates + PCs matrix for analysis
# orthonormalized unless normalize == FALSE
.getCovariates = function(
        param,
        rowsubset = NULL,
        normalize = TRUE,
        modelhasconstant){
    # Named covariates
    cvrt = param$covariates[ param$modelcovariates ];

    ### Add PCs as covariates
    if( param$modelPCs > 0 ){
        # e = readRDS(paste0(param$dirpca,"/eigen.rds"));
        eigenvectors = fm.open(
                filenamebase = paste0(param$dirpca, "/eigenvectors"),
                readonly = TRUE);
        PCs = eigenvectors[, seq_len(param$modelPCs)]; # , drop=FALSE
        close(eigenvectors);
        
        if(!is.null( rowsubset ))
            PCs = PCs[rowsubset,];
        mwascvrtqr = cbind(cvrt, PCs);
        # rm(e);
    } else {
        mwascvrtqr = cvrt;
    }
    # stopifnot( all.equal( tcrossprod(mwascvrtqr), diag(nrow(mwascvrtqr))) );
    if(normalize){
        rez = t(orthonormalizeCovariates(mwascvrtqr, modelhasconstant));
    } else {
        if(modelhasconstant){
            rez = t(cbind(rep(1, nrow(mwascvrtqr)),mwascvrtqr));
        } else {
            rez = t(mwascvrtqr); #;
        }
    }
    return(rez);
}


.set1MLKthread = "if(\"package:RevoUtilsMath\" %in% search())
    if(exists(\"setMKLthreads\", 
    where = \"package:RevoUtilsMath\"))
    RevoUtilsMath::setMKLthreads(1);";
# .set1MLKthread = function(){
#     if("package:RevoUtilsMath" %in% search())
#         if(exists("setMKLthreads", where = "package:RevoUtilsMath"))
#             RevoUtilsMath::setMKLthreads(1);
# }

# The logging function
.log = function(ld, fmt, ..., append = TRUE){
    msg = sprintf(fmt, ...);
    cat(file = paste0(ld, "/Log.txt"),
        msg, "\n", sep = "", append = append);
    message(msg);
    return(invisible(msg));
}

.logErrors = function(ld, fun){
    rez = function(...){
        withCallingHandlers(
            tryCatch(fun(...),
                error = function(e){
                    .log(ld, "%s, Process %06d\nError: %s\nIn: %s",
                        date(), Sys.getpid(),
                        conditionMessage(e), deparse(conditionCall(e)));
                }
            ),
            warning = function(w){
                .log(ld, "%s, Process %06d\nWarning: %s\nIn: %s",
                    date(), Sys.getpid(),
                    conditionMessage(w), deparse(conditionCall(w)));
                invokeRestart("muffleWarning");
            }
        )
    };
    env = new.env(parent = baseenv());
    assign("ld", ld, envir = env);
    assign(".log", .log, envir = env);
    assign("fun", fun, envir = env);
    environment(rez) = env;
    return(rez);
}

.showErrors = function(z){
    for( x in z )
        if(is.character(x))
            message(x);
}

# .logErrors = function(ld, fun)fun;

mat2cols = function(x){
    x = as.matrix(x);
    return(lapply(seq_len(ncol(x)), function(i)x[,i]));
};
        
