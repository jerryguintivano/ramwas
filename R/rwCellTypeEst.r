estimateCellTypes = function(param, ncpgs = NULL){
    pf = parameterPreprocess(param);
    
    if(is.null(ncpgs))
        ncpgs = pf$mmncpgs[1];
    
    message("Working in: ", pf$dirmwas);
    
    # Indices of reference samples
    ctindices = split.default(
        x = seq_len(nrow(pf$covariates)), 
        f = pf$covariates[ pf$modeloutcome ]);
    mycelltypes = names(ctindices);
    
    # Select CpGs to use
    fm = fm.open(paste0(pf$dirmwas, "/Stats_and_pvalues"));
    pv = fm[,3];
    close(fm)
    idset = findBestNpvs(pv, ncpgs);
    
    # Get covariates
    cvrtqr = t(.getCovariates(pf, modelhasconstant = TRUE));
    
    # Initialize loop
    fm = fm.open(paste0(pf$dircoveragenorm, "/Coverage"), readonly = TRUE)
    XTX = 0;
    XY = 0;
    # YY = 0;
    
    # main loop
    step1 = 1000;
    runto = ncpgs;
    nsteps = ceiling(runto/step1);
    for( part in seq_len(nsteps) ) { # part = 1
        fr = (part-1)*step1 + 1;
        to = min(part*step1, runto);
        message("Processing CpGs ", fr, " to ", to, " of ", runto);
        
        slice = fm[,idset[fr:to]];
        slice = slice - tcrossprod(cvrtqr, crossprod(slice, cvrtqr));
        
        mean1 = matrix(0, ncol(slice), length(mycelltypes));
        # mean2 = matrix(0, ncol(slice), length(mycelltypes));
        # colnames(mean1) = mycelltypes;
        # colnames(mean2) = mycelltypes;
        for( cti in seq_along(mycelltypes) ){ # cti = 1
            subdata = slice[ctindices[[cti]], , drop = FALSE];
            mean1[, cti] = colMeans(subdata);
            # mean2[, cti] = colSumsSq(subdata)/nrow(subdata);
        }
        
        lastCT =  mean1[,length(mycelltypes)];
        X = mean1[,-length(mycelltypes)] - lastCT;
        XTX = XTX + crossprod(X);
        XY = XY + (t(slice %*% X) - as.vector(crossprod(X, lastCT)));
        # YY = YY + rowSumsSq(slice - rep(lastCT, each = nrow(slice)));
        # YY = YY + (rowSumsSq(slice) - 
        #            2 * (slice %*% lastCT) + as.vector(crossprod(lastCT)));
        
    }
    rm(X, lastCT, cti, slice, subdata)
    rm(part, step1, runto, nsteps, fr, to);
    
    beta = solve(XTX, XY);
    fullbeta = cbind(t(beta), 1 - colSums(beta))
    
    colnames(fullbeta) = mycelltypes;
    rownames(fullbeta) = rownames(fm);
    
    saveRDS(
        file = sprintf("%s/CellTypeEstimatesR_%07d_CpGs.rds", pf$dirmwas,ncpgs), 
        object = fullbeta);
    write.table(
        file = sprintf("%s/CellTypeEstimatesT_%07d_CpGs.txt", pf$dirmwas,ncpgs), 
        x = data.frame(samples = rownames(fullbeta), fullbeta),
        row.names = FALSE,
        sep = "\t");
    
    close(fm);
    
    return(invisible(fullbeta));
}
