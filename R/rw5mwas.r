# test a phenotype against the data
# accounting for covariates
# cvrtqr - must be orthonormalized
# Supports categorical outcomes (text of factor)
testPhenotype = function(phenotype, data, cvrtqr){
    # Preserve original variables for debugging
    pheno = matrix(phenotype, ncol = 1);
    slice = data;
    cvqr0 = cvrtqr;

    # In rare case of missing phenotype values
    # Subset and re-orthonormalize covariates
    if( any(is.na(pheno)) ){
        keep = which(rowSums(is.na(pheno))==0);

        pheno = pheno[keep, , drop = FALSE];
        slice = slice[keep, , drop = FALSE];
        cvqr0 = cvqr0[keep, , drop = FALSE];
        cvqr0 = qr.Q(qr(cvqr0));
    }

    # Process factor phenotypes
    # Othonormalize them
    # Exclude redundant
    if( is.character(pheno) || is.factor(pheno) ){
        fctr = as.factor(pheno);
        # Single factor phenotype
        # Occur by mistake and
        # after heavy subsetting
        if( length(levels(fctr)) <= 1 ){
            return(list(
                    Rsquared = 0,
                    Fstat = 0,
                    pvalue = 1,
                    nVarTested = 0,
                    dfFull = ncol(cvqr0) - nrow(cvqr0),
                    statname = paste0("-F_",0)));
        }
        dummy = model.matrix(~fctr)[,-1];
        dummy = dummy - cvqr0 %*% crossprod(cvqr0, dummy);

        q = qr(dummy);
        keep = abs(diag(qr.R(q))) > (.Machine$double.eps * length(pheno) * 10);

        pheno = qr.Q(qr(dummy));
        if( !all(keep) )
            pheno = pheno[, keep, drop = FALSE];
    } else {
        cvsumsq1 = sum( pheno^2 );
        pheno = pheno - cvqr0 %*% crossprod(cvqr0, pheno);
        cvsumsq2 = sum( pheno^2 );
        if( cvsumsq2 <= (cvsumsq1 * .Machine$double.eps * ncol(pheno)) ){
            return(list(
                    correlation = 0,
                    tstat = 0,
                    pvalue = 1,
                    nVarTested = 0,
                    dfFull = 0,
                    statname = ""));
        } else {
            pheno = pheno / sqrt(sum(pheno^2));
        }
    }

    ###
    nVarTested = ncol(pheno);
    dfFull = nrow(pheno) - ncol(cvqr0) - nVarTested;
    if( nVarTested == 1 ){
        if( dfFull <= 0 )
            return(list(
                    correlation = 0,
                    tstat = 0,
                    pvalue = 1,
                    nVarTested = nVarTested,
                    dfFull = dfFull,
                    statname = ""));

        
        # SST - Total variation
        # cvD^2 - Variation explained by covariates
        # cvC^2 - Variation explained by phenotype
        # cr = cvD / sqrt(SST - cvC^2) with precautions
        
        # If data is residualized
        # cr = cvD / sqrt(SST)
        
        # SST = rowSums(slice^2);
        SST = colSumsSq(slice);

        cvD = crossprod(pheno, slice);
        # cvD2 = colSums(cvD^2);
        # cvD2 = cvD2^2;

        cvC = crossprod(cvqr0, slice);
        cvC2 = colSumsSq( cvC );
        # cvC2 = colSums( cvC^2 );

        # SSR = colSums( cvD^2 );
        cr = cvD / sqrt(pmax(SST - cvC2, 1e-50, SST/1e15));

        cor2tt = function(x){
            return( x * sqrt(dfFull / (1 - pmin(x^2,1))));
        }
        tt2pv = function(x){
            return( (pt(-abs(x), dfFull)*2) );
        }
        tt = cor2tt(cr);
        pv = tt2pv(tt);

        ### Check:
        # c(tt[1], pv[1])
        # summary(lm( as.vector(covariate) ~ 0 + data[1,] +
        #         t(cvrtqr)))$coefficients[1,]
        return(list(
                correlation = cr,
                tstat = tt,
                pvalue = pv,
                nVarTested = nVarTested,
                dfFull = dfFull,
                statname = ""));
    } else {
        if(dfFull <= 0)
            return(list(
                    Rsquared = 0,
                    Fstat = 0,
                    pvalue = 1,
                    nVarTested = nVarTested,
                    dfFull = dfFull,
                    statname = paste0("-F_",nVarTested)) );

        # SST - Total variation
        # cvD2 - Variation explained by covariates
        # cvC2 - Variation explained by phenotype
        # cr  = cvD / sqrt(SST - cvC^2) with precautions
        # rsq = cvD2 / (SST - cvC2)
        # SST = rowSums(slice^2);
        SST = colSumsSq(slice);

        cvD = crossprod(pheno, slice);
        # cvD2 = colSums(cvD^2);
        cvD2 = colSumsSq(cvD);

        cvC = crossprod(cvqr0, slice);
        # cvC2 = colSums( cvC^2 );
        cvC2 = colSumsSq( cvC );

        # SSR = colSums( cvD^2 );
        rsq = cvD2 / pmax(SST - cvC2, 1e-50, SST/1e15);

        # rsq = colSums(cr^2);
        rsq2F = function(x){
            return( x / (1 - pmin(x,1)) * (dfFull/nVarTested) );
        }
        F2pv = function(x){
            return( pf(x, nVarTested, dfFull, lower.tail = FALSE) );
        }
        ff = rsq2F(rsq);
        pv = F2pv(ff);

        ### Check:
        # c(ff[1], pv[1])
        # anova(lm( data[1,] ~ 0 + t(cvrtqr) +
        #     as.factor(as.vector(as.character(round(covariate))))))
        return(list(
                Rsquared = rsq,
                Fstat = ff,
                pvalue = pv,
                nVarTested = nVarTested,
                dfFull = dfFull,
                statname = paste0("-F_",nVarTested)));
    }
}



# Save top findings in a text file
# with annotation
ramwas5saveTopFindings = function(param){
    # library(filematrix)
    param = parameterPreprocess(param);

    message("Working in: ", param$dirmwas);

    message("Loading MWAS results");
    mwas = getMWASandLocations(param);
    
    message("Finding top MWAS hits");
    keep = findBestNpvs(mwas$`p-value`, param$toppvthreshold);
    # keep = which(mwas[,3] < param$toppvthreshold);
    ord = keep[sort.list(abs(mwas[[5]][keep]), decreasing = TRUE)];

    toptable = mwas[ord,];

    message("Saving top MWAS hits");
    write.table(
        file = paste0(param$dirmwas, "/Top_tests.txt"),
        x = toptable,
        sep = "\t", quote = FALSE, row.names = FALSE
    );
    return(invisible(NULL));
}

# Job function for MWAS
.ramwas5MWASjob = function(rng, param){
    # rng = rangeset[[1]];
    # library(filematrix);
    ld = param$dirmwas;

    .log(ld, "%s, Process %06d, Job %02d, Start MWAS, CpG range %d-%d",
        date(), Sys.getpid(), rng[3], rng[1], rng[2]);

    # Get data access
    data = new("rwDataClass", param = param, lockfile = param$lockfile2);

    outmat = double(3*(rng[2]-rng[1]+1));
    dim(outmat) = c((rng[2]-rng[1]+1),3);

    step1 = ceiling( 128*1024*1024 / data$ndatarows / 8);
    mm = rng[2]-rng[1]+1;
    nsteps = ceiling(mm/step1);
    for( part in seq_len(nsteps) ){ # part = 1
        .log(ld, "%s, Process %06d, Job %02d, Processing slice: %03d of %d",
            date(), Sys.getpid(), rng[3], part, nsteps);
        fr = (part-1)*step1 + rng[1];
        to = min(part*step1, mm) + rng[1] - 1;

        slice = data$getDataRez(fr:to, resid = FALSE);
        
        rez = testPhenotype(
                    phenotype = param$covariates[[param$modeloutcome]],
                    data = slice,
                    cvrtqr = data$cvrtqr)

        outmat[(fr:to) - (rng[1] - 1), ] = cbind(rez[[1]], rez[[2]], rez[[3]]);

        rm(slice);
    }
    data$close();

    if( rng[1] == 1 ){
        writeLines( 
                con = paste0(param$dirmwas, "/DegreesOfFreedom.txt"),
                text = as.character(c(rez$nVarTested, rez$dfFull)));
    }

    fmout = fm.open(
                filenamebase = paste0(param$dirmwas, "/Stats_and_pvalues"),
                lockfile = param$lockfile2);
    fmout[rng[1]:rng[2], 1:3] = outmat;
    close(fmout);
    .log(ld, "%s, Process %06d, Job %02d, Done MWAS, CpG range %d-%d",
        date(), Sys.getpid(), rng[3], rng[1], rng[2]);

    return(invisible(NULL));
}

# Setp 5 of RaMWAS
ramwas5MWAS = function( param ){
    # library(filematrix)
    param = parameterPreprocess(param);
    ld = param$dirmwas;
    dir.create(param$dirmwas, showWarnings = FALSE, recursive = TRUE);
    .log(ld, "%s, Start ramwas5MWAS() call", date(), append = FALSE);

    if(is.null(param$modeloutcome))
        stop(
            "Model outcome variable not defined\n",
            "See \"modeloutcome\" parameter");
    if( !any(names(param$covariates) == param$modeloutcome) )
        stop(
            "Model outcome is not found among covariates.\n",
            "See \"modeloutcome\" parameter");
    
    parameterDump(dir = param$dirmwas, param = param,
        toplines = c(   "dirmwas", "dirpca", "dircoveragenorm",
                        "filecovariates", "covariates",
                        "modeloutcome", "modelcovariates",
                        "modelPCs", "modelhasconstant",
                        "qqplottitle",
                        "diskthreads"));


    ### Kill NA outcomes
    killset = is.na(param$covariates[[ param$modeloutcome ]]);
    if(any(killset)){
        .log(ld, "%s, Removing observations with missing outcome", date());
        param$covariates = data.frame(lapply(
            param$covariates,
            `[`,
            !killset), stringsAsFactors = FALSE);
    }
    
    outcome = param$covariates[[ param$modeloutcome ]];
    
    # Get data access
    data = new("rwDataClass", param = param);

    ### Outpout matrix. Cor / t-test / p-value / q-value
    ### Outpout matrix. R2  / F-test / p-value / q-value
    {
        .log(ld, "%s, Creating output matrix", date());
        fm = fm.create( 
                    filenamebase = paste0(param$dirmwas, "/Stats_and_pvalues"),
                    nrow = data$ncpgs,
                    ncol = 4);
        if( is.numeric( outcome ) ){
            colnames(fm) = c("cor","t-test","p-value","q-value");
        } else {
            colnames(fm) = c("R-squared","F-test","p-value","q-value");
        }
        close(fm);
    }

    ### Running MWAS in parallel
    {
        .log(ld, "%s, Start Association Analysis", date());

        step1 = ceiling( 128*1024*1024 / data$ndatarows / 8);
        mm = data$ncpgs;
        nsteps = ceiling(mm/step1);

        nthreads = min(param$cputhreads, nsteps);
        rm(step1, mm, nsteps);
        if( nthreads > 1 ){
            rng = round(seq(1, data$ncpgs+1, length.out = nthreads+1));
            rangeset = rbind( 
                            rng[-length(rng)],
                            rng[-1] - 1,
                            seq_len(nthreads));
            rangeset = mat2cols(rangeset);

            if(param$usefilelock) param$lockfile2 = tempfile();
            # library(parallel);
            cl = makeCluster(nthreads);
            on.exit({
                stopCluster(cl);
                .file.remove(param$lockfile2);
            });
            clusterExport(cl, "testPhenotype");
            logfun = .logErrors(ld, .ramwas5MWASjob);
            clusterExport(  
                        cl = cl, 
                        varlist = ".set1MLKthread", 
                        envir = asNamespace("ramwas"));
            clusterEvalQ(cl, eval(parse(text = .set1MLKthread)));
            z = clusterApplyLB(
                        cl = cl,
                        x = rangeset,
                        fun = logfun,
                        param = param);
            .showErrors(z);
            tmp = sys.on.exit();
            eval(tmp);
            rm(tmp);
            on.exit();
            rm(cl, rng, rangeset);

        } else {
            z = .ramwas5MWASjob(
                        rng = c(1, data$ncpgs, 0),
                        param = param);
        }
        
        .log(ld, "%s, Done Association Analysis", date());
    }
    data$close();

    ### Fill in FDR column
    {
        .log(ld, "%s, Calculating FDR (q-values)", date());
        fm = fm.open(paste0(param$dirmwas, "/Stats_and_pvalues"));
        pvalues = fm[,3];
        pvalues[pvalues==0] = .Machine$double.xmin;
        fm[,4] = pvalue2qvalue( pvalues );
        close(fm);
    } # sortedpv

    ### QQ-plot
    {
        .log(ld, "%s, Creating QQ-plot", date());
        pdf(paste0(param$dirmwas, "/QQ_plot.pdf"), 7, 7);
        qq = qqPlotPrepare(pvalues);
        qqPlotFast(qq, lwd = 1);
        title(param$qqplottitle);
        dev.off();
        saveRDS(file = paste0(param$dirmwas,"/z_QQinfo.rds"), object = qq);
    }
    
    ### Creating QQ-plot and Manhattan plot pair
    locs = getLocations(param);
    if(!is.null(locs)){
        .log(ld, "%s, Creating QQ-plot and Manhattan plot pair", date());
        
        png(
            filename = paste0(param$dirmwas, "/ManPlot.png"), 
            width = 420*9*1.5, 
            height = 420*3.8*1.5, 
            pointsize = 16*4);

        ylim = c(0, max(qq$xpvs[1], qq$ypvs[1])*1.05);
                
        layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(1,2.2));

        qqPlotFast(qq, lwd = 7, ylim = ylim);
        
        man = manPlotPrepare(
                pvalues = pvalues,
                chr = locs$chr, 
                pos = (locs$start + locs$end) %/% 2L);
        
        manPlotFast(man, lwd = 7, ylim = ylim);
        
        dev.off();
        saveRDS(file = paste0(param$dirmwas,"/z_MANinfo.rds"), object = man);
    }
    
    if(!is.null(locs))
        ramwas5saveTopFindings(param);

    .log(ld, "%s, Done ramwas5MWAS() call", date());
    
    return(invisible(NULL));
}
