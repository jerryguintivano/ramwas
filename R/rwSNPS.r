# find how samples in CpG score matrix
# match those in "covariates" parameter
# get the total number of CpGs along the way
.matchCovmatCovar = function( param ){

    # Sample names in covariates
    cvsamples = param$covariates[[1]];

    # Grab info form coverage matrix
    fm = fm.open( paste0(param$dircoveragenorm, "/Coverage"), readonly = TRUE);
    fmsamples = rownames(fm);
    ncpgs = ncol(fm);
    nsamplesall = nrow(fm);
    close(fm);

    # Match samples in covariates with those in coverage matrix
    if(is.null(cvsamples)){
        # if covariates are not set, assume they match.
        rowsubset = NULL;
    } else {
        rowsubset = match(cvsamples, fmsamples, nomatch = 0L);
        if( any(rowsubset==0) )
            stop( paste("Unknown samples in covariate file:",
                        cvsamples[head(which(rowsubset==0))]) );

        # if no reordering is required, set rowsubset=NULL
        if( length(cvsamples) == length(fmsamples) ){
            if( all(rowsubset == seq_along(rowsubset)) ){
                rowsubset = NULL;
            }
        }
    }
    return(list(
        rowsubset = rowsubset,
        ncpgs = ncpgs,
        cvsamples = cvsamples, nsamplesall = nsamplesall));
}

# Save top findings in a text file
# with annotation
ramwas5saveTopFindingsSNPs = function(param){
    # library(filematrix)
    param = parameterPreprocess(param);

    message("Working in: ", param$dirSNPs);

    message("Loading MWAS results");
    mwas = fm.load( paste0(param$dirSNPs, "/Stats_and_pvalues") );

    message("Loading CpG locations");
    cpgloc = fm.load(
        filenamebase = paste0(param$dircoveragenorm, "/CpG_locations") );
    chrnames = readLines(
        con = paste0(param$dircoveragenorm, "/CpG_chromosome_names.txt") );

    message("Finding top MWAS hits");
    keep = findBestNpvs(mwas[,3], param$toppvthreshold);
    # keep = which(mwas[,3] < param$toppvthreshold);
    ord = keep[sort.list(abs(mwas[keep,2]),decreasing = TRUE)];

    toptable = data.frame( 
                    chr = chrnames[cpgloc[ord,1]],
                    position =     cpgloc[ord,2],
                    Ftest  = mwas[ord,2],
                    pvalue = mwas[ord,3],
                    qvalue = mwas[ord,4]);

    # saveRDS(file = paste0(param$dirSNPs,"/Top_tests.rds"), object = toptable);


    message("Saving top MWAS hits");
    write.table(
        file = paste0(param$dirSNPs,"/Top_tests.txt"),
        sep = "\t", quote = FALSE, row.names = FALSE,
        x = toptable
    );
    return(invisible(NULL));
}

# Job function for MWAS
.ramwasSNPsJob = function(rng, param, mwascvrtqr, rowsubset){

    testPhenotypeSNPs = function(phenotype, slice0, cvrtqr0, snps0){
        mycov = matrix(phenotype, nrow = 1);
        slice = slice0;
        snpss = snps0;
        cvrtq = cvrtqr0;

        nsamples = nrow(slice);

        if( any(is.na(mycov)) ){
            keep = which(colSums(is.na(mycov))==0);

            mycov = mycov[, keep, drop=FALSE];
            slice = slice[keep, , drop=FALSE];
            snpss = snpss[keep, , drop=FALSE];
            cvrtq = cvrtq[, keep, drop=FALSE];
            cvrtq = t( qr.Q(qr(t(cvrtq))) );
            rm(keep);
        }


        # methlation
        # slice

        # outcome variable
        varot = rep(mycov, ncol(slice));
        dim(varot) = dim(slice);

        # outcome*SNP matrix
        matos = snpss * varot;

        # pure SNPs
        # snpss


        # Rezidualize w.r.t. covariates
        slice = slice - crossprod(cvrtq, cvrtq %*% slice);
        varot = varot - crossprod(cvrtq, cvrtq %*% varot);
        matos = matos - crossprod(cvrtq, cvrtq %*% matos);
        snpss = snpss - crossprod(cvrtq, cvrtq %*% snpss);

        # Rezidualize w.r.t. SNPs
        snpss = snpss / rep(pmax(sqrt(colSums(snpss^2)), 1e-16),
                            each = nsamples);

        slice = slice - snpss * rep(colSums(slice*snpss), each = nsamples);
        varot = varot - snpss * rep(colSums(varot*snpss), each = nsamples);
        matos = matos - snpss * rep(colSums(matos*snpss), each = nsamples);

        # Normalize

        slice = slice / rep(pmax(sqrt(colSums(slice^2)), 1e-16),
                            each = nsamples);
        varot = varot / rep(pmax(sqrt(colSums(varot^2)), 1e-16),
                            each = nsamples);
        matos = matos / rep(pmax(sqrt(colSums(matos^2)), 1e-16),
                            each = nsamples);

        cr10 = colSums(slice*varot);
        cr20 = colSums(slice*matos);
        cr12 = colSums(varot*matos);

        R2 = (cr10^2 + cr20^2 - 2*cr12*cr10*cr20) / (1 - cr12^2);

        # hist(R2,100);

        ###
        nVarTested = 2;
        dfFull = ncol(cvrtq) - nrow(cvrtq) - nVarTested - 1; # One for SNPs

        if(dfFull <= 0)
            return(list(
                    Rsquared = 0,
                    Fstat = 0,
                    pvalue = 1,
                    nVarTested = nVarTested,
                    dfFull = dfFull,
                    statname = paste0("-F_",nVarTested)));

        rsq2F = function(x){
            return( x / (1 - pmin(x,1)) * (dfFull/nVarTested) );
        }
        F2pv = function(x){
            return( pf(x, nVarTested, dfFull, lower.tail = FALSE) );
        }
        ff = rsq2F(R2);
        pv = F2pv(ff);

        ### Check
        # lm1 = lm(slice[,1] ~ phenotype + phenotype*snps0[,1] +
        #                            snps0[,1] + t(cvrtq)[,-1])
        # lm0 = lm(slice[,1] ~
        #                            snps0[,1] + t(cvrtq)[,-1])
        # anova(lm0,lm1);
        # pv[1]
        return(list(
                Rsquared = R2,
                Fstat = ff,
                pvalue = pv,
                nVarTested = nVarTested,
                dfFull = dfFull,
                statname = paste0("-F_",nVarTested)) );
    }


    # rng = rangeset[[1]];
    # library(filematrix);
    fmm = fm.open( 
                filenamebase = paste0(param$dircoveragenorm, "/Coverage"),
                readonly = TRUE,
                lockfile = param$lockfile2);

    fms = fm.open(
                filenamebase = paste0(param$fileSNPs),
                readonly = TRUE,
                lockfile = param$lockfile2);

    outmat = double(3*(rng[2]-rng[1]+1));
    dim(outmat) = c((rng[2]-rng[1]+1),3);

    step1 = ceiling( 512*1024*1024 / nrow(fmm) / 8);
    mm = rng[2]-rng[1]+1;
    nsteps = ceiling(mm/step1);
    for( part in seq_len(nsteps) ){ # part = 1
        message("Slice ", part, " of ", nsteps);
        fr = (part-1)*step1 + rng[1];
        to = min(part*step1, mm) + rng[1] - 1;

        slice = fmm[,fr:to];
        snps =  fms[,fr:to];
        if( !is.null(rowsubset) ){
            slice = slice[rowsubset,];
            snps  = snps[ rowsubset,];
        }

        rez = testPhenotypeSNPs(
            phenotype = param$covariates[[param$modeloutcome]],
            slice0 = slice,
            cvrtqr0 = mwascvrtqr,
            snps0 = snps)

        outmat[(fr:to) - (rng[1] - 1),] = cbind(rez[[1]], rez[[2]], rez[[3]]);

        fmm$filelock$lockedrun( {
            cat(file = paste0(param$dirSNPs,"/Log.txt"),
                date(), ", Process ", Sys.getpid(), ", Job ", rng[3],
                ", processing slice ", part, " of ", nsteps, "\n",
                sep = "", append = TRUE);
        });
        rm(slice, snps);
    }
    close(fmm);
    close(fms);

    if(rng[1]==1){
        writeLines( 
            con = paste0(param$dirSNPs, "/DegreesOfFreedom.txt"),
            text = as.character(c(rez$nVarTested, rez$dfFull)))
    }

    fmout = fm.open(
                filenamebase = paste0(param$dirSNPs, "/Stats_and_pvalues"),
                lockfile = param$lockfile2);
    fmout[rng[1]:rng[2],seq_len(3)] = outmat;
    close(fmout);

    return("OK");
}

# Setp 5 of RaMWAS
ramwasSNPs = function( param ){
    # library(filematrix)
    param = parameterPreprocess(param);

    # Fix parameters

    dir.create(param$dirSNPs, showWarnings = FALSE, recursive = TRUE);

    parameterDump(dir = param$dirSNPs, param = param,
        toplines = c(   "fileSNPs", "dirSNPs",
                        "dirpca", "dircoveragenorm",
                        "filecovariates", "covariates",
                        "modeloutcome", "modelcovariates",
                        "modelPCs",
                        "qqplottitle",
                        "diskthreads"));


    message("Preparing for MWAS with SNPs");

    ### Kill NA outcomes
    if(any(is.na(param$covariates[[ param$modeloutcome ]]))){
        param$covariates = data.frame(lapply(
            X = param$covariates,
            FUN = `[`,
            !is.na(param$covariates[[ param$modeloutcome ]])));
    }
    ### Get and match sample names
    {
        message("Matching samples in covariates and data matrix");
        rez = .matchCovmatCovar( param );
        # rez = ramwas:::.matchCovmatCovar( param );
        rowsubset = rez$rowsubset;
        ncpgs     = rez$ncpgs;
        cvsamples = rez$cvsamples;
        rm(rez);
    } # rowsubset, ncpgs

    ### Prepare covariates, defactor
    {
        message("Preparing covariates (splitting dummies, orthonormalizing)");
        mwascvrtqr = .getCovariates(param = param,
                                    rowsubset = rowsubset,
                                    modelhasconstant = param$modelhasconstant);
        # mwascvrtqr = ramwas:::.getCovariates(param, rowsubset);
    } # mwascvrtqr


    ### Outpout matrix. Cor / t-test / p-value / q-value
    ### Outpout matrix. R2  / F-test / p-value / q-value
    {
        message("Creating output matrix");
        fm = fm.create( 
                    filenamebase = paste0(param$dirSNPs, "/Stats_and_pvalues"),
                    nrow = ncpgs,
                    ncol = 4);
        if( !is.character( param$covariates[[param$modeloutcome]] ) ){
            # colnames(fm) = c("cor","t-test","p-value","q-value");
            colnames(fm) = c("R-squared","F-test","p-value","q-value");
        } else {
            stop("Categorical outcome in not supported for analysis with SNPs")
        }
        close(fm);
    }

    ### Running MWAS in parallel
    {
        message("Running MWAS");
        cat(file = paste0(param$dirSNPs,"/Log.txt"),
            date(), ", Running methylome-wide association study with SNPs.",
            "\n",
            sep = "", append = FALSE);

        step1 = ceiling( 512*1024*1024 / length(cvsamples) / 8);
        mm = ncpgs;
        nsteps = ceiling(mm/step1);

        nthreads = min(param$diskthreads, nsteps);
        rm(step1, mm, nsteps);
        if( nthreads > 1 ){
            rng = round(seq(1, ncpgs+1, length.out = nthreads+1));
            rangeset = rbind( 
                            rng[-length(rng)],
                            rng[-1]-1,
                            seq_len(nthreads));
            rangeset = lapply(
                            X = seq_len(ncol(rangeset)),
                            FUN = function(i)rangeset[,i]);

            if(param$usefilelock) param$lockfile2 = tempfile();
            # library(parallel);
            cl = makeCluster(nthreads);
            on.exit({
                stopCluster(cl);
                .file.remove(param$lockfile2);});
            # clusterExport(cl, "testPhenotypeSNPs", envir = )
            clusterApplyLB(
                        cl = cl, 
                        x = rangeset,
                        fun = .ramwasSNPsJob,
                        param = param,
                        mwascvrtqr = mwascvrtqr,
                        rowsubset = rowsubset);
            tmp = sys.on.exit();
            eval(tmp);
            rm(tmp);
            on.exit();
            rm(cl, rng, rangeset);
        } else {
            covmat = .ramwasSNPsJob(
                        rng = c(1, ncpgs, 0),
                        param = param, 
                        mwascvrtqr = mwascvrtqr, 
                        rowsubset = rowsubset);
        }
        cat(file = paste0(param$dirSNPs,"/Log.txt"),
            date(), ", Done running methylome-wide association study.", "\n",
            sep = "", append = TRUE);
    }

    ### Fill in FDR column
    {
        message("Calculating FDR (q-values)");
        fm = fm.open( paste0(param$dirSNPs, "/Stats_and_pvalues"));
        pvalues = fm[,3];
        pvalues[pvalues==0] = .Machine$double.xmin;
        fm[,4] = pvalue2qvalue( pvalues );
        close(fm);

        rm(fm);
    } # sortedpv

    ### QQ-plot
    {
        if( is.null(param$qqplottitleSNPs) ){
            qqplottitleSNPs = paste0(
                    "Testing ", param$modeloutcome,
                    " with SNPs\n",
                    param$modelPCs,
                    " PC",if(param$modelPCs!=1)"s"else"");
            if(length(param$modelcovariates)>0)
                qqplottitleSNPs = paste0(
                    qqplottitleSNPs, " and ",
                    length(param$modelcovariates)," covariate",
                    if(length(param$modelcovariates)!=1)
                        "s:\n"else": ",
                    paste0(param$modelcovariates,collapse = ", "))
            param$qqplottitleSNPs = qqplottitleSNPs;
            rm(qqplottitleSNPs);
        }

        message("Creating QQ-plot");
        pdf(paste0(param$dirSNPs, "/QQ_plot.pdf"),7,7);
        qqPlotFast(pvalues);
        title(param$qqplottitleSNPs);
        dev.off();
    }
    ramwas5saveTopFindingsSNPs(param);
    return(invisible(NULL));
}
