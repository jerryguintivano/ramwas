# test every covariate agains the data
.testCovariates = function(covariates1, data, cvrtqr){
    # covariates1 = param$covariates[-1]
    # data = t(e$vectors[,seq_len(nonzeroPCs)])
    crF = vector("list", length(covariates1));
    pv  = vector("list", length(covariates1));
    nms = character(length(covariates1));
    for( i in seq_along(covariates1) ){ # i=1
        rez = testPhenotype(covariates1[[i]], data, cvrtqr);
        pv[[i]] = as.vector(rez[[3]]);
        # nms[i] = rez$statname;
        if(nchar(rez$statname)==0){
            nms[i] = "";
            crF[[i]] = as.vector(rez$correlation);
        } else {
            nms[i] = "_R2";
            crF[[i]] = as.vector(rez$Rsquared);
        }
    }
    crF = data.frame(crF);
    names(crF) = paste0(names(covariates1),nms);

    pv = data.frame(pv);
    names(pv) = paste0(names(covariates1),nms);
    return(list(crF = crF, pv = pv));
}

# Job function for PCA analysis
# (covariance matrix calculation)
.ramwas4PCAjob = function(rng, param){
    # library(filematrix);
    ld = param$dirpca;
    
    .log(ld, "%s, Process %06d, Job %02d, Start PCA, CpG range %d-%d",
        date(), Sys.getpid(), rng[3], rng[1], rng[2]);

    # Get data access
    data = new("rwDataClass", param = param, getPCs = FALSE);

    covmat = 0;

    step1 = ceiling( 128*1024*1024 / data$ndatarows / 8);
    mm = rng[2] - rng[1] + 1;
    nsteps = ceiling(mm/step1);
    for( part in seq_len(nsteps) ){ # part = 1
        .log(ld, "%s, Process %06d, Job %02d, Processing slice: %03d of %d",
            date(), Sys.getpid(), rng[3], part, nsteps);
        fr = (part-1)*step1 + rng[1];
        to = min(part*step1, mm) + rng[1] - 1;

        slice = data$getDataRez(fr:to);

        slice = slice /
            rep( pmax(sqrt(colSums(slice^2)), 1e-3), each = data$nsamples);
        
        covmat = covmat + tcrossprod(slice);
        
        rm(slice);
    }
    data$close();
    .log(ld, "%s, Process %06d, Job %02d, Done PCA, CpG range %d-%d",
        date(), Sys.getpid(), rng[3], rng[1], rng[2]);

    gc();
    return(covmat);
}

plotPCvalues = function(values, n = 40, ylim = NULL, col = "blue"){
    pc100 = head(values,n)/sum(values)*100;
    if(is.null(ylim))
        ylim = c(0, pc100[1]*1.05);
    plot(
        x = pc100,
        pch = 19,
        col = col,
        ylim = ylim,
        xlim = c(0, length(pc100)+0.5),
        xlab = "Principal components (PCs)",
        ylab = "Variation Explained, %",
        yaxs = "i",
        xaxs = "i",
        axes = FALSE);
    axis(1);
    axis(2);
}

plotPCvectors = function(e, i, col = "blue1"){
    plot(
        x = e$vectors[,i],
        main = paste("PC",i),
        xlab = "Samples",
        ylab = "PC components",
        pch = 19,
        col = col,
        xlim = c(0.001, length(e$values)+0.999),
        xaxs = "i");
    abline(h = 0, col = "grey");
}

postPCAprocessing = function(param, e = NULL, plotPCs = 20){
    param = parameterPreprocess(param);
    ld = param$dirpca;
    .log(ld, "%s, Start postPCAprocessing() call", date());

    # Get data access
    data = new("rwDataClass", param = param, getPCs = FALSE);

    if(is.null(e))
        e = readRDS(file = paste0(param$dirpca, "/eigen.rds"))

    nonzeroPCs = sum(   
                    abs(e$values/e$values[1]) >
                    length(e$values)*.Machine$double.eps );

    # PCA plots
    {
        plotfilename = paste0(param$dirpca, "/PC_plot_covariates_removed.pdf")
        .log(ld, "%s, Saving PCA plots in: %s", date(), plotfilename);
        pdf(plotfilename, 7, 7);
        plotPCvalues(e$values, n = 40);
        title("Principal components");
        for( i in seq_len(min(plotPCs,nonzeroPCs)))
            plotPCvectors(e,i);
        dev.off();
    }

    # Save PCs and loadings
    {
        .log(ld, "%s, Saving PC values and vectors", date());
        PC_loads = e$vectors[,seq_len(min(20,nonzeroPCs))];
        rownames(PC_loads) = data$samplenames;
        colnames(PC_loads) = paste0("PC", seq_len(ncol(PC_loads)));
        write.table(file = paste0(param$dirpca, "/PC_loadings.txt"),
                    x = data.frame(name=rownames(PC_loads), PC_loads),
                    sep="\t",
                    row.names = FALSE);
        PC_values = data.frame(
                    PC_num = paste0("PC", seq_len(length(e$values))),
                    VarianceExplained = e$values/sum(e$values));
        write.table(file = paste0(param$dirpca, "/PC_values.txt"),
                    x = PC_values,
                    sep = "\t",
                    row.names = FALSE,
                    quote = FALSE);
    }

    # Saving PC vs. covariates association
    if( NCOL(param$covariates) > 1 ){
        .log(ld, "%s, Saving PC vs. covariates associations", date());
        testcov = .testCovariates(
                        covariates1 = param$covariates[-1],
                        data = e$vectors[, seq_len(nonzeroPCs)],
                        cvrtqr = data$cvrtqr);
        write.table(file = paste0(param$dirpca, "/PC_vs_covs_corr.txt"),
                    x = data.frame(
                        name = paste0("PC", seq_len(nonzeroPCs)),
                        testcov$crF,
                        check.names = FALSE),
                    sep="\t", 
                    row.names = FALSE);
        write.table(file = paste0(param$dirpca, "/PC_vs_covs_pvalue.txt"),
                    x = data.frame(
                        name = paste0("PC", seq_len(nonzeroPCs)),
                        testcov$pv,
                        check.names = FALSE),
                    sep="\t", 
                    row.names = FALSE);
    }
    data$close();
    .log(ld, "%s, Done postPCAprocessing() call", date());
    return(invisible(NULL));
}

# Step 4 of RaMWAS
ramwas4PCA = function( param ){
    # library(filematrix)
    param = parameterPreprocess(param);
    ld = param$dirpca;
    dir.create(param$dirpca, showWarnings = FALSE, recursive = TRUE);
    .log(ld, "%s, Start ramwas4PCA() call", date(), append = FALSE);
    

    parameterDump(dir = param$dirpca, param = param,
        toplines = c(   "dirpca", "dircoveragenorm",
                        "filecovariates", "covariates",
                        "modelcovariates",
                        "cputhreads", "diskthreads"));

    # Get data access
    data = new("rwDataClass", param = param, getPCs = FALSE);

    ### PCA part
    {
        ### Calculate covmat from the data matrix
        {
            .log(ld, "%s, Calculating covariance matrix", date());

            step1 = ceiling( 128*1024*1024 / data$ndatarows / 8);
            mm = data$ncpgs;
            nsteps = ceiling(mm/step1);

            # Memory concerns
            totmem = memory.limit() * 1048576;
            thrmem = (data$nsamples^2 * 8) * 3;
            maxthr = max(totmem / thrmem, 1)
            thrmem = max(thrmem, 1);
            
            nthreads = min(param$cputhreads, nsteps, thrmem);
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
                logfun = .logErrors(ld, .ramwas4PCAjob);
                clusterExport(  
                            cl = cl,
                            varlist = ".set1MLKthread", 
                            envir = asNamespace("ramwas"));
                clusterEvalQ(cl, eval(parse(text = .set1MLKthread)));
                # clusterCall(cl, function(){RevoUtilsMath::setMKLthreads()});
                
                covlist = clusterApplyLB(
                            cl = cl,
                            x = rangeset,
                            fun = logfun,
                            param = param);
                .showErrors(covlist);
                covmat = Reduce(f = `+`, x = covlist);
                tmp = sys.on.exit();
                eval(tmp);
                rm(tmp);
                on.exit();
                rm(cl, rng, rangeset, covlist);
                gc();
            } else {
                covmat = .ramwas4PCAjob( 
                                rng = c(1, data$ncpgs, 0),
                                param = param);
                if(is.character(covmat))
                    stop(covmat);
            }
            .log(ld, "%s, Done calculating covariance matrix", date());

            .log(ld, "%s, Saving covariance matrix", date());
            saveRDS(
                file = paste0(param$dirpca, "/covmat.rds"),
                object = covmat,
                compress = FALSE);
        } # covmat

        ### Eigenvalue decomposition
        {
            .log(ld, "%s, Performing Eigenvalue Decomposition", date());
            e = eigen(covmat, symmetric = TRUE);
            .log(ld, "%s, Saving Eigenvalue Decomposition", date());
            saveRDS(
                file = paste0(param$dirpca, "/eigen.rds"),
                object = e,
                compress = FALSE);
            # e = readRDS(paste0(param$dirpca, "/eigen.rds"));
        } # e
    }
    data$close();
    postPCAprocessing(param, e);
    .log(ld, "%s, Done ramwas4PCA() call", date());
    return(invisible(NULL));
}

# Interactive covariate selection via
# correlations with PCs
ramwasPCsCovariateSelection = function(param){
    message("Processing parameters");
    param = parameterPreprocess(param);

    if( is.null(param$covselpcs) )
        param$covselpcs = 3;

    # covariates
    ann = param$covariates;

    # covariance matrix
    filecovmat = paste0(param$dirpca, "/covmat.rds")
    message("Loading sample covariance matrix\n",filecovmat);
    covmat = readRDS(filecovmat);

    covset = param$modelcovariates;

    repeat {

        # Orthogonalize covariates

        message("Removing covariates from covariance matrix");

        cvtrqr = orthonormalizeCovariates(cvrt = ann[covset]);
        covmat1 = covmat;
        covmat1 = covmat1 - tcrossprod(covmat1 %*% cvtrqr, cvtrqr);
        covmat1 = covmat1 - cvtrqr %*% crossprod(cvtrqr, covmat1);

        ### Eigenvalue decomposition
        message("Performing eigenvalue decomposition");
        e = eigen(covmat1, symmetric = TRUE);

        ### First PC
        testslist = vector("list", param$covselpcs);
        for( i in seq_len(param$covselpcs) ){ # i=1
            pc = e$vectors[,i, drop=FALSE];

            message("Testing PC",i," vs. covariates");
            # testPhenotype(phenotype = ann[[2]], data = pc, cvrtqr = cvtrqr)

            # Run testPhenotype for each covariate
            temp1 = lapply( 
                            X = ann[-1],
                            FUN = testPhenotype,
                            data = pc,
                            cvrtqr = cvtrqr);
            
            # Extract first three elements
            temp2 = lapply( X = temp1, FUN = `[`, 1:3);
            
            testslist[[i]] = t(sapply(temp2, unlist));
        }
        tstatsabs = lapply(testslist, function(x)abs(x[,3]))
        tstatsabs = do.call(pmin, tstatsabs);

        ord = sort.list(tstatsabs, decreasing = FALSE);

        pvalues = lapply(testslist, function(x)abs(x[ord,3]));
        names(pvalues) = paste0("PC",seq_len(param$covselpcs), "_pvalues");

        tests = data.frame( covariates = rownames(testslist[[1]])[ord],
                            pvalues);
        # tests = tests[order(tests[,3], -abs(tests[,2])), ];
        # tests = data.frame(covariates = rownames(tests), tests);
        rownames(tests) = NULL;
        cat("\n",
            paste0(
                "Covariates Included: \n ",
                paste(covset, collapse = ",")),
            "\n");
        print(head(tests,15));

        cat("\n");
        cat("Enter the line number for the new covariate, (0 to stop):","\n");
        if(interactive()){
            newcov = readline();
            if(nchar(newcov) == 0){
                newcov = 0;
            } else {
                newcov = as.integer(newcov);
            }
        } else {
            newcov = scan("stdin", integer(), n = 1);
        }
        if( newcov == 0 )
            break;
        covset = c(covset, as.character(tests$covariates[as.integer(newcov)]) );
    }
    cat("New covariate line for parameter file:","\n","\n");
    cat(paste0(
        "modelcovariates = c(\n",
        paste0("  ",sapply(covset, deparse), collapse = ",\n"),
        ")"), "\n","\n");
    return(invisible(newcov));
}
