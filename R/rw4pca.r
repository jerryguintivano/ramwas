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
    return(list(crF=crF, pv=pv));
}

# Job function for PCA analysis
# (covariance matrix calculation)
.ramwas4PCAjob = function(rng, param, cvrtqr, rowsubset){
    # rng = rangeset[[1]];
    # library(filematrix);
    fm = fm.open( paste0(param$dircoveragenorm, "/Coverage"),
                  readonly = TRUE,
                  lockfile = param$lockfile2);

    covmat = 0;

    step1 = ceiling( 128*1024*1024 / nrow(fm) / 8);
    mm = rng[2]-rng[1]+1;
    nsteps = ceiling(mm/step1);
    for( part in 1:nsteps ){ # part = 1
        message("Slice ", part, " of ", nsteps);
        fr = (part-1)*step1 + rng[1];
        to = min(part*step1, mm) + rng[1] - 1;

        slice = fm[,fr:to];
        if( !is.null(rowsubset) )
            slice = slice[rowsubset,];
        slice = t(slice);

        slice = slice - tcrossprod(slice, cvrtqr) %*% cvrtqr;
        ### rowMeans(slice) == 0
        slice = slice / pmax(sqrt(rowSums(slice^2)), 1e-3);
        ### rowSums(slice^2) == 1

        covmat = covmat + crossprod(slice);
        fm$filelock$lockedrun( {
            cat(file = paste0(param$dirpca,"/Log.txt"),
                 date(), ", Process ", Sys.getpid(), ", Job ", rng[3],
                 ", processing slice ", part, " of ", nsteps, "\n",
                 sep = "", append = TRUE);
        });
        rm(slice);
    }
    close(fm)
    return(covmat);
}

plotPCvalues = function(values, n = 40){
    pc100 = head(values,n)/sum(values)*100;
    plot(pc100,
         pch = 19,
         col="blue",
         ylim = c(0, pc100[1]*1.05),
         xlim = c(0, length(pc100)+0.5),
         main = "Principal components",
         xlab = "PCs",
         ylab = "Variation Explained (%)",
         yaxs = "i",
         xaxs = "i")
}

plotPCvectors = function(e, i){
     plot(e$vectors[,i],
                 main = paste("PC",i),
                 xlab = "Samples",
                 ylab = "PC components",
                 pch = 19,
                 col = "blue1",
                 xlim = c(0, length(e$values)+0.5),
                 xaxs = "i");
    abline(h = 0, col = "grey");
}

postPCAprocessing = function(param, e = NULL, plotPCs = 20){
    param = parameterPreprocess(param);

    ### Get and match sample names
    {
        message("Matching samples in covariates and data matrix");
        rez = .matchCovmatCovar( param );
        # rez = ramwas:::.matchCovmatCovar( param );
        rowsubset = rez$rowsubset;
        ncpgs     = rez$ncpgs;
        cvsamples = rez$cvsamples;
        rm(rez);
    } # rowsubset, ncpgs, cvsamples

    ### Prepare covariates, defactor
    {
        message("Preparing covariates (splitting dummies, orthonormalizing)");
        cvrtqr = .getCovariates(param = param,
                                rowsubset = rowsubset,
                                modelhasconstant = param$modelhasconstant);
    } # cvrtqr

    if(is.null(e))
        e = readRDS(file = paste0(param$dirpca,"/eigen.rds"))

    nonzeroPCs = sum(   abs(e$values/e$values[1]) >
                        length(e$values)*.Machine$double.eps );

    # PCA plots
    {
        message("Saving PCA plots");
        pdf(paste0(param$dirpca, "/PC_plot_covariates_removed.pdf"),7,7);
        plotPCvalues(e$values,40);
        for( i in seq_len(min(plotPCs,nonzeroPCs)))
            plotPCvectors(e,i);
        dev.off();
    }

    # Save PCs and loadings
    {
        message("Saving PC values and vectors");
        PC_loads = e$vectors[,seq_len(min(20,nonzeroPCs))];
        rownames(PC_loads) = cvsamples;
        colnames(PC_loads) = paste0("PC",seq_len(ncol(PC_loads)));
        write.table(file = paste0(param$dirpca, "/PC_loadings.txt"),
                    x = data.frame(name=rownames(PC_loads),PC_loads),
                    sep="\t",
                    row.names = FALSE);
        PC_values = data.frame(
            PC_num = paste0("PC",seq_len(length(e$values))),
            VarianceExplained = e$values/sum(e$values));
        write.table(file = paste0(param$dirpca, "/PC_values.txt"),
                    x = PC_values,
                    sep = "\t",
                    row.names = FALSE,
                    quote = FALSE);
    }

    # Saving PC vs. covariates association
    if(NCOL(param$covariates) > 1){
        message("Saving PC vs. covariates associations");
        testcov = .testCovariates(covariates1 = param$covariates[-1],
                                  data = e$vectors[,seq_len(nonzeroPCs)],
                                  cvrtqr = cvrtqr);

        write.table(file = paste0(param$dirpca, "/PC_vs_covs_corr.txt"),
                    x = data.frame(
                        name=paste0("PC",seq_len(nonzeroPCs)),
                        testcov$crF,
                        check.names = FALSE),
                    sep="\t", row.names = FALSE);
        write.table(file = paste0(param$dirpca, "/PC_vs_covs_pvalue.txt"),
                    x = data.frame(
                        name=paste0("PC",seq_len(nonzeroPCs)),
                        testcov$pv,
                        check.names = FALSE),
                    sep="\t", row.names = FALSE);
    }
}

# Step 4 of RaMWAS
ramwas4PCA = function( param ){
    # library(filematrix)
    param = parameterPreprocess(param);
    dir.create(param$dirpca,  showWarnings = FALSE, recursive = TRUE);

    parameterDump(dir = param$dirpca, param = param,
                      toplines = c("dirpca", "dircoveragenorm",
                                       "filecovariates", "covariates",
                                       "modelcovariates",
                                       "cputhreads", "diskthreads"));

    ### Get and match sample names
    {
        message("Matching samples in covariates and data matrix");
        rez = .matchCovmatCovar( param );
        # rez = ramwas:::.matchCovmatCovar( param );
        rowsubset = rez$rowsubset;
        ncpgs     = rez$ncpgs;
        cvsamples = rez$cvsamples;
        rm(rez);
    } # rowsubset, ncpgs, cvsamples

    ### Prepare covariates, defactor
    {
        message("Preparing covariates (splitting dummies, orthonormalizing)");
        param$modelPCs = 0;
        mwascvrtqr = .getCovariates(param = param,
                                    rowsubset = rowsubset,
                                    modelhasconstant = param$modelhasconstant);
    } # mwascvrtqr

    ### PCA part
    {
        ### Calculate covmat from the data matrix
        {
            message("Calculating Principal Components");
            cat(file = paste0(param$dirpca,"/Log.txt"),
                 date(), ", Running Principal Component Analysis.", "\n",
                 sep = "", append = FALSE);

            step1 = ceiling( 128*1024*1024 / length(cvsamples) / 8);
            mm = ncpgs;
            nsteps = ceiling(mm/step1);

            nthreads = min(param$diskthreads, nsteps);
            rm(step1, mm, nsteps);
            if( nthreads > 1 ){
                rng = round(seq(1, ncpgs+1, length.out = nthreads+1));
                rangeset = rbind( rng[-length(rng)],
                                  rng[-1]-1,
                                  seq_len(nthreads));
                rangeset = lapply(seq_len(ncol(rangeset)),
                                  function(i) rangeset[,i])

                if(param$usefilelock) param$lockfile2 = tempfile();
                # library(parallel);
                cl = makeCluster(nthreads);
                on.exit({stopCluster(cl);});
                covlist = clusterApplyLB(cl,
                                         rangeset,
                                         .ramwas4PCAjob,
                                         param = param,
                                         cvrtqr = mwascvrtqr,
                                         rowsubset = rowsubset);
                covmat = Reduce(f = `+`, x = covlist);
                tmp = sys.on.exit();
                eval(tmp);
                rm(tmp);
                on.exit();
                rm(cl, rng, rangeset, covlist);
                .file.remove(param$lockfile2);
            } else {
                covmat = .ramwas4PCAjob( rng = c(1, ncpgs, 0),
                                         param = param,
                                         cvrtqr = mwascvrtqr,
                                         rowsubset = rowsubset);
            }
            cat(file = paste0(param$dirpca,"/Log.txt"),
                 date(), ", Done running Principal Component Analysis.", "\n",
                 sep = "", append = TRUE);

            saveRDS(file = paste0(param$dirpca,"/covmat.rds"),
                    object = covmat,
                    compress = FALSE);
        } # covmat

        ### Eigenvalue decomposition
        {
            message("Performing Eigenvalue Decomposition");
            e = eigen(covmat, symmetric=TRUE);
            saveRDS(file = paste0(param$dirpca,"/eigen.rds"),
                    object = e,
                    compress = FALSE);
            # e = readRDS(paste0(param$dirpca,"/eigen.rds"));
        } # e
    }

    postPCAprocessing(param, e);
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
        e = eigen(covmat1, symmetric=TRUE);

        ### First PC
        testslist = vector("list", param$covselpcs);
        for( i in seq_len(param$covselpcs) ){ # i=1
            pc = e$vectors[,i, drop=FALSE];

            message("Testing PC",i," vs. covariates");
            # testPhenotype(phenotype = ann[[2]], data = pc, cvrtqr = t(cvtrqr))

            testslist[[i]] = t(sapply( lapply( lapply( ann[-1],
                                                       testPhenotype,
                                                       data=pc,
                                                       cvrtqr=t(cvtrqr)),
                                               `[`,
                                               1:3),
                                       unlist));
        }
        tstatsabs = lapply(testslist, function(x)abs(x[,3]))
        tstatsabs = do.call(pmin, tstatsabs);

        ord = sort.list(tstatsabs, decreasing = FALSE);

        pvalues = lapply(testslist, function(x)abs(x[ord,3]))
        names(pvalues) = paste0("PC",seq_len(param$covselpcs), "_pvalues");

        tests = data.frame( covariates = rownames(testslist[[1]])[ord],
                            pvalues);
        # tests = tests[order(tests[,3], -abs(tests[,2])), ];
        # tests = data.frame(covariates = rownames(tests), tests);
        rownames(tests) = NULL;
        cat("\n",
            paste0("Covariates Included: \n ",
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
            newcov = scan("stdin", integer(), n=1);
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
