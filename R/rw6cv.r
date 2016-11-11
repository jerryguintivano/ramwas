# Run 10 (cvnfolds) MWASes 
# each on 90% of the data.
ramwas6crossValidation = function(param) {
	param = parameterPreprocess(param);
	param$toppvthreshold = 1e-300;
	dir.create(param$dircv, showWarnings = FALSE, recursive = TRUE);
	parameterDump(dir = param$dircv, param = param,
					  toplines = c("dircv", "mmncpgs", "mmalpha", "cvnfolds","randseed",
					  				 "dirmwas", "dirpca", "dircoveragenorm",
					  				 "filecovariates", "covariates",
					  				 "modeloutcome", "modelcovariates", "modelPCs",
					  				 "qqplottitle",
					  				 "cputhreads"));
	
	if(any(is.na(param$covariates[[ param$modeloutcome ]]))){
		param$covariates = data.frame(lapply( param$covariates, `[`, !is.na(param$covariates[[ param$modeloutcome ]])));
	}
	nms = param$covariates[[1]];
	nsamples = length(nms);
	
	set.seed(param$randseed);
	shuffle = sample.int(nsamples);
	starts = floor(seq(1, nsamples+1, length.out = param$cvnfolds+1));
	
	
	for( fold in seq_len(param$cvnfolds) ) { # fold = 9
		
		message("Running MWAS for fold ",fold," of ",param$cvnfolds);
		
		exclude = logical(nsamples);
		exclude[ shuffle[starts[fold]:(starts[fold+1]-1)] ] = TRUE;
		names(exclude) = nms;
		
		outcome = param$covariates[[ param$modeloutcome ]];
		
		param2 = param;
		param2$dirmwas = sprintf("%s/fold_%02d", param$dircv, fold);
		param2$covariates[[ param$modeloutcome ]][exclude] = NA;
		
		ramwas5MWAS(param2);
		saveRDS( file = paste0(param2$dirmwas, "/exclude.rds"), object = exclude);
	}
}

# Plot true outcome vs. prediction
# with correlations and p-value in the title
plotPrediction = function(param, outcome, forecast, main, dfFull = NULL){
    rng = range(c(outcome, forecast));
    c1 = cor(outcome, forecast, use = "complete.obs", method = "pearson");
    c2 = cor(outcome, forecast, use = "complete.obs", method = "spearman");
    c1p = max(c1, 0);
    c2p = max(c2, 0);
    if( is.null(dfFull) )
        dfFull = length(forecast) - 1 - 1;
    cor2tt = function(x) { return( x * sqrt( dfFull / (1 - pmin(x^2,1))));	}
    tt2pv = function(x) { return( (pt(-x,dfFull))); }
    
    MSE = sqrt(mean( (outcome - forecast)^2 ))
    MAD = median( abs(outcome - forecast) )
    # z = summary(lm(outcome ~ forecast)); z$coefficients[2,4]
    
    plot( outcome, forecast, pch = 19, col = "blue", xlab = param$modeloutcome, ylab = "CV prediction",
          xlim = rng, ylim = rng,
          main = sprintf("%s\nRMSE = %.3f, MAD = %.3f, cor = %.3f / %.3f (P/S)\nR2 = %.3f / %.3f, p-value = %.1e / %.1e",
                         main, MSE, MAD, c1, c2, c1p^2, c2p^2, tt2pv(cor2tt(c1p)), tt2pv(cor2tt(c2p)))
    );
    legend(x = "bottomright", 
           legend = c(paste0("# CpGs = ",   param$mmncpgs), 
                      paste0("EN alpha = ", param$mmalpha))
    );
    abline(a = 0, b = 1, col = "gray")
}

# Apply Elastic Net 10 times and collect
# out of sample predictions
ramwas7multiMarker = function(param){
    # library(glmnet);library(filematrix)
    # library(ramwas);
    param = parameterPreprocess(param);
    message("Working in: ",param$dircv);
    {
        if(any(is.na(param$covariates[[ param$modeloutcome ]]))){
            param$covariates = data.frame(lapply( param$covariates, `[`, !is.na(param$covariates[[ param$modeloutcome ]])));
        } 
    } # remove samples with NA's in outcome
    {
        message("Matching samples in covariates and data matrix");
        rez = .matchCovmatCovar( param );
        rowsubset = rez$rowsubset;
        ncpgs     = rez$ncpgs;
        rm(rez);
    } # rowsubset, ncpgs
    {
        outcome = param$covariates[[ param$modeloutcome ]];
    } # outcome
    forecast0 = NULL;
    for( fold in seq_len(param$cvnfolds) ) { # fold = 1
        message("Processing fold ", fold, " of ", param$cvnfolds);
        {
            dircvmwas = sprintf("%s/fold_%02d", param$dircv, fold);
            rdsfile = paste0(dircvmwas, "/exclude.rds");
            if( !file.exists( rdsfile ) ) {
                message("Missing CV MWAS: ", rdsfile);
                stop("Missing CV MWAS: ", rdsfile);
                next;
            }
            exclude = readRDS( rdsfile );
            rm(rdsfile);
        } # dircvmwas, exclude		
        {
            if( is.null(forecast0) ) {
                forecast0 = double(length(exclude)*2);
                dim(forecast0) = c(2, length(exclude));
                colnames(forecast0) = names(exclude);
            }
        } # forecast0
        {
            # get p-values
            fm = fm.open(paste0(dircvmwas, "/Stats_and_pvalues"))
            # colnames(fm)
            pv = fm[,3];
            close(fm);
            rm(fm);
        } # pv
        {
            # Find top param$mmncpgs CpGs
            # Faster than:
            # 	cpgset = sort.list(pv)[seq_len(param$mmncpgs)];
            # 	cpgset = sort.int(cpgset);
            # tic = proc.time();
            pvthr = 10^((-100):0);
            fi = findInterval( pv, pvthr);
            tab = cumsum(tabulate(fi));
            upperfi = which(tab > param$mmncpgs)[1];
            set1 = which(fi <= upperfi);
            cpgsetraw = set1[sort.list(pv[set1])[seq_len(param$mmncpgs)]];
            cpgset = sort.int(cpgsetraw);
            rm(pvthr, fi, tab, upperfi, set1, cpgsetraw, pv);
            # toc = proc.time();
            # show(toc-tic);
        } # cpgset, -pv
        {
            # get raw coverage
            fm = fm.open( paste0(param$dircoveragenorm,"/Coverage") );
            coverage = fm[, cpgset];
            # rownames(coverage) = rownames(fm);
            if( !is.null(rowsubset) )
                coverage = coverage[rowsubset,];
            close(fm);
            rm(fm);
        } # coverage
        {
            coverageTRAIN = coverage[!exclude,];
            coverageTEST  = coverage[ exclude,];
            rm(coverage);
            gc();
        } # coverageTRAIN, coverageTEST, -coverage
        {
            z1 = cv.glmnet(x = coverageTRAIN, 
                           y = as.vector(outcome[!exclude]), 
                           nfolds = param$cvnfolds, 
                           # keep = TRUE, 
                           parallel = FALSE, 
                           alpha = param$mmalpha);
            z2 = predict.cv.glmnet(object = z1, 
                                   newx = coverageTEST, 
                                   type = "response", 
                                   s = "lambda.min", 
                                   alpha = param$mmalpha);
            forecast0[1,exclude] = forecast0[1,exclude] + z2;
            forecast0[2,exclude] = forecast0[2,exclude] + 1;
            rm(z1, z2);
        } # forecast0
        rm(cpgset, coverageTRAIN, coverageTEST);
    }
    {
        rez = data.frame(samples = colnames(forecast0),
                         outcome = outcome,
                         forecast = forecast0[1,]/forecast0[2,]
        )
        write.table( file = sprintf("%s/MMCVN_prediction_folds=%02d_CpGs=%06d_alpha=%s.txt", 
                                    param$dircv, param$cvnfolds, param$mmncpgs, param$mmalpha),
                     x = rez,
                     sep = "\t", row.names = FALSE);
    } # rez
    {
        pdf( sprintf("%s/MMCVN_prediction_folds=%02d_CpGs=%06d_alpha=%s.pdf", 
                     param$dircv, param$cvnfolds, param$mmncpgs, param$mmalpha) );
        plotPrediction(param, outcome, rez$forecast, 
                       main = "Prediction success (EN on coverage)");
        dev.off();
    } # pdf()
    return( invisible(rez) );
}
