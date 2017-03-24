# Ploc ROC curve for an outcome/forecast pair
plotROC = function(outcome, forecast){
    b = outcome[order(forecast)];
    b = (b == max(b));
    
    countLeftOf = 0:length(b);
    countHereOrRight = length(b) - countLeftOf;
    onesLeftOf = c(0L,cumsum(b));
    totalOnes = tail(onesLeftOf,1);
    
    onesHereOrRight = totalOnes - onesLeftOf
    zeroesHereOrRight = countHereOrRight - onesHereOrRight;
    
    fpr = zeroesHereOrRight / zeroesHereOrRight[1]
    tpr = onesHereOrRight / onesHereOrRight[1]
    
    plot(x = fpr,
         y = tpr, 
         xlim = c(0,1), 
         ylim = c(0,1), 
         xaxs = 'i', 
         yaxs = 'i', 
         type = 'l', 
         col = 'blue', 
         lwd = 2,
         xlab = 'False positive rate', 
         ylab = 'True positive rate');
    abline(a = 0, b = 1);
}

# Run 10 (cvnfolds) MWASes
# each on 90% of the data.
ramwas7ArunMWASes = function(param){
    param = parameterPreprocess(param);
    if( param$toppvthreshold < 1 )
    	param$toppvthreshold = 50;
    dir.create(param$dircv, showWarnings = FALSE, recursive = TRUE);
    parameterDump(dir = param$dircv, param = param,
                    toplines = c("dircv", "mmncpgs", "mmalpha",
                                "cvnfolds","randseed",
                                "dirmwas", "dirpca", "dircoveragenorm",
                                "filecovariates", "covariates",
                                "modeloutcome", "modelcovariates",
                                "modelPCs",
                                "qqplottitle",
                                "cputhreads"));

    if(any(is.na(param$covariates[[ param$modeloutcome ]]))){
        param$covariates = data.frame(lapply(
            param$covariates,
            `[`,
            !is.na(param$covariates[[ param$modeloutcome ]])));
    }
    nms = param$covariates[[1]];
    nsamples = length(nms);

    set.seed(param$randseed);
    shuffle = sample.int(nsamples);
    starts = floor(seq(1, nsamples+1, length.out = param$cvnfolds+1));


    for( fold in seq_len(param$cvnfolds) ){ # fold = 9

        message("Running MWAS for fold ",fold," of ",param$cvnfolds);

        exclude = logical(nsamples);
        exclude[ shuffle[starts[fold]:(starts[fold+1]-1)] ] = TRUE;
        names(exclude) = nms;

        outcome = param$covariates[[ param$modeloutcome ]];

        param2 = param;
        param2$dirmwas = sprintf("%s/fold_%02d", param$dircv, fold);
        param2$covariates[[ param$modeloutcome ]][exclude] = NA;

        ramwas5MWAS(param2);
        saveRDS( file = paste0(param2$dirmwas, "/exclude.rds"),
                 object = exclude);
    }
    return( invisible(NULL) );
}

predictionStats = function(outcome, forecast, dfFull = NULL){
    if( is.null(dfFull) )
        dfFull = length(forecast) - 1 - 1;
    cor2tt = function(x){ return( x * sqrt(dfFull / (1 - pmin(x^2,1)))) }
    tt2pv = function(x){ return( (pt(-x,dfFull))) }
    corp = cor(outcome, forecast, use = "complete.obs", method = "pearson");
    cors = cor(outcome, forecast, use = "complete.obs", method = "spearman");
    MSE = sqrt(mean( (outcome - forecast)^2 ));
    MAD = median( abs(outcome - forecast) );
    R2p = max(corp, 0)^2;
    R2s = max(cors, 0)^2;
    tp = cor2tt(corp);
    ts = cor2tt(cors);
    pvp = tt2pv(max(tp,0));
    pvs = tt2pv(max(ts,0));
    return(list(
        outcome = outcome,
        forecast = forecast,
        dfFull = dfFull,
        corp = corp,
        cors = cors,
        MSE = MSE,
        MAD = MAD,
        R2p = R2p,
        R2s = R2s,
        tp = tp,
        ts = ts,
        pvp = pvp,
        pvs = pvs));
}

# Plot true outcome vs. prediction
# with correlations and p-value in the title
plotPrediction = function(param, outcome, forecast, 
                          cpgs2use, main, dfFull = NULL){
    rng = range(c(outcome, forecast));
    stats = predictionStats(outcome, forecast, dfFull);

    plot( outcome,
          forecast,
          pch = 19,
          col = "blue",
          xlab = param$modeloutcome,
          ylab = "CV prediction",
          xlim = rng,
          ylim = rng,
          main = sprintf(paste0(
              "%s\n",
              "RMSE = %.3f, MAD = %.3f, cor = %.3f / %.3f (P/S)\n",
              "R2 = %.3f / %.3f, p-value = %.1e / %.1e"),
            main, stats$MSE, stats$MAD,
            stats$corp, stats$cors,
            stats$R2p, stats$R2s,
            stats$pvp, stats$pvs)
    );
    legend(x = "bottomright",
           legend = c(paste0("# CpGs = ",   cpgs2use),
                      paste0("EN alpha = ", param$mmalpha))
    );
    abline(a = 0, b = 1, col = "gray")
}

# Apply Elastic Net 10 times and collect
# out of sample predictions
ramwas7BrunElasticNet = function(param){
    # library(glmnet); library(filematrix); library(ramwas);
    param = parameterPreprocess(param);
    message("Working in: ",param$dircv);
    {
        if(any(is.na(param$covariates[[ param$modeloutcome ]]))){
            param$covariates = data.frame(lapply(
                param$covariates,
                `[`,
                !is.na(param$covariates[[ param$modeloutcome ]])));
        }
    } # remove samples with NAs in outcome
    {
        message("Matching samples in covariates and data matrix");
        rez = .matchCovmatCovar( param );
	    # rez = ramwas:::.matchCovmatCovar( param );
        rowsubset = rez$rowsubset;
        ncpgs     = rez$ncpgs;
        rm(rez);
    } # rowsubset, ncpgs
    {
        outcome = param$covariates[[ param$modeloutcome ]];
    } # outcome

    for( cpgs2use in param$mmncpgs ){ # cpgs2use = param$mmncpgs[1]
        message("Applying Elasting Net to ",cpgs2use," top CpGs");

        forecast0 = NULL;
        for( fold in seq_len(param$cvnfolds) ){ # fold = 1
            message("Processing fold ", fold, " of ", param$cvnfolds);
            {
                dircvmwas = sprintf("%s/fold_%02d", param$dircv, fold);
                rdsfile = paste0(dircvmwas, "/exclude.rds");
                if( !file.exists( rdsfile ) ){
                    message("Missing CV MWAS: ", rdsfile);
                    stop("Missing CV MWAS: ", rdsfile);
                    next;
                }
                exclude = readRDS( rdsfile );
                rm(rdsfile);
            } # dircvmwas, exclude
            {
                if( is.null(forecast0) ){
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
                cpgset = findBestNpvs(pv, n = cpgs2use);
                rm(pv);
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
                # library(glmnet)
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
            	# 
            	# z1a = cv.glmnet(x = coverageTRAIN,
            	# 					y = factor(as.vector(outcome[!exclude])),
            	# 					nfolds = param$cvnfolds,
            	# 					# keep = TRUE,
            	# 					parallel = FALSE,
            	# 					alpha = param$mmalpha,
            	# 					family = "binomial");
            	# z2a = predict.cv.glmnet(object = z1a,
            	# 							  newx = coverageTEST,
            	# 							  type = "response",
            	# 							  s = "lambda.min",
            	# 							  alpha = param$mmalpha);
            	
                forecast0[1,exclude] = forecast0[1,exclude] + z2;
                forecast0[2,exclude] = forecast0[2,exclude] + 1;
                rm(z1, z2);
            } # forecast0
            rm(cpgset, coverageTRAIN, coverageTEST);
        }
        forecast = forecast0[1,]/forecast0[2,];
        {
            rez = data.frame(samples = colnames(forecast0),
                             outcome = outcome,
                             forecast = forecast
            )
            write.table( file = sprintf(
                "%s/MMCVN_prediction_folds=%02d_CpGs=%06d_alpha=%f.txt",
                param$dircv,
                param$cvnfolds,
                cpgs2use,
                param$mmalpha),
                         x = rez,
                         sep = "\t", row.names = FALSE);
            dirr = paste0(param$dircv,"/rds");
            dir.create(dirr, showWarnings = FALSE);
            saveRDS(file = sprintf(
                "%s/CpGs=%06d_alpha=%f.rds",
                dirr,
                cpgs2use,
                param$mmalpha),
                object = list(outcome = outcome, forecast = forecast));
        } # rez
        {
            pdf(sprintf("%s/MMCVN_prediction_folds=%02d_CpGs=%06d_alpha=%f.pdf",
                         param$dircv,
                         param$cvnfolds,
                         cpgs2use,
                         param$mmalpha) );
            plotPrediction(param, outcome, rez$forecast, cpgs2use,
                           main = "Prediction success (EN on coverage)");
            dev.off();
            if( length(unique(outcome)) == 2 ){
                pdf( sprintf("%s/ROC_CpGs=%06d_alpha=%f.pdf",
                             param$dircv,
                             cpgs2use,
                             param$mmalpha) );
                plotROC(outcome = outcome,
                        forecast = forecast)
                legend(x = "bottomright",
                       legend = c(paste0("# CpGs = ",   cpgs2use),
                                  paste0("EN alpha = ", param$mmalpha)))
                title(paste0(
                    "ROC curve for prediction of \"", param$modeloutcome,"\"\n",
                    param$cvnfolds, "-fold cross validation"));
                dev.off();
            }
        } # pdf()
    }
    return( invisible(NULL) );
}

plotCVcors = function(cl, param){
    aymax = max(abs(cl$cors),abs(cl$corp));
    plot( x = cl$x,
          y = cl$cors,
          col = "red",
          pch = 19,
          ylim = c(-1,1)*aymax,
          log = "x",
          xlab = "Number of markers",
          ylab = "Correlation")
    points( cl$x, cl$corp, col = "cyan4", pch = 17)
    abline(h=0, col = "grey")
    legend(x = "bottomleft", legend = paste0("EN alpha = ", param$mmalpha))
    legend("bottomright",
           legend = c("Correlations:","Pearson", "Spearman"),
           pch = c(NA_integer_,19,17),
           col=c("red","red","cyan4"));
    title(paste0("Prediction of \"", param$modeloutcome,"\"\n",
                 param$cvnfolds,"-fold cross validation"));
}

ramwas7CplotByNCpGs = function(param){
    param = parameterPreprocess(param);
    if( length(param$mmncpgs) <= 1 )
        return();
    message("Working in: ",param$dircv);
    cors = double(length( param$mmncpgs ));
    corp = double(length( param$mmncpgs ));
    datalist = vector('list', length(param$mmncpgs));
    for( i in seq_along(param$mmncpgs) ){ # i = 1
        cpgs2use = param$mmncpgs[i];
        data = readRDS(sprintf(
            "%s/rds/CpGs=%06d_alpha=%f.rds",
            param$dircv,
            cpgs2use,
            param$mmalpha));
        datalist[[i]] = data;
        stats = predictionStats(outcome = data$outcome,
                                forecast = data$forecast);
        cors[i] = stats$cors;
        corp[i] = stats$corp;
    }

    cl = list(x = param$mmncpgs,
              cors = cors, 
              corp = corp);
    saveRDS(file = sprintf( "%s/rds/cor_data_alpha=%f.rds",
                            param$dircv,
                            param$mmalpha),
            object = cl);

    pdf( sprintf("%s/Prediction_alpha=%f.pdf",
                 param$dircv,
                 param$mmalpha) );
    plotCVcors(cl, param);
    dev.off();
    
    if( length(unique(datalist[[1]]$outcome)) == 2 ){
        pdf( sprintf("%s/ROC_alpha=%f.pdf",
                     param$dircv,
                     param$mmalpha) );
        for( i in seq_along(param$mmncpgs) ){ # i = 1
            cpgs2use = param$mmncpgs[i];
            plotROC(outcome = datalist[[i]]$outcome,
                    forecast = datalist[[i]]$forecast)
            legend(x = "bottomright",
                   legend = c(paste0("# CpGs = ",   cpgs2use),
                              paste0("EN alpha = ", param$mmalpha)))
            title(paste0(
                "ROC curve for prediction of \"", param$modeloutcome,"\"\n",
                param$cvnfolds,"-fold cross validation"));
        }
        dev.off();
    }
    return( invisible(NULL) );
}

ramwas7riskScoreCV = function(param){
    param = parameterPreprocess(param);
    ramwas7ArunMWASes(param);
    ramwas7BrunElasticNet(param);
    ramwas7CplotByNCpGs(param);
}
