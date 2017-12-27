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

    plot(
        x = fpr,
        y = tpr,
        xlim = c(0,1),
        ylim = c(0,1),
        xaxs = "i",
        yaxs = "i",
        type = "l",
        col = "blue",
        lwd = 4,
        xlab = "False positive rate",
        ylab = "True positive rate");
    abline(a = 0, b = 1, col = "grey");
    legend(
        "bottomright", 
        legend = c("ROC curve","diagonal"),
        lwd = c(4,1), 
        lty = 1, 
        col = c("blue","grey"),
        bg = "transparent");
    
    auc = -sum(tpr[-1] * diff(fpr));
    
    legend(
        "bottom",
        legend = sprintf("AUC = %.3f", auc),
        bg = "transparent",
        bty = "n");
    return(auc);
}

# Run 10 (cvnfolds) MWASes
# each on 90% of the data.
ramwas7ArunMWASes = function(param){
    param = parameterPreprocess(param);
    if( param$toppvthreshold < 1 )
        param$toppvthreshold = 50;
    ld = param$dircv;
    dir.create(param$dircv, showWarnings = FALSE, recursive = TRUE);
    .log(ld, "%s, Start ramwas7ArunMWASes() call", date(), append = FALSE);

    if(is.null(param$modeloutcome))
        stop(
            "Model outcome variable not defined\n",
            "See \"modeloutcome\" parameter");
    if( !any(names(param$covariates) == param$modeloutcome) )
        stop(
            "Model outcome is not found among covariates.\n",
            "See \"modeloutcome\" parameter");

    outcome = param$covariates[[ param$modeloutcome ]];
    
    if( !is.numeric(outcome) )
        stop("Categorical outcomes are not supported by Elastic Net");
    
    parameterDump(dir = param$dircv, param = param,
        toplines = c(   "dircv", "mmncpgs", "mmalpha",
                        "cvnfolds","randseed",
                        "dirmwas", "dirpca", "dircoveragenorm",
                        "filecovariates", "covariates",
                        "modeloutcome", "modelcovariates",
                        "modelPCs",
                        "qqplottitle",
                        "cputhreads"));


    
    # Get data dimensions
    data = new("rwDataClass", param = param, getPCs = FALSE);
    data$close();
    # data$samplenames
    
    exclude0 = is.na(outcome);
    names(exclude0) = data$samplenames; #param$covariates[[1]];
    goodids = which(!exclude0);
    nids = length(goodids);

    set.seed(param$randseed);
    shuffle = sample(x = goodids);
    starts = floor(seq(1, nids+1, length.out = param$cvnfolds+1));


    for( fold in seq_len(param$cvnfolds) ){ # fold = 1

        .log(ld, "%s, Start MWAS for fold %02d of %d",
            date(), fold, param$cvnfolds);

        exclude = exclude0;
        exclude[ shuffle[starts[fold]:(starts[fold+1]-1)] ] = TRUE;
        
        param2 = param;
        param2$dirmwas = sprintf("%s/fold_%02d", param$dircv, fold);
        param2$covariates[[ param$modeloutcome ]][exclude] = NA;

        ramwas5MWAS(param2);
        saveRDS(
            file = paste0(param2$dirmwas, "/exclude.rds"),
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
    MSE = sqrt(mean( (outcome - forecast)^2, na.rm = TRUE ));
    MAD = median( abs(outcome - forecast)  , na.rm = TRUE );
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
plotPrediction = function(
        param,
        outcome,
        forecast,
        cpgs2use,
        main,
        dfFull = NULL){
    rng = range(c(outcome, forecast), na.rm = TRUE);
    stats = predictionStats(outcome, forecast, dfFull);

    plot( 
        x = outcome,
        y = forecast,
        pch = 19,
        col = "blue",
        xlab = param$modeloutcome,
        ylab = "CV prediction",
        xlim = rng,
        ylim = rng,
        main = sprintf(paste0(
                    "%s\n",
                    "RMSE = %.3f, MAD = %.3f, cor = %.3f / %.3f (P/S)\n",
                    "R2 = %.3f / %.3f, p-value = %.3e / %.3e"),
                main, stats$MSE, stats$MAD,
                stats$corp, stats$cors,
                stats$R2p, stats$R2s,
                stats$pvp, stats$pvs));
    legend(
        x = "bottomright",
        legend = c(
                paste0("# CpGs = ",   cpgs2use),
                paste0("EN alpha = ", param$mmalpha)));
    abline(a = 0, b = 1, col = "gray")
}

# Apply Elastic Net 10 times and collect
# out of sample predictions
ramwas7BrunElasticNet = function(param){
    # library(glmnet); library(filematrix); library(ramwas);
    param = parameterPreprocess(param);
    ld = param$dircv;
    .log(ld, "%s, Start ramwas7BrunElasticNet() call", date());
    
    if(is.null(param$modeloutcome))
        stop(
            "Model outcome variable not defined\n",
            "See \"modeloutcome\" parameter");
    if( !any(names(param$covariates) == param$modeloutcome) )
        stop(
            "Model outcome is not found among covariates.\n",
            "See \"modeloutcome\" parameter");

    
    # Get data access
    data = new("rwDataClass");
    data$open(param, getPCs = FALSE);
    
    outcome = param$covariates[[ param$modeloutcome ]];

    for( cpgs2use in param$mmncpgs ){ # cpgs2use = param$mmncpgs[1]
        .log(ld, "%s, Applying Elasting Net to %d top CpGs", date(), cpgs2use);

        forcastlist = list();
        forecast0 = NULL;
        for( fold in seq_len(param$cvnfolds) ){ # fold = 1
            .log(ld, "%s, Processing fold %02d of %d", 
                date(), fold, param$cvnfolds);
            {
                dircvmwas = sprintf("%s/fold_%02d", param$dircv, fold);
                rdsfile = paste0(dircvmwas, "/exclude.rds");
                if( !file.exists( rdsfile ) ){
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
                mwas = getMWAS(dircvmwas);
                pv = mwas$`p-value`;
                rm(mwas);
            } # pv
            {
                cpgset = findBestNpvs(pv, n = cpgs2use);
                rm(pv);
            } # cpgset, -pv
            {
                coverage = data$getDataRez(cpgset, resid = FALSE)
            } # coverage
            {
                coverageTRAIN = coverage[!exclude,];
                coverageTEST  = coverage[ exclude,];
                rm(coverage);
                gc();
            } # coverageTRAIN, coverageTEST, -coverage
            {
                # library(glmnet)
                net = cv.glmnet(
                            x = coverageTRAIN,
                            y = as.vector(outcome[!exclude]),
                            nfolds = param$cvnfolds,
                            # keep = TRUE,
                            parallel = FALSE,
                            alpha = param$mmalpha);
                pred = predict.cv.glmnet(
                            object = net,
                            newx = coverageTEST,
                            type = "response",
                            s = "lambda.min",
                            alpha = param$mmalpha);

                forecast0[1,exclude] = forecast0[1,exclude] + pred;
                forecast0[2,exclude] = forecast0[2,exclude] + 1;

                adt = data.frame(
                            ind = which(exclude), 
                            outcome = outcome[exclude], 
                            forecast = pred,
                            row.names = names(exclude)[exclude]);
                forcastlist[[length(forcastlist)+1]] = adt;
                rm(net, pred, adt);
            } # forecast0, forcastlist
            rm(cpgset, coverageTRAIN, coverageTEST);
            gc();
        } # fold in seq_len(param$cvnfolds)
        .log(ld, "%s, Saving Cross Validation Predictions and Plots", date());
        forecast = forecast0[1,]/forecast0[2,];
        {
            rez = data.frame(
                        samples = colnames(forecast0),
                        outcome = outcome,
                        forecast = forecast)
            write.table( 
                file = sprintf(
                    "%s/MMCVN_prediction_folds=%02d_CpGs=%06d_alpha=%f.txt",
                    param$dircv, param$cvnfolds, cpgs2use, param$mmalpha),
                x = rez,
                sep = "\t", 
                row.names = FALSE);
            dirr = paste0(param$dircv,"/rds");
            dir.create(dirr, showWarnings = FALSE);
            saveRDS(file = sprintf(
                "%s/CpGs=%06d_alpha=%f.rds",
                dirr,
                cpgs2use,
                param$mmalpha),
                object = list(
                    outcome = outcome, 
                    forecast = forecast, 
                    forcastlist = forcastlist));
        } # rez
        {
            pdf(sprintf("%s/MMCVN_prediction_folds=%02d_CpGs=%06d_alpha=%f.pdf",
                    param$dircv, param$cvnfolds, cpgs2use, param$mmalpha) );
            plotPrediction(param, outcome, rez$forecast, cpgs2use,
                    main = "Prediction success (EN on coverage)");
            dev.off();
            if( length(unique(outcome)) == 2 ){
                pdf( sprintf("%s/ROC_CpGs=%06d_alpha=%f.pdf",
                        param$dircv, cpgs2use, param$mmalpha) );
                plotROC(outcome = outcome,
                        forecast = forecast)
                legend(
                    "right",
                    legend = c(
                        paste0("# CpGs = ",   cpgs2use),
                        paste0("EN alpha = ", param$mmalpha)),
                    bg = "transparent",
                    bty = "n");
                title(paste0(
                    "ROC curve for prediction of \"", param$modeloutcome,"\"\n",
                    param$cvnfolds, "-fold cross validation"));
                dev.off();
            }
        } # pdf()
    } # cpgs2use in param$mmncpgs
    data$close();
    .log(ld, "%s, Done ramwas7BrunElasticNet() call", date());
    return( invisible(NULL) );
}

plotCVcors = function(cl, param){
    aymax = max(abs(cl$cors),abs(cl$corp));
    plot( 
        x = cl$x,
        y = cl$corp,
        col = "red",
        pch = 19,
        ylim = c(-1,1)*aymax,
        log = "x",
        xlab = "Number of markers",
        ylab = "Correlation")
    points(cl$x, cl$cors, col = "cyan4", pch = 17)
    abline(h=0, col = "grey")
    legend(x = "bottomleft", legend = paste0("EN alpha = ", param$mmalpha))
    legend(
        "bottomright",
        legend = c("Correlations:", "Pearson", "Spearman"),
        pch = c(NA_integer_,19,17),
        col = c("red","red","cyan4"));
    title(paste0(
            "Prediction of \"", param$modeloutcome,"\"\n",
            param$cvnfolds,"-fold cross validation"));
    return(invisible(NULL));
}

ramwas7CplotByNCpGs = function(param){
    param = parameterPreprocess(param);
    if( length(param$mmncpgs) <= 1 )
        return();
    message("Working in: ",param$dircv);
    cors = double(length( param$mmncpgs ));
    corp = double(length( param$mmncpgs ));
    datalist = vector("list", length(param$mmncpgs));
    for( i in seq_along(param$mmncpgs) ){ # i = 1
        cpgs2use = param$mmncpgs[i];
        data = readRDS(sprintf(
            "%s/rds/CpGs=%06d_alpha=%f.rds",
            param$dircv, cpgs2use, param$mmalpha));
        datalist[[i]] = data;
        stats = predictionStats(
                    outcome = data$outcome,
                    forecast = data$forecast);
        cors[i] = stats$cors;
        corp[i] = stats$corp;
    }

    cl = list(
            x = param$mmncpgs,
            cors = cors,
            corp = corp);
    saveRDS(file = sprintf( "%s/rds/cor_data_alpha=%f.rds",
                            param$dircv, param$mmalpha),
            object = cl);

    pdf( sprintf("%s/Prediction_alpha=%f.pdf", param$dircv, param$mmalpha) );
    plotCVcors(cl, param);
    dev.off();

    if( length(unique(datalist[[1]]$outcome)) == 2 ){
        pdf( sprintf("%s/ROC_alpha=%f.pdf",
                    param$dircv, param$mmalpha) );
        for( i in seq_along(param$mmncpgs) ){ # i = 1
            cpgs2use = param$mmncpgs[i];
            plotROC(outcome = datalist[[i]]$outcome,
                    forecast = datalist[[i]]$forecast)
            legend(
                "right",
                legend = c(
                    paste0("# CpGs = ",   cpgs2use),
                    paste0("EN alpha = ", param$mmalpha)),
                bg = "transparent",
                bty = "n");
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
