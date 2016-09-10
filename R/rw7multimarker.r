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
ramwas7multiMarkerNoCvrt = function(param){
	# library(glmnet);library(filematrix)
	# library(ramwas);
	{
		if(any(is.na(param$covariates[[ param$modeloutcome ]]))){
			param$covariates = data.frame(lapply( param$covariates, `[`, !is.na(param$covariates[[ param$modeloutcome ]])));
		} 
	} # remove samples with NA's in outcome
	{
		message("Matching samples in covariates and data matrix");
		rez = ramwas:::.matchCovmatCovar( param );
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

ramwas7multiMarkerWithCvrt = function(param){
	# library(glmnet);library(filematrix)
	# library(ramwas);
	{
		if(any(is.na(param$covariates[[ param$modeloutcome ]]))){
			param$covariates = data.frame(lapply( param$covariates, `[`, !is.na(param$covariates[[ param$modeloutcome ]])));
		} 
	} # remove samples with NA's in outcome
	{
		message("Matching samples in covariates and data matrix");
		rez = ramwas:::.matchCovmatCovar( param );
		rowsubset = rez$rowsubset;
		ncpgs     = rez$ncpgs;
		rm(rez);
	} # rowsubset, ncpgs
	{
		outcome = param$covariates[[ param$modeloutcome ]];
		cvrt = ramwas:::.getCovariates(param, rowsubset, FALSE); # cvrtqr = ramwas:::.getCovariates(param)
		# cvrt = cvrt - rowMeans(cvrt);
		cvrt = rbind(rep(1,ncol(cvrt)), cvrt);
	} # outcome, cvrt (with constant)
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
				forecastLM = forecast0;
				forecastEN = forecast0;
			}
		} # forecast0, forecastLM, forecastEN
		{
			cvrtTRAIN = cvrt[ , !exclude, drop = FALSE];
			qq = qr(t(cvrtTRAIN));		
			cvrtQRTRAIN = t(qr.Q(qq));
			# range( solve(t(qr.R(qq))) %*% cvrtTRAIN - cvrtQRTRAIN )
			cvrtQRFULL = solve(t(qr.R(qq))) %*% cvrt;
			# range(cvrtQRFULL[,!exclude] - cvrtQRTRAIN)
			cvrtQRTEST  = solve(t(qr.R(qq))) %*% cvrt[, exclude];
			# range(cvrtQRFULL[, exclude] - cvrtQRTEST)
			rm(qq, cvrtTRAIN);
		} # cvrtQRTRAIN, cvrtQRFULL, cvrtQRTEST
		{
			# Predict with EM using just covariates
			z1 = cv.glmnet(x = t(cvrt[ ,!exclude, drop = FALSE]),
								y = as.vector(outcome[!exclude]),
								nfolds = param$cvnfolds,
								keep = TRUE,
								parallel = FALSE,
								alpha = param$mmalpha);
			z2 = predict.cv.glmnet(object = z1,
										  newx = t(cvrt[ ,exclude, drop = FALSE]), 
										  type = "response",
										  s = "lambda.min",
										  alpha = param$mmalpha);
			forecastEN[1,exclude] = forecastEN[1,exclude] + z2;
			forecastEN[2,exclude] = forecastEN[2,exclude] + 1;
			rm(z1,z2)
		} # forecastEN
		{
			z2 = tcrossprod(outcome[!exclude],cvrtQRTRAIN) %*% cvrtQRTEST;
			forecastLM[1,exclude] = forecastLM[1,exclude] + z2;
			forecastLM[2,exclude] = forecastLM[2,exclude] + 1;
			rm(z2);
		} # forecastLM
		{
			# get p-values
			fm = fm.open(paste0(dircvmwas, "/Stats_and_pvalues"))
			# colnames(fm)
			pv = fm[,3];
			close(fm)
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
			outcomeR = outcome - crossprod(cvrtQRFULL, cvrtQRTRAIN %*% outcome[!exclude ]);
			weights = cvrtQRTRAIN %*% coverage[!exclude,];
			coverageTRAIN = coverage[!exclude,] - crossprod(cvrtQRTRAIN, weights);
			coverageTEST  = coverage[ exclude,] - crossprod(cvrtQRTEST,  weights);
			rm(coverage, weights);
			gc();
		} # coverageTRAIN, coverageTEST, outcomeR, -coverage
		{
			z1 = cv.glmnet(x = coverageTRAIN, 
								y = as.vector(outcomeR[!exclude]), 
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
							  forecastRESID   = forecast0[1,]/forecast0[2,],
							  forecastEN      = forecast0[1,]/forecast0[2,] + forecastEN[1,]/forecastEN[2,], 
							  forecastLM      = forecast0[1,]/forecast0[2,] + forecastLM[1,]/forecastLM[2,],
							  forecastCVRT_EN = forecastEN[1,]/forecastEN[2,],
							  forecastCVRT_LM = forecastLM[1,]/forecastLM[2,]
		)
		write.table( file = sprintf("%s/MMCV4_prediction_folds=%02d_CpGs=%06d_alpha=%s.txt", 
											 param$dircv, param$cvnfolds, param$mmncpgs, param$mmalpha),
						 x = rez,
						 sep = "\t", row.names = FALSE);
	} # rez 
	{
		pdf( sprintf("%s/MMCV4_prediction_folds=%02d_CpGs=%06d_alpha=%s.pdf",
						 param$dircv, param$cvnfolds, param$mmncpgs, param$mmalpha) );
		plotPrediction(param, outcome, rez$forecastEN, 
							main = "Prediction success (EN on cvrt + EN on resid-d coverage)");
		plotPrediction(param, outcome, rez$forecastLM, 
							main = "Prediction success (LM on cvrt + EN on resid-d coverage)");
		plotPrediction(param, outcome, rez$forecastRESID + mean(outcome), 
							main = "Prediction success (EN on resid-d coverage only + mean)");
		plotPrediction(param, outcome, rez$forecastCVRT_EN, 
							main = "Prediction success (EN on cvrt only)");
		plotPrediction(param, outcome, rez$forecastCVRT_LM, 
							main = "Prediction success (LM on cvrt only)");
		dev.off();
	} # pdf()
	return( invisible(rez) );
}

ramwas7multiMarker = function(param){
	param = parameterPreprocess(param);
	message("Working in: ",param$dircv)
	if( length(param$modelcovariates) == 0){
		ramwas7multiMarkerNoCvrt(param);
	} else {
		ramwas7multiMarkerWithCvrt(param);
	}
}
