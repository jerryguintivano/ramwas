ramwas6crossValidation = function(param) {
	param = parameterPreprocess(param);
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

