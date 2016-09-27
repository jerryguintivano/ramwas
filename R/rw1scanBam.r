pipelineProcessBam = function(bamname, param){
	# Used parameters: scoretag, minscore, filecpgset, maxrepeats
	
	param = parameterPreprocess(param);
	
	if( !is.null(param$filecpgset) && is.null(param$maxfragmentsize) )
		return("Parameter not set: maxfragmentsize");
	
	bamname = gsub("\\.bam$", "", bamname, ignore.case = TRUE);
	bamfullname = .makefullpath(param$dirbam, paste0(bamname,".bam"))
	
	dir.create(param$dirrbam, showWarnings = FALSE, recursive = TRUE)
	dir.create(param$dirrqc, showWarnings = FALSE, recursive = TRUE)
	
	rdsbmfile = paste0( param$dirrbam, "/", basename(bamname), ".rbam.rds" );
	rdsqcfile = paste0( param$dirrqc, "/", basename(bamname), ".qc.rds" );
	
	savebam = TRUE;
	if( file.exists( rdsbmfile ) ) {
		if( param$recalculate.QCs ) {
			### Precache the input rds file
			{invisible(readBin( rdsbmfile, "raw", file.size(rdsbmfile)));}
			rbam = readRDS(rdsbmfile); 
			savebam = FALSE;
		} else {
			return(paste0("Rbam rds file already exists: ",rdsqcfile));
		}
	} else {
		if( !file.exists( bamfullname ) )
			return(paste0("Bam file does not exist: ",bamfullname));
		rbam = bam.scanBamFile(bamfilename = bamfullname, scoretag = param$scoretag, minscore = param$minscore);
	}
	
	rbam2 = bam.removeRepeats(rbam, param$maxrepeats);
	rbam2 = bam.chrXY.qc(rbam2); 
	rbam2$qc$nbams = 1L;
	
	if( !is.null(param$filecpgset) ) {
		cpgset = cachedRDSload(param$filecpgset);
		noncpgset = cachedRDSload(param$filenoncpgset);
		isocpgset = isocpgSitesFromCpGset(cpgset = cpgset, distance = param$maxfragmentsize);
		rbam3 = bam.hist.isolated.distances(rbam = rbam2, isocpgset = isocpgset, distance = param$maxfragmentsize);
		rbam4 = bam.coverage.by.density(rbam = rbam3, cpgset = cpgset, noncpgset = noncpgset,
												  minfragmentsize = param$minfragmentsize, maxfragmentsize = param$maxfragmentsize);
		rbam5 = bam.count.nonCpG.reads(rbam = rbam4, cpgset = cpgset, distance = param$maxfragmentsize);
		
		### QC plots
		pipelineSaveQCplots(param, rbam5, bamname);
		
	} else {
		rbam5 = rbam2;
	}
	# .qc qcmean(rbam5$qc$chrX.count)  rbam5$qc$chrX.count[1]/rbam5$qc$chrX.count[2]
	# message(.qcTextLine(rbam5$qc, bamname))
	
	if(savebam)
		saveRDS( object = rbam5, file = rdsbmfile, compress = "xz");
	rbam6 = rbam5;
	rbam6$startsfwd=NULL;
	rbam6$startsrev=NULL;
	saveRDS( object = rbam6, file = rdsqcfile, compress = "xz");
	
	return(paste0("OK. ", bamname));
}

.ramwas1scanBamJob = function(bamname, param){
	cat(file = paste0(param$dirfilter,"/Log.txt"),
		 date(), ", Process ", Sys.getpid(),", Processing BAM: ", bamname, "\n", sep = "", append = TRUE);
	pipelineProcessBam(bamname = bamname, param = param);
}

ramwas1scanBams = function( param ){
	param = parameterPreprocess(param);
	stopifnot( !is.null(param$bamnames));
	
	dir.create(param$dirfilter, showWarnings = FALSE, recursive = TRUE)
	cat(file = paste0(param$dirfilter,"/Log.txt"), 
		 date(), ", Scanning bams.", "\n", sep = "", append = FALSE);
	if( param$cputhreads > 1) {
		cl = makeCluster(param$cputhreads);
		z = clusterApplyLB(cl, param$bamnames, .ramwas1scanBamJob, param = param); #[1:64]
		stopCluster(cl);
	} else {
		z = character(length(param$bamnames));
		names(z) = param$bamnames;
		for(i in seq_along(param$bamnames)) { # i=1
			z[i] = .ramwas1scanBamJob(bamname = param$bamnames[i], param = param);
		}
	}
	cat(file = paste0(param$dirfilter,"/Log.txt"), 
		 date(), ", Done scanning bams.", "\n", sep = "", append = TRUE);
	return(z);
}
