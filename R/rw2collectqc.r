.combine.bams.qc = function( bamlist ){
	# bamlist = curbams
	if(length(bamlist)==1)
		return(bamlist[[1]]);
	
	### Deal with QCs
	qclist = lapply(bamlist, `[[`, "qc");
	
	qcnames = lapply(qclist, names);
	qcnames = unique(unlist(qcnames, use.names = FALSE))
	
	bigqc = vector("list", length(qcnames));
	names(bigqc) = qcnames;
	
	for( nm in qcnames) { # nm = qcnames[1]
		bigqc[[nm]] = Reduce(`%add%`, lapply(qclist, `[[`, nm));
	}
	return(list(qc = bigqc));
}

ramwas2collectqc = function( param ){
	param = parameterPreprocess(param);
	dir.create(param$dirqc, showWarnings = FALSE, recursive = TRUE);
	
	parameterDump(dir = param$dirqc, param = param,
					  toplines = c("dirrqc",
					  				 "filebamlist", "filebam2sample",
					  				 "bam2sample", "bamnames",
					  				 "scoretag", "minscore",
					  				 "minfragmentsize", "maxfragmentsize", "maxrepeats",
					  				 "filecpgset", "filenoncpgset"));
	
	bams = NULL;
	if( !is.null(param$bamnames) )
		bams = c(bams, param$bamnames);
	if( !is.null(param$bam2sample) )
		bams = c(bams, unlist(param$bam2sample, use.names = FALSE));
	
	bams = unique(basename(bams));
	
	message("Load BAM QC info");
	rbamlist = vector("list", length(bams));
	names(rbamlist) = bams;
	for( bamname in bams) {
		rdsqcfile = paste0( param$dirrqc, "/", bamname, ".qc.rds" );
		if(file.exists(rdsqcfile))
			rbamlist[[bamname]] = readRDS(rdsqcfile);
	}
	
	collect.qc.summary = function(bamset, dirname) {
		dirloc = paste0(param$dirqc, "/", dirname);
		dir.create(dirloc, showWarnings = FALSE, recursive = TRUE);
		
		message("Saving text summary");
		bigqc = vector("list", length(bamset));
		names(bigqc) = names(bamset);
		text = character(length(bamset));
		textR = character(length(bamset));
		for( ibam in seq_along(bamset) ) { # ibam=7
			curbams = rbamlist[bamset[[ibam]]];
			qc = .combine.bams.qc(curbams)$qc;
			if( length(qc) > 0 ) {
				bigqc[[ibam]] = qc;
			} else {
				qc = NULL;
			}
			text[ibam]  = .qcTextLine(  qc, names(bamset)[ibam] );
			textR[ibam] = .qcTextLineR( qc, names(bamset)[ibam] );
		}
		writeLines(con = paste0(dirloc, "/Summary_QC.txt"),   text = c(.qcTextHeader,  text));
		writeLines(con = paste0(dirloc, "/Summary_QC_R.txt"), text = c(.qcTextHeaderR, textR));
		
		figfun = function(qcname, plotname) {
			message("Saving plots ", plotname);
			pdf(paste0(dirloc,"/Fig_",plotname,".pdf"));
			for( ibam in seq_along(bamset) ) {
				plotinfo = bigqc[[ibam]][[qcname]];
				if( !is.null(bigqc[[ibam]][[qcname]]))
					plot(plotinfo, samplename = names(bamset)[ibam]);
				rm(plotinfo);
			}
			dev.off();
		}
		
		figfun(qcname = "hist.score1", plotname = "score");
		figfun(qcname = "bf.hist.score1", plotname = "score_before_filter");
		figfun(qcname = "hist.edit.dist1", plotname = "edit_distance");
		figfun(qcname = "bf.hist.edit.dist1", plotname = "edit_distance_before_filter");
		figfun(qcname = "hist.length.matched", plotname = "matched_length");
		figfun(qcname = "bf.hist.length.matched", plotname = "matched_length_before_filter");
		figfun(qcname = "hist.isolated.dist1", plotname = "isolated_distance");
		figfun(qcname = "avg.coverage.by.density", plotname = "coverage_by_density");
		return(invisible(bigqc));
	}
	
	# By bam
	bamset = bams; 
	names(bamset) = bams; 
	dirname = "summary_bams";
	message("Saving QC info by BAM");
	collect.qc.summary(bamset, dirname);
	rm(bamset, dirname);
	
	if( !is.null(param$bam2sample) ) {
		# by sample
		message("Saving QC info by BAMs in bam2sample");
		collect.qc.summary(bamset = unlist(param$bam2sample), dirname = "summary_bams_in_bam2sample");
		message("Saving QC info by SAMPLE");
		collect.qc.summary(bamset = param$bam2sample, dirname = "summary_by_sample");
		message("Saving QC info TOTAL (all BAMs in bam2sample)");
		bigqc = collect.qc.summary(bamset = list(total=unique(unlist(param$bam2sample, use.names = FALSE))), dirname = "summary_total");
	} else {
		message("Saving QC info TOTAL (all BAMs in bamnames/filebamlist)");
		bigqc = collect.qc.summary(bamset = list(total=bams), dirname = "summary_total");
	}
	
	### Fragment size
	frdata = bigqc$total$hist.isolated.dist1;
	estimate = estimateFragmentSizeDistribution(frdata, param$minfragmentsize);
	writeLines(con = paste0(param$dirfilter,"/Fragment_size_distribution.txt"), text = as.character(estimate));
	
	pdf(paste0(param$dirqc,"/Fragment_size_distribution_estimate.pdf"),8,8);
	lz = lm(frdata[seq_along(estimate)] ~ estimate)
	plot(as.vector(frdata)/1000, pch = 19, col="blue", main="Isolated CpG coverage vs.\nfragment size distribution estimate", ylab="count, thousands", xlab="Distance to isolated CpGs", xaxs="i");
	lines((estimate*lz$coefficients[2]+lz$coefficients[1])/1000, lwd = 4, col="red");
	dev.off();
	return(invisible(NULL));
}
