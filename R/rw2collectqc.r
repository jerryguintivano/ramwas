### Make text line for the QC set
.qcTextHeaderT = {paste(sep = "\t",
							  "Sample",
							  "# BAMs",
							  "Total reads",
							  "Reads aligned\tRA % of total",
							  "Reads after filter\tRAF % of aligned",
							  "Reads removed as duplicate\tRRAD % of aligned",
							  "Reads used for coverage\tRUFC % of aligned",
							  "Forward strand (%)",
							  "Avg alignment score",
							  "Avg aligned length",
							  "Avg edit distance",
							  "Non-CpG reads (%)",
							  "Avg non-CpG coverage",
							  "Avg CpG coverage",
							  "Non-Cpg/CpG coverage ratio",
							  "ChrX reads (%)",
							  "ChrY reads (%)",
							  "Peak SQRT")};
.qcTextHeaderR = {paste(sep = "\t",
								"Sample",
								"NBAMs",
								"TotalReads",
								"ReadsAligned\tReadsAlignedPct",
								"ReadsAfterFilter\tReadsAfterFilterPct",
								"ReadsRemovedAsDuplicate\tReadsRemovedAsDuplicatePct",
								"ReadsUsedForCoverage\tReadsUsedForCoveragePct",
								"ForwardStrandPct",
								"AvgAlignmentScore",
								"AvgAlignedLength",
								"AvgEditDistance",
								"NonCpGreadsPct",
								"AvgNonCpGcoverage",
								"AvgCpGcoverage",
								"NonCpg2CpGcoverageRatio",
								"ChrXreadsPct",
								"ChrYreadsPct",
								"PeakSQRT")};

.qccols = length(strsplit(.qcTextHeaderT,"\t",fixed = TRUE)[[1]])
.qcTextLineT = function(qc){
	name = qc$name;
	if( is.null(qc) )
		return(paste0(name,paste0(rep("\tNA",.qccols-1), collapse = "")));
	afracb = function(a,b) { sprintf("%s\t%.1f%%",s(a),100*a/b) };
	perc = function(x){ sprintf("%.2f%%",100*x) };
	twodig = function(x){ sprintf("%.2f",x) };
	s = function(x)formatC(x=x,digits=ceiling(log10(max(x)+1)),big.mark=",",big.interval=3);
	
	if( !is.null(qc$hist.score1))							class(qc$hist.score1) = "qcHistScore";
	if( !is.null(qc$bf.hist.score1))			 			class(qc$bf.hist.score1) = "qcHistScoreBF";
	if( !is.null(qc$hist.edit.dist1))					class(qc$hist.edit.dist1) = "qcEditDist";
	if( !is.null(qc$bf.hist.edit.dist1))				class(qc$bf.hist.edit.dist1) = "qcEditDistBF";
	if( !is.null(qc$hist.length.matched))				class(qc$hist.length.matched) = "qcLengthMatched";
	if( !is.null(qc$bf.hist.length.matched))			class(qc$bf.hist.length.matched) = "qcLengthMatchedBF";
	if( !is.null(qc$frwrev) ) 								class(qc$frwrev) = "qcFrwrev";
	
	rez = paste( sep = "\t",
					 name, # Sample
					 if(is.null(qc$nbams)){1}else{qc$nbams}, # Number of BAMs
					 s(qc$reads.total), # Total reads
					 afracb(qc$reads.aligned, qc$reads.total), # Reads aligned, % of total
					 afracb(qc$reads.recorded, qc$reads.aligned), # Reads after filter, % of aligned
					 afracb(qc$reads.recorded - qc$reads.recorded.no.repeats, qc$reads.aligned), # Reads removed as repeats\tRRAR % of aligned
					 afracb(qc$reads.recorded.no.repeats, qc$reads.aligned), # Reads used for coverage, % of aligned
					 perc(qcmean( qc$frwrev.no.repeats )), # Forward strand (%)
					 twodig(qcmean( qc$hist.score1 )), # Avg alignment score
					 twodig(qcmean( qc$hist.length.matched )), # Avg aligned length
					 twodig(qcmean( qc$hist.edit.dist1 )), # Avg edit distance
					 perc(qcmean( qc$cnt.nonCpG.reads )), # Non-CpG reads (%)
					 twodig( qc$avg.noncpg.coverage ), # Avg non-CpG coverage
					 twodig( qc$avg.cpg.coverage ), # Avg CpG coverage
					 perc( qc$avg.noncpg.coverage / qc$avg.cpg.coverage), # Non-Cpg/CpG coverage ratio
					 perc(qcmean( qc$chrX.count )), # ChrX reads (%)
					 perc(qcmean( qc$chrY.count )), # ChrY reads (%)
					 twodig(qcmean( qc$avg.coverage.by.density )) # Peak SQRT
	);
	# message(rez);
	return(rez);
}
.qcTextLineR = function(qc){
	name = qc$name;
	if( is.null(qc) )
		return(paste0(name,paste0(rep("\tNA",.qccols-1), collapse = "")));
	afracb = function(a,b) { paste0(a, "\t", a/b) };
	perc = identity;
	twodig = identity;
	s = identity;
	
	if( !is.null(qc$hist.score1))							class(qc$hist.score1) = "qcHistScore";
	if( !is.null(qc$bf.hist.score1))			 			class(qc$bf.hist.score1) = "qcHistScoreBF";
	if( !is.null(qc$hist.edit.dist1))					class(qc$hist.edit.dist1) = "qcEditDist";
	if( !is.null(qc$bf.hist.edit.dist1))				class(qc$bf.hist.edit.dist1) = "qcEditDistBF";
	if( !is.null(qc$hist.length.matched))				class(qc$hist.length.matched) = "qcLengthMatched";
	if( !is.null(qc$bf.hist.length.matched))			class(qc$bf.hist.length.matched) = "qcLengthMatchedBF";
	if( !is.null(qc$frwrev) ) 								class(qc$frwrev) = "qcFrwrev";
	
	rez = paste( sep = "\t",
					 name, # Sample
					 if(is.null(qc$nbams)){1}else{qc$nbams}, # Number of BAMs
					 s(qc$reads.total), # Total reads
					 afracb(qc$reads.aligned, qc$reads.total), # Reads aligned, % of total
					 afracb(qc$reads.recorded, qc$reads.aligned), # Reads after filter, % of aligned
					 afracb(qc$reads.recorded - qc$reads.recorded.no.repeats, qc$reads.aligned), # Reads removed as repeats\tRRAR % of aligned
					 afracb(qc$reads.recorded.no.repeats, qc$reads.aligned), # Reads used for coverage, % of aligned
					 perc(qcmean( qc$frwrev.no.repeats )), # Forward strand (%)
					 twodig(qcmean( qc$hist.score1 )), # Avg alignment score
					 twodig(qcmean( qc$hist.length.matched )), # Avg aligned length
					 twodig(qcmean( qc$hist.edit.dist1 )), # Avg edit distance
					 perc(qcmean( qc$cnt.nonCpG.reads )), # Non-CpG reads (%)
					 twodig( qc$avg.noncpg.coverage ), # Avg non-CpG coverage
					 twodig( qc$avg.cpg.coverage ), # Avg CpG coverage
					 perc( qc$avg.noncpg.coverage / qc$avg.cpg.coverage), # Non-Cpg/CpG coverage ratio
					 perc(qcmean( qc$chrX.count )), # ChrX reads (%)
					 perc(qcmean( qc$chrY.count )), # ChrY reads (%)
					 twodig(qcmean( qc$avg.coverage.by.density )) # Peak SQRT
	);
	return(rez);
}
if(FALSE){
	name = "My_bam";
	rbam = readRDS("D:/Cell_type/rds_rbam/150114_WBCS014_CD20_150.rbam.rds");
	rbam$qc$avg.cpg.coverage
	cpgset = cachedRDSload("C:/AllWorkFiles/Andrey/VCU/RaMWAS_2/code/Prepare_CpG_list/hg19/cpgset_hg19_SNPS_at_MAF_0.05.rds");
	minfragmentsize = 50;
	maxfragmentsize = 200;
	rbam = bam.coverage.by.density(rbam, cpgset, minfragmentsize, maxfragmentsize);
	rbam$qc$avg.cpg.coverage
	qc = rbam$qc;
	class(qc$frwrev) = "qcFrwrev";
	class(qc$frwrev.no.repeats) = "qcFrwrev";
	class(qc$cnt.nonCpG.reads) = "qcNonCpGreads"
	qc$nbams = 1;
	
	cat( qcTextLine(qc, "bam name"), "\n")
	
	as.matrix(names(qc))
	
	with(qc, paste0("avg.noncpg.coverage", avg.noncpg.coverage, "\n"));
	with(qc, is.null(qc$nbams))
}

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

loadBamQC = function(param, bams){
	rbamlist = vector("list", length(bams));
	names(rbamlist) = bams;
	for( bamname in bams) { # bamname = bams[1]
		rdsqcfile = paste0( param$dirrqc, "/", bamname, ".qc.rds" );
		if(file.exists(rdsqcfile)) {
			rbamlist[[bamname]] = readRDS(rdsqcfile);
		} else {
			message("QC file not found: ",rdsqcfile)
		}
	}
	return(rbamlist)
}
combineBamQcIntoSamples = function(rbamlist, bamset){
	bigqc = vector("list", length(bamset));
	names(bigqc) = names(bamset);
	for( ibam in seq_along(bamset) ) { # ibam=1
		curbams = rbamlist[bamset[[ibam]]];
		qc = .combine.bams.qc(curbams)$qc;
		if( length(qc) > 0 ) {
			bigqc[[ibam]] = qc;
		} else {
			bigqc[[ibam]] = list();
		}
		bigqc[[ibam]]$name = names(bamset)[ibam];
	}
	return(bigqc);
}

estimateFragmentSizeDistribution = function(hist.isolated.distances, seqLength){
	
	if( length(hist.isolated.distances) == seqLength )
		return( rep(1, seqLength) );
	
	### Point of crossing the middle
	ytop = median(hist.isolated.distances[1:seqLength]);
	ybottom = median(tail(hist.isolated.distances,seqLength));
	ymidpoint = ( ytop + ybottom )/2;
	yrange = ( ytop - ybottom );
	overymid = (hist.isolated.distances > ymidpoint)
	xmidpoint = which.max( cumsum( overymid - mean(overymid) ) );
	
	### interquartile range estimate
	xqrange = 
		which.max(cumsum( ((hist.isolated.distances > quantile(hist.isolated.distances,0.25))-0.75) )) -
		which.max(cumsum( ((hist.isolated.distances > quantile(hist.isolated.distances,0.75))-0.25) ));
	
	logitrange = diff(qlogis(c(0.25,0.75)));
	
	initparam = c(xmidpoint = xmidpoint, 
					  xdivider = (xqrange/logitrange)/2, 
					  ymultiplier = yrange, 
					  ymean = ybottom);
	
	fsPredict = function( x, param) {
		(plogis((param[1]-x)/param[2]))*param[3]+param[4]
	}
	
	x = seq_along(hist.isolated.distances);
	
	# plot( x, hist.isolated.distances)
	# lines( x, fsPredict(x, initparam), col="blue", lwd = 3)
	
	fmin = function(param) {
		fit2 = fsPredict(x, param); 
		# (plogis((param[1]-x)/param[2]))*param[3]+param[4];
		error = hist.isolated.distances - fit2;
		e2 = error^2;
		e2s = sort.int(e2,decreasing = TRUE);
		return(sum(e2s[-(1:10)]));
	}
	
	estimate = optim(par = initparam, fn = fmin, method = "BFGS");
	param = estimate$par;
	fit = fsPredict(x, param);
	
	rezfit = plogis((param[1]-x)/param[2]);
	keep = rezfit>0.05;
	keep[length(keep)] = FALSE;
	rezfit = rezfit - max(rezfit[!keep],0)
	rezfit[1:seqLength] = rezfit[seqLength];
	rezfit = rezfit[keep];
	rezfit = rezfit / rezfit[1];
	
	# lz = lm(hist.isolated.distances[seq_along(rezfit)] ~ rezfit)
	# lines(rezfit*lz$coefficients[2]+lz$coefficients[1], lwd = 4, col="red");
	
	return(rezfit);
}
if(FALSE){ # test code
	x = seq(0.01,0.99,0.01);
	y = sqrt(abs(x-0.5))*sign(x-0.5)
	plot(x,y)
	log.ss = nls(y ~ SSlogis(x, phi1, phi2, phi3))
	z = SSlogis(x, 0.59699, 0.61320, 0.04599)
	lines(x, z, col="blue")
	
	
	# setwd("D:/RW/NESDA/ramwas/AS120_sba/");
	setwd("D:/RW/RC2/ramwas/AS38_gap0/");
	# setwd("D:/RW/Celltype//ramwas/AS120_sba/");
	lst = list.files(pattern = "\\.qc\\.");
	qcs = lapply(lst, function(x){load(x);return(bam);})
	histinfo = Reduce( `+`, lapply( lapply(qcs, `[[`, "qcflt"), `[[`, "hist.iso.dist.250"), init = 0);
	rng = range(histinfo[-(1:10)]);
	plot(histinfo/1e3, ylim = rng/1e3, pch=19, col="blue")
	
	hist.isolated.distances = histinfo;
	seqLength = 50;
	
	fit = estimateFragmentSizeDistribution(hist.isolated.distances, seqLength)
	
	x = seq_along(hist.isolated.distances)
	plot( x, hist.isolated.distances)
	lz = lm(hist.isolated.distances[seq_along(fit)] ~ fit)
	lines(fit*lz$coefficients[2]+lz$coefficients[1], lwd = 4, col="red");
	
}


# pipelineEstimateFragmentSizeDistribution = function( param ){
# 	
# 	param = parameterPreprocess(param);
# 	
# 	if( !is.null(param$bam2sample) ) {
# 		bams = unlist( param$bam2sample, use.names = FALSE);
# 	} else if (!is.null(param$bamnames)) {
# 		bams = param$bamnames;
# 	} else {
# 		stop("Bams are not defined. Set filebam2sample, filebamlist, bam2sample or bamnames.","\n");
# 	}
# 	bams = unique(basename(bams));
# 	
# 	qclist = vector("list", length(bams));
# 	names(qclist) = bams;
# 	
# 	for( bamname in bams) {
# 		rdsqcfile = paste0( param$dirrqc, "/", bamname, ".qc.rds" );
# 		qclist[[bamname]] = readRDS(rdsqcfile);
# 	}
# 	
# 	qcset = lapply(lapply( qclist, `[[`, "qc"),`[[`,"hist.isolated.dist1")
# 	bighist = Reduce(`%add%`, qcset);
# 	estimate = estimateFragmentSizeDistribution(bighist, param$minfragmentsize);
# 	
# 	writeLines(con = paste0(param$dirfilter,"/Fragment_size_distribution.txt"), text = as.character(estimate));
# 	
# 	return(estimate);
# }

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
	{
		bams = NULL;
		if( !is.null(param$bamnames) )
			bams = c(bams, param$bamnames);
		if( !is.null(param$bam2sample) )
			bams = c(bams, unlist(param$bam2sample, use.names = FALSE));
		bams = unique(basename(bams));
	} # bams
	
	{
		message("Load BAM QC info");
		rbamlist = loadBamQC(param, bams);
	}
	
	collect.qc.summary = function(bamset, dirname) {
		dirloc = paste0(param$dirqc, "/", dirname);
		dir.create(dirloc, showWarnings = FALSE, recursive = TRUE);
		
		bigqc = combineBamQcIntoSamples(rbamlist = rbamlist, bamset = bamset);
		saveRDS(file = paste0(dirloc,"/qclist.rds"), object = bigqc);
		{
			textT = sapply(bigqc, .qcTextLineT)
			textR = sapply(bigqc, .qcTextLineR)
			writeLines(con = paste0(dirloc, "/Summary_QC.txt"),   text = c(.qcTextHeaderT, textT));
			writeLines(con = paste0(dirloc, "/Summary_QC_R.txt"), text = c(.qcTextHeaderR, textR));
			rm(textT, textR);
		} # text summary
		
		histqc = function(qcfun, plottitle, filename){
			vec = unlist(lapply(bigqc, qcfun));
			pdf(paste0(dirloc,"/Fig_hist_",filename,".pdf"));
			hist(vec, breaks = 3*round(sqrt(length(vec))), main = plottitle, col = 'lightblue', xlab = 'value', yaxs = 'i')
			dev.off()
		}
		
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
		# bigqc[[1]]$cnt.nonCpG.reads
		if(length(bigqc) >= 10) {
			histqc(qcfun = function(x){x$avg.cpg.coverage / x$avg.noncpg.coverage},
					 plottitle = "Enrichment lower bound\n(Avg CpG / avg non-CpG coverage)",
					 filename = "enrichment")
			histqc(qcfun = function(x){x$avg.noncpg.coverage / x$avg.cpg.coverage * 100},
					 plottitle = "Background noise level\n(Avg non-CpG / avg CpG coverage, %)",
					 filename = "noise")
			histqc(qcfun = function(x){x$reads.recorded.no.repeats/1e6},
					 plottitle = "Number of reads after filters, millions",
					 filename = "Nreads")
			histqc(qcfun = function(x){qcmean(x$hist.edit.dist1)},
					 plottitle = "Average edit distance of aligned reads",
					 filename = "edit_dist")
			histqc(qcfun = function(x){qcmean(x$avg.coverage.by.density)},
					 plottitle = "CpG density at the peak sensitivity (SQRT)\n(a value per BAM / sample)",
					 filename = "peak")
			histqc(qcfun = function(x){qcmean(x$cnt.nonCpG.reads)},
					 plottitle = "Fraction of reads not covering any CpGs",
					 filename = "noncpg_reads")
		}
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