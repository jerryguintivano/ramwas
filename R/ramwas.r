### Caching environment
.ramwasEnv = new.env()

`%add%` <- function(x, y) {
	if(is.null(x)) return(y);
	if(is.null(y)) return(x);
	l <- max(length(x), length(y))
	length(x) <- l
	length(y) <- l
	x[is.na(x)] <- 0
	y[is.na(y)] <- 0
	return(x + y)
}
.isAbsolutePath = function(pathname) {
	if( grepl("^~/", pathname) ) 
		return(TRUE)
	if( grepl("^.:(/|\\\\)", pathname) ) 
		return(TRUE)
	if( grepl("^(/|\\\\)", pathname) ) 
		return(TRUE)
	return(FALSE);
}
makefullpath = function(path, filename) {
	if( is.null(path) )
		return(filename);
	if(.isAbsolutePath(filename)) {
		return(filename) 
	} else {
		return( paste0(path, "/", filename) );
	}
}
if(FALSE) {
	.isAbsolutePath( "C:/123" );  # TRUE
	.isAbsolutePath( "~123" );    # FALSE
	.isAbsolutePath( "~/123" );   # TRUE
	.isAbsolutePath( "/123" );    # TRUE
	.isAbsolutePath( "\\123" );    # TRUE
	.isAbsolutePath( "asd\\123" ); # FALSE
	.isAbsolutePath( "a\\123" );   # FALSE
}

### Scan a file for parameters
parametersFromFile = function(.parameterfile){
	source(.parameterfile, local = TRUE);
	.nms = ls();
	return(mget(.nms));
}
if(FALSE) { # test code
	param = parametersFromFile(.parameterfile = "D:/RW/NESDA/ramwas/param_file.txt");
	param
}

parseBam2sample = function(lines) {
	# remove trailing commas
	lines = gsub(pattern = ",$", replacement = "", lines);
	
	lines = gsub(pattern = "\\.bam,", replacement = ",", lines);
	lines = gsub(pattern = "\\.bam$", replacement = "",  lines);
	
	split.eq = strsplit(lines, split = "=", fixed = TRUE);
	samplenames = sapply(split.eq, `[`, 1);
	bamlist = strsplit(sapply(split.eq,tail,1), split = ",", fixed = TRUE);
	names(bamlist) = samplenames;
	
	# bamvec = unlist(bamlist, use.names = FALSE)
	# bamlist = lapply(bamlist, basename);
	return(bamlist);
}

### Fill in gaps in the parameter list
parameterPreprocess = function(param) {
	### Get from a file if param is not a list
	if(is.character(param)) {
		param = parametersFromFile(param);
	}
	
	# Set up directories 
	if( is.null(param$dirproject)) {
		param$dirproject = ".";
	}
	if( is.null(param$dirfilter) ) {
		param$dirfilter = FALSE;
	}
	if( is.logical(param$dirfilter) ) {
		if( param$dirfilter ) {
			param$dirfilter = paste0( param$dirproject, "/Filter_", param$scoretag, "_", param$minscore);
		} else {
			param$dirfilter = param$dirproject;
		}
	} else {
		param$dirfilter = makefullpath(param$dirproject, param$dirfilter);
	}
	if( is.null(param$dirrbam) ) param$dirrbam = paste0( param$dirfilter, "/rds_rbam");
	if( is.null(param$dirrqc) ) param$dirrqc = paste0( param$dirfilter, "/rds_qc");
	if( is.null(param$dirqc) ) param$dirqc = paste0( param$dirfilter, "/qc");
	if( is.null(param$dirlog) ) param$dirlog = paste0( param$dirfilter, "/log");
	if( is.null(param$dirtemp) ) param$dirtemp  = "coverage_raw_temp";
	param$dirtemp  = makefullpath(param$dirproject, param$dirtemp );
	if( is.null(param$dircoveragenorm) ) param$dircoveragenorm = "coverage_norm";
	param$dircoveragenorm = makefullpath(param$dirproject, param$dircoveragenorm);
	
	
	
	### Filter parameters
	if( is.null(param$scoretag) ) param$scoretag = "mapq";
	if( is.null(param$minscore) ) param$minscore = 4;
	
	### More analysis parameters
	if( is.null(param$maxrepeats) ) param$maxrepeats = 0;
	if( is.null(param$cputhreads) ) param$cputhreads = detectCores();
	
	### BAM list processing
	if( is.null(param$bamnames) & !is.null(param$filebamlist)) {
		
		param$bamnames = readLines(makefullpath(param$dirproject,param$filebamlist));
		param$bamnames = gsub(pattern = "\\.bam$", replacement = "", param$bamnames);
	}
	### BAM2sample processing
	if( !is.null(param$filebam2sample) & is.null(param$bam2sample)) {
		filename = makefullpath(param$dirproject, param$filebam2sample);
		param$bam2sample = parseBam2sample( readLines(filename) );
		rm(filename);
	}
	### Analysis variables file
	if( !is.null(param$fileanalysis) ) {
		sep = "\t";
		if(grepl("\\.csv$",param$fileanalysis)) 
			sep = ",";
		filename = makefullpath(param$dirproject, param$fileanalysis);
		param$covariates = read.table(filename, header = TRUE, sep = sep, stringsAsFactors = FALSE);
		rm(filename);
	}

	### CpG set should exist
	if( !is.null(param$filecpgset) ) {
		stopifnot( file.exists(param$filecpgset) );
	}
	if( is.null(param$doublesize) ) param$doublesize = 4;
	if( is.null(param$recalculate.QCs) ) param$recalculate.QCs = FALSE;
	if( is.null(param$buffersize) ) param$buffersize = 8e9;

	if( is.null(param$minavgcpgcoverage) ) param$minavgcpgcoverage = 0.3;
	if( is.null(param$minnonzerosamples) ) param$minnonzerosamples = 0.3;

	return(param);
}

###
### BAM processing
###

bam.scanBamFile = function( bamfilename, scoretag = "mapq", minscore = 4){
	
	# header = scanBamHeader(bamfilename)
	# chrnames = names(header[[1]]$targets)
	
	### constants
	### More than 500 alignments per read is unlikely (although very possible)
	# max.alignments.per.read = 500; 
	
	fields = c("qname","rname","pos","cigar","flag") 		
	# "qname" is read name, "rname" is chromosome
	tags = "NM";# character();
	
	### Open the BAM file
	{
		if(scoretag == "mapq") {
			fields = c(fields,scoretag);
		} else {
			tags = c(tags, scoretag);
		}
		
		flag = scanBamFlag(isUnmappedQuery=NA, isSecondMateRead=FALSE);
		param = ScanBamParam(flag=flag, what=fields, tag=tags);
		bf <- BamFile(bamfilename, yieldSize=1e6) ## typically, yieldSize=1e6
		open(bf);	
		rm(fields, tags, flag);
	} # bf, param
	
	qc = list(nbams = 1L);
	# qc.frwrev = c(0,0);
	# qc.reads = 0;
	# qc.aligned = 0;
	# qc.aligned.and.used = 0;
	# qc.hist.length.matched = 0;
	# qc.hist.edit.dist = 0;
	
	startlistfwd = NULL;
	# keep the tail of previous yield.
	# NULL - just started
	# FALSE - last iteration
	# list - in the process
	# oldtail = NULL; 
	repeat{
		### break condition
		### is at the begining of the loop to support calls of "next"
		# if(is.logical(oldtail))
		# 	break;
		
		### Read "yieldSize" rows
		bb = scanBam(file=bf, param=param)[[1]];
		if( length(bb[[1]])==0 )
			break;
		
		### Put tags in the main list
		bb = c(bb[names(bb) != "tag"], bb$tag);
		# data.frame(lapply(bb,`[`, 1:60), check.rows = FALSE, stringsAsFactors = FALSE)
		
		# stopifnot( length(bb[[scoretag]]) == length(bb[[1]]) )
		
		### Create output lists
		if(is.null(startlistfwd)) {
			startlistfwd = vector("list",length(levels(bb$rname)));
			startlistrev = vector("list",length(levels(bb$rname)));
			names(startlistfwd) = levels(bb$rname);
			names(startlistrev) = levels(bb$rname);
			startlistfwd[] = list(list())
			startlistrev[] = list(list())
		} # startlistfwd, startlistrev 
		
		# ### Cut the new tail and append the previous one
		# {
		# 	bbsize = length(bb[[1]]);
		# 	if( bbsize == 0L ) {
		# 		if(is.null(oldtail))
		# 			stop("Empty BAM file (?).");
		# 		bb = oldtail;
		# 		oldtail = FALSE;
		# 	} else {
		# 		taillength = sum(tail(bb$qname,max.alignments.per.read) == tail(bb$qname,1));
		# 		newtail = lapply(bb, `[`, (bbsize-taillength+1):bbsize);
		# 		bb = lapply(bb, `[`, 1:(bbsize-taillength));
		# 		if( !is.null(oldtail) ) {
		# 			bb = combine.2.lists(oldtail,bb);
		# 		}
		# 		oldtail = newtail;
		# 		rm(newtail, taillength);
		# 	}
		# 	rm(bbsize);
		# } # oldtail
		
		### Keep only primary reads
		{
			keep = bitwAnd(bb$flag, 256L) == 0L;
			if(!all(keep))
				bb = lapply(bb,`[`,which(keep));
			rm(keep);
		}
		
		qc$reads.total = qc$reads.total %add% length(bb[[1]]);
		
		### Keep only aligned reads
		{
			keep = bitwAnd(bb$flag, 4L) == 0L;
			if(!all(keep))
				bb = lapply(bb,`[`,which(keep));
			rm(keep);
		}
		
		if(length(bb[[1]])==0) {
			cat(sprintf("Recorded %.f of %.f reads",qc$reads.recorded,qc$reads.total),"\n");
			next;
		}
		
		bb$matchedAlongQuerySpace = cigarWidthAlongQuerySpace(bb$cigar,after.soft.clipping = TRUE);
		
		qc$reads.aligned = qc$reads.aligned %add% length(bb[[1]]);
		qc$bf.hist.score1 = qc$bf.hist.score1 %add% tabulate(pmax(bb[[scoretag]]+1L,1L));
		qc$bf.hist.edit.dist1 = qc$bf.hist.edit.dist1 %add% tabulate(bb$NM+1L);
		qc$bf.hist.length.matched = qc$bf.hist.length.matched %add% tabulate(bb$matchedAlongQuerySpace);
		
		### Keep score >= minscore
		if( !is.null(minscore) ) {
			score = bb[[scoretag]];
			keep = score >= minscore;
			keep[is.na(keep)] = FALSE;
			if(!all(keep))
				bb = lapply(bb,`[`,which(keep));
			rm(keep);
		}
		
		qc$reads.recorded = qc$reads.recorded %add% length(bb[[1]]);
		qc$hist.score1 = qc$hist.score1 %add% tabulate(pmax(bb[[scoretag]]+1L,1L));
		qc$hist.edit.dist1 = qc$hist.edit.dist1 %add% tabulate(bb$NM+1L);
		qc$hist.length.matched = qc$hist.length.matched %add% tabulate(bb$matchedAlongQuerySpace);
		
		### Forward vs. Reverse strand
		bb$isReverse = bitwAnd(bb$flag, 0x10) > 0;
		qc$frwrev = qc$frwrev %add% tabulate(bb$isReverse + 1L)
		
		
		### Read start positions (accounting for direction)
		{
			bb$startpos = bb$pos;
			bb$startpos[bb$isReverse] = bb$startpos[bb$isReverse] + 
				(cigarWidthAlongReferenceSpace(bb$cigar[bb$isReverse])-1L) - 1L;
			# Last -1L is for shift from C on reverse strand to C on the forward
		} # startpos
		
		### Split and store the start locations
		{
			offset = length(startlistfwd);
			split.levels = as.integer(bb$rname) + offset*bb$isReverse;
			levels(split.levels) = c(names(startlistfwd),paste0(names(startlistfwd),"-"));
			class(split.levels) = "factor";
			splt = split( bb$startpos, split.levels, drop = FALSE);
			# print(sapply(splt,length))
			for( i in seq_along(startlistfwd) ) {
				if( length(splt[i]) > 0 ) {
					startlistfwd[[i]][[length(startlistfwd[[i]])+1L]] = splt[[i]];
				}
				if( length(splt[i+offset]) > 0 ) {
					startlistrev[[i]][[length(startlistrev[[i]])+1L]] = splt[[i+offset]];
				}
			}
			rm(offset, split.levels, splt);
		} # startlistfwd, startlistrev	
		cat(sprintf("Recorded %.f of %.f reads",qc$reads.recorded,qc$reads.total),"\n");
	}
	close(bf);
	rm(bf); # , oldtail
	
	startsfwd = startlistfwd;
	startsrev = startlistrev;
	
	### combine and sort lists in "outlist"
	for( i in seq_along(startlistfwd) ) {
		startsfwd[[i]] = sort.int(unlist(startlistfwd[[i]]));
		startsrev[[i]] = sort.int(unlist(startlistrev[[i]]));
	}		
	gc();
	
	if( !is.null(qc$hist.score1))							class(qc$hist.score1) = "qcHistScore";
	if( !is.null(qc$bf.hist.score1))			 			class(qc$bf.hist.score1) = "qcHistScoreBF";
	if( !is.null(qc$hist.edit.dist1))					class(qc$hist.edit.dist1) = "qcEditDist";
	if( !is.null(qc$bf.hist.edit.dist1))				class(qc$bf.hist.edit.dist1) = "qcEditDistBF";
	if( !is.null(qc$hist.length.matched))				class(qc$hist.length.matched) = "qcLengthMatched";
	if( !is.null(qc$bf.hist.length.matched))			class(qc$bf.hist.length.matched) = "qcLengthMatchedBF";
	if( !is.null(qc$frwrev) ) 								class(qc$frwrev) = "qcFrwrev";
	
	info = list(bamname = bamfilename, scoretag = scoretag, minscore = minscore, filesize = file.size(bamfilename));
	rbam = list(startsfwd = startsfwd, startsrev = startsrev, qc = qc, info = info);
	return( rbam );
}
if(FALSE) { # test code
	### Rsamtools to read BAMs
	library(Rsamtools)
	### GenomicAlignments to process CIGAR strings
	library(GenomicAlignments)
	
	bamfilename = "D:/NESDA_07D00232.bam"; scoretag = "AS"; minscore = 60;
	rbam = bam.scanBamFile(bamfilename = bamfilename, scoretag = scoretag, minscore = minscore);

	bamfilename = "D:/01A04SM1429Qa.bam"
	scoretag = "AS"
	minscore = 38
	rbam = bam.scanBamFile(bamfilename = bamfilename, scoretag = scoretag, minscore = minscore);
	plot(rbam$qc$hist.score1)
	plot(rbam$qc$bf.hist.score1)
	
	
	rbam$info
	
	plot(rbam$qc$hist.score1)
	plot(rbam$qc$hist.score1, col="red")
	plot(rbam$qc$hist.score1, xstep=15)
	plot(rbam$qc$hist.edit.dist1)
	plot(rbam$qc$hist.edit.dist1, xstep=2)
	plot(rbam$qc$hist.length.matched)
	plot(rbam$qc$hist.length.matched, xstep=5)
	
	qcmean(rbam$qc$hist.score1)
	qcmean(rbam$qc$hist.edit.dist1)
	qcmean(rbam$qc$hist.length.matched)
}

.my.hist.plot = function(values, main2, firstvalue=0, xstep = 10, ...) {
	maxval = max(values);
	thresholds = c(-Inf, 1e3, 1e6, 1e9)*1.5;
	bin = findInterval(maxval, thresholds)
	switch(bin,
			 {ylab = "count"},
			 {ylab = "count, thousands"; values=values/1e3;},
			 {ylab = "count, millions"; values=values/1e6;},
			 {ylab = "count, billions"; values=values/1e9;}
	)
	param = list(...);
	plotparam = list(height = as.vector(values), width = 1, space = 0, 
						  col = "royalblue", border = "blue", 
						  main = main2, xaxs="i", yaxs="i", ylab = ylab);
	plotparam[names(param)] = param;
	do.call(barplot, plotparam);
	# barplot(, ...);
	at = seq(0, length(values)+xstep, xstep);
	at[1] = firstvalue;
	axis(1,at = at+0.5-firstvalue, labels = at)
}
plot.qcHistScore = function(x, samplename="", xstep = 25, ...) {
	.my.hist.plot(as.vector(x), main2 = paste0("Distribution of read scores\n",samplename), firstvalue=0, xstep = xstep, ...);
}
plot.qcHistScoreBF = function(x, samplename="", xstep = 25, ...) {
	.my.hist.plot(as.vector(x), main2 = paste0("Distribution of read scores\n(including excluded reads)\n",samplename), firstvalue=0, xstep = xstep, ...);
}
plot.qcEditDist = function(x, samplename="", xstep = 5, ...) {
	.my.hist.plot(as.vector(x), main2 = paste0("Distribution of edit distance\n",samplename), firstvalue=0, xstep = xstep, ...);
}
plot.qcEditDistBF = function(x, samplename="", xstep = 5, ...) {
	.my.hist.plot(as.vector(x), main2 = paste0("Distribution of edit distance\n(including excluded reads)\n",samplename), firstvalue=0, xstep = xstep, ...);
}
plot.qcLengthMatched = function(x, samplename="", xstep = 25, ...) {
	.my.hist.plot(as.vector(x), main2 = paste0("Distribution of length of aligned part of read\n",samplename), firstvalue=1, xstep = xstep, ...);
}
plot.qcLengthMatchedBF = function(x, samplename="", xstep = 25, ...) {
	.my.hist.plot(as.vector(x), main2 = paste0("Distribution of length of aligned part of read\n(including excluded reads)\n",samplename), firstvalue=1, xstep = xstep, ...);
}
plot.qcIsoDist = function(x, samplename="", xstep = 25, ...) {
	.my.hist.plot(as.vector(x), main2 = paste0("Distribution of distances from read starts to isolated CpGs\n",samplename), firstvalue=0, xstep = xstep, ...);
}
plot.qcCoverageByDensity = function(y, samplename="", ...) {
	# y = rbam$qc$avg.coverage.by.density
	x = (seq_along(y)-1)/100;
	param = list(...);
	plotparam = list(x = x, y = y, type = "l", col = "magenta", 
						  lwd = 3, xaxs="i", yaxs="i", axes=FALSE,
						  ylim = c(0, max(y, na.rm = TRUE)*1.1), xlim = range(x), 
						  xlab = "CpG density", ylab = "Coverage", 
						  main = paste0("Average coverage by CpG density\n",samplename));
	plotparam[names(param)] = param;
	do.call(plot, plotparam);
	axis(1, at = seq(0,tail(x,1)+2,by = 1), labels = seq(0,tail(x,1)+2,by=1)^2)
	axis(2);
}
.histmean = function(x) {
	return( sum(x * seq_along(x)) / pmax(sum(x),.Machine$double.xmin) );
}
qcmean <- function(x) UseMethod("qcmean", x)
qcmean.qcHistScore = function(x) { .histmean(x)-1 }
qcmean.qcHistScoreBF = function(x) { .histmean(x)-1 }
qcmean.qcEditDist = function(x) { .histmean(x)-1 }
qcmean.qcEditDistBF = function(x) { .histmean(x)-1 }
qcmean.qcLengthMatched = function(x) { .histmean(x) }
qcmean.qcLengthMatchedBF = function(x) { .histmean(x) }
qcmean.qcIsoDist = function(x) { .histmean(x) }
qcmean.qcFrwrev = function(x){ x[1]/(x[1]+x[2]) }
qcmean.qcNonCpGreads = function(x){ x[1]/(x[1]+x[2]) }
qcmean.qcCoverageByDensity = function(x){ (which.max(x)-1)/100 }
qcmean.qcChrX = function(x){ x[1]/x[2] }

### Make text line for the QC set
.qcTextHeader = {paste(sep = "\t",
	"Sample",
	"# BAMs",
	"Total reads",
	"Reads aligned\tRA % of total",
	"Reads after filter\tRAF % of aligned",
	"Reads removed as repeats\tRRAR % of aligned",
	"Reads used for coverage\tRUFC % of aligned",
	"Forward strand (%)",
	"Avg alignment score",
	"Avg aligned length",
	"Avg edit distance",
	"Non-CpG reads (%)",
	"Avg non-CpG coverage",
	"Avg CpG coverage",
	"Non-Cpg/CpG coverage ratio",
	"Reads on chrX (%)",
	"Peak SQRT")};
.qccols = length(strsplit(.qcTextHeader,"\t",fixed = TRUE)[[1]])
.qcTextLine = function(qc, name) {
	if( is.null(qc) )
		return(paste0(name,paste0(rep("\tNA",.qccols-1), collapse = "")));
	afracb = function(a,b) { paste0(s(a), "\t", round(100*a/b,1),"%") };
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
	
	rez = with(qc, paste( sep = "\t",
		name, # Sample
		if(is.null(qc$nbams)){1}else{qc$nbams}, # Number of BAMs
		s(reads.total), # Total reads
		afracb(reads.aligned, reads.total), # Reads aligned, % of total
		afracb(reads.recorded, reads.aligned), # Reads after filter, % of aligned
		afracb(reads.recorded - reads.recorded.no.repeats, reads.aligned), # Reads removed as repeats\tRRAR % of aligned
		afracb(reads.recorded.no.repeats, reads.aligned), # Reads used for coverage, % of aligned
		perc(qcmean( frwrev.no.repeats )), # Forward strand (%)
		twodig(qcmean( hist.score1 )), # Avg alignment score
		twodig(qcmean( hist.length.matched )), # Avg aligned length
		twodig(qcmean( hist.edit.dist1 )), # Avg edit distance
		perc(qcmean( cnt.nonCpG.reads )), # Non-CpG reads (%)
		twodig( avg.noncpg.coverage ), # Avg non-CpG coverage
		twodig( avg.cpg.coverage ), # Avg CpG coverage
		perc( avg.noncpg.coverage / avg.cpg.coverage), # Non-Cpg/CpG coverage ratio
		perc(qcmean( chrX.count )), # Reads on chrX (%)
		twodig(qcmean( avg.coverage.by.density )) # Peak SQRT
	));
	# cat(rez,"\n");
	return(rez);
}
if(FALSE) {
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

###
### BAM QC / preprocessing
###

remove.repeats.over.maxrep = function(vec, maxrep){
	if( is.unsorted(vec) )
		vec = sort.int(vec);
	if( maxrep > 0 ) {
		kill = which(diff(vec, maxrep) == 0L);
		if(length(kill)>0) {
			vec[kill] = 0L;
			vec = vec[vec!=0L];
		}
	}
	return(vec);
}
if(FALSE) { # test code
	remove.repeats.over.maxrep(rep(1:10,1:10), 5L)
}
bam.removeRepeats = function(rbam, maxrep){
	if(maxrep>0) {
		newbam = list(
			startsfwd = lapply( rbam$startsfwd, remove.repeats.over.maxrep, maxrep),
			startsrev = lapply( rbam$startsrev, remove.repeats.over.maxrep, maxrep),
			qc = rbam$qc);
	} else {
		newbam = rbam;
	}
	newbam$qc$frwrev.no.repeats = c(
		sum(sapply(newbam$startsfwd,length)),
		sum(sapply(newbam$startsrev,length)));
	class(newbam$qc$frwrev.no.repeats) = "qcFrwrev";
	newbam$qc$reads.recorded.no.repeats = sum(newbam$qc$frwrev.no.repeats);
	
	return(newbam);
}

### Non-CpG set of locations
noncpgSitesFromCpGset = function(cpgset, distance){
	noncpg = vector("list", length(cpgset));
	names(noncpg) = names(cpgset);
	for( i in seq_along(cpgset) ) { # i=1;
		pos = cpgset[[i]];
		difpos = diff(pos);
		keep = which(difpos>=(distance*2L));
		newpos = (pos[keep+1L] + pos[keep]) %/% 2L;
		noncpg[[i]] = newpos;
	}
	return(noncpg);
}
### Find isolated CpGs among the given set of CpGs
isocpgSitesFromCpGset = function(cpgset, distance){
	isocpg = vector("list",length(cpgset));
	names(isocpg) = names(cpgset);
	for( i in seq_along(cpgset) ) {	
		distbig = diff(cpgset[[i]]) >= distance;
		isocpg[[i]] = cpgset[[i]][ which( c(distbig[1],distbig[-1] & distbig[-length(distbig)], distbig[length(distbig)]) ) ];
	}
	return(isocpg);
}
if(FALSE) { # test code
	cpgset = readRDS("C:/AllWorkFiles/Andrey/VCU/RaMWAS_2/code/Prepare_CpG_list/hg19/spgset_hg19_SNPS_at_MAF_0.05.rds")
	noncpg = noncpgSitesFromCpGset(cpgset, 200);
	sapply(cpgset, typeof)
	sapply(noncpg, typeof)
	sapply(cpgset, length)
	sapply(noncpg, length)
	
	cpgset = lapply(1:10, function(x){return(c(1,1+x,1+2*x))})
	names(cpgset) = paste0("chr",seq_along(cpgset))
	show(cpgset);
	noncpg = noncpgSitesFromCpGset(cpgset, 3);
	show(noncpg);
	isocpg = isocpgSitesFromCpGset(cpgset, 3);
	show(isocpg);
}

### Count reads away from CpGs
.count.nonCpG.reads.forward = function( starts, cpglocations, distance){
	### count CpGs before the read
	### count CpGs before and covered by the read
	ind = findInterval(c(starts-1L,starts+(distance-1L)), cpglocations);
	dim(ind)=c(length(starts),2);
	# cbind(ind, starts)
	return(c(sum(ind[,1] == ind[,2]),length(starts)));
}
.count.nonCpG.reads.reverse = function( starts, cpglocations, distance){
	### count CpGs left of read (+distance)
	### count CpGs left of read start or at start
	ind = findInterval(c(starts-distance,starts), cpglocations);
	dim(ind)=c(length(starts),2);
	# cbind(ind, starts)
	return(c(sum(ind[,1] == ind[,2]),length(starts)));
}
bam.count.nonCpG.reads = function(rbam, cpgset, distance){
	result = c(nonCpGreads = 0,totalreads = 0);
	for( chr in names(cpgset) ) { # chr = names(cpgset)[1]
		frwstarts = rbam$startsfwd[[chr]];
		if( length(frwstarts)>0 )
			result = result + .count.nonCpG.reads.forward( starts = frwstarts, cpglocations = cpgset[[chr]], distance);
		revstarts = rbam$startsrev[[chr]];
		if( length(revstarts)>0 )
			result = result + .count.nonCpG.reads.reverse( starts = revstarts, cpglocations = cpgset[[chr]], distance);
	}
	rbam$qc$cnt.nonCpG.reads = result;
	class(rbam$qc$cnt.nonCpG.reads) = "qcNonCpGreads"
	return(rbam);
}
if(FALSE) { # test code
	rbam = list( startsfwd = list( chr1 = 1:100, chr2 = 1:100 ), startsrev = list(chr1 = 100:200) )
	data(toycpgset);
	cpgset = toycpgset
	show(toycpgset)
	distance = 10;
	
	starts = rbam$startsfwd$chr1;
	cpglocations = cpgset$chr1
	
	starts = rbam$startsrev$chr1;
	cpglocations = cpgset$chr1
	
	starts = sample(1e9, 1e2);
	cpglocations = sort.int(sample(1e9, 1e7));
	
	system.time( .count.nonCpG.reads.forward( starts, cpglocations, distance=20) )
	system.time( .count.nonCpG.reads.forward2(starts, cpglocations, distance=20) )
	
	rbam2 = bam.count.nonCpG.reads(rbam, toycpgset, 50)
	cat( rbam2$qc$bam.count.nonCpG.reads[1], "of",rbam2$qc$bam.count.nonCpG.reads[2],"reads are not covering CpGs","\n" );
}

### Get distribution of distances to isolated CpGs
.hist.isodist.forward = function( starts, cpglocations, distance){
	### count CpGs before the read
	### count CpGs before and covered by the read
	ind = findInterval(c(starts-1L,starts+(distance-1L)), cpglocations);
	dim(ind)=c(length(starts),2);
	# cbind(ind[,1] != ind[,2], ind, starts)
	set = which(ind[,1] != ind[,2]);
	dists = cpglocations[ind[set,2]] - starts[set];
	counts = tabulate(dists+1L, distance);
	return(counts);
}
.hist.isodist.reverse = function( starts, cpglocations, distance){
	### count CpGs left of read (+distance)
	### count CpGs left of read start or at start
	ind = findInterval(c(starts-distance,starts), cpglocations);
	dim(ind)=c(length(starts),2);
	# cbind(ind, starts)
	set = which(ind[,1] != ind[,2]);
	dists = starts[set] - cpglocations[ind[set,2]];
	counts = tabulate(dists+1L, distance);
	return(counts);
}
bam.hist.isolated.distances = function(rbam, isocpgset, distance){
	result = 0;
	for( chr in names(isocpgset) ) { # chr = names(cpgset)[1]
		frwstarts = rbam$startsfwd[[chr]];
		if( length(frwstarts)>0 )
			result = result + .hist.isodist.forward( starts = frwstarts, cpglocations = isocpgset[[chr]], distance);
		revstarts = rbam$startsrev[[chr]];
		if( length(revstarts)>0 )
			result = result + .hist.isodist.reverse( starts = revstarts, cpglocations = isocpgset[[chr]], distance);
	}
	rbam$qc$hist.isolated.dist1 = result;
	class(rbam$qc$hist.isolated.dist1) = "qcIsoDist";
	return(rbam);
}
if(FALSE) { # test code
	rbam = list( startsfwd = list(chr1=100), startsrev = list(chr1 = 103) );
	isocpgset = list(chr1 = 101);
	distance = 100;
	
	rbam2 = bam.hist.isolated.distances(rbam, isocpgset, distance);
	which(rbam2$qc$hist.isolated.dist1>0)
}

### Get average coverage vs. CpG density
bam.coverage.by.density = function( rbam, cpgset, noncpgset, minfragmentsize, maxfragmentsize) {
	
	fragdistr = c(rep(1, minfragmentsize-1),seq(1,0,length.out = (maxfragmentsize-minfragmentsize)/1.5+1));
	fragdistr = fragdistr[fragdistr>0];

	if( is.null(noncpgset) ) {
		noncpgset = noncpgSitesFromCpGset(cpgset = cpgset, distance = maxfragmentsize)
	}
	# sum(sapply(noncpgset,length))
	# newcpgset = noncpgset;
	# for( chr in seq_along(noncpgset) ) {
	# 	newcpgset[[chr]] = sort.int( c(cpgset[[chr]], noncpgset[[chr]]) );
	# }
	# rm(noncpgset);
	
	cpgdensity1 = calc.coverage(rbam = list(startsfwd = cpgset), cpgset = cpgset, fragdistr = fragdistr);
	cpgdensity2 = calc.coverage(rbam = list(startsrev = lapply(cpgset,`-`,1L)), cpgset = cpgset, fragdistr = fragdistr[-1]);
	cpgdensity = unlist(cpgdensity1, recursive = FALSE, use.names = FALSE) + 
					 unlist(cpgdensity2, recursive = FALSE, use.names = FALSE);
	rm(cpgdensity1,cpgdensity2);
	
	### It's zero, not needed
	# noncpgdensity1 = calc.coverage(rbam = list(startsfwd = noncpgset), cpgset = cpgset, fragdistr = fragdistr);
	# noncpgdensity2 = calc.coverage(rbam = list(startsrev = lapply(noncpgset,`-`,1L)), cpgset = cpgset, fragdistr = fragdistr[-1]);
	# noncpgdensity = unlist(noncpgdensity1, recursive = FALSE, use.names = FALSE) + 
	# 					 unlist(noncpgdensity2, recursive = FALSE, use.names = FALSE);
	# rm(noncpgdensity1,noncpgdensity2);
	
	cpgcoverage = calc.coverage(rbam, cpgset,    fragdistr);
	cpgcoverage = unlist(cpgcoverage, recursive = FALSE, use.names = FALSE);
	
	noncoverage = calc.coverage(rbam, noncpgset, fragdistr);
	noncoverage = unlist(noncoverage, recursive = FALSE, use.names = FALSE);

	# sqrtcover = sqrt(coverage);
	sqrtcpgdensity = sqrt(cpgdensity);
	rm(cpgdensity);
	
	axmax = ceiling(quantile(sqrtcpgdensity,0.99)*100)/100;
	# axmaxsafe = ceiling(quantile(sqrtcpgdensity,0.9)*100)/100;

	library(KernSmooth);
	z = locpoly(x = c(sqrtcpgdensity, double(length(noncoverage))),
					y = c(cpgcoverage, noncoverage), 
					bandwidth = 0.5, gridsize = axmax*100+1, range.x = c(0,axmax));
	z$y[is.na(z$y)] = 0;
	# z = locpoly(cpgdensity, coverage, bandwidth = 0.2, gridsize = axmax*100+1, range.x = c(0,axmax))
	# plot(z$x, z$y, type="l", ylim = c(0,max(z$y, na.rm = TRUE)*1.1), yaxs="i", xaxs="i");
	# # sum(sapply(rbam$startsfwd[names(cpgset)], length)) + sum(sapply(rbam$startsrev[names(cpgset)], length))
	# reads.used = sum(sapply(rbam$startsfwd, length)) + sum(sapply(rbam$startsrev, length));
	# additive.vector = c(reads.used, z$y);
	
	# bins = hexbin(sqrtcpgdensity[sqrtcover<5], sqrtcover[sqrtcover<5],xbins = 100, ybnds = c(0,5))
	# plot(bins, style = "colorscale", colramp= function(n){magent(n,beg=200,end=1)}, trans = function(x)x^0.6);
	rbam$qc$avg.coverage.by.density = z$y;
	class(rbam$qc$avg.coverage.by.density) = "qcCoverageByDensity";
	rbam$qc$avg.noncpg.coverage = mean(noncoverage);
	rbam$qc$avg.cpg.coverage = mean(cpgcoverage);
	
	return(rbam);
}
if(FALSE) {
	rbam = readRDS("D:/RW/RC2/rds_rbam/01A01SM1429N.rbam.rds");
	cpgset = cachedRDSload("C:/AllWorkFiles/Andrey/VCU/RaMWAS_2/code/Prepare_CpG_list/hg19/hg19_1kG_MAF_0.01_chr1-22XY_cushaw3_50bp.rds");
	noncpgset = cachedRDSload("C:/AllWorkFiles/Andrey/VCU/RaMWAS_2/code/Prepare_CpG_list/hg19/hg19_1kG_MAF_0.01_chr1-22XY_cushaw3_50bp_nonCpG.rds");
	minfragmentsize = 50;
	maxfragmentsize = 200;
	rbam = bam.coverage.by.density( rbam, cpgset, noncpgset, minfragmentsize, maxfragmentsize);
	plot(rbam$qc$avg.coverage.by.density, "name", col="blue");
	
	
	param = list(
		dirbam = "D:/Cell_type/bams/",
		dirproject = "D:/Cell_type/",
		filebamlist = "D:/Cell_type/000_list_of_files.txt",
		scoretag = "AS",
		minscore = 100,
		cputhreads = 8,
		filecpgset = "C:/AllWorkFiles/Andrey/VCU/RaMWAS_2/code/Prepare_CpG_list/hg19/cpgset_hg19_SNPS_at_MAF_0.05.rds",
		filenoncpgset = NULL,
		maxrepeats = 3,
		maxfragmentsize=200,
		minfragmentsize=50,
		bamnames = NULL
	);
	param = parameterPreprocess(param);

	pipelineSaveQCplots(param, rbam, bamname="bamname")
}

### Estimate fragment size distribution
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
	rezfit = rezfit - max(rezfit[!keep],0)
	rezfit[1:seqLength] = rezfit[seqLength];
	rezfit = rezfit[keep];
	rezfit = rezfit / rezfit[1];
	
	# lz = lm(hist.isolated.distances[seq_along(rezfit)] ~ rezfit)
	# lines(rezfit*lz$coefficients[2]+lz$coefficients[1], lwd = 4, col="red");
	
	return(rezfit);
}
if(FALSE) { # test code
	x = seq(0.01,0.99,0.01);
	y = sqrt(abs(x-0.5))*sign(x-0.5)
	plot(x,y)
	log.ss <- nls(y ~ SSlogis(x, phi1, phi2, phi3))
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

### Cache CpG location files to avoid reloading.
cachedRDSload = function(rdsfilename){
	if(is.null(rdsfilename))
		return(NULL);
	cachename = rdsfilename; #paste0(".ramwas.",rdsfilename);
	if( exists(x = cachename, envir = .ramwasEnv) ) {
		# cat("Using cache","\n");
		return(base::get(x = cachename, envir = .ramwasEnv));
	} else {
		# cat("Loading","\n");
		data = readRDS(rdsfilename);
		base::assign(x = cachename, value = data, envir = .ramwasEnv);
		return(data);
	}
}
if(FALSE) { # test code
	rdsfilename = "C:/AllWorkFiles/Andrey/VCU/RaMWAS_2/code/Prepare_CpG_list/hg19/cpgset_hg19_SNPS_at_MAF_0.05.rds";
	system.time({z = cachedRDSload(rdsfilename)});
	system.time({z = cachedRDSload(rdsfilename)});
	system.time({z = cachedRDSload(rdsfilename)});
}

### Coverage calculation
.calc.coverage.chr = function(startfrw, startrev, cpgs, fragdistr) {
	maxfragmentsize = length(fragdistr);
	# if(is.null(cover)) {
		cover = double(length(cpgs));
	# }
	
	if(length(startfrw) > 0) {
		ind1 = findInterval(cpgs - maxfragmentsize, startfrw);  # CpGs left of start
		ind2 = findInterval(cpgs,                   startfrw);  # CpGs left of start+250L
		# coverage of CpGs 
		# which(ind2>ind1) 
		# are covered by fragments ind1[which(ind2>ind1)]+1 .. ind2[which(ind2>ind1)]
		.Call("cover_frw_c", startfrw, cpgs, fragdistr, ind1, ind2, cover, PACKAGE = "ramwas");
	}
	
	if(length(startrev) > 0) {
		ind1 = findInterval(cpgs - 1L,                 startrev);  # CpGs left of start
		ind2 = findInterval(cpgs + maxfragmentsize-1L, startrev);  # CpGs left of start+250L
		# coverage of CpGs 
		# which(ind2>ind1) 
		# are covered by fragments ind1[which(ind2>ind1)]+1 .. ind2[which(ind2>ind1)]
		.Call("cover_rev_c", startrev, cpgs, fragdistr, ind1, ind2, cover, PACKAGE = "ramwas");
	}
	return( cover );
}
calc.coverage = function(rbam, cpgset, fragdistr) {
	# if( is.null(coveragelist) ) {
		coveragelist = vector("list", length(cpgset));
		names(coveragelist) = names(cpgset);
	# }
	
	for( chr in names(coveragelist) ) { # chr = names(coveragelist)[1]
		coveragelist[[chr]] = 
			.calc.coverage.chr(rbam$startsfwd[[chr]], rbam$startsrev[[chr]], cpgset[[chr]], fragdistr); 
	}
	return(coveragelist);
}
if(FALSE) {
	# testing calc.coverage
	cpgset = list(chr1 = 1:100);
	rbam = list(startsfwd = list(chr1  = c(10L, 20L)), startsrev = list(chr1  = c(80L, 90L)));
	fragdistr = c(4,3,2,1);
	
	cvl = calc.coverage(rbam, cpgset, fragdistr)
	cv = cvl$chr1;
	
	# testing .calc.coverage.chr
	chr = "chr1"
	# startfrw = rbam$startsfwd[[chr]]; startrev = rbam$startsrev[[chr]]; cpgs = cpgset[[chr]];
	cv = .calc.coverage.chr(rbam$startsfwd[[chr]], rbam$startsrev[[chr]], cpgset[[chr]], fragdistr)
	
	cv[rbam$startsfwd[[chr]]]
	cv[rbam$startsrev[[chr]]]
	cv[rbam$startsfwd[[chr]]+1]
	cv[rbam$startsrev[[chr]]+1]
	cv[rbam$startsfwd[[chr]]-1]
	cv[rbam$startsrev[[chr]]-1]
	
	
	# Timing CpG density calculation
	
	rdsfilename = "C:/AllWorkFiles/Andrey/VCU/RaMWAS_2/code/Prepare_CpG_list/hg19/cpgset_hg19_SNPS_at_MAF_0.05.rds";
	
	cpgset = cachedRDSload(rdsfilename);
	
	fragdistr = c(rep(1,75), seq(1,0,length.out = 76))
	fragdistr = fragdistr[fragdistr>0];
	
	system.time({ covlist1 = calc.coverage( rbam = list( startsfwd = cpgset, startsrev = cpgset), cpgset = cpgset, fragdistr = fragdistr) });
	# 4.72
	
	system.time({ covlist2 = calc.coverage.simple( bam = list( startlistfwd = cpgset, startlistrev = cpgset), cpgsloc = cpgset, fragdistr = fragdistr) });
	# 31.40
	
	range(covlist1$chr1 - covlist2$chr1)
	range(covlist1$chr2 - covlist2$chr2)
}

pipelineSaveQCplots = function(param, rbam, bamname) {
	filename = paste0(param$dirqc,"/score/hs_",bamname,".pdf");
	dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
	pdf(filename);
	plot(rbam$qc$hist.score1, samplename = bamname);
	plot(rbam$qc$bf.hist.score1, samplename = bamname);
	dev.off();
	rm(filename);
	
	filename = paste0(param$dirqc,"/edit_distance/ed_",bamname,".pdf");
	dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
	pdf(filename);
	plot(rbam$qc$hist.edit.dist1, samplename = bamname);
	plot(rbam$qc$bf.hist.edit.dist1, samplename = bamname);
	dev.off();
	rm(filename);
	
	filename = paste0(param$dirqc,"/matched_length/ml_",bamname,".pdf");
	dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
	pdf(filename);
	plot(rbam$qc$hist.length.matched, samplename = bamname);
	plot(rbam$qc$bf.hist.length.matched, samplename = bamname);
	dev.off();
	rm(filename);
	
	if( !is.null(rbam$qc$hist.isolated.dist1) ) {
		filename = paste0(param$dirqc,"/isolated_distance/id_",bamname,".pdf");
		dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
		pdf(filename);
		plot(rbam$qc$hist.isolated.dist1, samplename = bamname);
		dev.off();
		rm(filename);
	}
	if( !is.null(rbam$qc$avg.coverage.by.density) ) {
		filename = paste0(param$dirqc,"/coverage_by_density/cbd_",bamname,".pdf");
		dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
		pdf(filename);
		plot(rbam$qc$avg.coverage.by.density, samplename = bamname);
		dev.off();
		rm(filename);
	}	
}
if(FALSE) {
	param = list(
		dirbam = "D:/Cell_type/bams/",
		dirproject = "D:/Cell_type/",
		filebamlist = "D:/Cell_type/000_list_of_files.txt",
		scoretag = "AS",
		minscore = 100,
		cputhreads = 8,
		filecpgset = "C:/AllWorkFiles/Andrey/VCU/RaMWAS_2/code/Prepare_CpG_list/hg19/cpgset_hg19_SNPS_at_MAF_0.05.rds",
		filenoncpgset = NULL,
		maxrepeats = 3,
		maxfragmentsize=200,
		minfragmentsize=50,
		bamnames = NULL
	);
	param = parameterPreprocess(param);
	rbam = readRDS("D:/Cell_type/rds_rbam/150114_WBCS014_CD20_150.rbam.rds");
	pipelineSaveQCplots(param, rbam, bamname="150114_WBCS014_CD20_150");
}


bam.chrX.qc = function(rbam) {
	strandfun = function(st){c(length(st$chrX), sum(sapply(st,length)))};
	rbam$qc$chrX.count =  strandfun(rbam$startsfwd) + strandfun(rbam$startsfwd);
	class(rbam$qc$chrX.count) = "qcChrX"
	return(rbam);
}

### Pipeline parts
pipelineProcessBam = function(bamname, param) {
	# Used parameters: scoretag, minscore, filecpgset, maxrepeats
	
	param = parameterPreprocess(param);
	
	if( !is.null(param$filecpgset) && is.null(param$maxfragmentsize) )
		return("Parameter not set: maxfragmentsize");
	
	bamname = gsub("\\.bam$","",bamname);
	bamfullname = makefullpath(param$dirbam, paste0(bamname,".bam"))

	dir.create(param$dirrbam, showWarnings = FALSE, recursive = TRUE)
	dir.create(param$dirrqc, showWarnings = FALSE, recursive = TRUE)
	
	rdsbmfile = paste0( param$dirrbam, "/", basename(bamname), ".rbam.rds" );
	rdsqcfile = paste0( param$dirrqc, "/", basename(bamname), ".qc.rds" );
	
	savebam = TRUE;
	if( file.exists( rdsbmfile ) ) {
		if( param$recalculate.QCs ) {
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
	
	if( !is.null(param$filecpgset) ) {
		cpgset = cachedRDSload(param$filecpgset);
		noncpgset = cachedRDSload(param$filenoncpgset);
		isocpgset = isocpgSitesFromCpGset(cpgset = cpgset, distance = param$maxfragmentsize);
		rbam3 = bam.hist.isolated.distances(rbam = rbam2, isocpgset = isocpgset, distance = param$maxfragmentsize);
		rbam4 = bam.coverage.by.density(rbam = rbam3, cpgset = cpgset, noncpgset = noncpgset,
			minfragmentsize = param$minfragmentsize, maxfragmentsize = param$maxfragmentsize);
		
		# if( !is.null(param$noncpgfile)) {
		# 	noncpgset = cachedRDSload(param$noncpgfile);
		# } else {
		# 	noncpgset = noncpgSitesFromCpGset(cpgset = cpgset, distance = param$maxfragmentsize);
		# }
		rbam5 = bam.count.nonCpG.reads(rbam = rbam4, cpgset = cpgset, distance = param$maxfragmentsize);
		
		### QC plots
		pipelineSaveQCplots(param, rbam4, bamname);
		
	} else {
		rbam5 = rbam2;
	}
	if(savebam)
		saveRDS( object = rbam5, file = rdsbmfile, compress = "xz");
	rbam6 = rbam5;
	rbam6$startsfwd=NULL;
	rbam6$startsrev=NULL;
	rbam6$qc$nbams = 1L;
	saveRDS( object = rbam6, file = rdsqcfile, compress = "xz");
	
	return(paste0("OK. ", bamname));
}

### Fragment size distribution estimation
pipelineEstimateFragmentSizeDistribution = function(param) {
	
	param = parameterPreprocess(param);
	
	if( !is.null(param$bam2sample) ) {
		bams = unlist( param$bam2sample, use.names = FALSE);
	} else if (!is.null(param$bamnames)) {
		bams = param$bamnames;
	} else {
		stop("Bams are not defined. Set filebam2sample, filebamlist, bam2sample or bamnames.","\n");
	}
	bams = unique(basename(bams));

	qclist = vector("list", length(bams));
	names(qclist) = bams;
	
	for( bamname in bams) {
		rdsqcfile = paste0( param$dirrqc, "/", bamname, ".qc.rds" );
		qclist[[bamname]] = readRDS(rdsqcfile);
	}
	
	qcset = lapply(lapply( qclist, `[[`, "qc"),`[[`,"hist.isolated.dist1")
	bighist = Reduce(`%add%`, qcset);
	estimate = estimateFragmentSizeDistribution(bighist, param$minfragmentsize);
	
	writeLines(con = paste0(param$dirfilter,"/Fragment_size_distribution.txt"), text = as.character(estimate));
	
	return(estimate);
}

### Get coverage for 1 or more BAMs
pipelineCoverage1Sample = function(colnum, param){
	# library(ramwas);
	# library(filematrix);
	
	cpgset = cachedRDSload(param$filecpgset);

	bams = param$bam2sample[[colnum]];
	
	# bams = c(param$bam2sample[1], param$bam2sample[2]);

	if( param$maxrepeats == 0 ) {
		coverage = NULL;
		for( j in seq_along(bams)) { # j=1
			rbam = readRDS( paste0( param$dirrbam, "/", bams[j], ".rbam.rds" ) );
			cov = calc.coverage(rbam = rbam, cpgset = cpgset, fragdistr = param$fragdistr)
			if(is.null(coverage)) {
				coverage = cov;
			} else {
				for( i in seq_along(coverage) )
					coverage[[i]] = coverage[[i]] + cov[[i]]
			}
			rm(cov);
		}
	} else {
		rbams = vector("list",length(bams));
		for( j in seq_along(bams)) { # j=1
			rbams[[j]] = readRDS( paste0( param$dirrbam, "/", bams[j], ".rbam.rds" ) );
		}
		if(length(bams) > 1) {
			rbam = list(startsfwd = list(), startsrev = list());
			for( i in seq_along(cpgset) ) { # i=1
				nm = names(cpgset)[i];
				
				fwd = lapply(rbams, function(x,y){x$startsfwd[[y]]}, nm);
				fwd = sort.int( unlist(fwd, use.names = FALSE) );
				rbam$startsfwd[[nm]] = remove.repeats.over.maxrep(fwd, param$maxrepeats);
				rm(fwd);
				
				rev = lapply(rbams, function(x,y){x$startsrev[[y]]}, nm);
				rev = sort.int( unlist(rev, use.names = FALSE) );
				rbam$startsrev[[nm]] = remove.repeats.over.maxrep(rev, param$maxrepeats);
				rm(rev);
			}
		} else {
			rbam = bam.removeRepeats( rbams[[1]], param$maxrepeats );
		}
		rm(rbams);
		coverage = calc.coverage(rbam = rbam, cpgset = cpgset, fragdistr = param$fragdistr)
	}
	return(coverage);
}

### RaMWAS pipeline
.ramwas1scanBam = function(bamname, param){
	cat(file = paste0(param$dirlog,'/011_log_scanning_bams.txt'), bamname, '\n', append = TRUE);
	pipelineProcessBam(bamname = bamname, param = param);
}
ramwas1scanBams = function( param ){
	param = parameterPreprocess(param);
	stopifnot( !is.null(param$bamnames));
	
	if( param$cputhreads > 1) {
		dir.create(param$dirlog, showWarnings = FALSE, recursive = TRUE);
		cat(file = paste0(param$dirlog,'/011_log_scanning_bams.txt'), 'Log, scanning bams:', date(), '\n', append = FALSE);
		cl <- makeCluster(param$cputhreads);
		z = clusterApplyLB(cl, param$bamnames, .ramwas1scanBam, param = param);
		stopCluster(cl);
	} else {
		z = character(length(param$bamnames));
		names(z) = param$bamnames;
		for(i in seq_along(param$bamnames)) { # i=1
			z[i] = .ramwas1scanBam(bamname = param$bamnames[i], param = param);
		}
	}
	return(z);
}

ramwas2collectqc = function( param ){
	param = parameterPreprocess(param);
	
	bams = NULL;
	if( !is.null(param$bamnames) )
		bams = c(bams, param$bamnames);
	if( !is.null(param$bam2sample) )
		bams = c(bams, unlist(param$bamnames, use.names = FALSE));
	
	bams = unique(basename(bams));
	
	message('Load BAM QC info','\n');
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
		
		cat('Saving text summary','\n');
		bigqc = vector("list", length(bamset));
		names(bigqc) = names(bamset);
		text = character(length(bamset));
		for( ibam in seq_along(bamset) ) { # ibam=7
			curbams = rbamlist[bamset[[ibam]]];
			qc = .combine.bams.qc(curbams)$qc;
			if( length(qc) > 0 ) {
				bigqc[[ibam]] = qc;
			} else {
				qc = NULL;
			}
			text[ibam] = .qcTextLine( qc, names(bamset)[ibam] );
		}
		writeLines( con = paste0(dirloc, "/Summary_QC.txt"), text = c(.qcTextHeader, text));
		
		figfun = function(qcname, plotname) {
			cat('Saving plots',plotname,'\n');
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
	dirname = "summary_by_bam";
	cat('Saving QC info by BAM','\n');
	collect.qc.summary(bamset, dirname);
	rm(bamset, dirname);
	
	if( !is.null(param$bam2sample) ) {
		# by sample
		cat('Saving QC info by SAMPLE','\n');
		collect.qc.summary(bamset = param$bam2sample, dirname = "summary_by_sample");
		cat('Saving QC info TOTAL (all BAMs in bam2sample)','\n');
		bigqc = collect.qc.summary(bamset = list(total=unique(unlist(param$bam2sample, use.names = FALSE))), dirname = "summary_total");
	} else {
		cat('Saving QC info TOTAL (all BAMs in bamnames/filebamlist)','\n');
		bigqc = collect.qc.summary(bamset = list(total=bams), dirname = "summary_total");
	}
	
	### Fragment size
	frdata = bigqc$total$hist.isolated.dist1;
	estimate = estimateFragmentSizeDistribution(frdata, param$minfragmentsize);
	writeLines(con = paste0(param$dirfilter,"/Fragment_size_distribution.txt"), text = as.character(estimate));
	
	pdf(paste0(param$dirqc,'/Fragment_size_distribution_estimate.pdf'),8,8);
	lz = lm(frdata[seq_along(estimate)] ~ estimate)
	plot(as.vector(frdata)/1000, pch = 19, col='blue', main='Isolated CpG coverage vs.\nfragment size distribution estimate', ylab='count, thousands', xlab='Distance to isolated CpGs', xaxs="i");
	lines((estimate*lz$coefficients[2]+lz$coefficients[1])/1000, lwd = 4, col='red');
	dev.off();
	return(invisible(NULL));
}

.ramwas3coverageJob = function(colnum, param, nslices){
	# colnum = 3
	library(ramwas);
	library(filematrix)
	coverage = pipelineCoverage1Sample(colnum, param);
	coverage = unlist(coverage, use.names = FALSE);

	start = 1;
	for( part in seq_len(nslices) ) {
		cat('colnum =',colnum,'part =', part,'\n');
		fmname = paste0(param$dirtemp,"/RawCoverage_part",part);
		fm = fm.open(fmname, lockfile = param$lockfile);
		ntowrite = nrow(fm);
		fm$writeCols(colnum, coverage[start:(start+ntowrite-1)]);
		close(fm);
		start = start + ntowrite;
	}
	
	fmname = paste0(param$dirtemp,"/RawCoverage_part",1);
	fm = fm.open(fmname, lockfile = param$lockfile);
	fm$filelock$lockedrun( {
		cat(file = paste0(param$dirlog,'/000_log_raw_coverage_calc.txt'), colnum, '\n', append = TRUE);
	});
	close(fm);
	
	return("OK");
}
.ramwas3transposeFilterJob = function(fmpart, param) {
	library(filematrix);
	# fmpart = 1
	filename = paste0(param$dirtemp,"/RawCoverage_part",fmpart);
	if( !file.exists(paste0(filename,'.bmat')) || !file.exists(paste0(filename,'.desc.txt')) )
		return(paste0('Raw coverage slice filematrix not found: ', filename));
	fmraw = fm.open(filename, lockfile = param$lockfile2);
	mat = fmraw[];
	# mat = as.matrix(fmraw);

	fmout = fm.create( paste0(param$dirtemp,"/TrCoverage_part",fmpart), 
							 nrow = ncol(mat), ncol = 0, size = param$doublesize, lockfile = param$lockfile2);
	fmpos = fm.create( paste0(param$dirtemp,"/TrCoverage_loc",fmpart), 
							 nrow = 1, ncol = 0, type = 'integer', lockfile = param$lockfile2);
	
	samplesums = rep(0, ncol(mat));
	
	### Sliced loop
	step1 = max(floor(32*1024*1024 / 8 / ncol(mat)),1);
	mm = nrow(mat);
	nsteps = ceiling(mm/step1);
	for( part in 1:nsteps ) { # part = 1
		cat(fmpart, part, 'of', nsteps, '\n');
		fr = (part-1)*step1 + 1;
		to = min(part*step1, mm);
		
		subslice = mat[fr:to,];
		
		### Filtering criteria
		cpgmean = rowMeans( subslice );
		cpgnonz = rowMeans( subslice>0 );
		keep = (cpgmean >= param$minavgcpgcoverage) & (cpgnonz >= param$minnonzerosamples);
		if( !any(keep) )
			next;
		
		slloc = fr:to;
		
		if( !all(keep) ) {
			keep = which(keep);
			subslice = subslice[keep,,drop=FALSE];
			slloc = slloc[keep];
		}
		
		subslice = t(subslice);

		samplesums = samplesums + rowSums(subslice);
		
		fmout$appendColumns(subslice);
		fmpos$appendColumns(slloc);
		rm(subslice, slloc, keep, cpgmean, cpgnonz)
	}
	rm(part, step1, mm, nsteps, fr, to, mat);
	gc();
	
	close(fmout);
	close(fmpos);
	
	fmss = fm.open( paste0(param$dirtemp,"/0_sample_sums"), lockfile = param$lockfile2);
	fmss[,fmpart] = samplesums;
	close(fmss);
	
	fmraw$filelock$lockedrun( {
		cat(file = paste0(param$dirlog,'/000_log_transpose_and_filter.txt'), fmpart, '\n', append = TRUE);
	});
	closeAndDeleteFiles(fmraw);
	return("OK.");
}
.ramwas3normalizeJob = function(fmpart_offset, param, samplesums){
	# fmpart_offset = fmpart_offset_list[[2]]
	scale = as.vector(samplesums) / mean(samplesums);
	
	library(filematrix);
	
	filename = paste0(param$dirtemp, "/TrCoverage_part", fmpart_offset[1]);
	mat = fm.load(filename, param$lockfile1);
	mat = mat / scale;
	
	filename = paste0(param$dircoveragenorm, "/Coverage");
	fm = fm.open(filename, lockfile = param$lockfile2);
	fm$writeCols( start = fmpart_offset[2]+1L, mat);
	close(fm);
	
	fm$filelock$lockedrun( {
		cat(file = paste0(param$dirlog,'/000_log_normalize.txt'), fmpart_offset[1], '\n', append = TRUE);
	});

	
	rm(mat);
	gc();
	return("OK.");
}
ramwas3NormalizedCoverage = function( param ) {
	# Prepare
	param = parameterPreprocess(param);
	param$fragdistr = as.double( readLines(con = paste0(param$dirfilter,"/Fragment_size_distribution.txt")));
	dir.create(param$dirtemp, showWarnings = FALSE, recursive = TRUE);
	dir.create(param$dirlog, showWarnings = FALSE, recursive = TRUE);
	
	### data dimensions
	cpgset = cachedRDSload(param$filecpgset);
	ncpgs = sum(sapply(cpgset, length));
	nsamples = length(param$bam2sample);
	
	### Check is all rbams are in place
	{
		message('Checking if all required Rbam files present');
		bams = unlist(param$bam2sample);
		for( bname in bams) {
			filename = paste0( param$dirrbam, '/', bname);
			if( file.exists(filename) ) {
				stop(paste0('Rbam file from bam2sample does not exist: ', filename));
			}
		}
		rm(bams, bname, filename);
	}

	### Create raw coverage matrix slices
	{
		library(filematrix)
		# Sliced loop
		kbblock = (128*1024)/8;
		step1 = max(floor(param$buffersize / (8 * nsamples)/kbblock),1)*kbblock;
		mm = ncpgs;
		nslices = ceiling(mm/step1);
		message('Creating ',nslices,' file matrices for raw coverage at: ', param$dirtemp);
		for( part in 1:nslices ) { # part = 1
			# cat('Creating raw coverage matrix slices', part, 'of', nslices, '\n');
			fr = (part-1)*step1 + 1;
			to = min(part*step1, mm);
			fmname = paste0(param$dirtemp,"/RawCoverage_part",part);
			fm = fm.create(fmname, nrow = to-fr+1, ncol = nsamples, size = param$doublesize)
			close(fm);
		}
		rm(part, step1, mm, fr, to, fmname);
	} # nslices

	### Fill in the raw coverage files
	{
		message('Calculating and saving raw coverage');
		param$lockfile = tempfile();
		cat(file = paste0(param$dirlog,'/000_log_raw_coverage_calc.txt'), 'Log, raw coverage:', date(), '\n', append = FALSE);
		library(parallel)
		if( param$cputhreads > 1) {
			cl = makeCluster(param$cputhreads);
			# cl = makePSOCKcluster(rep("localhost", param$cputhreads))
			z = clusterApplyLB(cl, seq_len(nsamples), .ramwas3coverageJob, param = param, nslices = nslices);
			stopCluster(cl);
		} else {
			z = character(nsamples);
			names(z) = names(param$bam2sample);
			for(i in seq_along(param$bam2sample)) { # i=1
				z[i] = .ramwas3coverageJob(colnum = i, param = param, nslices = nslices);
				cat(i,z[i],"\n");
			}
		}
		file.remove(param$lockfile);
	}
	
	### Transpose the slices, filter by average and fraction of non-zeroes
	{
		message('Transposing coverage matrices and filtering CpGs by coverage');
		
		fm = fm.create( paste0(param$dirtemp,"/0_sample_sums"), nrow = nsamples, ncol = nslices);
		close(fm);
		
		cat(file = paste0(param$dirlog,'/000_log_transpose_and_filter.txt'), 'Log, transpose and filter:', date(), '\n', append = FALSE);
		if( param$cputhreads > 1 ) {
			param$lockfile2 = tempfile();
			library(parallel);
			cl = makePSOCKcluster(rep("localhost", param$cputhreads))
			z = clusterApplyLB(cl, 1:nslices, .ramwas3transposeFilterJob, param = param);
			stopCluster(cl);
			file.remove(param$lockfile2);
		} else {
			for( fmpart in seq_len(nslices) ) { # fmpart = 5
				.ramwas3transposeFilterJob( fmpart, param);
			}
		}
	}
	
	### Prepare CpG set for filtered CpGs	
	{
		message('Saving locations for CpGs which passed the filter');
		
		cpgsloc1e9 = cpgset;
		for( i in seq_along(cpgsloc1e9) ) {
			cpgsloc1e9[[i]] = cpgset[[i]] + i*1e9;
		}
		cpgsloc1e9 = unlist(cpgsloc1e9, recursive = FALSE, use.names = FALSE);
		
		kbblock = (128*1024)/8;
		step1 = max(floor(param$buffersize / (8 * nsamples)/kbblock),1)*kbblock;
		mm = ncpgs;
		nsteps = ceiling(mm/step1);
		cpgsloclist = vector('list',nsteps);
		for( part in 1:nsteps ) { # part = 1
			# cat( part, 'of', nsteps, '\n');
			fr = (part-1)*step1 + 1;
			to = min(part*step1, mm);
			
			indx = as.vector(fm.load( paste0(param$dirtemp,"/TrCoverage_loc",part) ));
			cpgsloclist[[part]] = cpgsloc1e9[fr:to][indx];
		}
		rm(part, step1, mm, nsteps, fr, to, kbblock, indx);
		sliceoffsets = c(0L, cumsum(sapply(cpgsloclist, length)));

		cpgslocvec = unlist(cpgsloclist, use.names = FALSE);
		cpgslocmat = cbind( chr = as.integer(cpgslocvec %/% 1e9), position = as.integer(cpgslocvec %% 1e9));
		
		fm = fm.create.from.matrix( filenamebase = paste0(param$dircoveragenorm, '/CpG_locations'), mat = cpgslocmat);
		close(fm);
		writeLines(con = paste0(param$dircoveragenorm, '/CpG_chromosome_names.txt'), text = names(cpgset));
		rm(cpgsloc1e9, cpgsloclist, cpgslocvec, cpgslocmat);
	} # /CpG_locations, sliceoffsets
	
	### Sample sums
	{
		message('Gathering sample sums from ', nslices,' slices');
		
		mat = fm.load( paste0(param$dirtemp,"/0_sample_sums") );
		samplesums = rowSums(mat);
		rm(mat);
		fm = fm.create.from.matrix( paste0(param$dircoveragenorm, '/raw_sample_sums'), samplesums);
		close(fm);
	}
	
	### Normalize and combine in one matrix
	{
		message('Normalizing coverage and saving in one matrix');
		
		fmpart_offset_list = as.list(data.frame(rbind( seq_len(nslices), sliceoffsets[-length(sliceoffsets)])));
		cat(file = paste0(param$dirlog,'/000_log_normalize.txt'), 'Log, normalize:', date(), '\n');
		
		### Create big matrix for normalized coverage
		fm = fm.create(paste0(param$dircoveragenorm, '/Coverage'), 
							nrow = nsamples, ncol = tail(sliceoffsets,1), size = param$doublesize);
		rownames(fm) = names(param$bam2sample);
		close(fm);
		
		### normalize and fill in
		if( param$cputhreads > 1 ) {
			
			param$lockfile1 = tempfile();
			param$lockfile2 = tempfile();
			library(parallel);
			cl = makePSOCKcluster(rep("localhost", param$cputhreads))
			z = clusterApplyLB(cl, fmpart_offset_list, .ramwas3normalizeJob, param = param, samplesums = samplesums);
			stopCluster(cl);
			file.remove(param$lockfile1);
			file.remove(param$lockfile2);
			
		} else {
			for( fmpart in seq_len(nslices) ) { # fmpart = 5
				.ramwas3normalizeJob( fmpart_offset_list[[fmpart]], param, samplesums);
			}
		}
		
		
		
	}
	
	### Cleanup
	{
		message('Removing temporary files');
		for( part in 1:nslices ) {
			fm = fm.open( paste0(param$dirtemp,"/TrCoverage_loc",part) );
			closeAndDeleteFiles(fm);
			fm = fm.open( paste0(param$dirtemp,"/TrCoverage_part",part) );
			closeAndDeleteFiles(fm);
		}
		fm = fm.open( paste0(param$dirtemp,"/0_sample_sums") );
		closeAndDeleteFiles(fm);
		
	}
	
}


if(FALSE){
	param = parameterPreprocess(param);
	cl = makePSOCKcluster(rep("localhost", param$cputhreads))
	clusterEvalQ(cl, library(filematrix));
	z = clusterApplyLB(cl, 1, ramwas:::.ramwas3transposeFilterJob, param = param);
	stopCluster(cl);
}

if(FALSE) { # Cell Type
	library(ramwas);
	param = list(
		dirbam = "D:/RW/CELL/bams/",
		dirproject = "D:/RW/CELL/",
		filebamlist = "000_list_of_files.txt",
		dirtemp = NULL,
		scoretag = "AS",
		minscore = 100,
		cputhreads = 8,
		filecpgset =    "C:/AllWorkFiles/Andrey/VCU/RaMWAS3/cpgset/hg19_1kG_MAF_0.01_chr1-22XY_bowtie2_75bp.rds",
		filenoncpgset = "C:/AllWorkFiles/Andrey/VCU/RaMWAS3/cpgset/hg19_1kG_MAF_0.01_chr1-22XY_bowtie2_75bp_nonCpG.rds",
		maxrepeats = 3,
		maxfragmentsize=200,
		minfragmentsize=50,
		filebam2sample = "bam2sample2.txt",
		recalculate.QCs = TRUE,
		buffersize = 1e9
	);
	if(Sys.getenv('computername') == "HP") {
		param$dirtemp = "E:/RaMWAS_temp_CELL/"
	} else {
		param$dirtemp = "C:/AllWorkFiles/temp/"
	}
	{
		tic = proc.time();
		rbam = pipelineProcessBam(bamname = '150114_WBCS014_CD20_150', param );
		toc = proc.time();
		show(toc-tic);
		# HP: 61.37   18.86  184.10 
	}
	{
		tic = proc.time();
		z = ramwas1scanBams( param );
		toc = proc.time();
		show(toc-tic);
		# HP: 61.37   18.86  184.10 
	}
	{
		tic = proc.time();
		z = ramwas2collectqc( param );
		toc = proc.time();
		show(toc-tic);
		# HP: 61.37   18.86  184.10 
	}
	{
		tic = proc.time();
		z = ramwas2collectqc( param );
		toc = proc.time();
		show(toc-tic);
		# DL: 33.42   10.06   80.37
	}
	{
		param$cputhreads = 8;
		tic = proc.time();
		ramwas3NormalizedCoverage(param);
		toc = proc.time();
		show(toc-tic);
	}
} # Cell Type
if(FALSE) { # NESDA
	library(ramwas);
	param = list(
		dirbam = "D:/RW/NESDA/bams/",
		dirproject = "D:/RW/NESDA/",
		dirtemp = NULL,
		filebamlist = "000_list_of_files.txt",
		scoretag = "AS",
		minscore = 120,
		cputhreads = 8,
		filecpgset =    "C:/AllWorkFiles/Andrey/VCU/RaMWAS3/cpgset/hg19_1kG_MAF_0.01_chr1-22XY_bowtie2_75bp.rds",
		filenoncpgset = "C:/AllWorkFiles/Andrey/VCU/RaMWAS3/cpgset/hg19_1kG_MAF_0.01_chr1-22XY_bowtie2_75bp_nonCpG.rds",
		maxrepeats = 3,
		maxfragmentsize = 200,
		minfragmentsize = 50,
		filebam2sample = "bam2samples_QC_04-15-2016.csv",
		recalculate.QCs = TRUE,
		buffersize = 1e9,
		chrnorm = 1:22
	);
	if(Sys.getenv('computername') == "HP") {
		param$dirtemp = "E:/RaMWAS_temp_NESDA/"
	} else {
		param$dirtemp = "C:/AllWorkFiles/temp/"
	}
	{
		tic = proc.time();
		rez = pipelineProcessBam(bamname = '150210_NS500538_0012_AH2T57BGXX_NESDA_04D02311', param );
		qc = readRDS('D:/RW/NESDA/rds_qc/150210_NS500538_0012_AH2T57BGXX_NESDA_04D02311.qc.rds')
		toc = proc.time();
		show(toc-tic);
		# HP: 61.37   18.86  184.10 
	}
	{
		tic = proc.time();
		z = ramwas1scanBams( param );
		toc = proc.time();
		show(toc-tic);
		# DL, rescan: 7867.78
	} 
	{
		tic = proc.time();
		z = ramwas2collectqc( param );
		toc = proc.time();
		show(toc-tic);
		# DL: 33.42   10.06   80.37
	}
	{
		tic = proc.time();
		z = ramwas3NormalizedCoverage( param );
		toc = proc.time();
		show(toc-tic);
		# DL, 2e9, 8thr:   
		# DL, 1e9, 8thr:  
		# HP, 2e9, 8thr:  
		# HP, 1e9, 8thr:  
	}
	{
		tic = proc.time();
		z = ramwas3NormalizedCoverage1( param );
		toc = proc.time();
		show(toc-tic);
		# DL, 2e9, 8thr:   13.92    4.79 6110.45
		# DL, 1e9, 8thr:   15.54    6.45 7125.52 
		# HP, 2e9, 8thr:   21.79   14.44 4946.97
		# HP, 1e9, 8thr:   23.14   16.61 5599.71
	}
	{
		tic = proc.time();
		z = ramwas3NormalizedCoverage2( param );
		toc = proc.time();
		show(toc-tic);
		# 2961.89  657.59 5140.73
		# 2966.73  650.12 5372.79
		# 2965.59  642.71 5434.28
	}
	{
		param$cputhreads = 2;
		tic = proc.time();
		ramwas3NormalizedCoverage3(param);
		toc = proc.time();
		show(toc-tic);
		# 1 thread: 1767.03  477.41 4122.09
		# 2 threads:   1.00    0.99 4067.15 
	}
} # NESDA
if(FALSE) { # RC2
	library(ramwas);
	param = list(
		dirbam = "D:/RW/RC2/bams/",
		dirproject = "D:/RW/RC2/",
		# dirtemp = "C:/AllWorkFiles/temp/",
		dirtemp = "E:/RaMWAS_temp/",
		filebamlist = "000_list_of_files.txt",
		scoretag = "AS",
		minscore = 100,
		cputhreads = 8,
		filecpgset = "C:/AllWorkFiles/Andrey/VCU/RaMWAS3/cpgset/hg19_1kG_MAF_0.01_chr1-22XY_cushaw3_50bp.rds",
		filenoncpgset = "C:/AllWorkFiles/Andrey/VCU/RaMWAS3/cpgset/hg19_1kG_MAF_0.01_chr1-22XY_cushaw3_50bp_nonCpG.rds",
		maxrepeats = 3,
		maxfragmentsize = 200,
		minfragmentsize = 50,
		filebam2sample = "bam2sample.txt",
		recalculate.QCs = TRUE,
		buffersize = 2e9
	);
	
	{
		tic = proc.time();
		rbam = pipelineProcessBam(bamname = '01A01SM1429N', param );
		readRDS('D:/RW/RC2/rds_qc/01A01SM1429N.qc.rds')
		toc = proc.time();
		show(toc-tic);
		# HP: 61.37   18.86  184.10 
	}
	{
		tic = proc.time();
		z = ramwas1scanBams( param );
		toc = proc.time();
		show(toc-tic);
		# HP: 61.37   18.86  184.10 
	}
	{
		tic = proc.time();
		z = ramwas2collectqc( param );
		toc = proc.time();
		show(toc-tic);
		# HP: 61.37   18.86  184.10 
	}
	{
		tic = proc.time();
		z = ramwas3NormalizedCoverage1( param );
		toc = proc.time();
		show(toc-tic);
		# HP, 4e8: 
	}
	{
		tic = proc.time();
		samplesums = ramwas3NormalizedCoverage2( param );
		toc = proc.time();
		show(toc-tic);
	}
	{
		param$cputhreads = 2;
		tic = proc.time();
		ramwas3NormalizedCoverage3(param);
		toc = proc.time();
		show(toc-tic);
	}
} # RC2


.combine.bams.qc = function( bamlist ) {
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


if(FALSE) { # test code
	### Test process single bam
	library(ramwas)
	param = list(
		dirbam = "D:/Cell_type/bams/",
		dirproject = "D:/Cell_type/",
		filebamlist = "D:/Cell_type/000_list_of_files.txt",
		scoretag = "AS",
		minscore = 100,
		cputhreads = 8,
		filecpgset = "C:/AllWorkFiles/Andrey/VCU/RaMWAS_2/code/Prepare_CpG_list/hg19/cpgset_hg19_SNPS_at_MAF_0.05.rds",
		filenoncpgset = NULL,
		maxrepeats = 3,
		maxfragmentsize=200,
		minfragmentsize=50,
		filebam2sample = "bam2sample.txt",
		recalculate.QCs = TRUE
	);
	
	{
		tic = proc.time();
		pipelineProcessBam(bamname="150114_WBCS014_CD20_150.bam", param);
		toc = proc.time();
		show(toc-tic);
	}
	{
		tic = proc.time();
		ramwas1scanBams(param);
		toc = proc.time();
		show(toc-tic);
	}
	
	
}
