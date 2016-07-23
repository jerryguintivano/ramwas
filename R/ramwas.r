### Caching environment
.ramwasEnv = new.env()

`%add%` <- function(x, y){
	if(is.null(x)) return(y);
	if(is.null(y)) return(x);
	l <- max(length(x), length(y))
	length(x) <- l
	length(y) <- l
	x[is.na(x)] <- 0
	y[is.na(y)] <- 0
	return(x + y)
}
.isAbsolutePath = function( pathname ){
	if( grepl("^~/", pathname) ) 
		return(TRUE)
	if( grepl("^.:(/|\\\\)", pathname) ) 
		return(TRUE)
	if( grepl("^(/|\\\\)", pathname) ) 
		return(TRUE)
	return(FALSE);
}
.makefullpath = function(path, filename){
	if( is.null(path) )
		return(filename);
	if(.isAbsolutePath(filename)) {
		return(filename) 
	} else {
		return( paste0(path, "/", filename) );
	}
}
if(FALSE){
	.isAbsolutePath( "C:/123" );  # TRUE
	.isAbsolutePath( "~123" );    # FALSE
	.isAbsolutePath( "~/123" );   # TRUE
	.isAbsolutePath( "/123" );    # TRUE
	.isAbsolutePath( "\\123" );    # TRUE
	.isAbsolutePath( "asd\\123" ); # FALSE
	.isAbsolutePath( "a\\123" );   # FALSE
	
	.makefullpath( "/dir", "C:/file");  # TRUE
	.makefullpath( "/dir", "~file");    # FALSE
	.makefullpath( "/dir", "~/file" );   # TRUE
	.makefullpath( "/dir", "/file" );    # TRUE
	.makefullpath( "\\dir", "\\file" );    # TRUE
	.makefullpath( "\\dir", "dir2\\file" ); # FALSE
	.makefullpath( "/dir", "dir2/file" );   # FALSE
}
.file.remove = function(x) {
	if( !is.null(x) )
		if( file.exists(x) )
			file.remove(x);
}

### Scan a file for parameters
parametersFromFile = function( .parameterfile ){
	source(.parameterfile, local = TRUE);
	.nms = ls();
	return(mget(.nms));
}
if(FALSE){ # test code
	param = parametersFromFile(.parameterfile = "D:/RW/NESDA/ramwas/param_file.txt");
	param
}

parseBam2sample = function( lines ){
	# remove trailing commas
	lines = gsub(pattern = ",$", replacement = "", lines);
	
	lines = gsub(pattern = "\\.bam,", replacement = ",", lines);
	lines = gsub(pattern = "\\.bam$", replacement = "",  lines);
	lines = gsub(pattern = " $", replacement = "",  lines);
	
	split.eq = strsplit(lines, split = "=", fixed = TRUE);
	samplenames = sapply(split.eq, `[`, 1);
	bamlist = strsplit(sapply(split.eq,tail,1), split = ",", fixed = TRUE);
	names(bamlist) = samplenames;
	
	# bamvec = unlist(bamlist, use.names = FALSE)
	# bamlist = lapply(bamlist, basename);
	return(bamlist);
}
if(FALSE) {
	lines = readLines( .makefullpath(param$dirproject, param$filebam2sample) );
}

processCommandLine = function(.arg = NULL){
	if( is.null(.arg))
		.arg=commandArgs(TRUE);
	# .arg = c("fileparam=\"D:/RW/CELL/param_file.txt\"","extraparameter=123123123")
	if( length(.arg)==0 ) {
		message("No arguments supplied"); 
	} else {
		for (.i in seq_along(.arg)) { # .i=1
			message("Input parameter: ", .arg[.i]);
			eval(parse(text=.arg[.i]));
			if(exists("fileparam")) {
				source(fileparam, local = TRUE);
				rm(fileparam);
			}
		}
	}
	return(mget(ls()));
}

### Fill in gaps in the parameter list
parameterPreprocess = function( param ){
	### Get from a file if param is not a list
	if(is.character(param)) {
		param = parametersFromFile(param);
	}
	
	# Set up directories 
	if( is.null(param$dirproject) ) param$dirproject = ".";
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
		param$dirfilter = .makefullpath(param$dirproject, param$dirfilter);
	}
	if( is.null(param$dirrbam) ) 
		param$dirrbam = "rds_rbam";
	param$dirrbam = .makefullpath( param$dirfilter, param$dirrbam);
	if( is.null(param$dirrqc) ) param$dirrqc = "rds_qc";
	param$dirrqc = .makefullpath( param$dirfilter, param$dirrqc);
	if( is.null(param$dirqc) ) param$dirqc = "qc";
	param$dirqc = .makefullpath( param$dirfilter, param$dirqc);
	
	### Filter parameters
	if( is.null(param$scoretag) ) param$scoretag = "mapq";
	if( is.null(param$minscore) ) param$minscore = 4;
	
	### More analysis parameters
	if( is.null(param$maxrepeats) ) param$maxrepeats = 0;
	if( is.null(param$cputhreads) ) param$cputhreads = detectCores();
	if( is.null(param$diskthreads) ) param$diskthreads = min(param$cputhreads,2);
	
	### BAM list processing
	if( is.null(param$bamnames) & !is.null(param$filebamlist)) {
		param$bamnames = readLines(.makefullpath(param$dirproject,param$filebamlist));
	}
	if( !is.null(param$bamnames)) {
		param$bamnames = gsub(pattern = "\\.bam$", replacement = "", param$bamnames);
	}
	
	### CV and MM
	if( is.null(param$cvnfolds) ) param$cvnfolds = 10;
	if( is.null(param$mmalpha) ) param$mmalpha = 0;
	if( is.null(param$mmncpgs) ) param$mmncpgs = 1000;
	
	
	### BAM2sample processing
	if( !is.null(param$filebam2sample) & is.null(param$bam2sample)) {
		filename = .makefullpath(param$dirproject, param$filebam2sample);
		param$bam2sample = parseBam2sample( readLines(filename) );
		rm(filename);
	}
	### Covariate file
	if( !is.null(param$filecovariates) & is.null(param$covariates)) {
		
		sep = "\t";
		if(grepl("\\.csv$",param$filecovariates)) 
			sep = ",";
		filename = .makefullpath(param$dirproject, param$filecovariates);
		param$covariates = read.table(filename, header = TRUE, sep = sep, stringsAsFactors = FALSE, check.names = FALSE);
		rm(filename);
	}
	if( !is.null(param$covariates)) {
		if( is.null(param$dircoveragenorm) ) param$dircoveragenorm = paste0("coverage_norm_",nrow(param$covariates));
		param$dircoveragenorm = .makefullpath(param$dirfilter, param$dircoveragenorm);	
		

		if( any(duplicated(param$covariates[[1]])) )
			stop("Repeated samples in the covariate file");
		
		if( !all(param$modelcovariates %in% names(param$covariates) ) )
			stop( paste("Covariates (modelcovariates) missing in covariates data frame:",
			 param$modelcovariates[ !(param$modelcovariates %in% names(param$covariates)) ]));
		if( is.null(param$modelPCs) ) 
			param$modelPCs = 0;
		if( !is.null(param$modeloutcome) )
			if( !( param$modeloutcome %in% names(param$covariates)) )
				stop( paste("Model outcome not present in covariate file:", param$modeloutcome));
		
		
		if( is.null(param$dirpca) ) {
			if( length(param$modelcovariates) > 0 ) {
				library(digest);
				hash = digest( object = paste(sort(param$modelcovariates), collapse = "\t"), algo = "crc32", serialize = FALSE);
				param$dirpca = sprintf("PCA_%02d_cvrts_%s",length(param$modelcovariates), hash);
			} else {
				param$dirpca = "PCA_00_cvrts";
			}
		}
		param$dirpca = .makefullpath(param$dircoveragenorm, param$dirpca);
		
		if( is.null(param$dirmwas) ) 
			param$dirmwas = paste0("Testing_",param$modeloutcome,"_",param$modelPCs,"_PCs");
		param$dirmwas = .makefullpath(param$dirpca, param$dirmwas);
		
		if( is.null(param$qqplottitle) ) {
			qqplottitle = paste0("Testing ",param$modeloutcome,"\n",param$modelPCs," PC",if(param$modelPCs!=1)"s"else"");
			if(length(param$modelcovariates)>0)
				qqplottitle = paste0(qqplottitle, " and ",length(param$modelcovariates)," covariate",if(length(param$modelcovariates)!=1)"s:\n"else": ",paste0(param$modelcovariates,collapse = ", "))
			param$qqplottitle = qqplottitle;
			rm(qqplottitle);
		}
		if( is.null(param$dircv) ) param$dircv = sprintf("%s/CV_%02d_folds", param$dirmwas, param$cvnfolds);
	} else if( !is.null(param$bam2sample) ) {
		if( is.null(param$dircoveragenorm) ) param$dircoveragenorm = paste0("coverage_norm_",length(param$bam2sample));
		param$dircoveragenorm = .makefullpath(param$dirfilter, param$dircoveragenorm);	
	} else {
		if( is.null(param$dircoveragenorm) ) param$dircoveragenorm = "coverage_norm";
		param$dircoveragenorm = .makefullpath(param$dirfilter, param$dircoveragenorm);	
	}

	if( is.null(param$dirtemp) ) param$dirtemp = "temp";
	param$dirtemp = .makefullpath(param$dircoveragenorm, param$dirtemp );

	### CpG set should exist
	if( !is.null(param$filecpgset) ) {
		stopifnot( file.exists(param$filecpgset) );
	}
	if( is.null(param$doublesize) ) param$doublesize = 4;
	if( is.null(param$recalculate.QCs) ) param$recalculate.QCs = FALSE;
	if( is.null(param$buffersize) ) param$buffersize = 1e9;

	if( is.null(param$minavgcpgcoverage) ) param$minavgcpgcoverage = 0.3;
	if( is.null(param$minnonzerosamples) ) param$minnonzerosamples = 0.3;

	if( is.null(param$usefilelock) ) param$usefilelock = FALSE;
	
	return(param);
}

### Save parameters to a file in output directory
parameterDump = function( dir, param, toplines = NULL) {
	message("Working in: ",dir);
	.dump = function(fid, param) {
		for( nm in names(param) ) { # nm = "modelcovariates"
			value = param[[nm]];
			if( is.data.frame(value) ) {
				txt = paste0("<Data frame ",nrow(value)," x ",ncol(value),">");
			} else if( is.list(value) ) {
				txt = paste0("<List of length ", length(value), ">");
			} else if( length(value) > 1 ) {
				if(nm == "bamnames") {
					txt = paste0("<Total ",length(value)," BAM names>");
				} else {
					txt = paste0(
						"c(\n", 
						paste0("  ",sapply(value, deparse),collapse = ",\n"),
						")");
				}
			} else {
				txt = deparse(value);
			}
			cat(file = fid, nm, "=", txt, "\n");
			# dput(param[[nm]], file = fid)
		}		
	}
	
	fid = file( paste0(dir, "/UsedSettings.txt"), "wt");
	writeLines(con = fid, c("### Parameters used to create the files in this directory",""));
	if( !is.null(toplines)) {
		.dump(fid, param[toplines]);
		writeLines(con = fid, text = "");
		.dump(fid, param[!(names(param) %in% toplines)]);
	} else {
		.dump(fid, param);
	}
	close(fid);
}
if(FALSE) {
	dir = "D:/RW";
	param = list(
		dirbam = "D:/RW/CELL/bams/",
		dirproject = "D:/RW/CELL/",
		filebamlist = "000_list_of_files.txt",
		dirtemp = NULL,
		scoretag = "AS",
		minscore = 100,
		cputhreads = 8,
		filecpgset    = "C:/AllWorkFiles/Andrey/VCU/RaMWAS3/cpgset/hg19_1kG_MAF_0.01_chr1-22_bowtie2_75bp.rds",
		filenoncpgset = "C:/AllWorkFiles/Andrey/VCU/RaMWAS3/cpgset/hg19_1kG_MAF_0.01_chr1-22_bowtie2_75bp_nonCpG.rds",
		maxrepeats = 3,
		maxfragmentsize=200,
		minfragmentsize=50,
		filebam2sample = "bam2sample1.txt",
		filecovariates = "Covariates.txt",
		# modelcovariates = NULL,
		# modelcovariates = "Peak SQRT",
		# modelcovariates = "Subject",
		# modelcovariates = c("RUFC % of aligned", "Peak SQRT"),
		# modelcovariates = c("RUFC % of aligned", "Peak SQRT", "Subject"),
		# modelcovariates = c("RUFC % of aligned", "Peak SQRT", "Subject", "Reads used for coverage"),
		# modelcovariates = c("Subject", "CellType", "Peak SQRT"),
		modelcovariates = c("Subject", "CellType", "Peak SQRT", "ChrY reads (%)"),
		modeloutcome = "CellType",
		modelPCs = 1,
		recalculate.QCs = FALSE,
		buffersize = 1e9,
		dirfilter = TRUE
	);
	param = parameterPreprocess(param);
	
	toplines = c("scoretag", "minscore")
	
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
			message(sprintf("Recorded %.f of %.f reads", qc$reads.recorded,qc$reads.total));
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
		message(sprintf("Recorded %.f of %.f reads", qc$reads.recorded,qc$reads.total));
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
if(FALSE){ # test code
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

.my.hist.plot = function(values, main2, firstvalue=0, xstep = 10, ...){
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
plot.qcHistScore = function(x, samplename="", xstep = 25, ...){
	.my.hist.plot(as.vector(x), main2 = paste0("Distribution of read scores\n",samplename), firstvalue=0, xstep = xstep, ...);
}
plot.qcHistScoreBF = function(x, samplename="", xstep = 25, ...){
	.my.hist.plot(as.vector(x), main2 = paste0("Distribution of read scores\n(including excluded reads)\n",samplename), firstvalue=0, xstep = xstep, ...);
}
plot.qcEditDist = function(x, samplename="", xstep = 5, ...){
	.my.hist.plot(as.vector(x), main2 = paste0("Distribution of edit distance\n",samplename), firstvalue=0, xstep = xstep, ...);
}
plot.qcEditDistBF = function(x, samplename="", xstep = 5, ...) {
	.my.hist.plot(as.vector(x), main2 = paste0("Distribution of edit distance\n(including excluded reads)\n",samplename), firstvalue=0, xstep = xstep, ...);
}
plot.qcLengthMatched = function(x, samplename="", xstep = 25, ...){
	.my.hist.plot(as.vector(x), main2 = paste0("Distribution of length of aligned part of read\n",samplename), firstvalue=1, xstep = xstep, ...);
}
plot.qcLengthMatchedBF = function(x, samplename="", xstep = 25, ...){
	.my.hist.plot(as.vector(x), main2 = paste0("Distribution of length of aligned part of read\n(including excluded reads)\n",samplename), firstvalue=1, xstep = xstep, ...);
}
plot.qcIsoDist = function(x, samplename="", xstep = 25, ...){
	.my.hist.plot(as.vector(x), main2 = paste0("Distribution of distances from read starts to isolated CpGs\n",samplename), firstvalue=0, xstep = xstep, ...);
}
plot.qcCoverageByDensity = function(y, samplename="", ...){
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
.histmean = function(x){
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
qcmean.qcChrY = function(x){ x[1]/x[2] }
qcmean.NULL = function(x){ NA }

### Make text line for the QC set
.qcTextHeader = {paste(sep = "\t",
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

.qccols = length(strsplit(.qcTextHeader,"\t",fixed = TRUE)[[1]])
.qcTextLine = function(qc, name){
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
.qcTextLineR = function(qc, name){
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
if(FALSE){ # test code
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
	message(rbam2$qc$bam.count.nonCpG.reads[1], "of", rbam2$qc$bam.count.nonCpG.reads[2], "reads are not covering CpGs");
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
if(FALSE){ # test code
	rbam = list( startsfwd = list(chr1=100), startsrev = list(chr1 = 103) );
	isocpgset = list(chr1 = 101);
	distance = 100;
	
	rbam2 = bam.hist.isolated.distances(rbam, isocpgset, distance);
	which(rbam2$qc$hist.isolated.dist1>0)
}

### Get average coverage vs. CpG density
bam.coverage.by.density = function( rbam, cpgset, noncpgset, minfragmentsize, maxfragmentsize){
	
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
	
	### CpG density of Non-CpG is always zero, not needed
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
if(FALSE){
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

### Fraction of reads on ChrX/Y
bam.chrXY.qc = function(rbam){
	strandfunX = function(st){c(length(st$chrX), sum(sapply(st,length)))};
	rbam$qc$chrX.count =  strandfunX(rbam$startsfwd) + strandfunX(rbam$startsfwd);
	class(rbam$qc$chrX.count) = "qcChrX"
	
	strandfunY = function(st){c(length(st$chrY), sum(sapply(st,length)))};
	rbam$qc$chrY.count =  strandfunY(rbam$startsfwd) + strandfunY(rbam$startsfwd);
	class(rbam$qc$chrY.count) = "qcChrY"
	
	return(rbam);
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
		# message("cachedRDSload: Using cache for: ", rdsfilename);
		return(base::get(x = cachename, envir = .ramwasEnv));
	} else {
		# message("cachedRDSload: Loading to cache: ", rdsfilename);
		data = readRDS(rdsfilename);
		base::assign(x = cachename, value = data, envir = .ramwasEnv);
		return(data);
	}
}
if(FALSE){ # test code
	rdsfilename = "C:/AllWorkFiles/Andrey/VCU/RaMWAS3/cpgset/hg19_1kG_MAF_0.01_chr1-22.rds";
	system.time({z = cachedRDSload(rdsfilename)});
	system.time({z = cachedRDSload(rdsfilename)});
	system.time({z = cachedRDSload(rdsfilename)});
}

### Coverage calculation
.calc.coverage.chr = function(startfrw, startrev, cpgs, fragdistr){
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
calc.coverage = function(rbam, cpgset, fragdistr){
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
if(FALSE){
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

pipelineSaveQCplots = function(param, rbam, bamname){
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
if(FALSE){
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

### Pipeline parts
pipelineProcessBam = function(bamname, param){
	# Used parameters: scoretag, minscore, filecpgset, maxrepeats
	
	param = parameterPreprocess(param);
	
	if( !is.null(param$filecpgset) && is.null(param$maxfragmentsize) )
		return("Parameter not set: maxfragmentsize");
	
	bamname = gsub("\\.bam$","",bamname);
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

### Fragment size distribution estimation
pipelineEstimateFragmentSizeDistribution = function( param ){
	
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

	cpgset = cachedRDSload(param$filecpgset);

	bams = param$bam2sample[[colnum]];
	
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
.ramwas1scanBamJob = function(bamname, param){
	cat(file = paste0(param$dirfilter,"/Log.txt"),
		 date(), ", Process ", Sys.getpid(),", Processing BAM: ", bamname, "\n", sep = "", append = TRUE);
	pipelineProcessBam(bamname = bamname, param = param);
}
ramwas1scanBams = function( param ){
	param = parameterPreprocess(param);
	stopifnot( !is.null(param$bamnames));

	cat(file = paste0(param$dirfilter,"/Log.txt"), 
		 date(), ", Scanning bams.", "\n", sep = "", append = FALSE);
	if( param$cputhreads > 1) {
		cl <- makeCluster(param$cputhreads);
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

.ramwas3coverageJob = function(colnum, param, nslices){
	library(ramwas);
	library(filematrix)
	coverage = pipelineCoverage1Sample(colnum, param);
	coverage = unlist(coverage, use.names = FALSE);

	start = 1;
	for( part in seq_len(nslices) ) {
		message("colnum =",colnum,"part =", part);
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
		cat(file = paste0(param$dircoveragenorm,"/Log.txt"),
			 date(), ", Process ", Sys.getpid(), 
			 ", Processing sample ", colnum, " ", names(param$bam2sample)[colnum], "\n",
			 sep = "", append = TRUE);
	});
	close(fm);
	
	return("OK");
}
.ramwas3transposeFilterJob = function(fmpart, param){
	library(filematrix);
	# fmpart = 1
	filename = paste0(param$dirtemp,"/RawCoverage_part",fmpart);
	if( !file.exists(paste0(filename,".bmat")) || !file.exists(paste0(filename,".desc.txt")) )
		return(paste0("Raw coverage slice filematrix not found: ", filename));
	fmraw = fm.open(filename, lockfile = param$lockfile2);
	mat = fmraw[];
	# mat = as.matrix(fmraw);

	fmout = fm.create( paste0(param$dirtemp,"/TrCoverage_part",fmpart), 
							 nrow = ncol(mat), ncol = 0, size = param$doublesize, lockfile = param$lockfile2);
	fmpos = fm.create( paste0(param$dirtemp,"/TrCoverage_loc",fmpart), 
							 nrow = 1, ncol = 0, type = "integer", lockfile = param$lockfile2);
	
	samplesums = rep(0, ncol(mat));
	
	### Sliced loop
	step1 = max(floor(32*1024*1024 / 8 / ncol(mat)),1);
	mm = nrow(mat);
	nsteps = ceiling(mm/step1);
	for( part in 1:nsteps ) { # part = 1
		message(fmpart, part, "of", nsteps);
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
		cat(file = paste0(param$dircoveragenorm,"/Log.txt"), 
			 date(), ", Process ", Sys.getpid(), 
			 ", Processing slice ", fmpart, "\n",  
			 sep = "", append = TRUE);
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
		cat(file = paste0(param$dircoveragenorm,"/Log.txt"),
			 date(), ", Process ", Sys.getpid(), ", Processing slice ", fmpart_offset[1], "\n", 
			 sep = "", append = TRUE);
	});

	rm(mat);
	gc();
	return("OK.");
}
ramwas3NormalizedCoverage = function( param ){
	# Prepare
	param = parameterPreprocess(param);
	param$fragdistr = as.double( readLines(con = paste0(param$dirfilter,"/Fragment_size_distribution.txt")));
	dir.create(param$dirtemp, showWarnings = FALSE, recursive = TRUE);
	dir.create(param$dircoveragenorm, showWarnings = FALSE, recursive = TRUE);

	if( !is.null(param$covariates) )
		param$bam2sample = param$bam2sample[param$covariates[[1]]];
	
	parameterDump(dir = param$dircoveragenorm, param = param,
					  toplines = c("dircoveragenorm", "dirtemp", "dirrbam",
					  				 "filebam2sample", "bam2sample",
					  				 "maxrepeats",
					  				 "minavgcpgcoverage", "minnonzerosamples",
					  				 "filecpgset",
					  				 "buffersize", "doublesize",
					  				 "cputhreads", "diskthreads"));
	
	### data dimensions
	cpgset = cachedRDSload(param$filecpgset);
	ncpgs = sum(sapply(cpgset, length));
	nsamples = length(param$bam2sample);
	
	### Check is all rbams are in place
	{
		message("Checking if all required Rbam files present");
		bams = unlist(param$bam2sample);
		for( bname in bams) {
			filename = paste0( param$dirrbam, "/", bname);
			if( file.exists(filename) ) {
				stop(paste0("Rbam file from bam2sample does not exist: ", filename));
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
		message("Creating ", nslices, " file matrices for raw coverage at: ", param$dirtemp);
		for( part in 1:nslices ) { # part = 1
			# cat("Creating raw coverage matrix slices", part, "of", nslices, "\n");
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
		message("Calculating and saving raw coverage");
		if(param$usefilelock) param$lockfile = tempfile();
		cat(file = paste0(param$dircoveragenorm,"/Log.txt"), 
			 date(), ", Calculating raw coverage.", "\n", sep = "", append = FALSE);
		library(parallel)
		if( param$cputhreads > 1) {
			cl = makeCluster(param$cputhreads);
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
		cat(file = paste0(param$dircoveragenorm,"/Log.txt"), 
			 date(), ", Done calculating raw coverage.", "\n", sep = "", append = TRUE);
		.file.remove(param$lockfile);
	}
	
	### Transpose the slices, filter by average and fraction of non-zeroes
	{
		message("Transposing coverage matrices and filtering CpGs by coverage");
		
		fm = fm.create( paste0(param$dirtemp,"/0_sample_sums"), nrow = nsamples, ncol = nslices);
		close(fm);
		
		cat(file = paste0(param$dircoveragenorm,"/Log.txt"), 
			 date(), ", Transposing coverage matrix, filtering CpGs.", "\n", sep = "", append = TRUE);
		if( param$diskthreads > 1 ) {
			if(param$usefilelock) param$lockfile2 = tempfile();
			library(parallel);
			cl = makeCluster(param$diskthreads);
			# cl = makePSOCKcluster(rep("localhost", param$diskthreads))
			z = clusterApplyLB(cl, 1:nslices, .ramwas3transposeFilterJob, param = param);
			stopCluster(cl);
			.file.remove(param$lockfile2);
		} else {
			for( fmpart in seq_len(nslices) ) { # fmpart = 5
				.ramwas3transposeFilterJob( fmpart, param);
			}
		}
		cat(file = paste0(param$dircoveragenorm,"/Log.txt"), 
			 date(), ", Done transposing coverage matrix, filtering CpGs.", "\n", sep = "", append = TRUE);
	}
	
	### Prepare CpG set for filtered CpGs	
	{
		message("Saving locations for CpGs which passed the filter");
		
		cpgsloc1e9 = cpgset;
		for( i in seq_along(cpgsloc1e9) ) {
			cpgsloc1e9[[i]] = cpgset[[i]] + i*1e9;
		}
		cpgsloc1e9 = unlist(cpgsloc1e9, recursive = FALSE, use.names = FALSE);
		
		kbblock = (128*1024)/8;
		step1 = max(floor(param$buffersize / (8 * nsamples)/kbblock),1)*kbblock;
		mm = ncpgs;
		nsteps = ceiling(mm/step1);
		cpgsloclist = vector("list",nsteps);
		for( part in 1:nsteps ) { # part = 1
			# cat( part, "of", nsteps, "\n");
			fr = (part-1)*step1 + 1;
			to = min(part*step1, mm);
			
			indx = as.vector(fm.load( paste0(param$dirtemp,"/TrCoverage_loc",part) ));
			cpgsloclist[[part]] = cpgsloc1e9[fr:to][indx];
		}
		rm(part, step1, mm, nsteps, fr, to, kbblock, indx);
		sliceoffsets = c(0L, cumsum(sapply(cpgsloclist, length)));

		cpgslocvec = unlist(cpgsloclist, use.names = FALSE);
		cpgslocmat = cbind( chr = as.integer(cpgslocvec %/% 1e9), position = as.integer(cpgslocvec %% 1e9));
		
		fm = fm.create.from.matrix( filenamebase = paste0(param$dircoveragenorm, "/CpG_locations"), mat = cpgslocmat);
		close(fm);
		writeLines(con = paste0(param$dircoveragenorm, "/CpG_chromosome_names.txt"), text = names(cpgset));
		rm(cpgsloc1e9, cpgsloclist, cpgslocvec, cpgslocmat);
	} # /CpG_locations, sliceoffsets
	
	### Sample sums
	{
		message("Gathering sample sums from ", nslices, " slices");
		
		mat = fm.load( paste0(param$dirtemp,"/0_sample_sums") );
		samplesums = rowSums(mat);
		rm(mat);
		fm = fm.create.from.matrix( paste0(param$dircoveragenorm,"/raw_sample_sums"), samplesums);
		close(fm);
	}
	
	### Normalize and combine in one matrix
	{
		message("Normalizing coverage and saving in one matrix");
		
		fmpart_offset_list = as.list(data.frame(rbind( seq_len(nslices), sliceoffsets[-length(sliceoffsets)])));

		### Create big matrix for normalized coverage
		fm = fm.create(paste0(param$dircoveragenorm, "/Coverage"), 
							nrow = nsamples, ncol = tail(sliceoffsets,1), size = param$doublesize);
		rownames(fm) = names(param$bam2sample);
		close(fm);
		
		cat(file = paste0(param$dircoveragenorm,"/Log.txt"),
			 date(), ", Normalizing coverage matrix.", "\n", sep = "", append = TRUE);
		### normalize and fill in
		if( param$diskthreads > 1 ) {
			
			if(param$usefilelock) param$lockfile1 = tempfile();
			if(param$usefilelock) param$lockfile2 = tempfile();
			library(parallel);
			cl = makeCluster(param$diskthreads);
			# cl = makePSOCKcluster(rep("localhost", param$diskthreads))
			z = clusterApplyLB(cl, fmpart_offset_list, .ramwas3normalizeJob, param = param, samplesums = samplesums);
			stopCluster(cl);
			.file.remove(param$lockfile1);
			.file.remove(param$lockfile2);
			
		} else {
			for( fmpart in seq_len(nslices) ) { # fmpart = 5
				.ramwas3normalizeJob( fmpart_offset_list[[fmpart]], param, samplesums);
			}
		}
		cat(file = paste0(param$dircoveragenorm,"/Log.txt"),
			 date(), ", Done normalizing coverage matrix.", "\n", sep = "", append = TRUE);
		
	}
	
	### Cleanup
	{
		message("Removing temporary files");
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

test1VariableOld = function(covariate, data, cvrtqr){
	# covariate = covariates1[[1]]
	mycov = matrix(covariate, nrow = 1);
	slice = data;
	cvqr0 = cvrtqr;
	
	# mycov[1000:1050] = NA;
	if( any(is.na(mycov)) ){
		keep = which(colSums(is.na(mycov))==0);
		
		mycov = mycov[,keep,drop=FALSE];
		slice = slice[,keep,drop=FALSE];
		cvqr0 = cvqr0[,keep,drop=FALSE];
		cvqr0 = t( qr.Q(qr(t(cvqr0))) );
	}
	
	# mycov = as.character(round(mycov));
	if( is.character(mycov) || is.factor(mycov)) {
		fctr = as.factor(mycov)
		dummy = t(model.matrix(~fctr)[,-1]);
		dummy = dummy - tcrossprod(dummy,cvqr0) %*% cvqr0;

		q = qr(t(dummy));
		keep = abs(diag(qr.R(q))) > .Machine$double.eps*ncol(mycov);
		
		mycov = t( qr.Q(qr(t(dummy))) );
		mycov[!keep,] = 0;
	} else {
		cvsumsq1 = sum( mycov^2 );
		mycov = mycov - tcrossprod(mycov,cvqr0) %*% cvqr0;
		cvsumsq2 = sum( mycov^2 );
		if( cvsumsq2 <= cvsumsq1 * .Machine$double.eps*ncol(mycov) ) {
			mycov[] = 0; 
		} else {
			mycov = mycov / sqrt(rowSums(mycov^2));
		}
	}
	
	# if( !is.character(mycov) & !any(is.na(mycov)) ){
	{
		slsumsq1 = rowSums(slice^2)
		slice = slice - tcrossprod(slice,cvqr0) %*% cvqr0;
		slsumsq2 = rowSums(slice^2)
		slice = slice / pmax(sqrt(rowSums(slice^2)), 1e-10);
		
		cr = tcrossprod(mycov, slice);
		
		nVarTested = nrow(cr)
		dfFull = ncol(slice) - nrow(cvqr0) - nVarTested;
		if( nVarTested == 1) {
			cor2tt = function(x) { return( x * sqrt( dfFull / (1 - pmin(x^2,1))));	}
			tt2pv = function(x) { return( (pt(-abs(x),dfFull)*2)); }
			tt = cor2tt(cr);
			pv = tt2pv(tt);
			
			### Check: 
			# c(tt[1], pv[1])
			# summary(lm( as.vector(covariate) ~ 0 + data[1,] + t(cvrtqr)))$coefficients[1,]
			return( list(correlation = cr, tstat = tt, pvalue = pv, nVarTested = nVarTested, dfFull = dfFull, statname = "") );
		} else {
			rsq = colSums(cr^2);
			rsq2F = function(x) { return( x / (1 - pmin(x,1)) * (dfFull/nVarTested) ); }
			F2pv = function(x) { return( pf(x, nVarTested, dfFull, lower.tail = FALSE) ); }
			ff = rsq2F(rsq);
			pv = F2pv(ff);
			
			### Check: 
			# c(ff[1], pv[1])
			# anova(lm( data[1,] ~ 0 + t(cvrtqr) + as.factor(as.vector(as.character(round(covariate))))))
			return( list(Rsquared = rsq, Fstat = ff, pvalue = pv, nVarTested = nVarTested, dfFull = dfFull, statname = paste0("-F_",nVarTested)) );
		}
	}
}

test1Variable = function(covariate, data, cvrtqr){
	# covariate = covariates1[[1]]
	mycov = matrix(covariate, nrow = 1);
	slice = data;
	cvqr0 = cvrtqr;
	
	# mycov[1000:1050] = NA;
	if( any(is.na(mycov)) ){
		keep = which(colSums(is.na(mycov))==0);
		
		mycov = mycov[,keep,drop=FALSE];
		slice = slice[keep,,drop=FALSE];
		cvqr0 = cvqr0[,keep,drop=FALSE];
		cvqr0 = t( qr.Q(qr(t(cvqr0))) );
	}
	
	# mycov = as.character(round(mycov));
	if( is.character(mycov) || is.factor(mycov)) {
		fctr = as.factor(mycov)
		dummy = t(model.matrix(~fctr)[,-1]);
		dummy = dummy - tcrossprod(dummy,cvqr0) %*% cvqr0;
		
		q = qr(t(dummy));
		keep = abs(diag(qr.R(q))) > .Machine$double.eps*ncol(mycov);
		
		mycov = t( qr.Q(qr(t(dummy))) );
		mycov[!keep,] = 0;
	} else {
		cvsumsq1 = sum( mycov^2 );
		mycov = mycov - tcrossprod(mycov,cvqr0) %*% cvqr0;
		cvsumsq2 = sum( mycov^2 );
		if( cvsumsq2 <= cvsumsq1 * .Machine$double.eps*ncol(mycov) ) {
			mycov[] = 0; 
		} else {
			mycov = mycov / sqrt(rowSums(mycov^2));
		}
	}
	
	### 
	nVarTested = nrow(mycov)
	dfFull = ncol(cvqr0) - nrow(cvqr0) - nVarTested;
	if(nVarTested == 1) {
		# SST = rowSums(slice^2);
		SST = colSumsSq(slice);
		
		cvD = (mycov %*% slice);
		# cvD2 = colSums(cvD^2);
		# cvD2 = cvD2^2;
		
		cvC = (cvqr0 %*% slice);
		cvC2 = colSumsSq( cvC );
		# cvC2 = colSums( cvC^2 );
		
		# SSR = colSums( cvD^2 );
		cr = cvD / sqrt(pmax(SST - cvC2, SST/1e16));
		
		cor2tt = function(x) { return( x * sqrt( dfFull / (1 - pmin(x^2,1))));	}
		tt2pv = function(x) { return( (pt(-abs(x),dfFull)*2)); }
		tt = cor2tt(cr);
		pv = tt2pv(tt);
		
		### Check: 
		# c(tt[1], pv[1])
		# summary(lm( as.vector(covariate) ~ 0 + data[1,] + t(cvrtqr)))$coefficients[1,]
		return( list(correlation = cr, tstat = tt, pvalue = pv, nVarTested = nVarTested, dfFull = dfFull, statname = "") );
		
	} else {
		# SST = rowSums(slice^2);
		SST = colSumsSq(slice);
		
		cvD = (mycov %*% slice);
		# cvD2 = colSums(cvD^2);
		cvD2 = colSumsSq(cvD);
		
		cvC = (cvqr0 %*% slice);
		# cvC2 = colSums( cvC^2 );
		cvC2 = colSumsSq( cvC );
		
		# SSR = colSums( cvD^2 );
		rsq = cvD2 / pmax(SST - cvC2, SST/1e16);
		
		# rsq = colSums(cr^2);
		rsq2F = function(x) { return( x / (1 - pmin(x,1)) * (dfFull/nVarTested) ); }
		F2pv = function(x) { return( pf(x, nVarTested, dfFull, lower.tail = FALSE) ); }
		ff = rsq2F(rsq);
		pv = F2pv(ff);
		
		### Check: 
		# c(ff[1], pv[1])
		# anova(lm( data[1,] ~ 0 + t(cvrtqr) + as.factor(as.vector(as.character(round(covariate))))))
		return( list(Rsquared = rsq, Fstat = ff, pvalue = pv, nVarTested = nVarTested, dfFull = dfFull, statname = paste0("-F_",nVarTested)) );
	}
}

if(FALSE){
	data = matrix(runif(1e6),2000,50000);
	cvrt = matrix(runif(1e6),5,2000);
	cvrt[1,] = 1;
	cvrtqr = t( qr.Q(qr(t(cvrt))) );
	
	### Full numerical
	covariate = runif(2000);
	# rez = test1Variable(covariate, data, cvrtqr)
	tm1 = system.time( {rez1 = test1Variable(covariate, data, cvrtqr)} )
	tm2 = system.time( {rez2 = test1VariableOld(covariate, t(data), cvrtqr)} )
	stopifnot(all.equal(rez1$pvalue,rez2$pvalue))
	sapply(rez1, `[`, 1)
	summary(lm( covariate ~ 0 + data[,1] + t(cvrtqr)))$coefficients[1,]
	show(cbind(tm1,tm2));
	# 6.14    1.42    4.56
	# 1.87    0.14    1.03
	
	### Full numerical with missing values
	covariate = runif(2000);
	covariate[1:100] = NA;
	tm1 = system.time( {rez1 = test1Variable(covariate, data, cvrtqr)} )
	tm2 = system.time( {rez2 = test1VariableOld(covariate, t(data), cvrtqr)} )
	stopifnot(all.equal(rez1$pvalue,rez2$pvalue))
	sapply(rez1, `[`, 1)
	summary(lm( covariate ~ 0 + data[,1] + t(cvrtqr)))$coefficients[1,]
	show(cbind(tm1,tm2));
	
	### Categorical
	covariate = as.character( round(runif(2000), 1) )
	tm1 = system.time( {rez1 = test1Variable(covariate, data, cvrtqr)} )
	tm2 = system.time( {rez2 = test1VariableOld(covariate, t(data), cvrtqr)} )
	stopifnot(all.equal(rez1$pvalue,rez2$pvalue))
	as.matrix(anova(lm( data[,1] ~ 0 + t(cvrtqr) + covariate)))
	sapply(rez1, `[`, 1)
	show(cbind(tm1,tm2));
	
	### Categorical, with missing values
	covariate = runif(2000);
	covariate[1:100] = NA;
	covariate = as.character( round(covariate, 1) )
	tm1 = system.time( {rez1 = test1Variable(covariate, data, cvrtqr)} )
	tm2 = system.time( {rez2 = test1VariableOld(covariate, t(data), cvrtqr)} )
	stopifnot(all.equal(rez1$pvalue,rez2$pvalue))
	as.matrix(anova(lm( data[,1] ~ 0 + t(cvrtqr) + covariate)))
	sapply(rez1, `[`, 1)
	show(cbind(tm1,tm2));
}
.testCovariates = function(covariates1, data, cvrtqr){
	# covariates1 = param$covariates[-1]
	# data = t(e$vectors[,seq_len(nonzeroPCs)])
	crF = vector("list", length(covariates1));
	pv  = vector("list", length(covariates1));
	nms = character(length(covariates1));
	for( i in seq_along(covariates1) ) { # i=1
		rez = test1Variable(covariates1[[i]], data, cvrtqr);
		pv[[i]] = as.vector(rez[[3]]);
		nms[i] = rez$statname;
		if(nchar(rez$statname)==0) {
			crF[[i]] = as.vector(rez$correlation);
		} else {
			crF[[i]] = as.vector(rez$Fstat);
		}
	}
	crF = data.frame(crF);
	names(crF) = paste0(names(covariates1),nms);

	pv = data.frame(pv);
	names(pv) = paste0(names(covariates1),nms);
	return(list(crF=crF, pv=pv));
	
}
.ramwas4PCAjob = function(rng, param, cvrtqr, rowsubset){
	# rng = rangeset[[1]];
	library(filematrix);
	fm = fm.open( paste0(param$dircoveragenorm, "/Coverage"), readonly = TRUE, lockfile = param$lockfile2);

	covmat = 0;

	step1 = ceiling( 128*1024*1024 / nrow(fm) / 8);
	mm = rng[2]-rng[1]+1;
	nsteps = ceiling(mm/step1);
	for( part in 1:nsteps ) { # part = 1
		cat( part, "of", nsteps, "\n");
		fr = (part-1)*step1 + rng[1];
		to = min(part*step1, mm) + rng[1] - 1;
		
		slice = fm[,fr:to];
		if( !is.null(rowsubset) )
			slice = slice[rowsubset,];
		slice = t(slice);
		
		slice = slice - tcrossprod(slice, cvrtqr) %*% cvrtqr; ### rowMeans(slice) == 0
		slice = slice / pmax(sqrt(rowSums(slice^2)), 1e-3);   ### rowSums(slice^2) == 1
		
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
orthoCovariates = function(covariates) {
	cvrtset = c(const = list(rep(1, nrow(covariates))), covariates);
	for( ind in which(sapply(cvrtset, typeof)=="character")) { # ind = 3
		fctr = factor(cvrtset[[ind]]);
		cvrtset[[ind]] = model.matrix(~fctr)[,-1];
		rm(fctr);
	}
	cvrtmat = matrix(unlist(cvrtset), nrow(covariates));
	cvrtqr = qr.Q(qr(cvrtmat));  ### tcrossprod(cvrtqr) - diag(nrow(cvrtqr))
	return(cvrtqr)
}

.ramwas45matchSamples = function( param ) {

	cvsamples = param$covariates[[1]];
	
	fm = fm.open( paste0(param$dircoveragenorm, "/Coverage"), readonly = TRUE);
	fmsamples = rownames(fm);
	ncpgs = ncol(fm);
	close(fm);
	
	rowsubset = match(cvsamples, fmsamples, nomatch = 0L);
	if( any(rowsubset==0) )
		stop( paste("Unknown samples in covariate file:", cvsamples[head(which(mch==0))]) );
	
	if( length(cvsamples) == length(fmsamples) )
	{
		if( all(rowsubset == seq_along(rowsubset)) )
		{
			rowsubset = NULL;
		}
	}
	return(list(rowsubset = rowsubset, ncpgs = ncpgs));
}

ramwas4PCA = function( param ){
	library(filematrix)
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
		rez = .ramwas45matchSamples( param );
		rowsubset = rez$rowsubset;
		ncpgs     = rez$ncpgs;
		cvsamples = param$covariates[[1]];
		rm(rez);
	} # rowsubset, ncpgs, cvsamples
	
	### Prepare covariates, defactor, 
	{
		message("Preparing covariates (splitting dummy variables, orthonormalizing)");
		cvrtqr = t(orthoCovariates( param$covariates[ param$modelcovariates ] ));
	} # cvrtqr
	
	### PCA part
	{
		### Calculate covmat from the data matrix
		{
			message("Calculating Principal Components");
			cat(file = paste0(param$dirpca,"/Log.txt"), 
				 date(), ", Running Principal Component Analysis.", "\n", sep = "", append = FALSE);

			if( param$diskthreads > 1 ) {
				rng = round(seq(1, ncpgs+1, length.out = param$diskthreads+1));
				rangeset = rbind( rng[-length(rng)], rng[-1]-1, seq_len(param$diskthreads));
				rangeset = lapply(seq_len(ncol(rangeset)), function(i) rangeset[,i])
				
				if(param$usefilelock) param$lockfile2 = tempfile();
				library(parallel);
				cl <- makeCluster(param$diskthreads);
				# cl = makePSOCKcluster(rep("localhost", param$diskthreads))
				covlist = clusterApplyLB(cl, rangeset, .ramwas4PCAjob, param = param, cvrtqr = cvrtqr, rowsubset = rowsubset);
				covmat = Reduce(f = `+`, x = covlist);
				stopCluster(cl);
				rm(cl, rng, rangeset, covlist);
				.file.remove(param$lockfile2);
			} else {
				covmat = .ramwas4PCAjob( rng = c(1, ncpgs, 0), param, cvrtqr, rowsubset);
			}
			cat(file = paste0(param$dirpca,"/Log.txt"), 
				 date(), ", Done running Principal Component Analysis.", "\n", sep = "", append = TRUE);
			
			saveRDS(file = paste0(param$dirpca,"/covmat.rds"), object = covmat, compress = FALSE);
			# covmat = readRDS(paste0(param$dirpca,"/covmat.rds"));
		} # covmat
		
		### Eigenvalue decomposition
		{
			message("Performing Eigenvalue Decomposition");
			e = eigen(covmat, symmetric=TRUE);
			nonzeroPCs = sum(abs(e$values/e$values[1]) > length(e$values)*.Machine$double.eps);
			saveRDS(file = paste0(param$dirpca,"/eigen.rds"), object = e, compress = FALSE);
			# e = readRDS(paste0(param$dirpca,"/eigen.rds"));
		} # e, nonzeroPCs
		
		# PCA plots
		{		
			message("Saving PCA plots");
			pdf(paste0(param$dirpca, "/PC_plot_covariates_removed.pdf"),7,7);
			pc100 = head(e$values,40)/sum(e$values)*100;
			plot(pc100, pch = 19, col="blue", main = "Principal components",
				  xlab = "PCs", ylab = "Variation Explained (%)",
				  yaxs = "i", ylim = c(0,pc100[1]*1.05), 
				  xaxs = "i", xlim = c(0, length(pc100)+0.5))
			for( i in 1:min(20,nonzeroPCs) ) { # i=1
				plot(e$vectors[,i], main=paste("PC",i), xlab = "Samples", ylab = "PC components", pch=19, col="blue1",
					  xlim = c(0, length(e$values)+0.5), xaxs="i");
				abline(h = 0, col = "grey")
			}
			dev.off();
		}
		
		# Save PCs and loadings
		{
			message("Saving PC values and vectors");
			PC_loads = e$vectors[,seq_len(min(20,nonzeroPCs))];
			rownames(PC_loads) = cvsamples;
			colnames(PC_loads) = paste0("PC",seq_len(ncol(PC_loads)));
			write.table(file = paste0(param$dirpca, "/PC_loadings.txt"), 
							x = data.frame(name=rownames(PC_loads),PC_loads), sep="\t", row.names = FALSE);
			PC_values = data.frame(PC_num = paste0("PC",seq_len(length(e$values))), e$values/sum(e$values))
			write.table(file = paste0(param$dirpca, "/PC_values.txt"),
						  x = data.frame(PC_num = paste0("PC",seq_len(length(e$values))),
						  					  var_explained = e$values/sum(e$values)),
						  sep="\t", row.names = FALSE, quote = FALSE);
		}
		
		# Saving PC vs. covariates association
		if(ncol(param$covariates) > 1) {
			message("Saving PC vs. covariates associations");
			cvrtqrconst = matrix(1/sqrt(length(e$values)),nrow = 1, ncol = length(e$values));
			testcov = .testCovariates(covariates1 = param$covariates[-1], data = e$vectors[,seq_len(nonzeroPCs)], cvrtqr = cvrtqrconst);
			write.table(file = paste0(param$dirpca, "/PC_vs_covariates_direct_corr.txt"), 
							x = data.frame(name=paste0("PC",seq_len(nonzeroPCs)), testcov$crF, check.names = FALSE),
							sep="\t", row.names = FALSE);
			write.table(file = paste0(param$dirpca, "/PC_vs_covariates_direct_pvalue.txt"), 
							x = data.frame(name=paste0("PC",seq_len(nonzeroPCs)), testcov$pv, check.names = FALSE),
							sep="\t", row.names = FALSE);
			if( length(param$modelcovariates) > 0) {
				testcov = .testCovariates(covariates1 = param$covariates[-1], data = e$vectors[,seq_len(nonzeroPCs)], cvrtqr = cvrtqr);
				write.table(file = paste0(param$dirpca, "/PC_vs_covariates_fixed_corr.txt"), 
								x = data.frame(name=paste0("PC",seq_len(nonzeroPCs)), testcov$crF, check.names = FALSE),
								sep="\t", row.names = FALSE);
				write.table(file = paste0(param$dirpca, "/PC_vs_covariates_fixed_pvalue.txt"), 
								x = data.frame(name=paste0("PC",seq_len(nonzeroPCs)), testcov$pv, check.names = FALSE),
								sep="\t", row.names = FALSE);
			}
		}
	}
}

.ramwas4MWASjob = function(rng, param, mwascvrtqr, rowsubset){
	# rng = rangeset[[1]];
	library(filematrix);
	fm = fm.open( paste0(param$dircoveragenorm, "/Coverage"), readonly = TRUE, lockfile = param$lockfile2);
	
	outmat = double(3*(rng[2]-rng[1]+1));
	dim(outmat) = c((rng[2]-rng[1]+1),3);
	
	# fmout = fm.open( paste0(param$dirmwas, "/Stats_and_pvalues"), lockfile = param$lockfile2);
	# covmat = 0;
	
	step1 = ceiling( 256*1024*1024 / nrow(fm) / 8);
	mm = rng[2]-rng[1]+1;
	nsteps = ceiling(mm/step1);
	for( part in 1:nsteps ) { # part = 1
		cat( part, "of", nsteps, "\n");
		fr = (part-1)*step1 + rng[1];
		to = min(part*step1, mm) + rng[1] - 1;
		
		slice = fm[,fr:to];
		if( !is.null(rowsubset) )
			slice = slice[rowsubset,];
		
		rez = test1Variable(covariate = param$covariates[[param$modeloutcome]],
								  data = slice, cvrtqr = mwascvrtqr)
		
		outmat[(fr:to) - (rng[1] - 1),] = cbind(rez[[1]], rez[[2]], rez[[3]]);
		
		fm$filelock$lockedrun( {
			cat(file = paste0(param$dirmwas,"/Log.txt"),
				 date(), ", Process ", Sys.getpid(), ", Job ", rng[3],
				 ", processing slice ", part, " of ", nsteps, "\n",
				 sep = "", append = TRUE);
		});
		rm(slice);
	}
	close(fm)
	
	fmout = fm.open(paste0(param$dirmwas, "/Stats_and_pvalues"), lockfile = param$lockfile2);
	fmout[rng[1]:rng[2],1:3] = outmat;
	close(fmout);
	
	return("OK");
}
qqplotFast = function(pvalues, ntests=NULL, ci.level=0.05) {
	
	if(is.null(ntests))
		ntests = length(pvalues);
	
	if(is.unsorted(pvalues))
		pvalues = sort.int(pvalues);
	
	ypvs = -log10(pvalues);
	xpvs = -log10(seq_along(ypvs) / ntests);
	
	if(length(ypvs) > 1000) {
		# need to filter a bit, make the plotting faster
		levels = as.integer( (xpvs - xpvs[1])/(tail(xpvs,1) - xpvs[1]) * 2000);
		keep = c(TRUE, diff(levels)!=0);
		levels = as.integer( (ypvs - ypvs[1])/(tail(ypvs,1) - ypvs[1]) * 2000);
		keep = keep | c(TRUE, diff(levels)!=0);
		keep = which(keep);
		ypvs = ypvs[keep];
		xpvs = xpvs[keep];
		# 		rm(keep)
	} else {
		keep = seq_along(ypvs)
	}
	mx = head(xpvs,1)*1.05;
	my = max(mx*1.15,head(ypvs,1))*1.05;
	plot(NA,NA, ylim = c(0,my), xlim = c(0,mx), xaxs="i", yaxs="i", 
		  xlab = expression("\u2013 log"[10]*"(p-value), expected under null"),
		  ylab = expression("\u2013 log"[10]*"(p-value), observed"));
	# xlab = "-Log10(p-value), expected under null", ylab = "-Log10(p-value), observed");
	lines(c(0,mx),c(0,mx),col="grey")
	points(xpvs, ypvs, col = "red", pch = 19, cex = 0.25);
	
	if(!is.null(ci.level)) {
		if((ci.level>0)&(ci.level<1)) {
			quantiles <- qbeta(p = rep(c(ci.level/2,1-ci.level/2),each=length(xpvs)), shape1 = keep, shape2 = ntests - keep + 1)
			quantiles = matrix(quantiles, ncol=2);
			
			lines( xpvs, -log10(quantiles[,1]), col="cyan4")
			lines( xpvs, -log10(quantiles[,2]), col="cyan4")
		}
	}
	legend("bottomright", c("P-values",sprintf("%.0f %% Confidence band",100-ci.level*100)),lwd = c(0,1), pch = c(19,NA_integer_), lty = c(0,1), col=c("red","cyan4"))
	if(length(pvalues)*2>ntests) {
		lambda = sprintf("%.3f",log10(pvalues[ntests/2]) / log10(0.5));
		legend("bottom", legend = bquote(lambda == .(lambda)), bty = "n")
		# 		text(mx, mx/2, bquote(lambda == .(lambda)), pos=2)
	}
}
.getCovariates = function(param){
	cvrtqr = t(orthoCovariates( param$covariates[ param$modelcovariates ] ));
	### Reading PCs, add as coveriates
	if( param$modelPCs > 0 ) {
		e = readRDS(paste0(param$dirpca,"/eigen.rds"));
		mwascvrtqr = rbind(cvrtqr, t(e$vectors[,seq_len(param$modelPCs)]));
		rm(e);
	} else {
		mwascvrtqr = cvrtqr;
	}
	stopifnot( all.equal( tcrossprod(mwascvrtqr), diag(nrow(mwascvrtqr))) );
	return(mwascvrtqr);
}
ramwas5MWAS = function( param ){
	library(filematrix)
	param = parameterPreprocess(param);
	dir.create(param$dirmwas, showWarnings = FALSE, recursive = TRUE);

	parameterDump(dir = param$dirmwas, param = param,
					  toplines = c("dirmwas", "dirpca", "dircoveragenorm",
					  				 "filecovariates", "covariates",
					  				 "modeloutcome", "modelcovariates", "modelPCs",
					  				 "qqplottitle",
					  				 "cputhreads"));
	
	
	message("Preparing for MWAS");
	
	
	### Get and match sample names
	{
		message("Matching samples in covariates and data matrix");
		rez = .ramwas45matchSamples( param );
		rowsubset = rez$rowsubset;
		ncpgs     = rez$ncpgs;
		cvsamples = param$covariates[[1]];
		rm(rez);
	} # rowsubset, ncpgs, cvsamples
	
	### Prepare covariates, defactor, 
	{
		message("Preparing covariates (splitting dummy variables, orthonormalizing)");
		mwascvrtqr = .getCovariates(param);
		# cvrtqr = t(orthoCovariates( param$covariates[ param$modelcovariates ] ));
	} # cvrtqr
	
	
	### Outpout matrix. Cor / t-test / p-value / q-value
	### Outpout matrix. R2  / F-test / p-value / q-value
	{
		message("Creating output matrix");
		fm = fm.create( paste0(param$dirmwas, "/Stats_and_pvalues"), nrow = ncpgs, ncol = 4);
		if( !is.character( param$covariates[[param$modeloutcome]] ) ) {
			colnames(fm) = c("cor","t-test","p-value","q-value");
		} else {
			colnames(fm) = c("R-squared","F-test","p-value","q-value");
		}
		close(fm);
	}
	
	### Running MWAS in parallel
	{
		message("Running MWAS");
		cat(file = paste0(param$dirmwas,"/Log.txt"), 
			 date(), ", Running methylome-wide association study.", "\n", sep = "", append = FALSE);
		if( param$cputhreads > 1 ) {
			rng = round(seq(1, ncpgs+1, length.out = param$cputhreads+1));
			rangeset = rbind( rng[-length(rng)], rng[-1]-1, seq_len(param$cputhreads));
			rangeset = lapply(seq_len(ncol(rangeset)), function(i) rangeset[,i])
			
			if(param$usefilelock) param$lockfile2 = tempfile();
			library(parallel);
			cl <- makeCluster(param$cputhreads);
			# cl = makePSOCKcluster(rep("localhost", param$cputhreads))
			clusterExport(cl, "test1Variable")
			clusterApplyLB(cl, rangeset, .ramwas4MWASjob, 
								param = param, mwascvrtqr = mwascvrtqr, rowsubset = rowsubset);
			stopCluster(cl);
			rm(cl, rng, rangeset);
			.file.remove(param$lockfile2);
		} else {
			covmat = .ramwas4MWASjob( rng = c(1, ncpgs, 0), param, cvrtqr, rowsubset);
		}
		cat(file = paste0(param$dirmwas,"/Log.txt"), 
			 date(), ", Done running methylome-wide association study.", "\n", sep = "", append = TRUE);
	}	
	
	### Fill in FDR column
	{
		message("Calculating FDR (q-values)");
		fm = fm.open( paste0(param$dirmwas, "/Stats_and_pvalues"));
		pvalues = fm[,3];
		pvalues[pvalues==0] = .Machine$double.xmin;
		ord = sort.list(pvalues);
		FDR = pvalues[ord] * length(pvalues) / seq_along(pvalues);
		FDR[length(FDR)] = min(FDR[length(FDR)], 1);
		FDR = rev(cummin(rev(FDR)));
		
		savevec = pvalues;
		savevec[ord] = FDR;
		fm[,4] = savevec;
		close(fm);
		
		sortedpv = pvalues[ord];
		rm(fm, pvalues, ord, FDR, savevec);
	} # sortedpv
		
	### QQ-plot
	{
		message("Creating QQ-plot");
		pdf(paste0(param$dirmwas, "/QQ_plot.pdf"),7,7);
		qqplotFast(sortedpv);
		title(param$qqplottitle);
		dev.off();
	}
}

if(FALSE){ # cluster
	param = parameterPreprocess(param);
	cl = makePSOCKcluster(rep("localhost", param$cputhreads))
	clusterEvalQ(cl, library(filematrix));
	z = clusterApplyLB(cl, 1, ramwas:::.ramwas3transposeFilterJob, param = param);
	stopCluster(cl);
} # cluster



ramwas6crossValidation = function(param) {
	param = parameterPreprocess(param);
	dir.create(param$dircv, showWarnings = FALSE, recursive = TRUE);
	parameterDump(dir = param$dircv, param = param,
					  toplines = c("dircv", "mmncpgs", "mmalpha", "cvnfolds",
					  				 "dirmwas", "dirpca", "dircoveragenorm",
					  				 "filecovariates", "covariates",
					  				 "modeloutcome", "modelcovariates", "modelPCs",
					  				 "qqplottitle",
					  				 "cputhreads"));
	
	nms = param$covariates[[1]];
	nsamples = length(nms);
	
	starts = floor(seq(1, nsamples+1, length.out = param$cvnfolds+1));
	shuffle = sample(nsamples);
	
	
	for( fold in seq_len(param$cvnfolds) ) { # fold = 1
		
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

ramwas7multiMarker = function(param) {
	library(glmnet)
	# library(filematrix);
	# library(ramwas);
	
	param = parameterPreprocess(param);
	
	forecastS = NULL;
	forecastC = NULL;
	
	cvrtqr = .getCovariates(param); # cvrtqr = ramwas:::.getCovariates(param)
	outcome = param$covariates[[ param$modeloutcome ]];
	outcomeR = outcome - crossprod(cvrtqr, cvrtqr %*% outcome);
	
	for( fold in seq_len(param$cvnfolds) ) { # fold = 1
		
		message("Processing fold ", fold, " of ", param$cvnfolds);
		dirmwas = sprintf("%s/fold_%02d", param$dircv, fold);
		rdsfile = paste0(dirmwas, "/exclude.rds");
		if( !file.exists( rdsfile ) )
			next;
		exclude = readRDS( rdsfile )
		
		if( is.null(forecastS) ) {
			forecastS = double(length(exclude));
			names(forecastS) = names(exclude);
			forecastC = integer(length(exclude));
		}
		
		# get p-values
		{
			fm = fm.open(paste0(dirmwas, "/Stats_and_pvalues"))
			colnames(fm)
			pv = fm[,3];
			close(fm)
		}
		
		# Find top param$mmncpgs CpGs
		{
			# tic = proc.time();
			pvthr = 10^((-300):0);
			fi = findInterval( pv, pvthr);
			tab = cumsum(tabulate(fi));
			upperfi = which(tab > param$mmncpgs)[1];
			set1 = which(fi <= upperfi);
			cpgset = set1[sort.list(pv[set1])[seq_len(param$mmncpgs)]];
			cpgset = sort.int(cpgset);
			rm(pvthr, fi, tab, upperfi, set1);
			# toc = proc.time();
			# show(toc-tic);
		} # 1.49			
		# {
		# 	tic = proc.time();
		# 	cpgset = sort.list(pv)[seq_len(param$mmncpgs)];
		# 	cpgset = sort.int(cpgset);
		# 	toc = proc.time();
		# 	show(toc-tic);
		# } # 38.92
		
		# get coverage
		{
			fm = fm.open( paste0(param$dircoveragenorm,"/Coverage"));
			coverage = fm[, cpgset];
			rownames(coverage) = rownames(fm);
			close(fm);
		}
		
		# Residualize
		resids = coverage - crossprod(cvrtqr, cvrtqr %*% coverage);
		
		if(param$mmncpgs == 1)
			resids = cbind(resids,resids);
		
		z = cv.glmnet(x = resids[!exclude,], y = outcome[!exclude], nfolds = param$cvnfolds, keep = TRUE, parallel = FALSE, alpha = param$mmalpha);
		z2 = predict(z, newx=resids[exclude,], type="response", s="lambda.min", alpha = param$mmalpha);
		
		forecastS[exclude] = forecastS[exclude] + z2;
		forecastC[exclude] = forecastC[exclude] + 1;
	}
	
	forecast = forecastS/forecastC;
	
	pdf( sprintf("%s/prediction_folds=%02d_CpGs=%d_alpha=%s.pdf", param$dircv, param$cvnfolds, param$mmncpgs, param$mmalpha) );
	plot( c(outcomeR), forecast, pch = 19, col = "blue", xlab = param$modeloutcome, ylab = "CV prediction",
			main = sprintf("Prediction success (residualized outcome)\n cor = %.3f / %.3f (Pearson / Spearman)", 
								cor(outcomeR, forecast, use = "complete.obs", method = "pearson"),
								cor(outcomeR, forecast, use = "complete.obs", method = "spearman")));
	legend(x = "bottomright", legend = c(paste0("# CpGs = ", param$mmncpgs), paste0("EN alpha = ", param$mmalpha)));
	plot( outcome, forecast, pch = 19, col = "blue", xlab = param$modeloutcome, ylab = "CV prediction",
			main = sprintf("Prediction success\n cor = %.3f / %.3f (Pearson / Spearman)", 
								cor(outcome, forecast, use = "complete.obs", method = "pearson"),
								cor(outcome, forecast, use = "complete.obs", method = "spearman")))
	legend(x = "bottomright", legend = c(paste0("# CpGs = ", param$mmncpgs), paste0("EN alpha = ", param$mmalpha)));
	dev.off();
	
	write.table( file = sprintf("%s/prediction_folds=%02d_CpGs=%d_alpha=%s.txt", param$dircv, param$cvnfolds, param$mmncpgs, param$mmalpha),
					 x = data.frame(samples = names(forecastS), outcome, outcomeR, forecast),
					 sep = "\t", row.names = FALSE);
	
	return( list(forecast = forecastS/forecastC, outcome = outcome, outcomeR = outcomeR) );
}



