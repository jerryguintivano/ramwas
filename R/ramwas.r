### Caching environment
.ramwasEnv = new.env()

# library(BiocCheck); setwd('C:/AllWorkFiles/Andrey/R/git/'); BiocCheck("ramwas_0.99.0.tar.gz")
# library(devtools); devtools::build_vignettes()
# browseVignettes(package = 'ramwas')
# BiocCheck Checking native routine registration
# getDLLRegisteredRoutines('ramwas')
# getDLLRegisteredRoutines.character
# getLoadedDLLs()
# package.skeleton(name = 'ramwas', path = 'C:/AllWorkFiles/Andrey/R/git/Skel/', environment = as.environment("package:ramwas"))

.notnull = function(x,replacement){if(is.null(x)){replacement}else{x}}
`%add%` = function(x, y){
	if(is.null(x)) return(y);
	if(is.null(y)) return(x);
	l = max(length(x), length(y))
	length(x) = l
	length(y) = l
	x[is.na(x)] = 0
	y[is.na(y)] = 0
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
	
	lines = gsub(pattern = "\\.bam,", replacement = ",", lines, ignore.case = TRUE);
	lines = gsub(pattern = "\\.bam$", replacement = "",  lines, ignore.case = TRUE);
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
			fileparam = NULL;
			message("Input parameter: ", .arg[.i]);
			eval(parse(text=.arg[.i]));
			if(!is.null(fileparam)) {
				source(fileparam, local = TRUE);
			}
		}
	}
	rm(fileparam);
	return(mget(ls()));
}

### Fill in gaps in the parameter list
parameterPreprocess = function( param ){
	### Get from a file if param is not a list
	if(is.character(param)) {
		param = parametersFromFile(param);
	}
	
	# Set up directories 
	if( is.null(param$dirproject) ) param$dirproject = getwd();
	param$dirbam = .makefullpath(param$dirproject, param$dirbam);
	
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
		param$bamnames = gsub(pattern = "\\.bam$", replacement = "", param$bamnames, ignore.case = TRUE);
	}
	
	### CV and MM
	if( is.null(param$cvnfolds) ) param$cvnfolds = 10;
	if( is.null(param$mmalpha) ) param$mmalpha = 0;
	if( is.null(param$mmncpgs) ) param$mmncpgs = 1000;
	stopifnot( param$mmncpgs > 1)
	
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
				# library(digest);
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
	    param$filecpgset = .makefullpath(param$dirproject, param$filecpgset);
		stopifnot( file.exists(param$filecpgset) );
	}
	if( !is.null(param$filenoncpgset) ) {
	    param$filenoncpgset = .makefullpath(param$dirproject, param$filenoncpgset);
		stopifnot( file.exists(param$filenoncpgset) );
    }
	
	if( is.null(param$doublesize) ) param$doublesize = 4;
	if( is.null(param$recalculate.QCs) ) param$recalculate.QCs = FALSE;
	if( is.null(param$buffersize) ) param$buffersize = 1e9;

	if( is.null(param$minavgcpgcoverage) ) param$minavgcpgcoverage = 0.3;
	if( is.null(param$minnonzerosamples) ) param$minnonzerosamples = 0.3;

	if( is.null(param$usefilelock) ) param$usefilelock = FALSE;
	
	if( is.null(param$randseed) ) param$randseed = 18090212; #Charles Darwin Date of birth: February 12, 1809
	
	if( is.null(param$toppvthreshold) ) param$toppvthreshold = 1e-6;
	
	# BioInformatics paramters
	
	if( is.null(param$bihost) ) param$bihost = "grch37.ensembl.org";
	if( is.null(param$bimart) ) param$bimart = "ENSEMBL_MART_ENSEMBL"; 
	
	# listDatasets(useMart(param$bimart))
	if( is.null(param$bidataset) ) {
		param$bidataset = "hsapiens_gene_ensembl"; 
	
		# listAttributes(useMart(biomart=param$bimart, dataset=param$bidataset))
		if( is.null(param$biattributes) ) param$biattributes = c("hgnc_symbol","entrezgene","strand");
		
		if( is.null(param$bifilters) ) param$bifilters = list(with_hgnc_transcript_name=TRUE);
		
		if( is.null(param$biflank) ) param$biflank = 0;
	}	
	
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
		.dump(fid, param[toplines[toplines %in% names(param)]]);
		writeLines(con = fid, text = "");
		.dump(fid, param[!(names(param) %in% toplines)]);
	} else {
		.dump(fid, param);
	}
	close(fid);
	return(invisible(NULL));
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
		bf = BamFile(bamfilename, yieldSize=1e6) ## typically, yieldSize=1e6
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
	.my.hist.plot(as.vector(x), main2 = paste0("Distribution of aligned length\n",samplename), firstvalue=1, xstep = xstep, ...);
}
plot.qcLengthMatchedBF = function(x, samplename="", xstep = 25, ...){
	.my.hist.plot(as.vector(x), main2 = paste0("Distribution of aligned length\n(including excluded reads)\n",samplename), firstvalue=1, xstep = xstep, ...);
}
plot.qcIsoDist = function(x, samplename="", xstep = 25, ...){
	.my.hist.plot(as.vector(x), main2 = paste0("Distances from read starts to isolated CpGs\n",samplename), firstvalue=0, xstep = xstep, ...);
}
plot.qcCoverageByDensity = function(x, samplename="", ...){
	# y = rbam$qc$avg.coverage.by.density
	y = x;
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

qcmean = function(x) UseMethod("qcmean", x)
qcmean.qcHistScore = function(x) { pmax(.histmean(x)-1,0) }
qcmean.qcHistScoreBF = function(x) { pmax(.histmean(x)-1,0) }
qcmean.qcEditDist = function(x) { pmax(.histmean(x)-1,0) }
qcmean.qcEditDistBF = function(x) { pmax(.histmean(x)-1,0) }
qcmean.qcLengthMatched = function(x) { .histmean(x) }
qcmean.qcLengthMatchedBF = function(x) { .histmean(x) }
qcmean.qcIsoDist = function(x) { .histmean(x) }
qcmean.qcFrwrev = function(x){ x[1]/(x[1]+x[2]) }
qcmean.qcNonCpGreads = function(x){ x[1]/(x[1]+x[2]) }
qcmean.qcCoverageByDensity = function(x){ (which.max(x)-1)/100 }
qcmean.qcChrX = function(x){ x[1]/x[2] }
qcmean.qcChrY = function(x){ x[1]/x[2] }
qcmean.NULL = function(x){ NA }



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

	# library(KernSmooth);
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



### RaMWAS pipeline functions

testPhenotype = function(phenotype, data, cvrtqr){
	# covariate = covariates1[[1]]
	mycov = matrix(phenotype, nrow = 1);
	slice = data;
	cvqr0 = cvrtqr;
	
	# mycov[1000:1050] = NA;
	if( any(is.na(mycov)) ){
		keep = which(colSums(is.na(mycov))==0);
		
		mycov = mycov[, keep, drop=FALSE];
		slice = slice[keep, , drop=FALSE];
		cvqr0 = cvqr0[, keep, drop=FALSE];
		cvqr0 = t( qr.Q(qr(t(cvqr0))) );
	}
	
	# mycov = as.character(round(mycov));
	if( is.character(mycov) || is.factor(mycov) ) {
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
	nVarTested = nrow(mycov);
	dfFull = ncol(cvqr0) - nrow(cvqr0) - nVarTested;
	if(nVarTested == 1) {
		if(dfFull <= 0)
			return(list(correlation = 0, tstat = 0, pvalue = 1, nVarTested = nVarTested, dfFull = dfFull, statname = ""));

		# SST = rowSums(slice^2);
		SST = colSumsSq(slice);
		
		cvD = (mycov %*% slice);
		# cvD2 = colSums(cvD^2);
		# cvD2 = cvD2^2;
		
		cvC = (cvqr0 %*% slice);
		cvC2 = colSumsSq( cvC );
		# cvC2 = colSums( cvC^2 );
		
		# SSR = colSums( cvD^2 );
		cr = cvD / sqrt(pmax(SST - cvC2, 1e-50, SST/1e15));
		
		cor2tt = function(x) { return( x * sqrt( dfFull / (1 - pmin(x^2,1))));	}
		tt2pv = function(x) { return( (pt(-abs(x),dfFull)*2)); }
		tt = cor2tt(cr);
		pv = tt2pv(tt);
		
		### Check: 
		# c(tt[1], pv[1])
		# summary(lm( as.vector(covariate) ~ 0 + data[1,] + t(cvrtqr)))$coefficients[1,]
		return( list(correlation = cr, tstat = tt, pvalue = pv, nVarTested = nVarTested, dfFull = dfFull, statname = "") );
		
	} else {
		if(dfFull <= 0)
			return( list(Rsquared = 0, Fstat = 0, pvalue = 1, nVarTested = nVarTested, dfFull = dfFull, statname = paste0("-F_",nVarTested)) );

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
	# rez = testPhenotype(covariate, data, cvrtqr)
	tm1 = system.time( {rez1 = testPhenotype(covariate, data, cvrtqr)} )
	tm2 = system.time( {rez2 = testPhenotypeOld(covariate, t(data), cvrtqr)} )
	stopifnot(all.equal(rez1$pvalue,rez2$pvalue))
	sapply(rez1, `[`, 1)
	summary(lm( covariate ~ 0 + data[,1] + t(cvrtqr)))$coefficients[1,]
	show(cbind(tm1,tm2));
	# 6.14    1.42    4.56
	# 1.87    0.14    1.03
	
	### Full numerical with missing values
	covariate = runif(2000);
	covariate[1:100] = NA;
	tm1 = system.time( {rez1 = testPhenotype(covariate, data, cvrtqr)} )
	tm2 = system.time( {rez2 = testPhenotypeOld(covariate, t(data), cvrtqr)} )
	stopifnot(all.equal(rez1$pvalue,rez2$pvalue))
	sapply(rez1, `[`, 1)
	summary(lm( covariate ~ 0 + data[,1] + t(cvrtqr)))$coefficients[1,]
	show(cbind(tm1,tm2));
	
	### Categorical
	covariate = as.character( round(runif(2000), 1) )
	tm1 = system.time( {rez1 = testPhenotype(covariate, data, cvrtqr)} )
	tm2 = system.time( {rez2 = testPhenotypeOld(covariate, t(data), cvrtqr)} )
	stopifnot(all.equal(rez1$pvalue,rez2$pvalue))
	as.matrix(anova(lm( data[,1] ~ 0 + t(cvrtqr) + covariate)))
	sapply(rez1, `[`, 1)
	show(cbind(tm1,tm2));
	
	### Categorical, with missing values
	covariate = runif(2000);
	covariate[1:100] = NA;
	covariate = as.character( round(covariate, 1) )
	tm1 = system.time( {rez1 = testPhenotype(covariate, data, cvrtqr)} )
	tm2 = system.time( {rez2 = testPhenotypeOld(covariate, t(data), cvrtqr)} )
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
		rez = testPhenotype(covariates1[[i]], data, cvrtqr);
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

orthonormalizeCovariates = function(covariates) {
	if(any(sapply(lapply(covariates, is.na),any)))
		stop("Missing values are not allowed in the covariates")
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

.matchCovmatCovar = function( param ) { # .ramwas45matchSamples

	cvsamples = param$covariates[[1]];
	
	fm = fm.open( paste0(param$dircoveragenorm, "/Coverage"), readonly = TRUE);
	fmsamples = rownames(fm);
	ncpgs = ncol(fm);
	close(fm);
	
	rowsubset = match(cvsamples, fmsamples, nomatch = 0L);
	if( any(rowsubset==0) )
		stop( paste("Unknown samples in covariate file:", cvsamples[head(which(rowsubset==0))]) );
	
	if( length(cvsamples) == length(fmsamples) )
	{
		if( all(rowsubset == seq_along(rowsubset)) )
		{
			rowsubset = NULL;
		}
	}
	return(list(rowsubset = rowsubset, ncpgs = ncpgs));
}


.getCovariates = function(param, rowsubset, normalize = TRUE){
	cvrtqr = param$covariates[ param$modelcovariates ];
	### Reading PCs, add as coveriates
	if( param$modelPCs > 0 ) {
		e = readRDS(paste0(param$dirpca,"/eigen.rds"));
		PCs = e$vectors[, seq_len(param$modelPCs), drop=FALSE];
		if(!is.null( rowsubset )) 
			PCs = PCs[rowsubset,];
		mwascvrtqr = cbind(cvrtqr, PCs);
		rm(e);
	} else {
		mwascvrtqr = cvrtqr;
	}
	# stopifnot( all.equal( tcrossprod(mwascvrtqr), diag(nrow(mwascvrtqr))) );
	if(normalize) {
		rez = t(orthonormalizeCovariates(mwascvrtqr));
	} else {
		rez = t(mwascvrtqr); #t(cbind(rep(1, nrow(mwascvrtqr)),mwascvrtqr));
	}
	return(rez);
}
pvalue2qvalue = function(pv, n = length(pv)){
	ord = sort.list(pv);
	FDR = pv[ord] * n / seq_along(pv);
	FDR[length(FDR)] = min(FDR[length(FDR)], 1);
	FDR = rev(cummin(rev(FDR)));
	
	rez = pv;
	rez[ord] = FDR;
	return(rez)
}

if(FALSE){ # cluster
	param = parameterPreprocess(param);
	cl = makePSOCKcluster(rep("localhost", param$cputhreads))
	clusterEvalQ(cl, library(filematrix));
	z = clusterApplyLB(cl, 1, ramwas:::.ramwas3transposeFilterJob, param = param);
	stopCluster(cl);
} # cluster

