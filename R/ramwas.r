### load required libraries
### Rsamtools to read BAMs
library(Rsamtools)
### GenomicAlignments to process CIGAR strings
library(GenomicAlignments)
library(snow)

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
if(FALSE) {
	.isAbsolutePath( 'C:/123' );  # TRUE
	.isAbsolutePath( '~123' );    # FALSE
	.isAbsolutePath( '~/123' );   # TRUE
	.isAbsolutePath( '/123' );    # TRUE
	.isAbsolutePath( '\\123' );    # TRUE
	.isAbsolutePath( 'asd\\123' ); # FALSE
	.isAbsolutePath( 'a\\123' );   # FALSE
	
}

### Scan a file for parameters

parametersFromFile = function(.parameterfile){
	source(.parameterfile, local = TRUE);
	.nms = ls();
	# .result = vector('list',length(.nms));
	# names(.result) = .nms;
	# for( .i in seq_along(.nms)) {
	# 	.result[[.i]] = get(.nms[.i]);
	# }
	# return(.result);
	return(mget(.nms));
}
if(FALSE) { # test code
	param = parametersFromFile(.parameterfile = 'D:/RW/NESDA/ramwas/param_file.txt');
	param
}

###
### BAM processing
###

bam.scanBamFile = function( bamfilename, scoretag = "mapq", minscore = 4){
	
	# header = scanBamHeader(bamfilename)
	# chrnames = names(header[[1]]$targets)
	
	### constants
	### More than 500 alignments per read is unlikely (although very possible)
	max.alignments.per.read = 500; 
	
	fields = c('qname','rname','pos','cigar','flag') 		
	# 'qname' is read name, 'rname' is chromosome
	tags = 'NM';# character();
	
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
	
	qc = list();
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
		### is at the begining of the loop to support calls of 'next'
		# if(is.logical(oldtail))
		# 	break;
		
		### Read 'yieldSize' rows
		bb = scanBam(file=bf, param=param)[[1]];
		if( length(bb[[1]])==0 )
			break;
		
		### Put tags in the main list
		bb = c(bb[names(bb) != 'tag'], bb$tag);
		# data.frame(lapply(bb,`[`, 1:60), check.rows = FALSE, stringsAsFactors = FALSE)

		stopifnot( length(bb[[scoretag]]) == length(bb[[1]]) )
		
		### Create output lists
		if(is.null(startlistfwd)) {
			startlistfwd = vector('list',length(levels(bb$rname)));
			startlistrev = vector('list',length(levels(bb$rname)));
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
		# 			stop('Empty BAM file (?).');
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
			keep = bitwAnd(bb$flag, 256L) == 0L 
			bb = lapply(bb,`[`,which(keep));
			rm(keep);
		}
		
		qc$reads.total = qc$reads.total %add% length(bb[[1]]);
		
		### Keep only aligned reads
		{
			keep = bitwAnd(bb$flag, 4L) == 0L;
			bb = lapply(bb,`[`,which(keep));
			rm(keep);
		}
		
		qc$reads.aligned = qc$reads.aligned %add% length(bb[[1]]);
			
		### Keep score >= minscore
		if( ! is.null(minscore) ) {
			score = bb[[scoretag]];
			keep = score >= minscore;
			keep[is.na(keep)] = FALSE
			bb = lapply(bb,`[`,which(keep));
			rm(keep);
		}
		
		qc$reads.recorded = qc$reads.recorded %add% length(bb[[1]]);
		qc$hist.score1 = qc$hist.score1 %add% tabulate(bb[[scoretag]]+1L);
		qc$hist.edit.dist1 = qc$hist.edit.dist1 %add% tabulate(bb$NM+1L);
		qc$hist.length.matched = qc$hist.length.matched %add% 
			tabulate(cigarWidthAlongQuerySpace(bb$cigar,after.soft.clipping = TRUE));
		
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
			levels(split.levels) = c(names(startlistfwd),paste0(names(startlistfwd),'-'));
			class(split.levels) = 'factor';
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
		cat(sprintf('Recorded %.f of %.f reads',qc$reads.recorded,qc$reads.total),'\n')
	}
	close(bf);
	rm(bf); # , oldtail
	
	startsfwd = startlistfwd;
	startsrev = startlistrev;
	
	### combine and sort lists in 'outlist'
	for( i in seq_along(startlistfwd) ) {
		startsfwd[[i]] = sort.int(unlist(startlistfwd[[i]]));
		startsrev[[i]] = sort.int(unlist(startlistrev[[i]]));
	}		

	bam = list(startsfwd = startsfwd, startsrev = startsrev, qc = qc);
	return( bam );
}
if(FALSE) { # test code
	bamfilename = 'D:/02H08SM142EZ.bam'; scoretag = "AS"; minscore = 10;
	rbam = bam.scanBamFile(bamfilename = 'D:/02H08SM142EZ.bam', scoretag = "AS", minscore = 10);
	sprintf('Recorded %.f of %.f reads',1e4,1e10)
}
###
### BAM QC / preprocessing
###

.remove.repeats.over.maxrep = function(vec, maxrep){
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
	.remove.repeats.over.maxrep(rep(1:10,1:10), 5L)
}
bam.removeRepeats = function(rbam, maxrep){
	if(maxrep<=0)
		return(rbam);
	# vec = c(floor(sqrt(0:99))); maxrep=5
	
	newbam = list(
		startsfwd = lapply( rbam$startsfwd, .remove.repeats.over.maxrep, maxrep),
		startsrev = lapply( rbam$startsrev, .remove.repeats.over.maxrep, maxrep),
		qc = rbam$qc);
	
	newbam$qc$frwrev.no.repeats = c(
		sum(sapply(newbam$startsfwd,length)),
		sum(sapply(newbam$startsrev,length)));
	
	newbam$qc$reads.recorded.no.repeats = sum(newbam$qc$frwrev.no.repeats);
	
	return(newbam);
}

### Non-CpG set of locations
noncpgSitesFromCpGset = function(cpgset, distance){
	noncpg = vector('list', length(cpgset));
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
	isocpg = vector('list',length(cpgset));
	names(isocpg) = names(cpgset);
	for( i in seq_along(cpgset) ) {	
		distbig = diff(cpgset[[i]]) >= distance;
		isocpg[[i]] = cpgset[[i]][ which( c(distbig[1],distbig[-1] & distbig[-length(distbig)], distbig[length(distbig)]) ) ];
	}
	return(isocpg);
}
if(FALSE) { # test code
	cpgset = readRDS('C:/AllWorkFiles/Andrey/VCU/RaMWAS_2/code/Prepare_CpG_list/hg19/spgset_hg19_SNPS_at_MAF_0.05.rds')
	noncpg = noncpgSitesFromCpGset(cpgset, 200);
	sapply(cpgset, typeof)
	sapply(noncpg, typeof)
	sapply(cpgset, length)
	sapply(noncpg, length)
	
	
	cpgset = lapply(1:10, function(x){return(c(1,1+x,1+2*x))})
	names(cpgset) = paste0('chr',seq_along(cpgset))
	show(cpgset);
	noncpg = noncpgSitesFromCpGset(cpgset, 3);
	show(noncpg);
	isocpg = isocpgSitesFromCpGset(cpgset, 3);
	show(isocpg);
}


### Count reads away from CpGs
.count.nonCpG.reads.forward = function( starts, cpglocations, distance ){
	### count CpGs before the read
	### count CpGs before and covered by the read
	ind = findInterval(c(starts-1L,starts+(distance-1L)), cpglocations);
	dim(ind)=c(length(starts),2);
	# cbind(ind, starts)
	return(c(sum(ind[,1] == ind[,2]),length(starts)));
}
.count.nonCpG.reads.reverse = function( starts, cpglocations, distance ){
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
	cat( rbam2$qc$bam.count.nonCpG.reads[1], 'of',rbam2$qc$bam.count.nonCpG.reads[2],'reads are not covering CpGs','\n' );
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
	return(rbam);
}
if(FALSE) { # test code
	rbam = list( startsfwd = list(chr1=100), startsrev = list(chr1 = 103) );
	isocpgset = list(chr1 = 101);
	distance = 100;
	
	rbam2 = bam.hist.isolated.distances(rbam, isocpgset, distance);
	which(rbam2$qc$hist.isolated.dist1>0)
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
	# lines( x, fsPredict(x, initparam), col='blue', lwd = 3)

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
	# lines(rezfit*lz$coefficients[2]+lz$coefficients[1], lwd = 4, col='red');
	
	return(rezfit);
}
if(FALSE) { # test code
	x = seq(0.01,0.99,0.01);
	y = sqrt(abs(x-0.5))*sign(x-0.5)
	plot(x,y)
	log.ss <- nls(y ~ SSlogis(x, phi1, phi2, phi3))
	z = SSlogis(x, 0.59699, 0.61320, 0.04599)
	lines(x, z, col='blue')
	
	
	# setwd('D:/RW/NESDA/ramwas/AS120_sba/');
	setwd('D:/RW/RC2/ramwas/AS38_gap0/');
	# setwd('D:/RW/Celltype//ramwas/AS120_sba/');
	lst = list.files(pattern = '\\.qc\\.');
	qcs = lapply(lst, function(x){load(x);return(bam);})
	histinfo = Reduce( `+`, lapply( lapply(qcs, `[[`, 'qcflt'), `[[`, 'hist.iso.dist.250'), init = 0);
	rng = range(histinfo[-(1:10)]);
	plot(histinfo/1e3, ylim = rng/1e3, pch=19, col='blue')
	
	hist.isolated.distances = histinfo;
	seqLength = 50;
	
	fit = estimateFragmentSizeDistribution(hist.isolated.distances, seqLength)
	
	x = seq_along(hist.isolated.distances)
	plot( x, hist.isolated.distances)
	lz = lm(hist.isolated.distances[seq_along(fit)] ~ fit)
	lines(fit*lz$coefficients[2]+lz$coefficients[1], lwd = 4, col='red');
	
}

### Cache CpG location files to avoid reloading.
cachedRDSload = function(rdsfilename){
	globalname = rdsfilename; #paste0('.ramwas.',rdsfilename);
	if( exists(x = globalname, envir = ramwas:::.ramwasEnv) ) {
		cat('Using cache','\n');
		return(base::mget(x = globalname, envir = ramwas:::.ramwasEnv));
	} else {
		cat('Loading','\n');
		data = readRDS(rdsfilename);
		base::assign(x = globalname, value = data, envir = ramwas:::.ramwasEnv);
		return(data);
	}
}
if(FALSE) { # test code
	rdsfilename = "C:/AllWorkFiles/Andrey/VCU/RaMWAS_2/code/Prepare_CpG_list/hg19//spgset_hg19_SNPS_at_MAF_0.05.rds";
	system.time({z = cachedRDSload(rdsfilename)});
	system.time({z = cachedRDSload(rdsfilename)});
	system.time({z = cachedRDSload(rdsfilename)});.
	
	globalname = paste0('.ramwas.',rdsfilename);
	if( exists(x = globalname, envir = ramwas:::.ramwasEnv) ) {
		return(base::mget(x = globalname, envir = ramwas:::.ramwasEnv));
	} else {
		data = readRDS(rdsfilename);
		base::assign(x = globalname, value = data, envir = ramwas:::.ramwasEnv);
		return(data);
	}
}

### Pipeline parts
pipelineProcessBam = function(bamname, param) {
	# Used parameters: scoretag, minscore, cpgfile, maxrepeats
	if( is.null(param$scoretag) )
		param$scoretag = "mapq";
	if( is.null(param$minscore) )
		param$minscore = 4;
	if( is.null(param$maxrepeats) )
		param$maxrepeats = 0;
	if( !is.null(param$cpgfile) && is.null(param$maxfragmentsize) )
		return('Parameter not set: maxfragmentsize');
	
	bamname = gsub('\\.bam$','',bamname);
	if( !.isAbsolutePath(bamname) && (length(param$bamdir)>0) ) {
		bamfullname = paste0(param$bamdir, '/', bamname, '.bam');
	} else {
		bamfullname = paste0(bamname, '.bam');
	}
	
	rdsbmdir = paste0( param$projectdir, "/", param$scoretag, "_", param$minscore, '/rbam_rds');
	rdsqcdir = paste0( param$projectdir, "/", param$scoretag, "_", param$minscore, '/rbam_qc_rds');
	dir.create(rdsbmdir, showWarnings = FALSE, recursive = TRUE)
	dir.create(rdsqcdir, showWarnings = FALSE, recursive = TRUE)
	
	rdsbmfile = paste0( rdsbmdir, '/', basename(bamname), '.rbam.rds' );
	rdsqcfile = paste0( rdsqcdir, '/', basename(bamname), '.qc.rds' );
	if( file.exists( rdsqcfile ) )
		return(paste0('Rbam qc rds file already exists: ',rdsqcfile));
	
	if( !file.exists( bamfullname ) )
		return(paste0('Bam file does not exist: ',bamfullname));
	
	rbam = bam.scanBamFile(bamfilename = bamfullname, scoretag = param$scoretag, minscore = param$minscore);
	
	rbam2 = bam.removeRepeats(rbam, param$maxrepeats);
	
	if( !is.null(param$cpgfile) ) {
		cpgset = cachedRDSload(param$cpgfile);
		isocpgset = isocpgSitesFromCpGset(cpgset = cpgset, distance = param$maxfragmentsize);
		rbam3 = bam.hist.isolated.distances(rbam = rbam2, isocpgset = isocpgset, distance = param$maxfragmentsize);
		
		# if( !is.null(param$noncpgfile)) {
		# 	noncpgset = cachedRDSload(param$noncpgfile);
		# } else {
		# 	noncpgset = noncpgSitesFromCpGset(cpgset = cpgset, distance = param$maxfragmentsize);
		# }
		rbam4 = bam.count.nonCpG.reads(rbam = rbam3, cpgset = cpgset, distance = param$maxfragmentsize);
	} else {
		rbam4 = rbam2;
	}
	
	saveRDS( object = rbam4, file = rdsbmfile, compress = 'xz');
	rbam5 = rbam4;
	rbam5$startsfwd=NULL;
	rbam5$startsrev=NULL;
	saveRDS( object = rbam5, file = rdsqcfile, compress = 'xz');
	return(paste0('OK. ', bamname));
}

### RaMWAS pipeline
ramwas1scanBams = function( param ){
	if(is.character(param)) {
		param = parametersFromFile(param);
	}
	
	bamnames = readLines(param$bamlistfile);
	
	if(is.null(param$cputhreads))
		param$cputhreads = 1;
	
	if( param$cputhreads > 1) {
		cl <- makeCluster(param$cputhreads)
		# clusterExport(cl, list = c("nms", "rvcfdir"))
		# nmslist = clusterSplit(cl, nms)
		# z = clusterApplyLB(cl, 1:8, function(i){ vcf = readRDS(paste0(rvcfdir,'/Rvcf_',nms[i],'.rds')); return(vcf$pos)})
		z = clusterApplyLB(cl, bamnames, pipelineProcessBam, param=param)
		stopCluster(cl)
	}
	return(z);
}
if(FALSE) { # test code
	param = list(
		bamdir = 'C:/Cell_type/bams/',
		projectdir = 'D:/Cell_type/',
		bamlistfile = 'c:/Cell_type/000_list_of_files.txt',
		scoretag = "AS",
		minscore = 100,
		cputhreads = 6,
		cpgfile = 'C:/AllWorkFiles/Andrey/VCU/RaMWAS_2/code/Prepare_CpG_list/hg19/spgset_hg19_SNPS_at_MAF_0.05.rds',
		noncpgfile = NULL,
		maxrepeats = 3,
		maxfragmentsize=200,
		minfragmentsize=50
	);
	
	param = list(
		bamdir = '.',
		projectdir = 'D:/Cell_type/',
		bamlistfile = 'c:/Cell_type/000_list_of_files.txt',
		scoretag = "AS",
		minscore = 100,
		cputhreads = 8,
		cpgfile = 'C:/AllWorkFiles/Andrey/VCU/RaMWAS_2/code/Prepare_CpG_list/hg19/spgset_hg19_SNPS_at_MAF_0.05.rds',
		noncpgfile = NULL,
		maxrepeats = 3,
		maxfragmentsize=200,
		minfragmentsize=50
	);
		
	library(ramwas)
	# ramwas1scanBams(param)
	
	bamname='150114_WBCS014_CD20_150.bam';
	{
		tic = proc.time();
		pipelineProcessBam(bamname, param);
		toc = proc.time();
		show(toc-tic);
	}
	
	bamnames = readLines(param$bamlistfile);
	bamnames[seq(1, length(bamnames), 2)] = paste0('C:/Cell_type/bams/', bamnames[seq(1, length(bamnames), 2)]);
	bamnames[seq(2, length(bamnames), 2)] = paste0('D:/Cell_type/bams/', bamnames[seq(2, length(bamnames), 2)]);
	
	library(snow)
	library(ramwas)
	# cl = makeCluster(param$cputhreads)
	cl = makeSOCKcluster(rep("localhost",param$cputhreads))
	# clusterCall(cl, function(){library('ramwas')});
	clusterEvalQ(cl, library(ramwas))
	
	# nmslist = clusterSplit(cl, nms)
	# z = clusterApplyLB(cl, 1:8, function(i){ vcf = readRDS(paste0(rvcfdir,'/Rvcf_',nms[i],'.rds')); return(vcf$pos)})
	report = clusterApplyLB(cl, bamnames, pipelineProcessBam, param=param)
	stopCluster(cl)
	
	file.exists('D:/Cell_type/bams/141106_WBCS011_BuCo_350.bam')
	file.exists('D:/Cell_type/bams//141201_WBCS011_BuCo_150.bam')
	
	### Cell type performance:
	# 15:35 - 17:40
	# 638 GB of BAM files
	# 2.7 GB of ramwas read start info
	# 60 KB of QC info
	
	cl = makeSOCKcluster(c("localhost","localhost"))
}

# snow makeCluster clusterEvalQ clusterApplyLB stopCluster

### Test C code wrapper
.conv <- function(a, b) .Call("convolve2", a, b)
