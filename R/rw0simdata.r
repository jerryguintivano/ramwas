ramwas0createArtificialData = function(dir, nsamples = 20, nreads = 1e6, chr = "chr22", randseed = 18090212, verbose = TRUE){

	# Create Directory
	dir.create(paste0(dir,"/bams"), showWarnings = FALSE, recursive = TRUE);
	
	# CpG locations
	{
		
		cpgset = readRDS(system.file("extdata", "hg19_1kG_MAF_0.01_chr1-22.rds", package="ramwas"))
		cpgset = cpgset[chr]
		locs = cpgset[[1]]
		saveRDS(object = cpgset, file = paste0(dir,"/Single_chromosome.rds"), compress = "xz")
	} # locs, cpgset, /Single_chromosome.rds
	
	# Fragment size distribution, with 150 median fragment size
	{
		x = 0:250;
		fragdistr = pmin(1.01*plogis((150-x)/20),1);
		# plot(x,fragdistr,pch=19)
		rm(x)
	} # fragdistr

	# Exclude CpGs with density over 15
	{
		coverage = calc.coverage(rbam = list(startsfwd = cpgset, startsrev = cpgset), 
													 cpgset = cpgset, 
													 fragdistr = fragdistr);
		locsgood = locs[ coverage[[1]] <= 15 ];
		rm(coverage, cpgset, locs)
		
	} # locsgood, -locs, -cpgset
	
	# Age covariate and effects
	{
		set.seed(randseed);
		age = sample(20:80, size = nsamples, replace = TRUE)
		cpgageset = rep(sort.int(sample(floor(length(locsgood)/6)-1, length(locsgood)/600)), each = 6)*6+(0:5);
		cpgage1 = cpgageset[  1:(length(cpgageset)/2) ];
		cpgage2 = cpgageset[-(1:(length(cpgageset)/2))];
		rm(cpgageset)
	} # age, cpgage1, cpgage2
	
	# casecontrol
	{
		ccs = seq_len(nsamples) %% 2L;
		cpgccset = rep(sort.int(sample(floor(length(locsgood)/6)-1, 100)), each = 6)*6+(0:5);
		cpgccc1 = cpgccset[  1:(length(cpgccset)/2) ];
		cpgccc2 = cpgccset[-(1:(length(cpgccset)/2))];
		rm(cpgccset)
	} # ccs, cpgccc1, cpgccc2
	
	# Covariate file, bam list
	{
		cvrt = data.frame(
			samples = sprintf("BAM%03d",1:nsamples),
			age = age,
			casecontrol = ccs,
			stringsAsFactors = FALSE);
		write.table(file = paste0(dir,"/covariates.txt"), x = cvrt, sep = "\t", row.names = FALSE);
		writeLines(con = paste0(dir,"/bam_list.txt"), text = cvrt$samples);
	} # cvrt, /covariates.txt, /bam_list.txt
	
	# Main loop
	cpgprob = rep(1,length(locsgood));
	for( bam in seq_len(nsamples)){ # bam = 1
		
		if(verbose)
			message("Creating BAM ", bam, " of ", nsamples);
		
		# Change CpG probabilities
		# by age and case-control status
		cpgprob[cpgage1] =     age[bam]/100;
		cpgprob[cpgage2] = 1 - age[bam]/100;
		cpgprob[cpgccc1] =     ccs[bam];
		cpgprob[cpgccc1] = 1 - ccs[bam];
		
		# Pick CpG locations by the probabilities above
		cpglocs = sample(locsgood, prob = cpgprob,
							  size = nreads, replace = TRUE);
		# read strand (0 - forward, 1 - reverse)
		readdir = sample(c(0L,1L), size = nreads, replace = TRUE);
		# distance to the CpG of interest
		roffset = sample(seq_along(fragdistr)-1L, prob = fragdistr,
							  size = nreads, replace = TRUE);
		# read start (and end) position, as reads are 1bp long
		readpos = cpglocs - (1L - 2L*readdir)*roffset;
		rm(roffset, cpglocs);
		
		# Sort reads by location
		ord = sort.list(readpos);
		readpos = readpos[ord];
		readdir = readdir[ord];
		rm(ord)
		
		# Write sam file
		filesam = sprintf("%s/bams/SAM%03d.sam", dir, bam);
		fid = file(description = filesam, open = "wt")
		writeLines(con = fid, sprintf("@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:%s\tLN:%d",chr,tail(locsgood,1)+1e5L));
		writeLines(con = fid, 
					  sprintf("%06d\t%d\tchr22\t%d\t65\t1M\t*\t0\t0\tA\t*", 
					  		  1:nreads, # 1 QNAME String [!-?A-~]{1,254} Query template NAME
					  		  readdir*16L, # FLAG Int [0,216-1] bitwise FLAG
					  		  # RNAME String \*|[!-()+-<>-~][!-~]* Reference sequence NAME
					  		  readpos # POS Int [0,231-1] 1-based leftmost mapping POSition
					  		  # mqscore #MAPQ Int [0,2^8-1] MAPping Quality
					  		  # 6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
					  		  # 7 RNEXT String \*|=|[!-()+-<>-~][!-~]* Ref. name of the mate/next read
					  		  # 8 PNEXT Int [0,231-1] Position of the mate/next read
					  		  # 9 TLEN Int [-231+1,231-1] observed Template LENgth
					  		  # 10 SEQ String \*|[A-Za-z=.]+ segment SEQuence
					  		  # 11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity+33
					  ));
		# flush(fid)
		close(fid);
		rm(readdir, readpos)
		
		asBam(file = filesam, destination = sprintf("%s/bams/BAM%03d", dir, bam),
				overwrite=TRUE, indexDestination=FALSE);
		file.remove(filesam)
	}
	return(invisible(NULL));
}