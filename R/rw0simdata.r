# Sample a subset in groups of 5 (gr)
groupSample = function(len, size, gr){
    # len = 1000; size = 36; gr = 6; set.seed(18090212);
    groupstarts = sample(floor(len/gr)-1, size/gr);
    rez = rep(groupstarts, each = gr)*gr + (seq_len(gr)-1L)
    return(rez);
}


# Create artificial data set for vignettes and examples
ramwas0createArtificialData = function(dir,
                                       nsamples = 20,
                                       nreads = 1e6,
                                       ncpgs = 500e3,
                                       randseed = 18090212,
                                       verbose = TRUE){

    # dir="D:/temp";nsamples=20;nreads=1e6;
    # ncpgs=500e3;randseed=18090212;verbose=TRUE
    set.seed(randseed);

    # Create Directory
    dir.create(paste0(dir,"/bams"), showWarnings = FALSE, recursive = TRUE);

    # CpG locations
    {
        # ncpgs = 500e3
        chrlen = ncpgs*50
        probs = seq(1,0,length.out = chrlen) * (2 * ncpgs / chrlen);
        locs = which(probs > runif(chrlen)) + 1e6L;
        locs[1] = 800e3; # for non-CpG count
        rm(probs);

        cpgset = list( chr1 = locs );
        saveRDS(object = cpgset,
                file = paste0(dir,"/Simulated_chromosome.rds"),
                compress = "xz");
    } # locs, cpgset, chrlen, /Single_chromosome.rds

    # Fragment size distribution, with 150 median fragment size
    {
        x = 0:250;
        fragdistr = pmin(1.01*plogis((150-x)/20),1);
        # plot(x,fragdistr,pch=19)
        rm(x)
    } # fragdistr

    # # Exclude CpGs with density over 15
    # {
    #     coverage = ramwas:::calc.coverage(
    #         rbam = list(startsfwd = cpgset, startsrev = cpgset),
    #         cpgset = cpgset,
    #         fragdistr = fragdistr);
    #     locsgood = locs[ coverage[[1]] <= 15 ];
    #     rm(coverage, cpgset, locs)
    # } # locsgood, -locs, -cpgset
    locsgood = locs;

    # Age covariate and effects
    {
        age = sample(20:80, size = nsamples, replace = TRUE)
        cpgageset = groupSample( len = length(locsgood),
                                 size = length(locsgood)/100,
                                 gr = 6);
        cpgage1 = cpgageset[  1:(length(cpgageset)/2) ];
        cpgage2 = cpgageset[-(1:(length(cpgageset)/2))];
        rm(cpgageset)
    } # age, cpgage1, cpgage2

    # sex
    {
        sex = seq_len(nsamples) %% 2L;
        cpgsexset = groupSample( len = length(locsgood), size = 600, gr = 6);
        cpgsex1 = cpgsexset[  1:(length(cpgsexset)/2) ];
        cpgsex2 = cpgsexset[-(1:(length(cpgsexset)/2))];
        rm(cpgsexset)
    } # sex, cpgsex1, cpgsex2

    # Covariate file, bam list
    {
        cvrt = data.frame(
            samples = sprintf("BAM%03d",1:nsamples),
            age = age,
            sex = sex,
            stringsAsFactors = FALSE);
        write.table(file = paste0(dir,"/covariates.txt"),
                    x = cvrt,
                    sep = "\t",
                    row.names = FALSE);
        writeLines(con = paste0(dir,"/bam_list.txt"),
                   text = cvrt$samples);
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
        cpgprob[cpgsex1] =     sex[bam];
        cpgprob[cpgsex2] = 1 - sex[bam];

        # Pick CpG locations by the probabilities above
        cpglocs = sample(locsgood, prob = cpgprob,
                              size = nreads, replace = TRUE);

        # Add non-CpG reads
        cpglocs[1:(length(cpglocs)/100)] = 900e3 + (1:(length(cpglocs)/100));

        # read strand (0 - forward, 1 - reverse)
        readdir = sample(c(0L,1L), size = nreads, replace = TRUE);
        readlen = sample(75:70, size = nreads, prob = (6:1)^4, replace = TRUE);
        # distance to the CpG of interest
        roffset = sample(seq_along(fragdistr)-1L, prob = fragdistr,
                              size = nreads, replace = TRUE);
        # read start (and end) position, as reads are 1bp long
        readpos = cpglocs - (1L - 2L*readdir)*roffset - readlen * readdir;
        rm(roffset, cpglocs);

        # Sort reads by location
        ord = sort.list(readpos);
        readpos = readpos[ord];
        readdir = readdir[ord];
        rm(ord)

        # Write sam file
        filesam = sprintf("%s/bams/SAM%03d.sam", dir, bam);
        fid = file(description = filesam, open = "wt")
        writeLines(con = fid,
                   sprintf("@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:%d",
                           chrlen+2e6));
        writeLines(con = fid,
            sprintf("%06d\t%d\tchr1\t%d\t%d\t%dM\t*\t0\t0\t*\t*\tNM:i:%d",
                1:nreads, # 1 QNAME String [!-?A-~]{1,254} Query template NAME
                readdir*16L, # FLAG Int [0,216-1] bitwise FLAG
                # RNAME String \*|[!-()+-<>-~][!-~]* Reference sequence NAME
                readpos, # POS Int [0,231-1] 1-based leftmost mapping POSition
                readlen, # mqscore #MAPQ Int [0,2^8-1] MAPping Quality
                readlen, # 6 CIGAR String  CIGAR string
                # 7 RNEXT String  Ref. name of the mate/next read
                # 8 PNEXT Int [0,231-1] Position of the mate/next read
                # 9 TLEN Int [-231+1,231-1] observed Template LENgth
                # 10 SEQ String \*|[A-Za-z=.]+ segment SEQuence
                # 11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity+33
                75-readlen
            ));
        # flush(fid)
        close(fid);
        rm(readdir, readpos)

        asBam(file = filesam,
              destination = sprintf("%s/bams/BAM%03d", dir, bam),
              overwrite=TRUE, indexDestination=FALSE);
        file.remove(filesam)
    }
    return(invisible(NULL));
}
