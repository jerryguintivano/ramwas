# Sample a subset in groups of 5 (gr)
groupSample = function(len, size, gr){
    # len = 1000; size = 36; gr = 6; set.seed(18090212);
    groupstarts = sample(floor(len/gr)-1, size/gr);
    rez = rep(groupstarts, each = gr)*gr + (seq_len(gr)-1L);
    return(rez);
}

.ramwas0cadJob = function(
    rng, ld, 
    locfile, nreads, randseed, fragdistr,
    age, cpgage1, cpgage2, 
    sex, cpgsex1, cpgsex2){
    
    .log(ld, "%s, Process %06d, Job %02d, Start Generating BAMs %03d-%03d",
        date(), Sys.getpid(), rng[3], rng[1], rng[2]);

    locs = readRDS(locfile);
    locs = locs$chr1;
    chrlen = tail(locs, 1);

    # First 10% CpGs are unmethylated
    cpgprob = rep(1, length(locs));
    cpgprob[1:(length(locs)/10)] = 0;
    
    for( bam in rng[1]:rng[2] ){ # bam = rng[1]

        .log(ld, "%s, Process %06d, Job %02d, BAM: %03d, Generating data",
            date(), Sys.getpid(), rng[3], bam);

        # Change CpG probabilities
        # by age and case-control status
        cpgprob[cpgage1] =     age[bam]/100;
        cpgprob[cpgage2] = 1 - age[bam]/100;
        cpgprob[cpgsex1] =     sex[bam];
        cpgprob[cpgsex2] = 1 - sex[bam];

        # Pick CpG locations by the probabilities above
        set.seed(randseed + bam);
        cpglocs = sample(
                        x = locs,
                        prob = cpgprob,
                        size = nreads,
                        replace = TRUE);
        # cpglocs = rep(locs, length.out = nreads);

        # Add non-CpG reads
        cpglocs[seq_len(length(cpglocs)/100)] = 
                seq_len(length(cpglocs)/100) + max(locs);

        # read strand (0 - forward, 1 - reverse)
        readdir = sample(
                        x = c(0L, 1L), 
                        size = nreads, 
                        replace = TRUE);
        readlen = sample(
                        x = 75:70, 
                        size = nreads, 
                        prob = (6:1)^4, 
                        replace = TRUE);
        # distance to the CpG of interest
        roffset = sample(
                        x = seq_along(fragdistr)-1L, 
                        prob = fragdistr,
                        size = nreads, 
                        replace = TRUE);
        # read start (and end) position, as reads are 1bp long
        readpos = cpglocs - (1L - 2L*readdir) * roffset - readlen * readdir;
        rm(roffset, cpglocs);

        # Sort reads by location
        ord = sort.list(readpos);
        readpos = readpos[ord];
        readdir = readdir[ord];
        rm(ord);

        # Write sam file
        .log(ld, "%s, Process %06d, Job %02d, BAM: %03d, Saving SAM",
            date(), Sys.getpid(), rng[3], bam);
        filesam = sprintf("%s/SAM%03d.sam", ld, bam);
        fid = file(description = filesam, open = "wt");
        writeLines(
            con = fid,
            sprintf("@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:%d", 
                chrlen + 2e6));
        
        step1 = 8192;
        runto = nreads;
        nsteps = ceiling(runto/step1);
        for( part in seq_len(nsteps) ) { # part = 1
            fr = (part-1)*step1 + 1;
            to = min(part*step1, runto);
            set = fr:to;
            
            writeLines(con = fid,
                sprintf("%06d\t%d\tchr1\t%d\t%d\t%dM\t*\t0\t0\t*\t*\tNM:i:%d",
                set, # 1 QNAME String [!-?A-~]{1,254} Query template NAME
                readdir[set]*16L, # FLAG Int [0,216-1] bitwise FLAG
                # RNAME String \*|[!-()+-<>-~][!-~]* Reference sequence NAME
                readpos[set], # POS Int [0,231-1] 1-based leftmost mapping POS
                readlen[set], # mqscore #MAPQ Int [0,2^8-1] MAPping Quality
                readlen[set], # 6 CIGAR String  CIGAR string
                # 7 RNEXT String  Ref. name of the mate/next read
                # 8 PNEXT Int [0,231-1] Position of the mate/next read
                # 9 TLEN Int [-231+1,231-1] observed Template LENgth
                # 10 SEQ String \*|[A-Za-z=.]+ segment SEQuence
                # 11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity+33
                75 - readlen[set]));
        }
        rm(part, step1, runto, nsteps, fr, to);
        close(fid);
        rm(readdir, readpos, readlen);

        .log(ld, "%s, Process %06d, Job %02d, BAM: %03d, Convert SAM into BAM",
            date(), Sys.getpid(), rng[3], bam);
        asBam(
            file = filesam,
            destination = sprintf("%s/BAM%03d", ld, bam),
            overwrite = TRUE,
            indexDestination = FALSE);
        file.remove(filesam);
        .log(ld, "%s, Process %06d, Job %02d, BAM: %03d, BAM created",
            date(), Sys.getpid(), rng[3], bam);
    }
    return(invisible(NULL));
}

# Create artificial data set for vignettes and examples
ramwas0createArtificialData = function(
            dir,
            nsamples = 20,
            nreads = 1e6,
            ncpgs = 500e3,
            randseed = 18090212,
            threads = 1){
    
    ld = paste0(dir, "/bams");
    dir.create(ld, showWarnings = FALSE, recursive = TRUE);
    .log(ld, "%s, Start ramwas0createArtificialData() call", date(), 
        append = FALSE);

    # CpG locations
    {
        .log(ld, "%s, Generating %d CpG locations", date(), ncpgs);

        set.seed(randseed);
        gaps = rnbinom(
                    n = ncpgs,
                    size = 1,
                    prob = seq(1/5, 1/1000, length.out = ncpgs));
        locs = cumsum(2L + gaps) + 1e3L;
        chrlen = tail(locs, 1);
        rm(gaps);
        
        .log(ld, "%s, Saving CpG locations in: %s", date(), 
            paste0(dir, "/Simulated_chromosome.rds"));
        cpgset = list(chr1 = locs);
        saveRDS(
            object = cpgset,
            file = paste0(dir, "/Simulated_chromosome.rds"),
            compress = "xz");
        
    } # locs, cpgset, chrlen, /Single_chromosome.rds

    # Fragment size distribution, with 150 median fragment size
    {
        .log(ld, "%s, Generating fragment size distribution", date());
        x = 0:250;
        fragdistr = pmin(1.01*plogis((150-x)/20),1);
        # plot(x, fragdistr, pch=19);
        rm(x);
    } # fragdistr

    # Generating covariates
    {
        # Age covariate and effects
        .log(ld, "%s, Generating covariates", date());

        age = sample(x = 20:80, size = nsamples, replace = TRUE)
        cpgageset = groupSample( 
                        len = length(locs),
                        size = length(locs)/100,
                        gr = 6);
        cpgage1 = cpgageset[  1:(length(cpgageset)/2) ];
        cpgage2 = cpgageset[-(1:(length(cpgageset)/2))];
        rm(cpgageset)

        sex = seq_len(nsamples) %% 2L;
        cpgsexset = groupSample( 
                        len = length(locs),
                        size = 600,
                        gr = 6);
        cpgsex1 = cpgsexset[  1:(length(cpgsexset)/2) ];
        cpgsex2 = cpgsexset[-(1:(length(cpgsexset)/2))];
        rm(cpgsexset)
    } # age, cpgage1, cpgage2, sex, cpgsex1, cpgsex2

    # Covariate file, bam list
    {
        .log(ld, "%s, Savings covariates in: %s", date(),
            paste0(dir, "/covariates.txt"));
        cvrt = data.frame(
            samples = sprintf("BAM%03d", seq_len(nsamples)),
            age = age,
            sex = sex,
            stringsAsFactors = FALSE);
        write.table(
            file = paste0(dir, "/covariates.txt"),
            x = cvrt,
            sep = "\t",
            row.names = FALSE);

        .log(ld, "%s, Savings BAM list in: %s", date(),
            paste0(dir, "/bam_list.txt"));
        writeLines(
            con = paste0(dir, "/bam_list.txt"),
            text = cvrt$samples);
    } # cvrt, /covariates.txt, /bam_list.txt

    {
        .log(ld, "%s, Start generating BAM files in: %s", date(), ld);

        nthreads = min(threads, nsamples);
        if( nthreads > 1 ){
            rng = round(seq(1, nsamples+1, length.out = nthreads+1));
            rangeset = rbind(
                            rng[-length(rng)],
                            rng[-1] - 1,
                            seq_len(nthreads));
            rangeset = mat2cols(rangeset);

            # library(parallel);
            cl = makeCluster(nthreads);
            on.exit({
                stopCluster(cl);
            });
            logfun = .logErrors(ld, .ramwas0cadJob);
            rez = clusterApplyLB(
                            cl = cl,
                            x = rangeset,
                            fun = logfun,
                            ld = ld,
                            locfile = paste0(dir, "/Simulated_chromosome.rds"),
                            nreads = nreads,
                            randseed = randseed,
                            fragdistr = fragdistr,
                            age = age,
                            cpgage1 = cpgage1,
                            cpgage2 = cpgage2,
                            sex = sex,
                            cpgsex1 = cpgsex1,
                            cpgsex2 = cpgsex2);
            .showErrors(rez);
            tmp = sys.on.exit();
            eval(tmp);
            rm(tmp);
            on.exit();
            rm(cl, rng, rangeset, rez);
            gc();
        } else {
            rez = .ramwas0cadJob( 
                            rng = c(1, nsamples, 0),
                            ld = ld,
                            locfile = paste0(dir, "/Simulated_chromosome.rds"),
                            nreads = nreads,
                            randseed = randseed,
                            fragdistr = fragdistr, 
                            age = age,
                            cpgage1 = cpgage1,
                            cpgage2 = cpgage2,
                            sex = sex,
                            cpgsex1 = cpgsex1,
                            cpgsex2 = cpgsex2);
            if(is.character(rez))
                stop(rez);
        }
        .log(ld, "%s, Done generating BAM files in: %s", date(), ld);
    }
    return(invisible(NULL));
}
