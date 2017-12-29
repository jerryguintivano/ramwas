# Make QC plots for an Rbam
pipelineSaveQCplots = function(param, rbam, bamname){
    filename = sprintf("%s/score/hs_%s.pdf", param$dirqc, bamname);
    dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE);
    pdf(filename);
    plot(rbam$qc$hist.score1, samplename = bamname);
    plot(rbam$qc$bf.hist.score1, samplename = bamname);
    dev.off();
    rm(filename);

    filename = sprintf("%s/edit_distance/ed_%s.pdf", param$dirqc, bamname);
    dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE);
    pdf(filename);
    plot(rbam$qc$hist.edit.dist1, samplename = bamname);
    plot(rbam$qc$bf.hist.edit.dist1, samplename = bamname);
    dev.off();
    rm(filename);

    filename = sprintf("%s/matched_length/ml_%s.pdf", param$dirqc, bamname);
    dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE);
    pdf(filename);
    plot(rbam$qc$hist.length.matched, samplename = bamname);
    plot(rbam$qc$bf.hist.length.matched, samplename = bamname);
    dev.off();
    rm(filename);

    if( !is.null(rbam$qc$hist.isolated.dist1) ){
        filename = sprintf("%s/isolated_distance/id_%s.pdf",
                        param$dirqc, bamname);
        dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE);
        pdf(filename);
        plot(rbam$qc$hist.isolated.dist1, samplename = bamname);
        dev.off();
        rm(filename);
    }
    if( !is.null(rbam$qc$avg.coverage.by.density) ){
        filename = sprintf("%s/isolated_distance/cbd_%s.pdf",
                        param$dirqc, bamname);
        dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE);
        pdf(filename);
        plot(rbam$qc$avg.coverage.by.density, samplename = bamname);
        dev.off();
        rm(filename);
    }
}

# Scan a BAM file into Rbam object
bam.scanBamFile = function(bamfilename, scoretag = "MAPQ", minscore = 4){
    fields = c("qname", "rname", "pos", "cigar", "flag");
    # "qname" is read name, "rname" is chromosome
    tags = "NM";# character();

    # Open the BAM file
    {
        if(nchar(scoretag) == 2){
            scoretag = toupper(scoretag);
            tags = c(tags, scoretag);
        } else {
            scoretag = tolower(scoretag);
            fields = c(fields,scoretag);
        }

        flag = scanBamFlag(isUnmappedQuery=NA, isSecondMateRead=FALSE);
        param = ScanBamParam(flag=flag, what=fields, tag=tags);
        bf = BamFile(bamfilename, yieldSize=1e6) ## typically, yieldSize=1e6
        open(bf);
        rm(fields, tags, flag);
    } # bf, param

    qc = list(nbams = 1L);

    startlistfwd = NULL;
    repeat{
        ### Read "yieldSize" rows
        bb = scanBam(file=bf, param=param)[[1]];
        if( length(bb[[1]])==0 )
            break;

        ### Put tags in the main list
        bb = c(bb[names(bb) != "tag"], bb$tag);
        # data.frame(lapply(bb,`[`, 1:60), check.rows = FALSE)

        # stopifnot( length(bb[[scoretag]]) == length(bb[[1]]) )

        ### Create output lists
        if(is.null(startlistfwd)){
            startlistfwd = vector("list",length(levels(bb$rname)));
            startlistrev = vector("list",length(levels(bb$rname)));
            names(startlistfwd) = levels(bb$rname);
            names(startlistrev) = levels(bb$rname);
            startlistfwd[] = list(list())
            startlistrev[] = list(list())
        } # startlistfwd, startlistrev

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

        if(length(bb[[1]])==0){
            message(sprintf("Recorded %.f of %.f reads",
                            qc$reads.recorded,qc$reads.total));
            next;
        }

        bb$matchedAlongQuerySpace =
            cigarWidthAlongQuerySpace(bb$cigar,after.soft.clipping = TRUE);

        qc$reads.aligned =
            qc$reads.aligned %add% length(bb[[1]]);
        qc$bf.hist.score1 =
            qc$bf.hist.score1 %add% tabulate(pmax(bb[[scoretag]]+1L,1L));
        qc$bf.hist.edit.dist1 =
            qc$bf.hist.edit.dist1 %add% tabulate(bb$NM+1L);
        qc$bf.hist.length.matched =
            qc$bf.hist.length.matched %add% tabulate(bb$matchedAlongQuerySpace);

        ### Keep score >= minscore
        if( !is.null(minscore) ){
            score = bb[[scoretag]];
            keep = score >= minscore;
            keep[is.na(keep)] = FALSE;
            if(!all(keep))
                bb = lapply(bb,`[`,which(keep));
            rm(keep);
        }

        qc$reads.recorded =
            qc$reads.recorded %add% length(bb[[1]]);
        qc$hist.score1 =
            qc$hist.score1 %add% tabulate(pmax(bb[[scoretag]]+1L,1L));
        qc$hist.edit.dist1 =
            qc$hist.edit.dist1 %add% tabulate(bb$NM+1L);
        qc$hist.length.matched =
            qc$hist.length.matched %add% tabulate(bb$matchedAlongQuerySpace);

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
            levels(split.levels) = c(
                                    names(startlistfwd),
                                    paste0(names(startlistfwd),"-"));
            class(split.levels) = "factor";
            splt = split( bb$startpos, split.levels, drop = FALSE);
            # print(sapply(splt,length))
            for( i in seq_along(startlistfwd) ){
                if( length(splt[i]) > 0 ){
                    startlistfwd[[i]][[length(startlistfwd[[i]])+1L]] =
                        splt[[i]];
                }
                if( length(splt[i+offset]) > 0 ){
                    startlistrev[[i]][[length(startlistrev[[i]])+1L]] =
                        splt[[i+offset]];
                }
            }
            rm(offset, split.levels, splt);
        } # startlistfwd, startlistrev
        message(sprintf("Recorded %.f of %.f reads",
                        qc$reads.recorded,qc$reads.total));
    }
    close(bf);
    rm(bf); # , oldtail

    startsfwd = startlistfwd;
    startsrev = startlistrev;

    ### combine and sort lists in "outlist"
    for( i in seq_along(startlistfwd) ){
        startsfwd[[i]] = sort.int(unlist(startlistfwd[[i]]));
        startsrev[[i]] = sort.int(unlist(startlistrev[[i]]));
    }
    gc();

    if( !is.null(qc$hist.score1))
        class(qc$hist.score1) = "qcHistScore";
    if( !is.null(qc$bf.hist.score1))
        class(qc$bf.hist.score1) = "qcHistScoreBF";
    if( !is.null(qc$hist.edit.dist1))
        class(qc$hist.edit.dist1) = "qcEditDist";
    if( !is.null(qc$bf.hist.edit.dist1))
        class(qc$bf.hist.edit.dist1) = "qcEditDistBF";
    if( !is.null(qc$hist.length.matched))
        class(qc$hist.length.matched) = "qcLengthMatched";
    if( !is.null(qc$bf.hist.length.matched))
        class(qc$bf.hist.length.matched) = "qcLengthMatchedBF";
    if( !is.null(qc$frwrev) )
        class(qc$frwrev) = "qcFrwrev";

    info = list(bamname = bamfilename,
                scoretag = scoretag,
                minscore = minscore,
                filesize = file.size(bamfilename));
    rbam = list(startsfwd = startsfwd,
                startsrev = startsrev,
                qc = qc,
                info = info);
    return( rbam );
}

# Scan a BAM, calculate QCs, generate plots
pipelineProcessBam = function(bamname, param){
    # Used parameters: scoretag, minscore, filecpgset, maxrepeats

    param = parameterPreprocess(param);

    if( !is.null(param$filecpgset) && is.null(param$maxfragmentsize) )
        return("Parameter not set: maxfragmentsize");

    bamname = gsub("\\.bam$", "", bamname, ignore.case = TRUE);
    bamfullname = makefullpath(param$dirbam, paste0(bamname,".bam"))

    dir.create(param$dirrbam, showWarnings = FALSE, recursive = TRUE)
    dir.create(param$dirrqc,  showWarnings = FALSE, recursive = TRUE)

    rdsbmfile = paste0(param$dirrbam, "/", basename(bamname), ".rbam.rds");
    rdsqcfile = paste0(param$dirrqc,  "/", basename(bamname), ".qc.rds");

    savebam = TRUE;
    rbam = NULL;
    if( file.exists(rdsqcfile) & file.exists(rdsbmfile) ){
        if( param$recalculate.QCs ){
            ### Precache the input rds file
            rbam = readRDS(rdsbmfile);
            savebam = FALSE;
        } else {
            if( !file.exists( bamfullname ) )
                return(paste0("Rbam rds file already exists",
                            " (no bam): ", rdsqcfile));
            if( file.mtime(rdsbmfile) > file.mtime(bamfullname) )
                return(paste0("Rbam rds file already exists",
                            " (newer than bam): ", rdsqcfile));  
        }
    }
    if( is.null(rbam) ){
        if( !file.exists( bamfullname ) )
            return(paste0("Bam file does not exist: ", bamfullname));
        rbam = bam.scanBamFile(
                    bamfilename = bamfullname,
                    scoretag = param$scoretag,
                    minscore = param$minscore);
    }

    rbam2 = bam.removeRepeats(rbam, param$maxrepeats);
    rbam2 = bam.chrXY.qc(rbam2);
    rbam2$qc$nbams = 1L;

    if( !is.null(param$filecpgset) ){
        cpgset = cachedRDSload(param$filecpgset);
        noncpgset = cachedRDSload(param$filenoncpgset);
        isocpgset = isocpgSitesFromCpGset(
                        cpgset = cpgset,
                        distance = param$maxfragmentsize);
        rbam3 = bam.hist.isolated.distances(
                        rbam = rbam2,
                        isocpgset = isocpgset,
                        distance = param$maxfragmentsize);
        rbam4 = bam.coverage.by.density(
                        rbam = rbam3,
                        cpgset = cpgset,
                        noncpgset = noncpgset,
                        minfragmentsize = param$minfragmentsize,
                        maxfragmentsize = param$maxfragmentsize)
        rbam5 = bam.count.nonCpG.reads(
                        rbam = rbam4,
                        cpgset = cpgset,
                        distance = param$maxfragmentsize);

        ### QC plots
        pipelineSaveQCplots(param, rbam5, basename(bamname));

    } else {
        rbam5 = rbam2;
    }
    # .qc qcmean(rbam5$qc$chrX.count)
    # rbam5$qc$chrX.count[1]/rbam5$qc$chrX.count[2]
    # message(.qcTextLine(rbam5$qc, bamname))

    if(savebam)
        saveRDS( object = rbam5, file = rdsbmfile, compress = "xz");
    rbam6 = rbam5;
    rbam6$startsfwd = NULL;
    rbam6$startsrev = NULL;
    saveRDS( object = rbam6, file = rdsqcfile, compress = "xz");

    return(NULL);
}

# Parallel job function
.ramwas1scanBamJob = function(bamname, param){
    ld = param$dirfilter;
    
    .log(ld, "%s, Process %06d, Starting BAM: %s",
        date(), Sys.getpid(), bamname);

    rez = pipelineProcessBam(bamname = bamname, param = param);
    
    .log(ld, "%s, Process %06d, Finished BAM: %s",
        date(), Sys.getpid(), bamname);
    return(rez);
}

# Step 1 of the pipeline
ramwas1scanBams = function( param ){
    param = parameterPreprocess(param);
    ld = param$dirfilter;
    
    # Parameter checks
    if( is.null(param$bamnames) )
        stop("BAM names must be specified. ",
            "See \"filebamlist\" or \"bamnames\" parameter.");

    if( !dir.exists(param$dirbam) )
        stop("Directory with BAM files not found: ",
            param$dirbam, "\n",
            "See \"dirbam\" parameter");
    
    for( nm in param$bamname ){ # nm = param$bamname[1]
        fn = makefullpath(param$dirbam, paste0(nm, ".bam"));
        if( !file.exists(fn) )
            stop("BAM file not found: ", fn);
    }
    
    if( !is.null(param$filecpgset) && is.null(param$maxfragmentsize) )
        stop("Parameter not set: maxfragmentsize");
    
    dir.create(param$dirfilter, showWarnings = FALSE, recursive = TRUE);

    .log(ld, "%s, Start ramwas1scanBams() call", date(), append = FALSE);
    
    nthreads = min(param$cputhreads, length(param$bamname));
    if( nthreads > 1 ){
        cl = makeCluster(nthreads);
        on.exit({stopCluster(cl);});
        logfun = .logErrors(ld, .ramwas1scanBamJob);
        z = clusterApplyLB(
                cl = cl,
                x = param$bamnames, 
                fun = logfun,
                param = param);
        tmp = sys.on.exit();
        eval(tmp);
        rm(tmp);
        on.exit();
    } else {
        z = vector("list", length(param$bamnames));
        names(z) = param$bamnames;
        for( i in seq_along(param$bamnames) ){ # i=1
            z[[i]] = .ramwas1scanBamJob(
                        bamname = param$bamnames[i],
                        param = param);
        }
    }
    .showErrors(z);
    .log(ld, "%s, Done ramwas1scanBams() call", date());
    return(invisible(z));
}
