# Scan a BAM file into Rbam object
bam.scanBamFile = function(bamfilename, scoretag = "MAPQ", minscore = 4){
    fields = c("qname","rname","pos","cigar","flag")
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
            levels(split.levels) = c(names(startlistfwd),
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
    bamfullname = .makefullpath(param$dirbam, paste0(bamname,".bam"))

    dir.create(param$dirrbam, showWarnings = FALSE, recursive = TRUE)
    dir.create(param$dirrqc, showWarnings = FALSE, recursive = TRUE)

    rdsbmfile = paste0( param$dirrbam, "/", basename(bamname), ".rbam.rds" );
    rdsqcfile = paste0( param$dirrqc, "/", basename(bamname), ".qc.rds" );

    savebam = TRUE;
    if( file.exists( rdsbmfile ) ){
        if( param$recalculate.QCs ){
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
        rbam = bam.scanBamFile(bamfilename = bamfullname,
                               scoretag = param$scoretag,
                               minscore = param$minscore);
    }

    rbam2 = bam.removeRepeats(rbam, param$maxrepeats);
    rbam2 = bam.chrXY.qc(rbam2);
    rbam2$qc$nbams = 1L;

    if( !is.null(param$filecpgset) ){
        cpgset = cachedRDSload(param$filecpgset);
        noncpgset = cachedRDSload(param$filenoncpgset);
        isocpgset = isocpgSitesFromCpGset(cpgset = cpgset,
                                          distance = param$maxfragmentsize);
        rbam3 = bam.hist.isolated.distances(rbam = rbam2,
                                            isocpgset = isocpgset,
                                            distance = param$maxfragmentsize);
        rbam4 = bam.coverage.by.density(rbam = rbam3,
                                        cpgset = cpgset,
                                        noncpgset = noncpgset,
                                        minfragmentsize = param$minfragmentsize,
                                        maxfragmentsize = param$maxfragmentsize)
        rbam5 = bam.count.nonCpG.reads(rbam = rbam4,
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
    rbam6$startsfwd=NULL;
    rbam6$startsrev=NULL;
    saveRDS( object = rbam6, file = rdsqcfile, compress = "xz");

    return(paste0("OK. ", bamname));
}

# Parallel job function
.ramwas1scanBamJob = function(bamname, param){
    cat(file = paste0(param$dirfilter,"/Log.txt"),
         date(), ", Process ", Sys.getpid(),
         ", Processing BAM: ", bamname, "\n",
         sep = "", append = TRUE);
    pipelineProcessBam(bamname = bamname, param = param);
}

# Step 1 of the pipeline
ramwas1scanBams = function( param ){
    param = parameterPreprocess(param);
    stopifnot( !is.null(param$bamnames));

    dir.create(param$dirfilter, showWarnings = FALSE, recursive = TRUE)
    cat(file = paste0(param$dirfilter,"/Log.txt"),
         date(), ", Scanning bams.", "\n", sep = "", append = FALSE);
    if( param$cputhreads > 1){
        cl = makeCluster(param$cputhreads);
        z = clusterApplyLB(cl,
                           param$bamnames,
                           .ramwas1scanBamJob,
                           param = param); #[1:64]
        stopCluster(cl);
    } else {
        z = character(length(param$bamnames));
        names(z) = param$bamnames;
        for(i in seq_along(param$bamnames)){ # i=1
            z[i] = .ramwas1scanBamJob(bamname = param$bamnames[i],
                                      param = param);
        }
    }
    cat(file = paste0(param$dirfilter,"/Log.txt"),
         date(), ", Done scanning bams.", "\n", sep = "", append = TRUE);
    return(invisible(z));
}
