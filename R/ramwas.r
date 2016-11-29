# Replace "x" with "replacement" if x is NULL
.notnull = function(x, replacement){
    if(is.null(x)){ replacement }else{ x }
}

# addition of vectors of mismatching length
`%add%` = function(x, y){
    if(is.null(x))
        return(y);
    if(is.null(y))
        return(x);
    if( length(x) > length(y) ){
        length(y) = length(x);
        y[is.na(y)] = 0;
    } else {
        length(x) = length(y);
        x[is.na(x)] = 0;
    }
    return(x + y);
}

# chech if path is absolute
# (common alternatives are flawed)
.isAbsolutePath = function( pathname ){
    if( grepl("^~/", pathname) )
        return(TRUE)
    if( grepl("^.:(/|\\\\)", pathname) )
        return(TRUE)
    if( grepl("^(/|\\\\)", pathname) )
        return(TRUE)
    return(FALSE);
}

# Get full path to the "filename" assuming current directory is "path"
.makefullpath = function(path, filename){
    if( is.null(path) )
        return(filename);
    if(.isAbsolutePath(filename)){
        return(filename)
    } else {
        return( paste0(path, "/", filename) );
    }
}

# Delete file, no warning on NULL or missing file
.file.remove = function(x){
    if( !is.null(x) )
        if( file.exists(x) )
            file.remove(x);
}

# Scan a file for parameters
parametersFromFile = function( .parameterfile ){
    source(.parameterfile, local = TRUE);
    .nms = ls();
    return(mget(.nms));
}

# Transform bam2sample file into a list
parseBam2sample = function( lines ){
    # remove trailing commas, .bam at name ends,
    # spaces around commas and "=", and trailing spaces
    lines = gsub(",$", "",       lines);
    lines = gsub("\\.bam,", ",", lines, ignore.case = TRUE);
    lines = gsub("\\.bam$", "",  lines, ignore.case = TRUE);
    lines = gsub(" $", "",       lines);
    lines = gsub(" ,", ",",      lines, fixed = TRUE);
    lines = gsub(", ", ",",      lines, fixed = TRUE);
    lines = gsub(" =", "=",      lines, fixed = TRUE);
    lines = gsub("= ", "=",      lines, fixed = TRUE);

    split.eq = strsplit(lines, split = "=", fixed = TRUE);
    samplenames = sapply(split.eq, `[`, 1);
    bamlist = strsplit(sapply(split.eq,tail,1), split = ",", fixed = TRUE);
    names(bamlist) = samplenames;

    return(bamlist);
}

# Get parameters from command line, return them in a list
# For "fileparam" parameter - source the file
processCommandLine = function(.arg = NULL){
    if( is.null(.arg))
        .arg=commandArgs(TRUE);
    # .arg = c("fileparam=\"D:/RW/CELL/param_file.txt\"","ext=123123123")
    if( length(.arg)==0 ){
        message("No arguments supplied");
    } else {
        for (.i in seq_along(.arg)){ # .i=1
            fileparam = NULL;
            message("Input parameter: ", .arg[.i]);
            eval(parse(text=.arg[.i]));
            if(!is.null(fileparam)){
                source(fileparam, local = TRUE);
            }
        }
    }
    rm(fileparam);
    return(mget(ls()));
}

# Fill in gaps in the parameter list
# Make paths absolute
parameterPreprocess = function( param ){
    ### Get from a file if param is not a list
    if(is.character(param)){
        param = parametersFromFile(param);
    }

    # Set up directories
    if( is.null(param$dirproject) ) param$dirproject = getwd();
    if(!is.null(param$dirbam))
        param$dirbam = .makefullpath(param$dirproject, param$dirbam);

    if( is.null(param$dirfilter) ){
        param$dirfilter = FALSE;
    }
    if( is.logical(param$dirfilter) ){
        if( param$dirfilter ){
            param$dirfilter = paste0( param$dirproject,
                            "/Filter_", param$scoretag, "_", param$minscore);
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
    if( is.null(param$maxrepeats) ) param$maxrepeats = 3;

    ### More analysis parameters
    if( is.null(param$cputhreads) ) param$cputhreads = detectCores();
    if( is.null(param$diskthreads)) param$diskthreads = min(param$cputhreads,2);

    ### BAM list processing
    if( is.null(param$bamnames) & !is.null(param$filebamlist)){
        param$filebamlist = .makefullpath(param$dirproject,param$filebamlist)
        param$bamnames = readLines(param$filebamlist);
    }
    if( !is.null(param$bamnames)){
        param$bamnames = gsub("\\.bam$", "", param$bamnames, ignore.case = TRUE)
    }

    ### CV and MM
    if( is.null(param$cvnfolds) ) param$cvnfolds = 10;
    if( is.null(param$mmalpha) ) param$mmalpha = 0;
    if( is.null(param$mmncpgs) ) param$mmncpgs = 1000;
    stopifnot(all( param$mmncpgs > 1 ))

    ### BAM2sample processing
    if( !is.null(param$filebam2sample) & is.null(param$bam2sample)){
        filename = .makefullpath(param$dirproject, param$filebam2sample);
        param$bam2sample = parseBam2sample( readLines(filename) );
        rm(filename);
    }
    if( is.null(param$bam2sample) & !is.null(param$bamnames) ){
        param$bam2sample = basename(param$bamnames);
        names(param$bam2sample) = basename(param$bamnames);
    }
    ### Covariate file
    if( !is.null(param$filecovariates) & is.null(param$covariates)){

        sep = "\t";
        if(grepl("\\.csv$",param$filecovariates))
            sep = ",";
        filename = .makefullpath(param$dirproject, param$filecovariates);
        param$covariates = read.table(filename, header = TRUE, sep = sep,
                              stringsAsFactors = FALSE, check.names = FALSE);
        rm(filename);
    }
    if( !is.null(param$covariates)){
        param$covariates[[1]] = as.character(param$covariates[[1]]);
        if( is.null(param$dircoveragenorm) )
            param$dircoveragenorm =
                paste0("coverage_norm_",nrow(param$covariates));
        param$dircoveragenorm =
            .makefullpath(param$dirfilter, param$dircoveragenorm);

        if( any(duplicated(param$covariates[[1]])) )
            stop("Repeated samples in the covariate file");

        if( !all(param$modelcovariates %in% names(param$covariates) ) )
            stop( paste("Covariates (modelcovariates) missing in covariates",
             param$modelcovariates[
                 !(param$modelcovariates %in% names(param$covariates)) ]));
        if( is.null(param$modelPCs) )
            param$modelPCs = 0;
        if( !is.null(param$modeloutcome) )
            if( !( param$modeloutcome %in% names(param$covariates)) )
                stop( paste("Model outcome not present in covariate file:",
                            param$modeloutcome));


        if( is.null(param$dirpca) ){
            if( length(param$modelcovariates) > 0 ){
                # library(digest);
                hash = digest(
                    object = paste(sort(param$modelcovariates),
                                   collapse = "\t"),
                    algo = "crc32", serialize = FALSE);
                param$dirpca = sprintf("PCA_%02d_cvrts_%s",
                                       length(param$modelcovariates), hash);
            } else {
                param$dirpca = "PCA_00_cvrts";
            }
        }
        param$dirpca = .makefullpath(param$dircoveragenorm, param$dirpca);

        if( is.null(param$dirmwas) )
            param$dirmwas = paste0("Testing_",param$modeloutcome,
                                   "_",param$modelPCs,"_PCs");
        param$dirmwas = .makefullpath(param$dirpca, param$dirmwas);

        if( is.null(param$qqplottitle) ){
            qqplottitle = paste0("Testing ",param$modeloutcome,"\n",
                         param$modelPCs," PC",if(param$modelPCs!=1)"s"else"");
            if(length(param$modelcovariates)>0)
                qqplottitle = paste0(qqplottitle, " and ",
                    length(param$modelcovariates)," covariate",
                    if(length(param$modelcovariates)!=1)"s:\n"else": ",
                    paste0(param$modelcovariates,collapse = ", "))
            param$qqplottitle = qqplottitle;
            rm(qqplottitle);
        }
        if( is.null(param$dircv) )
            param$dircv = sprintf("%s/CV_%02d_folds",
                                  param$dirmwas, param$cvnfolds);
    } else if( !is.null(param$bam2sample) ){
        if( is.null(param$dircoveragenorm) )
            param$dircoveragenorm =
                paste0("coverage_norm_",length(param$bam2sample));
        param$dircoveragenorm =
            .makefullpath(param$dirfilter, param$dircoveragenorm);
    } else {
        if( is.null(param$dircoveragenorm) )
            param$dircoveragenorm = "coverage_norm";
        param$dircoveragenorm =
            .makefullpath(param$dirfilter, param$dircoveragenorm);
    }

    if( is.null(param$dirtemp) ) param$dirtemp = "temp";
    param$dirtemp = .makefullpath(param$dircoveragenorm, param$dirtemp );

    ### CpG set should exist
    if( !is.null(param$filecpgset) ){
        param$filecpgset =
            .makefullpath(param$dirproject, param$filecpgset);
        stopifnot( file.exists(param$filecpgset) );
    }
    if( !is.null(param$filenoncpgset) ){
        param$filenoncpgset =
            .makefullpath(param$dirproject, param$filenoncpgset);
        stopifnot( file.exists(param$filenoncpgset) );
    }

    if( is.null(param$doublesize) ) param$doublesize = 4;
    if( is.null(param$recalculate.QCs) ) param$recalculate.QCs = FALSE;
    if( is.null(param$buffersize) ) param$buffersize = 1e9;

    if( is.null(param$minavgcpgcoverage) ) param$minavgcpgcoverage = 0.3;
    if( is.null(param$minnonzerosamples) ) param$minnonzerosamples = 0.3;

    if( is.null(param$usefilelock) ) param$usefilelock = FALSE;

    if( is.null(param$randseed) ) param$randseed = 18090212;
    # Famous person date of birth: February 12, 1809

    if( is.null(param$toppvthreshold) ) param$toppvthreshold = 50;

    # BioInformatics paramters

    if( is.null(param$bihost) ) param$bihost = "grch37.ensembl.org";
    if( is.null(param$bimart) ) param$bimart = "ENSEMBL_MART_ENSEMBL";

    # listDatasets(useMart(param$bimart))
    if( is.null(param$bidataset) ){
        param$bidataset = "hsapiens_gene_ensembl";

        # listAttributes(useMart(biomart=param$bimart, dataset=param$bidataset))
        if( is.null(param$biattributes) )
            param$biattributes = c("hgnc_symbol","entrezgene","strand");

        if( is.null(param$bifilters) )
            param$bifilters = list(with_hgnc_transcript_name=TRUE);

        if( is.null(param$biflank) )
            param$biflank = 0;
    }

    return(param);
}

# Save parameters "param" to a file in the "dir" directory
# Save "toplines" first, other parameters next
# Lists and "bamnames" are skipped
parameterDump = function(dir, param, toplines = NULL){
    message("Working in: ",dir);
    .dump = function(fid, param){
        for( nm in names(param) ){ # nm = "modelcovariates"
            pre = "";
            value = param[[nm]];
            if( is.data.frame(value) ){
                txt = paste0("<Data frame ",nrow(value)," x ",ncol(value),">");
                pre = "# ";
            } else if( is.list(value) ){
                txt = paste0("<List of length ", length(value), ">");
                pre = "# ";
            } else if( length(value) > 1 ){
                if(nm == "bamnames"){
                    txt = paste0("<Total ",length(value)," BAM names>");
                    pre = "# ";
                } else {
                    txt = paste0(
                        "c(\n",
                        paste0("  ",sapply(value, deparse),collapse = ",\n"),
                        ")");
                }
            } else {
                txt = deparse(value);
            }
            cat(file = fid, pre, nm, " = ", txt, "\n", sep = "");
            # dput(param[[nm]], file = fid)
        }
    }

    fid = file( paste0(dir, "/UsedSettings.txt"), "wt");
    writeLines(con = fid,
       text = c("## Parameters used to create the files in this directory",""));
    if( !is.null(toplines)){
        .dump(fid, param[toplines[toplines %in% names(param)]]);
        writeLines(con = fid, text = "");
        .dump(fid, param[!(names(param) %in% toplines)]);
    } else {
        .dump(fid, param);
    }
    close(fid);
    return(invisible(NULL));
}

# QC ploting functions
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
    .my.hist.plot(as.vector(x),
                  main2 = paste0("Distribution of read scores\n",samplename),
                  firstvalue=0,
                  xstep = xstep,
                  ...);
}
plot.qcHistScoreBF = function(x, samplename="", xstep = 25, ...){
    .my.hist.plot(as.vector(x),
        main2 = paste0("Distribution of read scores\n",
                       "(including excluded reads)\n",samplename),
        firstvalue=0,
        xstep = xstep,
        ...);
}
plot.qcEditDist = function(x, samplename="", xstep = 5, ...){
    .my.hist.plot(as.vector(x),
       main2 = paste0("Distribution of edit distance\n",samplename),
       firstvalue=0,
       xstep = xstep,
       ...);
}
plot.qcEditDistBF = function(x, samplename="", xstep = 5, ...){
    .my.hist.plot(as.vector(x),
        main2 = paste0("Distribution of edit distance\n",
                       "(including excluded reads)\n", samplename),
        firstvalue=0,
        xstep = xstep,
        ...);
}
plot.qcLengthMatched = function(x, samplename="", xstep = 25, ...){
    .my.hist.plot(as.vector(x),
        main2 = paste0("Distribution of aligned length\n", samplename),
        firstvalue=1,
        xstep = xstep,
        ...);
}
plot.qcLengthMatchedBF = function(x, samplename="", xstep = 25, ...){
    .my.hist.plot(as.vector(x),
        main2 = paste0("Distribution of aligned length\n",
                       "(including excluded reads)\n", samplename),
        firstvalue=1,
        xstep = xstep,
        ...);
}
plot.qcIsoDist = function(x, samplename="", xstep = 25, ...){
    .my.hist.plot(as.vector(x),
        main2 = paste0("Distances from read starts to isolated CpGs\n",
                       samplename),
        firstvalue=0,
        xstep = xstep,
        ...);
}
plot.qcCoverageByDensity = function(x, samplename="", ...){
    # y = rbam$qc$avg.coverage.by.density
    y = x;
    x = (seq_along(y)-1)/100;
    param = list(...);
    plotparam = list(
        x = x, y = y, type = "l", col = "magenta",
        lwd = 3, xaxs="i", yaxs="i", axes=FALSE,
        ylim = c(0, max(y, na.rm = TRUE)*1.1), xlim = range(x),
        xlab = "CpG density", ylab = "Coverage",
        main = paste0("Average coverage by CpG density\n", samplename));
    plotparam[names(param)] = param;
    do.call(plot, plotparam);
    axis(1, at = seq(0,tail(x,1)+2,by = 1), labels = seq(0,tail(x,1)+2,by=1)^2);
    axis(2);
}
.histmean = function(x){
    return( sum(x * seq_along(x)) / pmax(sum(x), .Machine$double.xmin) );
}

# QC single number summary functions
qcmean = function(x) UseMethod("qcmean", x)
qcmean.qcHistScore = function(x){ pmax(.histmean(x)-1,0) }
qcmean.qcHistScoreBF = function(x){ pmax(.histmean(x)-1,0) }
qcmean.qcEditDist = function(x){ pmax(.histmean(x)-1,0) }
qcmean.qcEditDistBF = function(x){ pmax(.histmean(x)-1,0) }
qcmean.qcLengthMatched = function(x){ .histmean(x) }
qcmean.qcLengthMatchedBF = function(x){ .histmean(x) }
qcmean.qcIsoDist = function(x){ .histmean(x) }
qcmean.qcFrwrev = function(x){ x[1]/(x[1]+x[2]) }
qcmean.qcNonCpGreads = function(x){ x[1]/(x[1]+x[2]) }
qcmean.qcCoverageByDensity = function(x){ (which.max(x)-1)/100 }
qcmean.qcChrX = function(x){ x[1]/x[2] }
qcmean.qcChrY = function(x){ x[1]/x[2] }
qcmean.NULL = function(x){ NA }

# Take a sorted vector "vec", remove repeated values
# repeating over "maxrep" times (keep first "maxrep")
remove.repeats.over.maxrep = function(vec, maxrep){
    if( is.unsorted(vec) )
        vec = sort.int(vec);
    if( maxrep > 0 ){
        kill = which(diff(vec, maxrep) == 0L);
        if(length(kill)>0){
            vec[kill] = 0L;
            vec = vec[vec!=0L];
        }
    }
    return(vec);
}

# Remove reads starting from the same position,
# on the same strand, and
# repeating over "maxrep" times (keep first "maxrep")
# Calculate some QC
bam.removeRepeats = function(rbam, maxrep){
    if(maxrep>0){
        newbam = list(
            startsfwd = lapply( rbam$startsfwd,
                                remove.repeats.over.maxrep,
                                maxrep),
            startsrev = lapply( rbam$startsrev,
                                remove.repeats.over.maxrep,
                                maxrep),
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

# Generate set of non-CpGs
noncpgSitesFromCpGset = function(cpgset, distance){
    noncpg = vector("list", length(cpgset));
    names(noncpg) = names(cpgset);
    for( i in seq_along(cpgset) ){ # i=1;
        pos = cpgset[[i]];
        difpos = diff(pos);
        keep = which( difpos >= (distance*2L) );
        newpos = (pos[keep+1L] + pos[keep]) %/% 2L;
        noncpg[[i]] = newpos;
    }
    return(noncpg);
}

# Find isolated CpGs among the given set of CpGs
isocpgSitesFromCpGset = function(cpgset, distance){
    isocpg = vector("list",length(cpgset));
    names(isocpg) = names(cpgset);
    for( i in seq_along(cpgset) ){
        distbig = diff(cpgset[[i]]) >= distance;
        isocpg[[i]] = cpgset[[i]][ which( c(
            distbig[1],
            distbig[-1] & distbig[-length(distbig)],
            distbig[length(distbig)]) ) ];
    }
    return(isocpg);
}

# Count reads away from all CpGs, forward looking reads
.count.nonCpG.reads.forward = function(starts, cpglocations, distance){
    ### count CpGs before the read
    ### count CpGs before and covered by the read
    ind = findInterval(c(starts-1L,starts+(distance-1L)), cpglocations);
    dim(ind)=c(length(starts),2);
    # cbind(ind, starts)
    return(c(sum(ind[,1] == ind[,2]),length(starts)));
}

# Count reads away from all CpGs, reverse looking reads
.count.nonCpG.reads.reverse = function(starts, cpglocations, distance){
    ### count CpGs left of read (+distance)
    ### count CpGs left of read start or at start
    ind = findInterval(c(starts-distance,starts), cpglocations);
    dim(ind)=c(length(starts),2);
    # cbind(ind, starts)
    return(c(sum(ind[,1] == ind[,2]),length(starts)));
}

# QC: Count reads away from all CpGs for an Rbam
bam.count.nonCpG.reads = function(rbam, cpgset, distance){
    result = c(nonCpGreads = 0,totalreads = 0);
    for( chr in names(cpgset) ){ # chr = names(cpgset)[1]
        frwstarts = rbam$startsfwd[[chr]];
        if( length(frwstarts)>0 )
            result = result + .count.nonCpG.reads.forward(
                starts = frwstarts, cpglocations = cpgset[[chr]], distance);
        revstarts = rbam$startsrev[[chr]];
        if( length(revstarts)>0 )
            result = result + .count.nonCpG.reads.reverse(
                starts = revstarts, cpglocations = cpgset[[chr]], distance);
    }
    rbam$qc$cnt.nonCpG.reads = result;
    class(rbam$qc$cnt.nonCpG.reads) = "qcNonCpGreads";
    return(rbam);
}

# Calculate distribution of distances to isolated CpGs, forward reads
.hist.isodist.forward = function(starts, cpglocations, distance){
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

# Calculate distribution of distances to isolated CpGs, reverse reads
.hist.isodist.reverse = function(starts, cpglocations, distance){
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

# QC: Calculate distribution of distances to isolated CpGs for an Rbam
bam.hist.isolated.distances = function(rbam, isocpgset, distance){
    result = 0;
    for( chr in names(isocpgset) ){ # chr = names(cpgset)[1]
        frwstarts = rbam$startsfwd[[chr]];
        if( length(frwstarts)>0 )
            result = result + .hist.isodist.forward(
                starts = frwstarts, cpglocations = isocpgset[[chr]], distance);
        revstarts = rbam$startsrev[[chr]];
        if( length(revstarts)>0 )
            result = result + .hist.isodist.reverse(
                starts = revstarts, cpglocations = isocpgset[[chr]], distance);
    }
    rbam$qc$hist.isolated.dist1 = result;
    class(rbam$qc$hist.isolated.dist1) = "qcIsoDist";
    return(rbam);
}

# QC: Calculate average coverage vs. CpG density
bam.coverage.by.density = function(rbam, cpgset, noncpgset,
                                   minfragmentsize, maxfragmentsize){

    fragdistr = c(
        rep(1, minfragmentsize-1),
        seq(1, 0, length.out = (maxfragmentsize-minfragmentsize)/1.5+1));
    fragdistr = fragdistr[fragdistr>0];

    if( is.null(noncpgset) ){
        noncpgset =
            noncpgSitesFromCpGset(cpgset = cpgset, distance = maxfragmentsize);
    }
    # sum(sapply(noncpgset,length))
    # newcpgset = noncpgset;
    # for( chr in seq_along(noncpgset) ){
    #     newcpgset[[chr]] = sort.int( c(cpgset[[chr]], noncpgset[[chr]]) );
    # }
    # rm(noncpgset);

    cpgdensity1 = calc.coverage(rbam = list(startsfwd = cpgset),
                                cpgset = cpgset,
                                fragdistr = fragdistr);
    cpgdensity2 = calc.coverage(rbam = list(startsrev = lapply(cpgset,`-`,1L)),
                                cpgset = cpgset,
                                fragdistr = fragdistr[-1]);
    cpgdensity = unlist(cpgdensity1, recursive = FALSE, use.names = FALSE) +
                 unlist(cpgdensity2, recursive = FALSE, use.names = FALSE);
    rm(cpgdensity1,cpgdensity2);

    cpgcoverage = calc.coverage(rbam, cpgset,    fragdistr);
    cpgcoverage = unlist(cpgcoverage, recursive = FALSE, use.names = FALSE);

    noncoverage = calc.coverage(rbam, noncpgset, fragdistr);
    noncoverage = unlist(noncoverage, recursive = FALSE, use.names = FALSE);

    # sqrtcover = sqrt(coverage);
    sqrtcpgdensity = sqrt(cpgdensity);
    rm(cpgdensity);

    axmax = ceiling(quantile(sqrtcpgdensity,0.99)*100)/100;

    # library(KernSmooth);
    z = locpoly(x = c(sqrtcpgdensity, double(length(noncoverage))),
                y = c(cpgcoverage, noncoverage),
                bandwidth = 0.5,
                gridsize = axmax*100+1,
                range.x = c(0,axmax));
    z$y[is.na(z$y)] = 0;

    rbam$qc$avg.coverage.by.density = z$y;
    class(rbam$qc$avg.coverage.by.density) = "qcCoverageByDensity";
    rbam$qc$avg.noncpg.coverage = mean(noncoverage);
    rbam$qc$avg.cpg.coverage = mean(cpgcoverage);

    return(rbam);
}

# QC: Fraction of reads on ChrX/Y
bam.chrXY.qc = function(rbam){
    strandfunX = function(st){c(length(st$chrX), sum(sapply(st,length)))};
    rbam$qc$chrX.count = strandfunX(rbam$startsfwd) +
                         strandfunX(rbam$startsfwd);
    class(rbam$qc$chrX.count) = "qcChrX"

    strandfunY = function(st){c(length(st$chrY), sum(sapply(st,length)))};
    rbam$qc$chrY.count =  strandfunY(rbam$startsfwd) +
                          strandfunY(rbam$startsfwd);
    class(rbam$qc$chrY.count) = "qcChrY"

    return(rbam);
}

# Caching environment
.ramwasEnv = new.env()

# Load an RDS file and cache it
# or load from cache
cachedRDSload = function(rdsfilename){
    if(is.null(rdsfilename))
        return(NULL);
    cachename = rdsfilename; #paste0(".ramwas.",rdsfilename);
    if( exists(x = cachename, envir = .ramwasEnv) ){
        # message("cachedRDSload: Using cache for: ", rdsfilename);
        return(base::get(x = cachename, envir = .ramwasEnv));
    } else {
        # message("cachedRDSload: Loading to cache: ", rdsfilename);
        data = readRDS(rdsfilename);
        base::assign(x = cachename, value = data, envir = .ramwasEnv);
        return(data);
    }
}

# CpG score calculation for a single chromosome
# Calling some C/C++ code
.calc.coverage.chr = function(startfrw, startrev, cpgs, fragdistr){
    maxfragmentsize = length(fragdistr);
    cover = double(length(cpgs));

    if(length(startfrw) > 0){
        # CpGs left of start
        ind1 = findInterval(cpgs - maxfragmentsize, startfrw);
        # CpGs left of start+250L
        ind2 = findInterval(cpgs,                   startfrw);
        # coverage of CpGs
        # which(ind2>ind1)
        # are covered by fragments
        # ind1[which(ind2>ind1)]+1 .. ind2[which(ind2>ind1)]
        .Call("cover_frw_c",
              startfrw, cpgs, fragdistr, ind1, ind2, cover, PACKAGE = "ramwas");
    }
    if(length(startrev) > 0){
        # CpGs left of start
        ind1 = findInterval(cpgs - 1L,                 startrev);
        # CpGs left of start+250L
        ind2 = findInterval(cpgs + maxfragmentsize-1L, startrev);
        # coverage of CpGs
        # which(ind2>ind1)
        # are covered by fragments
        # ind1[which(ind2>ind1)]+1 .. ind2[which(ind2>ind1)]
        .Call("cover_rev_c",
              startrev, cpgs, fragdistr, ind1, ind2, cover, PACKAGE = "ramwas");
    }
    return(cover);
}

# Calculate CpG scores for an Rbam
calc.coverage = function(rbam, cpgset, fragdistr){
    coveragelist = vector("list", length(cpgset));
    names(coveragelist) = names(cpgset);
    for( chr in names(coveragelist) ){ # chr = names(coveragelist)[1]
        coveragelist[[chr]] =
            .calc.coverage.chr(rbam$startsfwd[[chr]],
                               rbam$startsrev[[chr]],
                               cpgset[[chr]],
                               fragdistr);
    }
    return(coveragelist);
}

# Make QC plots for an Rbam
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

    if( !is.null(rbam$qc$hist.isolated.dist1) ){
        filename = paste0(param$dirqc,"/isolated_distance/id_",bamname,".pdf");
        dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
        pdf(filename);
        plot(rbam$qc$hist.isolated.dist1, samplename = bamname);
        dev.off();
        rm(filename);
    }
    if( !is.null(rbam$qc$avg.coverage.by.density) ){
        filename = paste0(param$dirqc,
                          "/coverage_by_density/cbd_",bamname,".pdf");
        dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
        pdf(filename);
        plot(rbam$qc$avg.coverage.by.density, samplename = bamname);
        dev.off();
        rm(filename);
    }
}

# test a phenotype against the data
# accounting for covariates
# cvrtqr - must be orthonormalized
# Supports categorical outcomes (text of factor)
testPhenotype = function(phenotype, data, cvrtqr){
    mycov = matrix(phenotype, nrow = 1);
    slice = data;
    cvqr0 = cvrtqr;

    if( any(is.na(mycov)) ){
        keep = which(colSums(is.na(mycov))==0);

        mycov = mycov[, keep, drop=FALSE];
        slice = slice[keep, , drop=FALSE];
        cvqr0 = cvqr0[, keep, drop=FALSE];
        cvqr0 = t( qr.Q(qr(t(cvqr0))) );
    }

    # mycov = as.character(round(mycov));
    if( is.character(mycov) || is.factor(mycov) ){
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
        if( cvsumsq2 <= cvsumsq1 * .Machine$double.eps*ncol(mycov) ){
            mycov[] = 0;
        } else {
            mycov = mycov / sqrt(rowSums(mycov^2));
        }
    }

    ###
    nVarTested = nrow(mycov);
    dfFull = ncol(cvqr0) - nrow(cvqr0) - nVarTested;
    if(nVarTested == 1){
        if(dfFull <= 0)
            return(list(correlation = 0,
                        tstat = 0,
                        pvalue = 1,
                        nVarTested = nVarTested,
                        dfFull = dfFull,
                        statname = ""));

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

        cor2tt = function(x){
            return( x * sqrt( dfFull / (1 - pmin(x^2,1))));
        }
        tt2pv = function(x){
            return( (pt(-abs(x),dfFull)*2));
        }
        tt = cor2tt(cr);
        pv = tt2pv(tt);

        ### Check:
        # c(tt[1], pv[1])
        # summary(lm( as.vector(covariate) ~ 0 + data[1,] +
        #         t(cvrtqr)))$coefficients[1,]
        return( list(correlation = cr,
                     tstat = tt,
                     pvalue = pv,
                     nVarTested = nVarTested,
                     dfFull = dfFull,
                     statname = "") );

    } else {
        if(dfFull <= 0)
            return( list(Rsquared = 0,
                         Fstat = 0,
                         pvalue = 1,
                         nVarTested = nVarTested,
                         dfFull = dfFull,
                         statname = paste0("-F_",nVarTested)) );

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
        rsq2F = function(x){
            return( x / (1 - pmin(x,1)) * (dfFull/nVarTested) );
        }
        F2pv = function(x){
            return( pf(x, nVarTested, dfFull, lower.tail = FALSE) );
        }
        ff = rsq2F(rsq);
        pv = F2pv(ff);

        ### Check:
        # c(ff[1], pv[1])
        # anova(lm( data[1,] ~ 0 + t(cvrtqr) +
        #     as.factor(as.vector(as.character(round(covariate))))))
        return( list(Rsquared = rsq,
                     Fstat = ff,
                     pvalue = pv,
                     nVarTested = nVarTested,
                     dfFull = dfFull,
                     statname = paste0("-F_",nVarTested)) );
    }
}

# Orthonormalize a set of covariates
orthonormalizeCovariates = function(cvrt){
    if(any(sapply(lapply(cvrt, is.na), any)))
        stop("Missing values are not allowed in the covariates")
    cvrtset = c(const = list(rep(1, nrow(cvrt))), cvrt);
    factorset = which(sapply(cvrtset, class) %in% c("character","factor"));
    for( ind in factorset ){ # ind = 3
        fctr = factor(cvrtset[[ind]]);
        cvrtset[[ind]] = model.matrix(~fctr)[,-1];
        rm(fctr);
    }
    cvrtmat = matrix(unlist(cvrtset), nrow(cvrt));
    cvrtqr = qr.Q(qr(cvrtmat));  ### tcrossprod(cvrtqr) - diag(nrow(cvrtqr))
    return(cvrtqr)
}

# find how samples in CpG score matrix
# match those in "covariates" parameter
# get the total number of CpGs along the way
.matchCovmatCovar = function( param ){
    cvsamples = param$covariates[[1]];

    fm = fm.open( paste0(param$dircoveragenorm, "/Coverage"), readonly = TRUE);
    fmsamples = rownames(fm);
    ncpgs = ncol(fm);
    close(fm);

    rowsubset = match(cvsamples, fmsamples, nomatch = 0L);
    if( any(rowsubset==0) )
        stop( paste("Unknown samples in covariate file:",
                    cvsamples[head(which(rowsubset==0))]) );

    if( length(cvsamples) == length(fmsamples) ){
        if( all(rowsubset == seq_along(rowsubset)) ){
            rowsubset = NULL;
        }
    }
    return(list(rowsubset = rowsubset, ncpgs = ncpgs));
}

# Get covariates + PCs matrix for analysis
# orthonormalized unless normalize == FALSE
.getCovariates = function(param, rowsubset, normalize = TRUE){
    cvrtqr = param$covariates[ param$modelcovariates ];
    ### Reading PCs, add as coveriates
    if( param$modelPCs > 0 ){
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
    if(normalize){
        rez = t(orthonormalizeCovariates(mwascvrtqr));
    } else {
        rez = t(mwascvrtqr); #t(cbind(rep(1, nrow(mwascvrtqr)),mwascvrtqr));
    }
    return(rez);
}

# Find best N p-values, in unsorted vector
findBestNpvs = function(pv, n){
    if(n < 1)
        return(which(pv <= n));

    pvthr = 10^((-100):0);
    fi = findInterval(pv, pvthr);
    tab = cumsum(tabulate(fi));
    upperfi = which(tab > n)[1];
    set1 = which(fi <= upperfi);
    cpgsetraw = set1[sort.list(pv[set1])[seq_len(n)]];
    cpgset = sort.int(cpgsetraw);
    return(cpgset);
}

# Standard BH p-value to q-value calculation
pvalue2qvalue = function(pv, n = length(pv)){
    ord = sort.list(pv);
    FDR = pv[ord] * n / seq_along(pv);
    FDR[length(FDR)] = min(FDR[length(FDR)], 1);
    FDR = rev(cummin(rev(FDR)));

    rez = double(length(pv));
    rez[ord] = FDR;
    return(rez)
}
