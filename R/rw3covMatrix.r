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
        .Call(
            "cover_frw_c",
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
        .Call(
            "cover_rev_c",
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
            .calc.coverage.chr(
                    startfrw = rbam$startsfwd[[chr]],
                    startrev = rbam$startsrev[[chr]],
                    cpgs = cpgset[[chr]],
                    fragdistr = fragdistr);
    }
    return(coveragelist);
}

# .readRDS = function(object){
#     # readBin(object, what = "raw", n = file.size(object));
#     return(readRDS(object));
# }

# calculate CpG score matrix for 1 sample
pipelineCoverage1Sample = function(colnum, param){

    cpgset = cachedRDSload(param$filecpgset);

    bams = param$bam2sample[[colnum]];

    if( param$maxrepeats == 0 ){
        coverage = NULL;
        for( j in seq_along(bams)){ # j=1
            rbam = readRDS(paste0(param$dirrbam, "/", bams[j], ".rbam.rds"));
            cov = calc.coverage(rbam = rbam,
                                cpgset = cpgset,
                                fragdistr = param$fragdistr)
            if(is.null(coverage)){
                coverage = cov;
            } else {
                for( i in seq_along(coverage) )
                    coverage[[i]] = coverage[[i]] + cov[[i]]
            }
            rm(cov);
        }
    } else {
        rbams = vector("list",length(bams));
        for( j in seq_along(bams)){ # j=1
            rbams[[j]] =
                readRDS(paste0( param$dirrbam,"/",bams[j],".rbam.rds"));
        }
        if(length(bams) > 1){
            rbam = list(startsfwd = list(), startsrev = list());
            for( i in seq_along(cpgset) ){ # i=1
                nm = names(cpgset)[i];

                fwd = lapply(rbams, function(x,y){x$startsfwd[[y]]}, nm);
                fwd = sort.int( unlist(fwd, use.names = FALSE) );
                rbam$startsfwd[[nm]] =
                    remove.repeats.over.maxrep(fwd, param$maxrepeats);
                rm(fwd);

                rev = lapply(rbams, function(x,y){x$startsrev[[y]]}, nm);
                rev = sort.int( unlist(rev, use.names = FALSE) );
                rbam$startsrev[[nm]] =
                    remove.repeats.over.maxrep(rev, param$maxrepeats);
                rm(rev);
            }
        } else {
            rbam = bam.removeRepeats( rbams[[1]], param$maxrepeats );
        }
        rm(rbams);
        coverage = calc.coverage(
                        rbam = rbam,
                        cpgset = cpgset,
                        fragdistr = param$fragdistr)
    }
    return(coverage);
}

# Job function for calculation of raw CpG score matrix
.ramwas3coverageJob = function(colnum, param, nslices){
    # library(ramwas);
    ld = param$dircoveragenorm;

    .log(ld, "%s, Process %06d, Start processing sample: %04d, %s",
        date(), Sys.getpid(), colnum, names(param$bam2sample)[colnum]);

    # Calculate coverage
    # coverage = ramwas:::pipelineCoverage1Sample(colnum, param);
    coverage = pipelineCoverage1Sample(colnum, param);
    coverage = unlist(coverage, use.names = FALSE);

    # Save it in file matrix slices
    start = 1;
    for( part in seq_len(nslices) ){ # part=1
        # message("colnum = ", colnum, " part = ", part);
        fmname = paste0(param$dirtemp1, "/RawCoverage_part", part);
        fm = fm.open(fmname, lockfile = param$lockfile);
        ntowrite = nrow(fm);
        fm$writeCols(colnum, coverage[start:(start+ntowrite-1)]);
        close(fm);
        start = start + ntowrite;
    }

    .log(ld, "%s, Process %06d,  Done processing sample: %04d, %s",
        date(), Sys.getpid(), colnum, names(param$bam2sample)[colnum]);

    return(NULL);
}

# Job function for filtering CpGs
.ramwas3transposeFilterJob = function(fmpart, param){
    ld = param$dircoveragenorm;
    # library(ramwas);
    
    # fmpart = 1
    # Open input file matrix slice
    filename = paste0(param$dirtemp1, "/RawCoverage_part", fmpart);
    fmraw = fm.open(
                filenamebase = filename,
                lockfile = param$lockfile2);

    .log(ld, "%s, Process %06d, Start transposing slice: %03d",
        date(), Sys.getpid(), fmpart);

    # Read the whole slice
    # mat = fmraw[]; # Not elegant
    mat = as.matrix(fmraw); # Failed with missing line in NAMESPACE

    # Create transposed+filtered output and the corresponding location files
    fmout = fm.create( 
                filenamebase = 
                    paste0(param$dirtemp2, "/TrCoverage_part", fmpart),
                nrow = ncol(mat),
                ncol = 0,
                size = param$doublesize,
                lockfile = param$lockfile2);
    fmpos = fm.create( 
                filenamebase = 
                    paste0(param$dirtemp2, "/TrCoverage_loc",  fmpart),
                nrow = 1,
                ncol = 0,
                type = "integer",
                lockfile = param$lockfile2);

    samplesums = rep(0, ncol(mat));

    ### Sliced loop
    step1 = max(floor(32*1024*1024 / 8 / ncol(mat)),1);
    mm = nrow(mat);
    nsteps = ceiling(mm/step1);
    for( part in seq_len(nsteps) ){ # part = 1
        # message(part, " of ", nsteps);
        fr = (part-1)*step1 + 1;
        to = min(part*step1, mm);

        subslice = mat[fr:to, , drop=FALSE];

        ### Filtering criteria
        cpgmean = rowMeans( subslice );
        cpgnonz = rowMeans( subslice>0 );
        keep = 
            (cpgmean >= param$minavgcpgcoverage) &
            (cpgnonz >= param$minnonzerosamples);
        if( !any(keep) )
            next;

        slloc = fr:to;

        if( !all(keep) ){
            keep = which(keep);
            subslice = subslice[keep, , drop=FALSE];
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

    fmss = fm.open( 
                filenamebase = paste0(param$dirtemp2, "/0_sample_sums"),
                lockfile = param$lockfile2);
    fmss[,fmpart] = samplesums;
    close(fmss);

    .log(ld, "%s, Process %06d,  Done transposing slice: %03d",
        date(), Sys.getpid(), fmpart);
    
    closeAndDeleteFiles(fmraw);
    return(NULL);
}

# Job function for sample normalization
.ramwas3normalizeJob = function(fmpart_offset, param, samplesums){
    # library(ramwas);
    ld = param$dircoveragenorm;
    
    scale = as.vector(samplesums) / mean(samplesums);

    # Open output file
    filename = paste0(param$dircoveragenorm, "/Coverage");
    fm = fm.open(filename, lockfile = param$lockfile2);
    
    .log(ld, "%s, Process %06d, Start normalizing slice: %03d",
        date(), Sys.getpid(), fmpart_offset[1]);

    # Read filtered+transposed 
    filename = paste0(param$dirtemp2, "/TrCoverage_part", fmpart_offset[1]);
    fmin = fm.open(filename, lockfile = param$lockfile1);

    step1 = max(floor(32*1024*1024 / 8 / nrow(fmin)),1);
    mm = ncol(fmin);
    nsteps = ceiling(mm/step1);
    for( part in seq_len(nsteps) ){ # part = 1L
        # message(part, " of ", nsteps);
        fr = (part-1)*step1 + 1L;
        to = min(part*step1, mm);

        fm$writeCols( 
                start = fmpart_offset[2]+fr,
                value = fmin[,fr:to] / scale);
    }
    rm(part, step1, mm, nsteps, fr, to);

    # Cleanup
    closeAndDeleteFiles(fmin);
    close(fm);

    .log(ld, "%s, Process %06d,  Done normalizing slice: %03d",
        date(), Sys.getpid(), fmpart_offset[1]);
    return(NULL);
}

# Step 3 of RaMWAS
ramwas3normalizedCoverage = function( param ){
    param = parameterPreprocess(param);
    ld = param$dircoveragenorm;
    
    # Fragment size estimate
    if( param$minfragmentsize < param$maxfragmentsize ){
        filename = paste0(param$dirfilter, "/Fragment_size_distribution.txt");
        if( !file.exists(filename) )
            stop(
                "Fragment size distribution estimate not found: ", filename,
                "\nRun ramwas2collectqc() function first")
        param$fragdistr = as.double(readLines(con = filename));
        rm(filename);
    } else {
        param$fragdistr = rep(1, param$maxfragmentsize);
    }
    
    # Is CpG set defined.
    if(is.null(param$filecpgset))
        stop("CpG set is not defined. See \"filecpgset\" parameter.")
    
    
    # Are samples in covariates are in bam2sample?
    if( !is.null(param$covariates) ){
        badset = !(names(param$bam2sample) %in% param$covariates[[1]]);
        if( any(badset) )
            stop("Covariate file has samples not present in \"bam2sample\" ",
                "parameter:\n ",
                paste0(head(param$covariates[[1]][badset]), collapse = "\n "));
        
        message("Using the ", length(param$covariates[[1]]),
                " samples in the covariate file.")
        param$bam2sample = param$bam2sample[param$covariates[[1]]];
    }
    # Check is all rbams are in place
    {
        message("Checking if all required Rbam files present");
        bams = unlist(param$bam2sample);
        for( bname in bams){
            filename = paste0( param$dirrbam, "/", bname, ".rbam.rds");
            if( !file.exists(filename) ){
                stop(paste0("Rbam file from bam2sample not found: ", filename));
            }
        }
        rm(bams, bname, filename);
        message("All required Rbam files present are present.");
    }
    
    dir.create(param$dircoveragenorm, showWarnings = FALSE, recursive = TRUE);
    dir.create(param$dirtemp, showWarnings = FALSE, recursive = TRUE);


    parameterDump(dir = param$dircoveragenorm, param = param,
        toplines = c(   "dircoveragenorm", "dirtemp", "dirrbam",
                        "dirtemp","dirtemp1","dirtemp2",
                        "filebam2sample", "bam2sample",
                        "maxrepeats",
                        "minavgcpgcoverage", "minnonzerosamples",
                        "filecpgset",
                        "buffersize", "doublesize",
                        "cputhreads", "diskthreads"));

    .log(ld, "%s, Start ramwas3normalizedCoverage()", date(), append = FALSE);
    
    ### data dimensions
    cpgset = cachedRDSload(param$filecpgset);
    ncpgs = sum(sapply(cpgset, length));
    nsamples = length(param$bam2sample);


    # Two temporary directories for faster processing
    if(is.null(param$dirtemp1))
        param$dirtemp1 = param$dirtemp;
    if(is.null(param$dirtemp2))
        param$dirtemp2 = param$dirtemp;
    
    stopifnot(dir.exists(param$dirtemp1));
    stopifnot(dir.exists(param$dirtemp2));
    
    ### Create raw coverage matrix slices
    {
        # library(filematrix)
        # Sliced loop
        kbblock = (128*1024)/8;
        step1 = max(floor(param$buffersize / (8 * nsamples)/kbblock),1)*kbblock;
        mm = ncpgs;
        nslices = ceiling(mm/step1);
        .log(ld, "%s, Creating %d file matrices for raw scores at: %s",
            date(), nslices, param$dirtemp1);
        
        for( part in seq_len(nslices) ){ # part = 1
            # cat("Creating raw  matrix slices", part, "of", nslices, "\n");
            fr = (part-1)*step1 + 1;
            to = min(part*step1, mm);
            fmname = paste0(param$dirtemp1, "/RawCoverage_part", part);
            fm = fm.create( 
                        filenamebase = fmname,
                        nrow = to-fr+1,
                        ncol = nsamples,
                        size = param$doublesize)
            close(fm);
        }
        rm(part, step1, mm, fr, to, fmname);
    } # nslices

    ### Fill in the raw coverage files
    {
        # library(parallel)
        nthreads = min(param$cputhreads, nsamples);
        .log(ld, "%s, Calculating raw coverage", date());
        if(param$usefilelock) param$lockfile = tempfile();
        if( nthreads > 1 ){
            cl = makeCluster(nthreads);
            on.exit({
                stopCluster(cl);
                .file.remove(param$lockfile);
            });
            logfun = .logErrors(ld, .ramwas3coverageJob);
            # clusterExport(cl, c(".log","ld",".ramwas3coverageJob"));
            # logfun = ramwas:::.ramwas3coverageJob
            z = clusterApplyLB(
                        cl = cl,
                        x = seq_len(nsamples),
                        fun = logfun,
                        param = param,
                        nslices = nslices);
            tmp = sys.on.exit();
            eval(tmp);
            rm(tmp);
            on.exit();
        } else {
            z = vector("list", nsamples);
            names(z) = names(param$bam2sample);
            for(i in seq_along(param$bam2sample)){ # i=1
                z[[i]] = .ramwas3coverageJob(  
                        colnum = i,
                        param = param,
                        nslices = nslices);
            }
        }
        .showErrors(z);
        .log(ld, "%s, Done calculating raw coverage", date());
    }
    
    ### Transpose the slices, filter by average and fraction of non-zeroes
    {
        .log(ld, "%s, Transposing score matrices and filtering CpGs by score",
            date());

        fm = fm.create( 
                    filenamebase = paste0(param$dirtemp2,"/0_sample_sums"),
                    nrow = nsamples,
                    ncol = nslices);
        close(fm);

        logfun = .logErrors(ld, .ramwas3transposeFilterJob);
        nthreads = min(param$diskthreads, nslices);
        if( nthreads > 1 ){
            if(param$usefilelock) param$lockfile2 = tempfile();
            # library(parallel);
            cl = makeCluster(nthreads);
            on.exit({
                stopCluster(cl);
                .file.remove(param$lockfile2);
            });
            z = clusterApplyLB(
                        cl = cl,
                        x = seq_len(nslices),
                        fun = logfun,
                        param = param);
            tmp = sys.on.exit();
            eval(tmp);
            rm(tmp);
            on.exit();
        } else {
            z = vector("list", nslices);
            for( fmpart in seq_len(nslices) ){ # fmpart = 1
                z[[fmpart]] = .ramwas3transposeFilterJob(fmpart, param);
            }
        }
        .showErrors(z);
        .log(ld, "%s, Done transposing score matrices, filtering CpGs", date());
    }

    ### Prepare CpG set for filtered CpGs
    {
        .log(ld, "%s, Saving locations for CpGs which passed the filter", 
            date());

        chr = rep(seq_along(cpgset), sapply(cpgset, length));
        pos = unlist(cpgset, recursive = FALSE, use.names = FALSE);
        keep = logical(ncpgs);
        
        kbblock = (128*1024)/8;
        step1 = max(floor(param$buffersize / (8*nsamples)/kbblock),1)*kbblock;
        mm = ncpgs;
        nsteps = ceiling(mm/step1);
        keeplist = vector("list",nsteps);
        sliceoffsets = integer(1+nsteps);
        for( part in seq_len(nsteps) ){ # part = 1
            # cat( part, "of", nsteps, "\n");
            fr = (part-1)*step1 + 1;
            to = min(part*step1, mm);

            indx = fm.load(paste0(param$dirtemp2, "/TrCoverage_loc", part));
            indx = as.vector(indx)
            keep[fr:to][indx] = TRUE;
            sliceoffsets[part+1] = sliceoffsets[part] + length(indx);
            
            # Cleanup
            fmlc = fm.open(paste0(param$dirtemp2, "/TrCoverage_loc", part));
            closeAndDeleteFiles(fmlc);
        }
        rm(part, step1, mm, nsteps, fr, to, kbblock, indx);

        cpgslocmat = cbind( 
                    chr = chr[keep],
                    position = pos[keep]);

        fm = fm.create.from.matrix(
            filenamebase = paste0(param$dircoveragenorm, "/CpG_locations"),
            mat = cpgslocmat);
        close(fm);
        
        writeLines(
            con = paste0(param$dircoveragenorm, "/CpG_chromosome_names.txt"),
            text = names(cpgset));
        rm(chr, pos, cpgslocmat);
    } # /CpG_locations, sliceoffsets

    ### Sample sums
    {
        .log(ld, "%s, Gathering sample sums from %d slices", date(), nslices);

        mat = fm.load(paste0(param$dirtemp2, "/0_sample_sums"));
        samplesums = rowSums(mat);
        rm(mat);
        fm = fm.create.from.matrix(
            filenamebase = paste0(param$dircoveragenorm, "/raw_sample_sums"),
            mat = samplesums);
        close(fm);
        
        # Cleanup
        fm = fm.open(paste0(param$dirtemp2, "/0_sample_sums"));
        closeAndDeleteFiles(fm);
    }

    ### Normalize and combine in one matrix
    {
        .log(ld, "%s, Normalizing coverage and saving in one matrix", date());

        fmpart_offset_list = mat2cols(
            rbind(
                seq_len(nslices),
                sliceoffsets[-length(sliceoffsets)]));

        ### Create big matrix for normalized coverage
        fm = fm.create(
                    filenamebase = paste0(param$dircoveragenorm, "/Coverage"),
                    nrow = nsamples,
                    ncol = tail(sliceoffsets,1),
                    size = param$doublesize);
        rownames(fm) = names(param$bam2sample);
        close(fm);

        
        ### normalize and fill in
        nthreads = min(param$diskthreads, nslices);
        if( nthreads > 1 ){
            if(param$usefilelock) param$lockfile1 = tempfile();
            if(param$usefilelock) param$lockfile2 = tempfile();
            # library(parallel);
            cl = makeCluster(nthreads);
            on.exit({
                stopCluster(cl);
                .file.remove(param$lockfile1);
                .file.remove(param$lockfile2);
            });
            logfun = .logErrors(ld, .ramwas3normalizeJob);
            z = clusterApplyLB(
                        cl = cl,
                        x = fmpart_offset_list,
                        fun = logfun,
                        param = param,
                        samplesums = samplesums);
            tmp = sys.on.exit();
            eval(tmp);
            rm(tmp);
            on.exit();
        } else {
            for( i in seq_len(nslices) ){ # i = 1
                z = vector("list", nslices);
                z[[i]] = .ramwas3normalizeJob( 
                        fmpart_offset = fmpart_offset_list[[i]],
                        param = param,
                        samplesums = samplesums);
            }
        }
        .showErrors(z);
        .log(ld, "%s, Done Normalizing coverage matrix.", date());
    }
        
    .log(ld, "%s, Done ramwas3normalizedCoverage()", date());
}
