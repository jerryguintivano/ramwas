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

.readRDS = function(object){
    # readBin(object, what = 'raw', n = file.size(object));
    return(readRDS(object));
}

# calculate CpG score matrix for 1 sample
pipelineCoverage1Sample = function(colnum, param){

    cpgset = cachedRDSload(param$filecpgset);

    bams = param$bam2sample[[colnum]];

    if( param$maxrepeats == 0 ){
        coverage = NULL;
        for( j in seq_along(bams)){ # j=1
            rbam = .readRDS(paste0(param$dirrbam, "/", bams[j], ".rbam.rds"));
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
                .readRDS(paste0( param$dirrbam,"/",bams[j],".rbam.rds"));
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
        coverage = calc.coverage(rbam = rbam,
                                 cpgset = cpgset,
                                 fragdistr = param$fragdistr)
    }
    return(coverage);
}

# Job function for calculation of raw CpG score matrix
.ramwas3coverageJob = function(colnum, param, nslices){
    # library(ramwas);
    coverage = pipelineCoverage1Sample(colnum, param);
    coverage = unlist(coverage, use.names = FALSE);

    start = 1;
    for( part in seq_len(nslices) ){ # part=1
        message("colnum = ", colnum, " part = ", part);
        fmname = paste0(param$dirtemp, "/RawCoverage_part", part);
        fm = fm.open(fmname, lockfile = param$lockfile);
        ntowrite = nrow(fm);
        fm$writeCols(colnum, coverage[start:(start+ntowrite-1)]);
        close(fm);
        start = start + ntowrite;
    }

    fmname = paste0(param$dirtemp, "/RawCoverage_part", 1);
    fm = fm.open(fmname, lockfile = param$lockfile);
    fm$filelock$lockedrun( {
        cat(file = paste0(param$dircoveragenorm,"/Log.txt"),
             date(), ", Process ", Sys.getpid(),
             ", Processing sample ", colnum, " ",
             names(param$bam2sample)[colnum], "\n",
             sep = "", append = TRUE);
    });
    close(fm);

    return("OK");
}

# Job function for filtering CpGs
.ramwas3transposeFilterJob = function(fmpart, param){
    # library(filematrix);
    # fmpart = 1
    filename = paste0(param$dirtemp, "/RawCoverage_part", fmpart);
    if( !file.exists(paste0(filename, ".bmat")) ||
        !file.exists(paste0(filename, ".desc.txt")) )
        return(paste0("Raw coverage slice filematrix not found: ", filename));
    fmraw = fm.open(filename, lockfile = param$lockfile2);
    mat = fmraw[];
    # mat = as.matrix(fmraw);

    fmout = fm.create( paste0(param$dirtemp,"/TrCoverage_part",fmpart),
                       nrow = ncol(mat),
                       ncol = 0,
                       size = param$doublesize,
                       lockfile = param$lockfile2);
    fmpos = fm.create( paste0(param$dirtemp,"/TrCoverage_loc",fmpart),
                       nrow = 1,
                       ncol = 0,
                       type = "integer",
                       lockfile = param$lockfile2);

    samplesums = rep(0, ncol(mat));

    ### Sliced loop
    step1 = max(floor(32*1024*1024 / 8 / ncol(mat)),1);
    mm = nrow(mat);
    nsteps = ceiling(mm/step1);
    for( part in 1:nsteps ){ # part = 1
        message(part, " of ", nsteps);
        fr = (part-1)*step1 + 1;
        to = min(part*step1, mm);

        subslice = mat[fr:to,,drop=FALSE];

        ### Filtering criteria
        cpgmean = rowMeans( subslice );
        cpgnonz = rowMeans( subslice>0 );
        keep = (cpgmean >= param$minavgcpgcoverage) &
               (cpgnonz >= param$minnonzerosamples);
        if( !any(keep) )
            next;

        slloc = fr:to;

        if( !all(keep) ){
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

    fmss = fm.open( paste0(param$dirtemp,"/0_sample_sums"),
                    lockfile = param$lockfile2);
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

# Job function for sample normalization
.ramwas3normalizeJob = function(fmpart_offset, param, samplesums){
    # fmpart_offset = fmpart_offset_list[[2]]
    scale = as.vector(samplesums) / mean(samplesums);

    # library(filematrix);

    filename = paste0(param$dirtemp, "/TrCoverage_part", fmpart_offset[1]);
    mat = fm.load(filename, param$lockfile1);
    mat = mat / scale;

    filename = paste0(param$dircoveragenorm, "/Coverage");
    fm = fm.open(filename, lockfile = param$lockfile2);
    fm$writeCols( start = fmpart_offset[2]+1L, mat);
    close(fm);

    fm$filelock$lockedrun( {
        cat(file = paste0(param$dircoveragenorm,"/Log.txt"),
             date(), ", Process ", Sys.getpid(),
             ", Processing slice ", fmpart_offset[1], "\n",
             sep = "", append = TRUE);
    });

    rm(mat);
    gc();
    return("OK.");
}

# Step 3 of RaMWAS
ramwas3normalizedCoverage = function( param ){
    # Prepare
    param = parameterPreprocess(param);
    if( param$minfragmentsize < param$maxfragmentsize ){
        filename = paste0(param$dirfilter,"/Fragment_size_distribution.txt");
        param$fragdistr = as.double( readLines(con = filename));
        rm(filename);
    } else {
        param$fragdistr = rep(1, param$maxfragmentsize);
    }
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
        for( bname in bams){
            filename = paste0( param$dirrbam, "/", bname);
            if( file.exists(filename) ){
                stop(paste0("Rbam file from bam2sample not found: ", filename));
            }
        }
        rm(bams, bname, filename);
    }

    ### Create raw coverage matrix slices
    {
        # library(filematrix)
        # Sliced loop
        kbblock = (128*1024)/8;
        step1 = max(floor(param$buffersize / (8 * nsamples)/kbblock),1)*kbblock;
        mm = ncpgs;
        nslices = ceiling(mm/step1);
        message("Creating ", nslices, " file matrices for raw scores at: ",
                param$dirtemp);
        for( part in 1:nslices ){ # part = 1
            # cat("Creating raw  matrix slices", part, "of", nslices, "\n");
            fr = (part-1)*step1 + 1;
            to = min(part*step1, mm);
            fmname = paste0(param$dirtemp,"/RawCoverage_part",part);
            fm = fm.create(fmname,
                           nrow = to-fr+1,
                           ncol = nsamples,
                           size = param$doublesize)
            close(fm);
        }
        rm(part, step1, mm, fr, to, fmname);
    } # nslices

    ### Fill in the raw coverage files
    {
        message("Calculating and saving raw coverage");
        if(param$usefilelock) param$lockfile = tempfile();
        cat(file = paste0(param$dircoveragenorm,"/Log.txt"),
             date(), ", Calculating raw coverage.", "\n",
             sep = "", append = FALSE);
        # library(parallel)
        nthreads = min(param$cputhreads, nsamples);
        if( nthreads > 1 ){
            cl = makeCluster(nthreads);
            on.exit({stopCluster(cl);});
            z = clusterApplyLB(cl,
                               seq_len(nsamples),
                               .ramwas3coverageJob,
                               param = param,
                               nslices = nslices);
            eval(sys.on.exit());
            on.exit();
        } else {
            z = character(nsamples);
            names(z) = names(param$bam2sample);
            for(i in seq_along(param$bam2sample)){ # i=1
                z[i] = .ramwas3coverageJob(colnum = i,
                                           param = param,
                                           nslices = nslices);
                message(i, z[i]);
            }
        }
        cat(file = paste0(param$dircoveragenorm,"/Log.txt"),
             date(), ", Done calculating raw coverage.", "\n",
             sep = "", append = TRUE);
        .file.remove(param$lockfile);
    }

    ### Transpose the slices, filter by average and fraction of non-zeroes
    {
        message("Transposing score matrices and filtering CpGs by score");

        fm = fm.create( paste0(param$dirtemp,"/0_sample_sums"),
                        nrow = nsamples,
                        ncol = nslices);
        close(fm);

        cat(file = paste0(param$dircoveragenorm,"/Log.txt"),
             date(), ", Transposing coverage matrix, filtering CpGs.", "\n",
             sep = "", append = TRUE);
        nthreads = min(param$diskthreads, nslices);
        if( nthreads > 1 ){
            if(param$usefilelock) param$lockfile2 = tempfile();
            # library(parallel);
            cl = makeCluster(nthreads);
            on.exit({
                stopCluster(cl);
                .file.remove(param$lockfile2);});
            z = clusterApplyLB(cl,
                               1:nslices,
                               .ramwas3transposeFilterJob,
                               param = param);
            eval(sys.on.exit());
            on.exit();
        } else {
            for( fmpart in seq_len(nslices) ){ # fmpart = 5
                .ramwas3transposeFilterJob(fmpart, param);
            }
        }
        cat(file = paste0(param$dircoveragenorm,"/Log.txt"),
             date(), ", Done transposing coverage matrix, filtering CpGs.",
             "\n", sep = "", append = TRUE);
    }

    ### Prepare CpG set for filtered CpGs
    {
        message("Saving locations for CpGs which passed the filter");

        cpgsloc1e9 = cpgset;
        for( i in seq_along(cpgsloc1e9) ){
            cpgsloc1e9[[i]] = cpgset[[i]] + i*1e9;
        }
        cpgsloc1e9 = unlist(cpgsloc1e9, recursive = FALSE, use.names = FALSE);

        kbblock = (128*1024)/8;
        step1 = max(floor(param$buffersize / (8 * nsamples)/kbblock),1)*kbblock;
        mm = ncpgs;
        nsteps = ceiling(mm/step1);
        cpgsloclist = vector("list",nsteps);
        for( part in 1:nsteps ){ # part = 1
            # cat( part, "of", nsteps, "\n");
            fr = (part-1)*step1 + 1;
            to = min(part*step1, mm);

            indx = fm.load( paste0(param$dirtemp,"/TrCoverage_loc",part) );
            indx = as.vector(indx)
            cpgsloclist[[part]] = cpgsloc1e9[fr:to][indx];
        }
        rm(part, step1, mm, nsteps, fr, to, kbblock, indx);
        sliceoffsets = c(0L, cumsum(sapply(cpgsloclist, length)));

        cpgslocvec = unlist(cpgsloclist, use.names = FALSE);
        cpgslocmat = cbind( chr = as.integer(cpgslocvec %/% 1e9),
                            position = as.integer(cpgslocvec %% 1e9));

        fm = fm.create.from.matrix(
            filenamebase = paste0(param$dircoveragenorm, "/CpG_locations"),
            mat = cpgslocmat);
        close(fm);
        writeLines(
            con = paste0(param$dircoveragenorm, "/CpG_chromosome_names.txt"),
            text = names(cpgset));
        rm(cpgsloc1e9, cpgsloclist, cpgslocvec, cpgslocmat);
    } # /CpG_locations, sliceoffsets

    ### Sample sums
    {
        message("Gathering sample sums from ", nslices, " slices");

        mat = fm.load( paste0(param$dirtemp, "/0_sample_sums") );
        samplesums = rowSums(mat);
        rm(mat);
        fm = fm.create.from.matrix(
            filenamebase = paste0(param$dircoveragenorm, "/raw_sample_sums"),
            mat = samplesums);
        close(fm);
    }

    ### Normalize and combine in one matrix
    {
        message("Normalizing coverage and saving in one matrix");

        fmpart_offset_list = as.list(data.frame(rbind(
            seq_len(nslices),
            sliceoffsets[-length(sliceoffsets)])));

        ### Create big matrix for normalized coverage
        fm = fm.create(paste0(param$dircoveragenorm, "/Coverage"),
                       nrow = nsamples,
                       ncol = tail(sliceoffsets,1),
                       size = param$doublesize);
        rownames(fm) = names(param$bam2sample);
        close(fm);

        cat(file = paste0(param$dircoveragenorm,"/Log.txt"),
             date(), ", Normalizing coverage matrix.", "\n",
             sep = "", append = TRUE);
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
                .file.remove(param$lockfile2);});
            z = clusterApplyLB(cl,
                               fmpart_offset_list,
                               .ramwas3normalizeJob,
                               param = param,
                               samplesums = samplesums);
            eval(sys.on.exit());
            on.exit();
        } else {
            for( fmpart in seq_len(nslices) ){ # fmpart = 5
                .ramwas3normalizeJob( fmpart_offset_list[[fmpart]],
                                      param,
                                      samplesums);
            }
        }
        cat(file = paste0(param$dircoveragenorm,"/Log.txt"),
             date(), ", Done normalizing coverage matrix.", "\n",
             sep = "", append = TRUE);

    }

    ### Cleanup
    {
        message("Removing temporary files");
        for( part in 1:nslices ){
            fm = fm.open( paste0(param$dirtemp,"/TrCoverage_loc",part) );
            closeAndDeleteFiles(fm);
            fm = fm.open( paste0(param$dirtemp,"/TrCoverage_part",part) );
            closeAndDeleteFiles(fm);
        }
        fm = fm.open( paste0(param$dirtemp,"/0_sample_sums") );
        closeAndDeleteFiles(fm);
    }
}
