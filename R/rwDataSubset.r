# Subset a coverage matrix and locations set.
subsetCoverageDirByLocation = function(x, chr, start, targetdir){
    if( is.list(x) ){
        param = parameterPreprocess(x);
        dircov = param$dircoveragenorm;
    } else {
        dircov = x;
    }
    chrnames = readLines(paste0(dircov, "/CpG_chromosome_names.txt"));
    if( is.factor(chr) )
        chr = as.character(chr);
    if( is.character(chr) ){
        chrn = match(chr, chrnames, nomatch = 0L);
    } else {
        chrn = chr;
    }
    
    locs = getLocations(dircov);
    maxmult = (max(locs$end)+2) * 2;

    locs$chr = as.integer(locs$chr);
    mch = match(
            x = chrn*maxmult + start,
            table = locs$chr * maxmult + locs$start,
            nomatch = 0L);
    if( any(mch == 0L) )
        stop("Some locations not found");

    mch = sort(mch);

    dir.create(targetdir, showWarnings = FALSE)

    file.copy(
            from = paste0(dircov,   "/CpG_chromosome_names.txt"),
            to =   paste0(targetdir, "/CpG_chromosome_names.txt"));
    
    fm = fm.create.from.matrix(
            filenamebase = paste0(targetdir, "/CpG_locations"),
            mat = as.matrix(locs[mch,]));
    close(fm);

    fi = fm.open(
                filenamebase = paste0(dircov, "/Coverage"),
                readonly = TRUE)
    fo = fm.create(
                filenamebase = paste0(targetdir, "/Coverage"),
                nrow = nrow(fi),
                ncol = length(mch),
                type = fi$type,
                size = fi$size);

    step1 = ceiling( 128*1024*1024 / nrow(fi) / 8);
    runto = length(mch);
    nsteps = ceiling(runto/step1);
    for( part in seq_len(nsteps) ){ # part = 1
        message("Loop filling coverage matrix, step ", part, " of ", nsteps);
        fr = (part-1)*step1 + 1;
        to = min(part*step1, runto);
        fo[, fr:to] = fi[, mch[fr:to]];
    }
    rm(part, step1, runto, nsteps, fr, to);

    if( !is.null(rownames(fi)))
        rownames(fo) = rownames(fi);
    if( !is.null(colnames(fi)))
        colnames(fo) = colnames(fi)[mch];

    close(fo);
    close(fi);
    return(invisible(NULL));
}

