# Subset a coverage matrix and locations set.
subsetCoverageDirByLocation = function(x, chr, position, targetdir){
    if( is.list(x) ){
        param = parameterPreprocess(x);
        dircov = param$dircoveragenorm;
    } else {
        dircov = x;
    }
    chrnames = readLines(paste0(dircov,"/CpG_chromosome_names.txt"));
    if( is.factor(chr) )
        chr = as.character(chr);
    if( is.character(chr) ){
        chrn = match(chr, chrnames, nomatch = 0L);
    } else {
        chrn = chr;
    }
    locmat = fm.load(paste0(dircov,"/CpG_locations"));
    maxmult = (max(locmat)+2)*2;

    mch = match(chrn*maxmult + position, locmat %*% c(maxmult,1), nomatch = 0L);
    if( any(mch == 0L) )
        stop("Some locations not found");

    mch = sort(mch);

    dir.create(targetdir, showWarnings = FALSE)
    setwd(targetdir)

    file.copy(
            from = paste0(dircov, "/CpG_chromosome_names.txt"),
            to = "CpG_chromosome_names.txt");
    fm = fm.create.from.matrix("CpG_locations", locmat[mch,]);
    close(fm);

    fi = fm.open(paste0(dircov, "/Coverage"), readonly = TRUE)
    fo = fm.create(
                filenamebase = "Coverage",
                nrow = nrow(fi),
                ncol = length(mch),
                type = fi$type,
                size = fi$size);

    step1 = 1024;
    runto = length(mch);
    nsteps = ceiling(runto/step1);
    for( part in seq_len(nsteps) ){ # part = 1
        message("Loop filling coverage matrix, step ", part, " of ", nsteps);
        fr = (part-1)*step1 + 1;
        to = min(part*step1, runto);
        fo[,fr:to] = fi[,mch[fr:to]];
    }
    rm(part, step1, runto, nsteps, fr, to);

    # all( fi[,mch] == as.matrix(fo) )

    if( !is.null(rownames(fi)))
        rownames(fo) = rownames(fi);
    if( !is.null(colnames(fi)))
        colnames(fo) = colnames(fi)[mch];

    close(fo);
    close(fi);
}

