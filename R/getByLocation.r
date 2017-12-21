# Extract the locations for the project
getLocations = function(x){
    if( is.list(x) ){
        param = parameterPreprocess(x);
        dircov = param$dircoveragenorm;
    } else {
        dircov = x;
    }
    
    chrfile = paste0(dircov, "/CpG_chromosome_names.txt");
    matfile = paste0(dircov, "/CpG_locations.bmat");
    
    if( !file.exists(chrfile) || !file.exists(matfile) )
        return(NULL);
    
    chrnames = readLines(chrfile);
    locmat = fm.load(paste0(dircov, "/CpG_locations"));
    chr = as.integer(locmat[,1]);
    levels(chr) = chrnames;
    class(chr) = "factor";
    locations = data.frame(
                    chr = chr,
                    start = locmat[,2]);
    if( ncol(locmat) >= 3 ){
        locations$end = locmat[,3];
    } else {
        locations$end = locations$start + 1L;
    }
    return(locations);
}

getMWAS = function(x){
    if( is.list(x) ){
        param = parameterPreprocess(x);
        dirmwas = param$dirmwas;
    } else {
        dirmwas = x;
        dircov = paste0(x, "/../..");
    }
    mwas = fm.load(filenamebase = paste0(dirmwas, "/Stats_and_pvalues"));
    return(data.frame(mwas, check.names = FALSE));
}

getMWASandLocations = function(x){
    if( is.list(x) ){
        param = parameterPreprocess(x);
        dirmwas = param$dirmwas;
        dircov = param$dircoveragenorm;
    } else {
        dirmwas = x;
        dircov = paste0(x, "/../..");
    }
    locations = getLocations(dircov);
    mwas = getMWAS(dirmwas);
    result = data.frame(
                locations,
                mwas,
                check.names = FALSE);
    return(result);
}

# Get MWAS results by locations
getMWASrange = function(x, chr, start, end){
    mwas = getMWASandLocations(x);
    chrnames = levels(mwas$chr);
    if( is.factor(chr) )
        chr = as.character(chr);
    if( is.character(chr) ){
        chrn = match(chr, chrnames, nomatch = 0L);
    } else {
        chrn = chr;
    }
    
    keep =  as.integer(mwas$chr) == chrn &
            (mwas$end >= start) & 
            (mwas$start <= end);
    
    return(mwas[keep,]);
}

# Get MWAS results by locations
getDataByLocation = function(x, chr, start, end){
    if( is.list(x) ){
        param = parameterPreprocess(x);
        dircov = param$dircoveragenorm;
    } else {
        dircov = x;
    }

    locs = getLocations(dircov);
    chrnames = levels(locs$chr);
    if( is.factor(chr) )
        chr = as.character(chr);
    if( is.character(chr) ){
        chrn = match(chr, chrnames, nomatch = 0L);
    } else {
        chrn = chr;
    }

    keep =  as.integer(locs$chr) == chrn &
            (locs$end >= start) & 
            (locs$start <= end);
    
    fm = fm.open(paste0(dircov, "/Coverage"));
    mat = fm[, keep];
    rownames(mat) = rownames(fm);
    close(fm);

    rez = list(
        locations = locs[keep,],
        matrix = mat);

    return(rez);
}

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

madeBED = function(x, filename){
    mwas = getMWASandLocations(x)
    bed = data.frame(
        chrom = mwas$chr,
        chromStart  = mwas$start - 1,
        chromEnd =  mwas$end,
        name = ".",
        score = mwas$`p-value`
    )
    write.table(
        file = filename,
        x = bed,
        quote = FALSE,
        sep = "\t",
        col.names = FALSE,
        row.names = FALSE);
    return(invisible(bed));
}

madeBEDrange = function(x, filename, chr, start, end){
    mwas = getMWASandLocations(x);
    chrnames = levels(mwas$chr);
    if( is.factor(chr) )
        chr = as.character(chr);
    if( is.character(chr) ){
        chrn = match(chr, chrnames, nomatch = 0L);
    } else {
        chrn = chr;
    }
    
    keep =  which(
                as.integer(mwas$chr) == chrn &
                (mwas$end >= start) & 
                (mwas$start <= end));
    
    bed = data.frame(
        chrom = mwas$chr[keep],
        chromStart  = mwas$start[keep] - 1,
        chromEnd =  mwas$end[keep],
        name = ".",
        score = mwas$`p-value`[keep]
    )
    write.table(
        file = filename,
        x = bed,
        quote = FALSE,
        sep = "\t",
        col.names = FALSE,
        row.names = FALSE);
    return(invisible(bed));
}

madeBEDgraph = function(x, filename){
    mwas = getMWASandLocations(x)
    bed = data.frame(
        chrom = mwas$chr,
        chromStart  = mwas$start - 1,
        chromEnd =  mwas$end,
        score = mwas$`p-value`
    )
    write.table(
        file = filename,
        x = bed,
        quote = FALSE,
        sep = "\t",
        col.names = FALSE,
        row.names = FALSE);
    return(invisible(bed));
}

madeBEDgraphGange = function(x, filename, chr, start, end){
    mwas = getMWASandLocations(x);
    chrnames = levels(mwas$chr);
    if( is.factor(chr) )
        chr = as.character(chr);
    if( is.character(chr) ){
        chrn = match(chr, chrnames, nomatch = 0L);
    } else {
        chrn = chr;
    }
    
    keep =  which(
                as.integer(mwas$chr) == chrn &
                (mwas$end >= start) & 
                (mwas$start <= end));
    
    bed = data.frame(
        chrom = mwas$chr[keep],
        chromStart  = mwas$start[keep] - 1,
        chromEnd =  mwas$end[keep],
        score = mwas$`p-value`[keep]
    )
    write.table(
        file = filename,
        x = bed,
        quote = FALSE,
        sep = "\t",
        col.names = FALSE,
        row.names = FALSE);
    return(invisible(bed));
}
