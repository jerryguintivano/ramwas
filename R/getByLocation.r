# Get MWAS results by locations
getTestsByLocation = function(x, chr, position){
    if(is.list(x)){
        param = parameterPreprocess(x);
        dirmwas = param$dirmwas;
        dircov = param$dircoveragenorm
    } else {
        dirmwas = x;
        dircov = paste0(x, "/../..");
    }
    chrnames = readLines(paste0(dircov,"/CpG_chromosome_names.txt"));
    if(is.factor(chr))
        chr = as.character(chr);
    if(is.character(chr)){
        chrn = match(chr, chrnames, nomatch = 0L);
    } else {
        chrn = chr;
    }
    locmat = fm.load(paste0(dircov,"/CpG_locations"));
    maxmult = (max(locmat)+2)*2;

    mch = match(chrn*maxmult + position, locmat %*% c(maxmult,1), nomatch = 0L);
    fm = fm.load(paste0(dirmwas, "/Stats_and_pvalues"));
    testsmch = fm[mch[mch>0L],];
    colnames(testsmch) = colnames(fm);
    # close(fm);

    rez = data.frame(
        chr = rep(chrnames[chrn], times = length(position)/length(chr))[mch>0L],
        position = position[mch>0L],
        testsmch);

    return(rez);
}

# Get MWAS results by locations
getDataByLocation = function(x, chr, position){
    if(is.list(x)){
        param = parameterPreprocess(x);
        dircov = param$dircoveragenorm
    } else {
        dircov = x;
    }
    chrnames = readLines(paste0(dircov,"/CpG_chromosome_names.txt"));
    if(is.factor(chr))
        chr = as.character(chr);
    if(is.character(chr)){
        chrn = match(chr, chrnames, nomatch = 0L);
    } else {
        chrn = chr;
    }
    locmat = fm.load(paste0(dircov,"/CpG_locations"));
    maxmult = (max(locmat)+2)*2;

    mch = match(chrn*maxmult + position, locmat %*% c(maxmult,1), nomatch = 0L);
    fm = fm.open(paste0(dircov, "/Coverage"));
    mat = fm[,mch[mch>0L]];
    rownames(mat) = rownames(fm);
    close(fm);

    rez = list(
        chr = rep(chrnames[chrn], times = length(position)/length(chr))[mch>0L],
        position = position[mch>0L],
        matrix = mat);

    return(rez);
}

getLocations = function(x){
    if(is.list(x)){
        param = parameterPreprocess(x);
        dircov = param$dircoveragenorm
    } else {
        dircov = x;
    }   
    chrnames = readLines(paste0(dircov,"/CpG_chromosome_names.txt"));
    locmat = fm.load(paste0(dircov,"/CpG_locations"));
    chr = locmat[,1];
    levels(chr) = chrnames;
    class(chr) = "factor";
    locations = data.frame(chr = chr, 
                           start = locmat[,2], 
                           end = locmat[,2] + 1L);
    if(ncol(locmat) >= 3)
        locations$end = locmat[,3];
    return(locations);
}

getMWASandLocations = function(x){
    if(is.list(x)){
        param = parameterPreprocess(x);
        dirmwas = param$dirmwas;
        dircov = param$dircoveragenorm
    } else {
        dirmwas = x;
        dircov = paste0(x, "/../..");
    }
    locations = getLocations(dircov);
    mwas = fm.load(paste0(dirmwas, "/Stats_and_pvalues"));
    result = data.frame(locations, ttest = mwas[,2], pvalue = mwas[,3], qvalue = mwas[,4]);
    return(result);
}

if(FALSE){
    library(filematrix)
    setwd('D:/')
    x = '/gpfs_fs/pharm/AUD_brain/RaMWAS/coverage_norm_50/PCA_05_cvrts_9a48f039/Testing_AUD_0_PCs'
    # locs = getLocations(x)
    mwas1 = ramwas:::getMWASandLocations(x);
    
    y = '/gpfs_fs/pharm/NESDA/RaMWAS/Alcohol/Split_cntrl_aaudit01_320/PCA_23_cvrts_81fb5393/Testing_cntrl_aaudit01_r_3_PCs'
    mwas2 = ramwas:::getMWASandLocations(y);
    
    head(mwas1)
    head(mwas2)
    mwas1$chr[1:10]
    mwas2$chr[1:10]
}

# setwd('D:/')
# dr = '/gpfs_fs/pharm/NESDA/RaMWAS/Alcohol/Split_cntrl_aaudit01_320'
# setwd(dr);
# tbl = read.table('AUDBrain_AUD_NESDA_aaudit01_overlap_all.txt', header = TRUE);
# head(tbl)
# 
subsetCoverageDirByLocation = function(x, chr, position, targetdir){
    if(is.list(x)){
        param = parameterPreprocess(x);
        dircov = param$dircoveragenorm
    } else {
        dircov = x;
    }
    chrnames = readLines(paste0(dircov,"/CpG_chromosome_names.txt"));
    if(is.factor(chr))
        chr = as.character(chr);
    if(is.character(chr)){
        chrn = match(chr, chrnames, nomatch = 0L);
    } else {
        chrn = chr;
    }
    locmat = fm.load(paste0(dircov,"/CpG_locations"));
    maxmult = (max(locmat)+2)*2;

    mch = match(chrn*maxmult + position, locmat %*% c(maxmult,1), nomatch = 0L);
    if(any(mch == 0L))
        stop('Some locations not found')

    mch = sort(mch);

    dir.create(targetdir, showWarnings = FALSE)
    setwd(targetdir)

    file.copy(from = paste0(dircov,'/CpG_chromosome_names.txt'),
              to = 'CpG_chromosome_names.txt');
    fm = fm.create.from.matrix('CpG_locations', locmat[mch,]);
    close(fm);

    fi = fm.open(paste0(dircov,'/Coverage'))
    fo = fm.create('Coverage', nrow = nrow(fi), ncol = length(mch), type = fi$type, size = fi$size);

    step1 = 1024;
    runto = length(mch);
    nsteps = ceiling(runto/step1);
    for( part in seq_len(nsteps) ) { # part = 1
        message('Loop filling coverage matrix, step ', part, ' of ', nsteps);
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

    close(fo)
    close(fi)
}
# 
# library(filematrix)
# subsetCoverageDirByLocation(
#     x = '/gpfs_fs/pharm/NESDA/RaMWAS/Alcohol/Split_cntrl_aaudit01_320', 
#     chr = tbl$chr, 
#     position = tbl$position, 
#     targetdir = '/gpfs_fs/pharm/NESDA/RaMWAS/Alcohol/Split_cntrl_aaudit01_320/subset')
# 


