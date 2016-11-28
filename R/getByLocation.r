# Get MWAS results by locations
getTestsByLocation = function(x, chr, position){
    if(is.list(x)){
        param = parameterPreprocess(x);
        dirmwas = param$dirmwas;
        dircov = param$dircoveragenorm
    } else {
        dirmwas = x;
        dircov = paste0(x, '/../..');
    }
    chrnames = readLines(paste0(dircov,'/CpG_chromosome_names.txt'));
    if(is.character(chr)) {
        chrn = match(chr, chrnames, nomatch = 0L);
    } else {
        chrn = chr;
    }
    locmat = fm.load(paste0(dircov,'/CpG_locations'));
    maxmult = (max(locmat)+2)*2;
    
    mch = match(chrn*maxmult + position, locmat %*% c(maxmult,1), nomatch = 0L);
    fm = fm.open(paste0(dirmwas, '/Stats_and_pvalues'));
    testsmch = fm[mch[mch>0L],];
    colnames(testsmch) = colnames(fm);
    close(fm);
        
    rez = data.frame(
        chr = rep(chr, times = length(position)/length(chr))[mch>0L],
        position = position[mch>0L],
        testsmch);
    
    return(rez);
}
