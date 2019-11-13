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
        sum(vapply(newbam$startsfwd,length,0)),
        sum(vapply(newbam$startsrev,length,0)));
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
bam.coverage.by.density = function(
        rbam,
        cpgset,
        noncpgset,
        minfragmentsize,
        maxfragmentsize){

    fragdistr = c(
        rep(1, minfragmentsize-1),
        seq(1, 0, length.out = (maxfragmentsize-minfragmentsize)/1.5+1));
    fragdistr = fragdistr[fragdistr>0];

    if( is.null(noncpgset) ){
        noncpgset =
            noncpgSitesFromCpGset(cpgset = cpgset, distance = maxfragmentsize);
    }
    # sum(vapply(noncpgset,length,0))
    # newcpgset = noncpgset;
    # for( chr in seq_along(noncpgset) ){
    #     newcpgset[[chr]] = sort.int( c(cpgset[[chr]], noncpgset[[chr]]) );
    # }
    # rm(noncpgset);

    cpgdensity1 = calc.coverage(
                        rbam = list(startsfwd = cpgset),
                        cpgset = cpgset,
                        fragdistr = fragdistr);
    cpgdensity2 = calc.coverage(
                        rbam = list(startsrev = lapply(cpgset,`-`,1L)),
                        cpgset = cpgset,
                        fragdistr = fragdistr[-1]);
    cpgdensity = 
            unlist(cpgdensity1, recursive = FALSE, use.names = FALSE) +
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
    z = locpoly(
            x = c(sqrtcpgdensity, double(length(noncoverage))),
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
    strandfunX = function(st){c(length(st$chrX), sum(vapply(st,length,0)))};
    rbam$qc$chrX.count =
                    strandfunX(rbam$startsfwd) +
                    strandfunX(rbam$startsfwd);
    class(rbam$qc$chrX.count) = "qcChrX"

    strandfunY = function(st){c(length(st$chrY), sum(vapply(st,length,0)))};
    rbam$qc$chrY.count =
                    strandfunY(rbam$startsfwd) +
                    strandfunY(rbam$startsfwd);
    class(rbam$qc$chrY.count) = "qcChrY"

    return(rbam);
}
