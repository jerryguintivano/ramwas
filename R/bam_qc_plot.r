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
