# Fast QQ-plot
# With confidence band and lambda
# Makes small PDFs even with billions of p-values
qqPlotFast = function(pvalues, ntests=NULL, ci.level=0.05){

    if(is.null(ntests))
        ntests = length(pvalues);

    if(is.unsorted(pvalues))
        pvalues = sort.int(pvalues);

    ypvs = -log10(pvalues);
    xpvs = -log10(seq_along(ypvs) / ntests);

    if(length(ypvs) > 1000){
        # need to filter a bit, make the plotting faster
        levels = as.integer( (xpvs - xpvs[1])/(tail(xpvs,1) - xpvs[1]) * 2000);
        keep = c(TRUE, diff(levels)!=0);
        levels = as.integer( (ypvs - ypvs[1])/(tail(ypvs,1) - ypvs[1]) * 2000);
        keep = keep | c(TRUE, diff(levels)!=0);
        keep = which(keep);
        ypvs = ypvs[keep];
        xpvs = xpvs[keep];
        #         rm(keep)
    } else {
        keep = seq_along(ypvs)
    }
    mx = head(xpvs,1)*1.05;
    my = max(mx*1.15,head(ypvs,1))*1.05;
    plot(NA,NA, ylim = c(0,my), xlim = c(0,mx), xaxs="i", yaxs="i",
          xlab = expression("- log"[10]*"(p-value), expected under null"),
          ylab = expression("- log"[10]*"(p-value), observed"));
    # xlab = "-Log10(p-value), expected under null",
    # ylab = "-Log10(p-value), observed");
    lines(c(0,mx),c(0,mx),col="grey")
    points(xpvs, ypvs, col = "red", pch = 19, cex = 0.25);

    if(!is.null(ci.level)){
        if((ci.level>0)&(ci.level<1)){
            quantiles = qbeta(p = rep(c(ci.level/2,1-ci.level/2),
                                      each=length(xpvs)),
                              shape1 = keep,
                              shape2 = ntests - keep + 1);
            quantiles = matrix(quantiles, ncol=2);

            lines( xpvs, -log10(quantiles[,1]), col="cyan4")
            lines( xpvs, -log10(quantiles[,2]), col="cyan4")
        }
    }
    legend("bottomright",
           c("P-values", sprintf("%.0f %% Confidence band",100-ci.level*100)),
           lwd = c(0,1),
           pch = c(19,NA_integer_),
           lty = c(0,1),
           col=c("red","cyan4"))
    if(length(pvalues)*2>ntests){
        lambda = sprintf("%.3f",
                         qchisq(pvalues[ntests/2],1,lower.tail = FALSE) /
                         qchisq(0.5,1,lower.tail = FALSE));
        legend("bottom", legend = bquote(lambda == .(lambda)), bty = "n")
        #         text(mx, mx/2, bquote(lambda == .(lambda)), pos=2)
    }
    return(invisible(NULL));
}
