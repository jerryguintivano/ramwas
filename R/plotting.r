# Prepare the compact summary of p-values for QQ-plot
qqPlotPrepare = function(pvalues, ntests = NULL, ismlog10 = FALSE){
    if(is.null(ntests))
        ntests = length(pvalues);
    
    
    if(ismlog10) {
        ypvs = pvalues;
    } else {
        ypvs = -log10(pvalues);
    }
    xpvs = -log10(seq_along(ypvs) / (ntests+1));
    
    if(is.unsorted(-ypvs))
        ypvs = sort.int(ypvs, decreasing = TRUE);
    
    if(length(ypvs)*2 > ntests) {
        lambda =
                qchisq(p = 0.1^ypvs[ntests/2], df = 1, lower.tail = FALSE) / 
                qchisq(p = 0.5, df = 1, lower.tail = FALSE);
        # lambda = sprintf("%.3f",log10(pvalues[ntests/2]) / log10(0.5));
        # legend("bottom", legend = bquote(lambda == .(lambda)), bty = "n")
        # 		text(mx, mx/2, bquote(lambda == .(lambda)), pos=2)
    } else {
        lambda = NULL;
    }
    
    if(length(ypvs) > 1000) {
        # need to filter a bit, make the plotting faster
        levels = as.integer( (xpvs - xpvs[1])/(tail(xpvs,1) - xpvs[1]) * 2000);
        keep = c(TRUE, diff(levels)!=0);
        levels = as.integer( (ypvs - ypvs[1])/(tail(ypvs,1) - ypvs[1]) * 2000);
        keep = keep | c(TRUE, diff(levels)!=0);
        keep = which(keep);
        ypvs = ypvs[keep];
        xpvs = xpvs[keep];
        # 		rm(keep)
    } else {
        keep = seq_along(ypvs)
    }
    
    qq = list(
        xpvs = xpvs, # p-values expected under null
        ypvs = ypvs, # observer p-values
        keep = keep, # indices of the preserved p-values
        ntests = ntests, # Number of tests
        lambda = lambda  # Estimate of inflation factor lambda
    );
    class(qq) = "qqPlotInfo";
    return(qq);
}

qqPlotFast = function(
            x, 
            ntests = NULL, 
            ismlog10 = FALSE, 
            ci.level = 0.05, 
            ylim = NULL, 
            newplot = TRUE, 
            col = "#D94D4C", 
            cex = 0.5, 
            yaxmax = NULL, 
            lwd = 3, 
            axistep = 2, 
            col.band = "#ECA538"){
    
    # Get compact summary of p-values for QQ-plot
    if( class(x) == "qqPlotInfo" ){
        qq = x;
    } else {
        qq = qqPlotPrepare(pvalues = x, ntests = ntests, ismlog10 = ismlog10);
    }

    # Axis ranges
    mx = head(qq$xpvs,1) * 1.05;
    if( is.null(ylim) ) {
        my = max(mx*1.15,head(qq$ypvs,1)) * 1.05;
        ylim = c(0,my);
    } else {
        my = ylim[2];
    }
    if(is.null(yaxmax))
        yaxmax = floor(my);
    
    if(newplot){
        plot(
            x = NA, 
            y = NA, 
            ylim = ylim, 
            xlim = c(0, mx), 
            xaxs = "i", 
            yaxs = "i", 
            xlab = expression(
                paste("\u2013", " log"[10]*"(", italic("P"), "), null")),
            ylab = expression(
                paste("\u2013", " log"[10]*"(", italic("P"), "), observed")),
            axes = FALSE);
        axis(1, seq(0, mx + 2, axistep), lwd = lwd);
        axis(2, seq(0, yaxmax, axistep), lwd = lwd);
    }
    abline(a = 0, b = 1, col = "grey", lwd = lwd)
    points(qq$xpvs, qq$ypvs, col = col, cex = cex, pch = 19);
    
    if( !is.null(ci.level) ){
        if( (ci.level>0)&(ci.level<1) ){
            quantiles = qbeta(
                p = rep(c(ci.level/2,1-ci.level/2), each=length(qq$xpvs)), 
                shape1 = qq$keep, 
                shape2 = qq$ntests - qq$keep + 1);
            quantiles = matrix(quantiles, ncol=2);
            
            lines( qq$xpvs, -log10(quantiles[,1]), col = col.band, lwd = lwd);
            lines( qq$xpvs, -log10(quantiles[,2]), col = col.band, lwd = lwd);
        }
    }
    
    if( !is.null(ci.level) ){
        legend(
            "topleft", 
            legend = c(
                    expression(paste(italic("P"), " value")),
                    sprintf("%.0f%% CI",100-ci.level*100)),
            lwd = c(0, lwd), 
            pch = c(19, NA_integer_), 
            lty = c(0, 1), 
            col = c(col, col.band),
            box.col = "transparent",
            bg = "transparent");
    } else {
        legend(
            "topleft", 
            legend = expression(paste(italic("P"), " value")),
            lwd = 0, 
            pch = 19, 
            lty = 0, 
            col = col,
            box.col = "transparent",
            bg = "transparent");
    }
    if( !is.null(qq$lambda) ){
        lastr = sprintf("%.3f", qq$lambda);
        legend("bottom", legend = bquote(lambda == .(lastr)), bty = "n")
    }
    return(invisible(qq));
}


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
