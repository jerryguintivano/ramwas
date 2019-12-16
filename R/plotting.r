# Prepare the compact summary of p-values for QQ-plot
qqPlotPrepare = function(pvalues, ntests = NULL, ismlog10 = FALSE){
    if( is.null(ntests) )
        ntests = length(pvalues);
    
    if( ismlog10 ){
        ypvs = pvalues;
    } else {
        ypvs = -log10(pvalues);
    }
    xpvs = -log10(seq_along(ypvs) / (ntests+1));
    
    if( is.unsorted(-ypvs) )
        ypvs = sort.int(ypvs, decreasing = TRUE);
    
    if( length(ypvs)*2 > ntests ){
        lambda =
            qchisq(p = 0.1^ypvs[ntests/2], df = 1, lower.tail = FALSE) /
            qchisq(p = 0.5, df = 1, lower.tail = FALSE);
    } else {
        lambda = NULL;
    }
    
    if( length(ypvs) > 1000 ){
        # need to filter a bit, make the plotting faster
        levels = as.integer( (xpvs - xpvs[1])/(tail(xpvs,1) - xpvs[1]) * 2000);
        keep = c(TRUE, diff(levels)!=0);
        levels = as.integer( (ypvs - ypvs[1])/(tail(ypvs,1) - ypvs[1]) * 2000);
        keep = keep | c(TRUE, diff(levels)!=0);
        keep = which(keep);
        ypvs = ypvs[keep];
        xpvs = xpvs[keep];
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

# Create QQ-plot from p-values or prepared summary
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
            col.band = "#ECA538",
            makelegend = TRUE,
            xlab = expression(
                paste("\u2013", " log"[10]*"(", italic("P"), "), null")),
            ylab = expression(
                paste("\u2013", " log"[10]*"(", italic("P"), "), observed"))
        ){
    
    # Get compact summary of p-values for QQ-plot
    if( methods::is(x, "qqPlotInfo") ){
        qq = x;
    } else {
        qq = qqPlotPrepare(pvalues = x, ntests = ntests, ismlog10 = ismlog10);
    }

    # Axis ranges
    mx = head(qq$xpvs,1) * 1.05;
    if( is.null(ylim) ){
        my = max(mx, head(qq$ypvs,1) * 1.05) ;
        ylim = c(0, my);
    } else {
        my = ylim[2];
    }
    if( is.null(yaxmax) )
        yaxmax = floor(my);
    
    if( newplot ){
        plot(
            x = NA,
            y = NA,
            ylim = ylim,
            xlim = c(0, mx),
            xaxs = "i",
            yaxs = "i",
            xlab = xlab,
            ylab = ylab,
            axes = FALSE);
        axis(1, seq(0, mx + axistep, axistep), lwd = lwd);
        axis(2, seq(0, yaxmax, axistep), lwd = lwd);
    }
    abline(a = 0, b = 1, col = "grey", lwd = lwd);
    points(qq$xpvs, qq$ypvs, col = col, cex = cex, pch = 19);
    
    if( !is.null(ci.level) ){
        if( (ci.level>0)&(ci.level<1) ){
            quantiles = qbeta(
                p = rep(c(ci.level/2,1-ci.level/2), each=length(qq$xpvs)), 
                shape1 = qq$keep, 
                shape2 = qq$ntests - qq$keep + 1);
            quantiles = matrix(quantiles, ncol = 2);
            
            lines( qq$xpvs, -log10(quantiles[,1]), col = col.band, lwd = lwd);
            lines( qq$xpvs, -log10(quantiles[,2]), col = col.band, lwd = lwd);
        }
    }
    if( makelegend ){
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
    }
    if( !is.null(qq$lambda) ){
        lastr = sprintf("%.3f", qq$lambda);
        legend("bottom", legend = bquote(lambda == .(lastr)), bty = "n")
    }
    return(invisible(qq));
}

# Prepare the compact summary of p-values for Manhattan plot
manPlotPrepare = function(
            pvalues,
            chr,
            pos,
            ismlog10 = FALSE,
            chrmargins = 5e6){
    # chr = locs[,1]; pos = locs[,2]; pvalues = mwas[,3]
    # z = getMWASandLocations(param)
    # chr = z$chr; pos = z$start; pvalues = z$`p-value`; 
    # ismlog10 = FALSE; chrmargins = 0;

    stopifnot( length(pvalues) == length(chr) );
    stopifnot( length(pvalues) == length(pos) );
        
    # Factorize chromosome
    if( is.double(chr) )
        chr = as.integer(chr);
    if( is.integer(chr) ){
        levels(chr) = as.character(seq_len(max(chr)));
        class(chr) = "factor";
    }
    if( is.character(chr) ){
        levels = str_sort(unique(chr), numeric = TRUE);
        chr = factor(chr, levels = levels);
    }

    # max of each chromosome
    poslist = split(pos, chr, drop = FALSE);
    poslist[vapply(poslist, length, 0)==0] = list(0);
    chrmax = as.double(vapply(poslist, max, 0));# + chrmargins;
    
    # chromosome starts on the plot
    names(chrmax) = NULL;
    offsets = c(0, cumsum(chrmax)) + chrmargins;
    names(offsets)[seq_along(poslist)] = levels(chr);
    
    # within plot coordinates
    x0 = offsets[unclass(chr)] + pos;
    if( ismlog10 ){
        y0 = pvalues;
    } else {
        y0 = -log10(pvalues);
    }
    
    # Prune the data
    yfac = as.integer(y0*100)+1L;
    yorder = sort.list(yfac);
    levels(yfac) = as.character(seq_len(max(yfac)));
    class(yfac) = "factor";
    
    ygroup = split(seq_along(yfac), yfac);
    for( i in seq_along(ygroup) ){ # i=1
        if( length(ygroup[[i]]) > 300 ){
            ygroup[[i]] = sample(ygroup[[i]], size = 300, replace = FALSE);
        }
    }
    # sum(vapply(ygroup, length, 0))
    keep = unlist(ygroup, use.names = FALSE);
    
    # Color code
    colindex = unclass(chr);
    
    # Chromosome names
    chrnames = gsub("chr", "", levels(chr));
    
    # Return minimum object;
    man = list(
        x = x0[keep],
        y = y0[keep],
        colindex = colindex[keep],
        offsets = offsets,
        chrnames = chrnames,
        chrmargins = chrmargins
    );
    class(man) = "manPlotInfo";
    return(man);
}

# Create Manhattan plot from prepared summary
manPlotFast = function(
            man,
            ylim = NULL,
            colorSet = c("steelblue4", "#2C82D1", "#4CB2D1"),
            yaxmax = NULL,
            lwd = 3,
            axistep = 2,
            cex = 1){
    
    if( !methods::is(man, "manPlotInfo") )
        stop("The \"man\" parameter is not produced by manPlotPrepare().");
    
    # Axis ranges
    if( is.null(ylim) ){
        my = max(man$y) * 1.05;
        ylim = c(0,my);
    } else {
        my = ylim[2];
    }
    if( is.null(yaxmax) )
        yaxmax = floor(my);
    
    # Plot frame
    plot(
        x = NA,
        y = NA,
        xlim = c(0, tail(man$offsets,1)),
        ylim = ylim,
        xaxs = "i",
        yaxs = "i",
        xlab = "Chromosome",
        ylab = expression(
            paste("\u2013", " log"[10]*"(", italic("P"), "), observed")),
        axes = FALSE);
    axis(
        side = 1,
        at = man$offsets,
        labels = rep("", length(man$offsets)),
        lwd = lwd);
    axis(
        side = 1,
        at = (man$offsets[-1] + man$offsets[-length(man$offsets)])/2,
        labels = man$chrnames,
        tick = FALSE,
        lwd = lwd);
    axis(
        side = 2,
        at = seq(0, yaxmax, axistep),
        lwd = lwd);
    
    # plot points in palette color
    oldPal = palette(colorSet);
    points(
        x = man$x,
        y = man$y,
        pch = 20,
        col = ((man$colindex-1L) %% length(colorSet)) + 1L,
        cex = cex);
    palette(oldPal);
}
