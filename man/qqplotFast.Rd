\name{qqPlotFast}
\alias{qqPlotFast}
\alias{qqPlotPrepare}
\title{
    Fast QQ-plot for Large Number of P-values
}
\description{
    Function \code{qqPlotFast} creates a QQ-plot
    with a confidence band and
    an estimate of inflation factor lambda.
    It optimized to work quickly even for tens of millions of p-values.
}
\usage{
qqPlotPrepare(
        pvalues,
        ntests = NULL,
        ismlog10 = FALSE)
qqPlotFast(
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
            paste("\u2013", " log"[10]*"(", italic("P"), "), observed")))
}
\arguments{
    \item{pvalues}{
        Vector of p-values. 
        As is (if \code{ismlog10 = FALSE}) or minus log10 transformed
        (if \code{ismlog10 = TRUE}).
    }
    \item{ntests}{
        If only significant p-values are provided,
        the total number of tests performed.\cr
        By default \code{ntests} is  equal to the length of \code{pvalues}.
    }
    \item{ismlog10}{
        Specifies whether the provides p-values (\code{pvalues} parameter)
        are minus log10 transformed (\code{- log10(pv)})
    }
    \item{x}{
        Either a vector of p-values, as in \code{qqPlotPrepare},
        or the object returned by \code{qqPlotPrepare}.
    }
    \item{ci.level}{
        Significance level of the confidence band.
        Set to \code{NULL} avoid plotting the confidence band.
    }
    \item{ylim}{
        Numeric vectors of length 2, giving the y coordinate range.
        Exactly as in \link[graphics:plot.window]{Plotting Parameters}.
    }
    \item{newplot}{
        If \code{TRUE}, the function creates a new plot window.
    }
    \item{col}{
        The QQ-plot curve color.
    }
    \item{col.band}{
        Confidence band curve color.
    }
    \item{cex}{
        The size of QQ-plot points.
        As in \link[graphics:par]{Graphics Parameters}.
    }
    \item{lwd}{
        The line width.\cr
        As in \link[graphics:par]{Graphics Parameters}.
    }
    \item{axistep}{
        Distance between axis label ticks for both axis.
    }
    \item{yaxmax}{
        Maximum reach of the y axis.
    }
    \item{makelegend}{
        If true, add legend to the plot.
    }
    \item{xlab, ylab}{
        Axis labels. As in \link[graphics:plot]{plot function}.
    }
}
\details{
    The function \code{qqPlotFast} creates a QQ-plot. \cr
    The function \code{qqPlotPrepare} extracts the necessary information
    from a vector of p-values sufficient for creating QQ-plot. \cr
    The resulting object is many times smaller than the vector of p-values.
}
\note{
    The plot has no title. To add a title use \code{\link[graphics]{title}}.
}
\value{
    The function \code{qqPlotPrepare} returns an object with
    the necessary information from a vector of p-values
    sufficient for creating QQ-plot.
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
\note{
    The function works faster if the p-values are sorted.
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{
# Million p-values
n = 1e6

# Null p-values
pv = runif(n)

# QQ-plot should be nearly diagonal
qqPlotFast(pv)
title("QQ-plot")

# Size of p-values before extraction of QQ-plot info
object.size(pv)

# Extract the QQ-plot info
qq = qqPlotPrepare(pv)

# Size of the QQ-plot info object
object.size(qq)

# Create QQ-plot, it is the same
qqPlotFast(qq)

# Create QQ-plot with plotting parameters
qqPlotFast(qq, ylim = c(0,10), yaxmax = 9, axistep = 3, lwd = 3, cex = 1)
}
