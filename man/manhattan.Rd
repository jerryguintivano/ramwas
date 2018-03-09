\name{manhattan}
\alias{manPlotPrepare}
\alias{manPlotFast}
\title{
    Fast Manhattan plot for Large Number of P-values
}
\description{
    The function \code{manPlotFast} creates a Manhattan plot.\cr
    The function \code{manPlotPrepare} extracts necessary information
    from a vector of p-values sufficient for creating a Manhattan plot.\cr
    It optimized to work quickly even for tens of millions of p-values.
}
\usage{
manPlotPrepare(
        pvalues,
        chr,
        pos,
        ismlog10 = FALSE,
        chrmargins = 5e6)

manPlotFast(
        man,
        ylim = NULL,
        colorSet = c('steelblue4',"#2C82D1","#4CB2D1"),
        yaxmax = NULL,
        lwd = 3,
        axistep = 2,
        cex = 1)
}
\arguments{
    \item{pvalues}{
        Vector of p-values. 
        As is (if \code{ismlog10 = FALSE}) or minus log10 transformed
        (if \code{ismlog10 = TRUE}).
    }
    \item{chr, pos}{
        Vectors indicating the chromosomes and genomic positions  (in basepairs)
        for each p-value in \code{pvalues}.
    }
    \item{ismlog10}{
        Specifies whether the provides p-values (\code{pvalues} parameter)
        are minus log10 transformed (\code{- log10(pv)})
    }
    \item{chrmargins}{
        The plot margins at the ends of chromosomes (in basepairs).
    }
    \item{man}{
        Object returned by \code{manPlotPrepare}.
    }
    \item{ylim}{
        Numeric vectors of length 2, giving the y coordinate range.
        Exactly as in \link[graphics:plot.window]{Plotting Parameters}.
    }
    \item{colorSet}{
        Colors of points, rotating over chromosomes.\cr
        Points for first chromosome have color \code{colorSet[1]},
        next \code{colorSet[2]}, etc. Once the colors are exhausted,
        the colors are reused from the beginning.
    }
    \item{yaxmax}{
        Maximum reach of the y axis.
    }
    \item{lwd}{
        The line width.\cr
        As in \link[graphics:par]{Graphics Parameters}.
    }
    \item{axistep}{
        Distance between axis label ticks for y axis.
    }
    \item{cex}{
        The size of Manhattan plot points.
        As in \link[graphics:par]{Graphics Parameters}.
    }
}
\details{
    The function \code{manPlotFast} creates Manhattan plot. \cr
    It requires the use of the function \code{manPlotPrepare}
    which extracts the necessary information
    from a vector of p-values sufficient for creating Manhattan plot. \cr
    The resulting object is many times smaller than the vector of p-values.
}
\note{
    The plot has no title. To add a title use \code{\link[graphics]{title}}.
}
\value{
    This function \code{manPlotPrepare} returns an object with
    information for creating Manhattan plot.
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{
# Simulate data (9 chromosomes, million tests each)
chr = rep(paste0('chr',1:9), each = 1e6)
pos = rep(1:1e6, 9)
pv = runif(9e6)^1.1

# Extract the Manhattan plot info
man = manPlotPrepare(pv, chr, pos, chrmargins = 1000)

# Create Manhattan plot
manPlotFast(man)
title("Manhattan plot")

# Size of p-values before extraction of Manhattan plot info
object.size(list(pv, chr, pos))

# Size of the Manhattan plot info object
object.size(man)
}
