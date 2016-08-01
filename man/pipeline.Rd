\name{pipeline}
\alias{ramwas1scanBams}
\alias{ramwas2collectqc}
\alias{ramwas3NormalizedCoverage}
\alias{ramwas4PCA}
\alias{ramwas5MWAS}
\alias{ramwas6crossValidation}
\alias{ramwas7multiMarker}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	High Level Pipeline Functions.
}
\description{
	These functions provide a simple way to run all steps of RaMWAS pipeline.
}
\usage{
ramwas1scanBams(param)
ramwas2collectqc(param)
ramwas3NormalizedCoverage(param)
ramwas4PCA(param)
ramwas5MWAS(param)
ramwas6crossValidation(param)
ramwas7multiMarker(param)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{
%%     ~~Describe \code{param} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (param) 
{
    param = parameterPreprocess(param)
    stopifnot(!is.null(param$bamnames))
    cat(file = paste0(param$dirfilter, "/Log.txt"), date(), ", Scanning bams.", 
        "\n", sep = "", append = FALSE)
    if (param$cputhreads > 1) {
        cl = makeCluster(param$cputhreads)
        z = clusterApplyLB(cl, param$bamnames, .ramwas1scanBamJob, 
            param = param)
        stopCluster(cl)
    }
    else {
        z = character(length(param$bamnames))
        names(z) = param$bamnames
        for (i in seq_along(param$bamnames)) {
            z[i] = .ramwas1scanBamJob(bamname = param$bamnames[i], 
                param = param)
        }
    }
    cat(file = paste0(param$dirfilter, "/Log.txt"), date(), ", Done scanning bams.", 
        "\n", sep = "", append = TRUE)
    return(z)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
