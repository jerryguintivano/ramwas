\name{plotCV}
\alias{plotPrediction}
\alias{plotCVcors}
\alias{plotROC}
\title{
    Plotting Functions used in Cross Validation Analysis
    (Methylation Risk Score).
}
\description{
    The function \code{plotPrediction} plots 
    cross validation predictions of a phenotype
    against true values of the phenotype
    with multiple summary stats in the title.
    
    The function \code{plotCVcors} plots
    the predictive power (correlations) across
    predictions using various numbers of markers.

    The function \code{plotROC} plots an ROC (Receiver operating characteristic)
    curve for predictions of a binary outcome.
}
\usage{
plotPrediction(
        param,
        outcome,
        forecast,
        cpgs2use,
        main,
        dfFull = NULL)
        
plotCVcors(cl, param)

plotROC(outcome, forecast)
}
\arguments{
    \item{param}{
        List of parameters as described in the "RW6_param.Rmd" vignette. \cr
        Try: \code{vignette("RW6_param","ramwas")}.\cr
        Only \code{modeloutcome}, \code{cvnfolds}
        and \code{mmalpha} elements are used.
    }
    \item{outcome}{
        Values of a phenotype. Must be binary for \code{plotROC}.
    }
    \item{forecast}{
        Predictions for the phenotype.
    }
    \item{cpgs2use}{
        Number of variables used for prediction (for the legend).
    }
    \item{main}{
        Part of the title (summary stats are added beneath).
    }
    \item{dfFull}{
        Number of degrees of freedom for the significance testing.\cr
        Default is: \code{length(forecast) - 2}
    }
    \item{cl}{
        List with three elements:
        \itemize{
            \item \code{x} - vector with the number of variables 
                used for prediction
            \item \code{corp} - Pearson correlations between 
                the predictions and the true value of the phenotype.
            \item \code{cors} - Spearman correlations between 
                the predictions and the true value of the phenotype.
        }
    }
    
}
\details{
    The \code{plotROC} and plot has no title.\cr 
    To add a title use \code{\link[graphics]{title}}.
}
\value{
    The \code{plotROC} returns the area under the curve (AUC) for the ROC.\cr
    The \code{plotPrediction} function returns the list of 
    calculated statistics printed in the title.\cr
    The \code{plotCVcors} returns nothing (\code{NULL}).
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{
# Sample data
n = 300
param = list(modeloutcome = "Age", mmalpha = 0, cvnfolds = 5)
outcome = rnorm(n, mean = 50, sd = 20)
forecast = outcome + rnorm(n, mean = 0, sd = 20)
cpgs2use = 1000
main = "Prediction success (simulated data)"

# Plot phenotype-prediction plot
plotPrediction(
        param,
        outcome,
        forecast,
        cpgs2use,
        main)

# Artificial data for plotCVcors()
cl = list(
    x =    c(50,  100, 200, 500,  1000),
    corp = c(0.1, 0.6, 0.7, 0.85, 0.8),
    cors = c(0.1, 0.6, 0.7, 0.85, 0.8) + rnorm(5)*0.1)

# Plot prediction performance by the number of markers    
plotCVcors(cl, param)  

# Make the outcome binary for ROC plot
outcome = (outcome > 50)

# Plot ROC curve and calculate the AUC
plotROC(outcome, forecast)
}
