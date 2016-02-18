\name{summary}
\alias{summary.lmp}
\alias{summary.mlmp}
\alias{print.summary.lmp}
\alias{summary.aovp}
\alias{anova.lmp}


\title{Summarizing functions for linear models}
\description{Replaces corresponding functions in base package.}
\usage{
\method{summary}{lmp}(object, correlation = FALSE, symbolic.cor = FALSE, ...)
\method{summary}{mlmp}(object, ...)
\method{print}{summary.lmp}(x, digits = max(3, getOption("digits") - 3),
              symbolic.cor = x$symbolic.cor,
	      signif.stars= getOption("show.signif.stars"),	...)
\method{summary}{aovp}(object, intercept = FALSE, split,
                        expand.split = TRUE, keep.zero.df = TRUE, ...)
\method{anova}{lmp}(object, ...)
}
\arguments{ Same as for the corresponding functions in base package:

  \item{object}{an object of class \code{"lm"}, usually, a result of a
    call to \code{\link{lm}}.}
  \item{x}{an object of class \code{"summary.lm"}, usually, a result of a
    call to \code{summary.lm}.}
  \item{correlation}{logical; if \code{TRUE}, the correlation matrix of
    the estimated parameters is returned and printed.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{symbolic.cor}{logical. If \code{TRUE}, print the correlations in
    a symbolic form (see \code{\link{symnum}}) rather than as numbers.}
  \item{signif.stars}{logical. If \code{TRUE}, \dQuote{significance stars}
    are printed for each coefficient.}
   \item{intercept}{logical: should intercept terms be included?}
  \item{split}{an optional named list, with names corresponding to terms
    in the model.  Each component is itself a list with integer
    components giving contrasts whose contributions are to be summed.}
  \item{expand.split}{logical: should the split apply also to
    interactions involving the factor?}
  \item{keep.zero.df}{logical: should terms with no degrees of freedom
    be included?}
 \item{\dots}{further arguments passed to or from other methods.}


}

\details{
These modified functions are needed because the perm values, which are 
attached to the object, replace the usual test columns in the output 
from these functions.
}

\author{
	Bob Wheeler \email{rwheeler@echip.com}
}

\keyword{htest}