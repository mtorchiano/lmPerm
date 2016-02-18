\name{summaryC}
\alias{summaryC}


\title{Coefficient summary from aovp}
\description{Summarize an analysis of variance model showing coefficients}
\usage{

## S3 method for class 'aovlist':
summaryC(object, ...)

}
\arguments{

	\item{object}{an object of class \code{aovlist}}
	\item{\dots}{further arguments passed on.}

}

\details{
Summarizes \code{aovp} output in terms of coefficients with p-values.
}
\examples{
	data(Alfalfa)
	summaryC(aovp(Yield~Variety*Date+Error(Block/Variety),Alfalfa))}

\author{
	Bob Wheeler \email{rwheeler@echip.com}
}

\keyword{htest}



