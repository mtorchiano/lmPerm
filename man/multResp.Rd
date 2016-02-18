\name{multResp}
\alias{multResp}
\title{Multple response creation}
\description{
	This function creates a multiple response matrix for its argument variables. When used 
	on the lhs of the formula in lmp() or aovp() it will create a matrix containing one or more
	response columns from variables defined in the data argument.
}
\usage{
	multResp(...)
}

\arguments{
	\item{...}{Variable names separated by commas.}
}


\value{
	A matrix with named columns
}

\examples{

	A<-1:5
	B<-1:5
	multResp(A,B)

	data(Plasma)
	anova(lmp(multResp(Amin,Pct,sinpoly)~.,Plasma))
}
\author{
	Bob Wheeler \email{bwheelerg@gmail.com}
	
	Please cite this program as follows:

Wheeler, R.E. (2010). multResp() \emph{lmPerm}. The R project for statistical computing \url{http://www.r-project.org/}
}


\keyword{design}

