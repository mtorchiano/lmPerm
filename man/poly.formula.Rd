\name{poly.formula}
\alias{poly.formula}
\title{Makes a polynomial model}
\description{
 Formulas are expanded to accommodate special functions for continuous and mixture variables.
}
\usage{
	poly.formula(frml)
}

\arguments{
	\item{frml}{A formula using ~ in the usual way.}
}
\details{
This function expands formulas to accommodate polynomial models for which R has minimal support.
Assuming for illustration that there are three variables, A, B, and C, the following expressions
may be used.  

All agruments to quad(), cubic(), and cubicS() must be numeric.

\itemize{
	\item quad(A,B,C) makes \eqn{(A+B+C)^2+I(A^2)+I(B^2)+I(C^2)}
	\item cubic(A,B,C) makes \eqn{(A+B+C)^3+I(A^2)+I(B^2)+I(C^2)+I(A^3)+I(B^3)+I(C^3)}
 	\item cubicS(A,B,C) makes \eqn{(A+B+C)^3+I(A*B*(A-B))+I(A*C*(A-C))+I(B*C*(B-C))}
}

The cubicS() function produces a non-singular representation of a cubic model, when the 
variables are mixture variables, that is when the rows of \code{data} sum to a constant
value, usually 1.0. Because of the mixture constraint, models containing mixture variables
should not have a constant term. The linear and quadratic models for mixture variables
A, B, and C are given by \eqn{-1+(A+B+C)} and \eqn{-1+(A+B+C)^2} respectively. See Gorman and Hinman [1962] for 
details.

}

\value{
	An expanded formula is returned.
}

\examples{

	poly.formula(y~quad(A,B,C)+Error(block))
}
\author{
	Bob Wheeler \email{bwheelerg@gmail.com}
	
	Please cite this program as follows:

Wheeler, R.E. (2010). poly.formula \emph{lmPerm}. The R project for statistical computing \url{http://www.r-project.org/}
}

\references{
  Gorman, J.W. and Hinman, J.E. (1962). Simplex lattice designs for 
  multicomponent systems. \emph{Technometrics}. 4-4. 463-487.
}
\keyword{design}

