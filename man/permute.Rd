\name{permute}
\alias{permute}
\title{Permutation Generation}
\description{
  Generates permutations a pair at a time.
}
\usage{
permute(N=4,K=1,initialize=0)
}
\arguments{
	\item{N}{The number of elements in the permutation vector. The elements are numbered 1:N}
	\item{K}{The number of exchange pairs to be returned at each access}
	\item{initialize}{Set to 1 to initialize. Set to 0 for additional pairs.}
}
\details{
On first call, set initialize to 1. On subsequent accesses, initialize should be 0. 
The function returns at most K exchange pairs at each access. The function should 
be repeatedly called until result is false.

The function returns a list.
}
\note{This is the minimum time and effort routine. It is used for Exact permutation
calculations. The present function is a wrapper for the C code and may be useful
for other purposes. The code follows Reingole and Nievergelt (1977).}

\value{
	\item{result}{Returns false when the last permutation has been generated or a failure
		has occurred. Returns true for success and at most K new pairs in vec.}
	\item{N}{The input parameter}
	\item{K}{The input parameter}
	\item{vec}{A 2*K vector of element numbers to be exchanged -- permutations may be
		obtained by successively applying these pairs.}
	\item{initialize}{Always returns 1}
	\item{count}{The number of pairs in vec. It may be less than K for the final set.}
	}

\references{
	\describe{
	\item{Reingold, E.M., Nievergelt, J., and Deo, N. (1977)}{Combinatorial 
		Algorithms Theory and Practice. Prentice Hall, New Jersey. p170.}
	}
}
\author{
	Bob Wheeler \email{rwheeler@echip.com}
}

\keyword{htest}