\name{aovp}
\alias{aovp}

\title{Fitting and testing ANOVA using permutation tests} 

\description{\code{aovp} is \code{\link{aov}} modified to use permutation tests instead of normal
theory tests. Like \code{aov}, the ANOVA model is fitted by a call to \code{\link{lmp}} for each
stratum. Timing differences between \code{aovp} and \code{aov} are negligible.} 

\usage{
aovp(formula, data = NULL, perm="Exact", seqs=FALSE, 
              center=TRUE, projections = FALSE, qr = TRUE,
              contrasts = NULL, ...)
}

\arguments{
	The arguments are the same as for \code{\link{aov}}. 
Additional parameters may be included. They are described in the 
"Additional Parameters" section below. These additional parameters
are the same as for \code{\link{lmp}}.

  \item{formula}{A formula specifying the model.}
  \item{data}{A data frame in which the variables specified in the
    formula will be found. If missing, the variables are searched for in
    the standard way.}
  \item{perm}{"Exact", "Prob", "SPR" will produce permutation probabilities. 
	Anything else, such as "", will produce F-test probabilities.}
  \item{seqs}{If TRUE, will calculate sequential SS. If FALSE, unique SS.}
  \item{center}{If TRUE, numerical variables will be centered.}
  \item{projections}{Logical flag: should the projections be returned?}
  \item{qr}{Logical flag: should the QR decomposition be returned?}
  \item{contrasts}{A list of contrasts to be used for some of the factors
    in the formula. These are not used for any \code{Error} term, and
    supplying contrasts for factors only in the \code{Error} term will give
    a warning.}
  \item{\dots}{Arguments to be passed to \code{aovp}, such as those listed
	below.}

}
 
\section{Additional Parameters}{
	These are the same as for \code{\link{lmp}}. 

	\describe{
	\item{settings}{If TRUE, settings such as sequential or unique will
		be printed. Default TRUE}
	\item{useF}{If TRUE, SS/Resid SS will be used, otherwise SS. The 
		default is TRUE}
	\item{maxIter}{For Prob and SPR: The maximum number of iterations. 
		Default 1000.}
	\item{Ca}{For Prob: Stop iterations when estimated standard error 
		of the estimated p is less than Ca*p. Default 0.1}
	\item{p0}{For SPR: Null hypothesis probability. Default 0.05}
	\item{p1}{For SPR: Alternative hypothesis probability. Default 0.06}
	\item{alpha}{For SPR: Size of SPR test. Default 0.01}
	\item{beta}{For SPR: Type II error for SPR test. Default 0.01}
	\item{maxExact}{For Exact: maximum number of observations allowed. 
		If data exceeds this, Prob is used. Default 10.}
	\item{nCycle}{For Prob and SPR: Performs a complete random permutation, 
		instead of pairwise exchanges, every nCycle cycles. Default 1000.}
	}
}

\details{

The model Y=Xb+Zg+e is assumed, where X is the incidence matrix for fixed effects,
and Z is an incidence matrix for random effects, with columns representing the several error strata.
The \code{aovp()} algorithm projects Y into strata such that each stratum has a single error term, such
as a stratum defined by whole blocks. X is also projected so that the model
in this stratum becomes P(Y)=P(X)bi+ei. 

The vector bi is divided into sources with dfj degrees of freedom for the 
jth source, and \code{summary(aovp())} will produce an ANOVA table for these sources
in the ith strata. See Venables and Ripley for details.

Either permutation test p-values or the usual F-test 
p-values will be output. Polynomial model terms are collected into 
sources, so that \code{Y~A+B+I(A^2)} will contain two sources, one for A with 2 df,
and one for B with 1 df. Sources for factors are treated as usual, and polynomial
terms and factors may be mixed in one model. The function \code{\link{poly.formula}} may 
be used to create polynomial models, and the function \code{\link{multResp}} may be used to
create multiresponse matrices for the lhs from variables defined in \code{data}.

The \code{Exact} method will permute the values exactly. The \code{Prob} and \code{SPR} methods will approximate
the permutation distribution by randomly exchanging pairs of Y elements. The \code{Exact}
method will be used by default when the number of observations is less than 
or equal to \code{maxExact}, otherwise \code{Prob} will be used.

Prob:  Iterations terminate when the estimated standard error of the estimated 
proportion p is less than p*Ca. The iteration continues until all sources and 
coefficients meet this criterion or until \code{maxIter} is reached. See Anscome(1953) 
for the origin of the criterion. 

SPR: This method uses sequential probability ratio tests to decide between 
the hypotheses \code{p0} and \code{p1} for a strength \code{(alpha, beta)} test. The test terminates 
upon the acceptance or rejection of p0 or if \code{maxIter} is reached. See Wald (1947).
The power of the SPR is beta at p0 and increases to 1-beta at p1. Placing \code{p0} and
\code{p1} close together makes the cut off sharp.  

Exact: This method generates all permutations of Y. It will generally be found 
too time consuming for more than 10 or 11 observations, but note that \code{\link{aovp}}
may be used to divide the data into small enough blocks for which exact 
permutation tests may be possible.

For \code{Prob} and \code{SPR}, one may set \code{nCycle} to unity to exchange all elements instead 
of just pairs at each iteration, but there seems to be no advantage to doing this
unless the number of iterations is small -- say less than 100.


The SS will be calculated \emph{sequentially}, just as \code{\link{lm}} does; or they may be
calculated \emph{uniquely}, which means that the SS for each source is calculated 
conditionally on all other sources. This is SAS type III, which is also what \code{drop1()} 
produces, except that \code{drop1()} will not drop main effects when interactions are present. 
The parameter \code{seqs} may be used to override the default unique calculation behavior.





}

\note{

There is a vignette with more details and an example. To access it, type 

vignette("lmPerm")

The default contrasts are set internally to \code{(contr.sum, contr.poly)}, which means
that coefficients are either pairwise contrasts with the last level or polynomial contrasts.

Numerical variables should be centered in order to make them orthogonal to the constant when
ANOVA is to be done.

Factors with many levels may cause problems for the methodology used by R. It is
designed to work well with both balanced and unbalanced data, but it is best used 
for factors with no more than four or five levels. If difficulties occur, degrees 
of freedom for high order sources will be reduced, and sometimes the sources will 
be omitted entirely. An error message usually accompanies such behavior. 

The variables inside the \code{Error()} term are treated in the order given. Thus 
\code{Error(A+B)} will usually produce components in the order A,B, with B orthogonalized
with respect to A. This may cause confusion for unbalanced block structures.

This function will behave identically to \code{aov()} if the following parameters are set:
\code{perm="", seq=TRUE, center=FALSE}.

}
\value{
	The usual output from \code{\link{aov}}, with permutation p-values instead of normal 
	theory p-values.
	\item{Iter}{For \code{Prob} and \code{SPR}: The number of iterations until the criterion 
		is met.}
	\item{Accept}{For \code{SPR}: 1 the p0 hypothesis is accepted, and 0 for rejection or when no decision before maxIter occurs.} 
}
\seealso{\code{\link{summary.aovp}}, \code{\link{lmp}}}

\examples{

## A simple randomized block example. 
# There are 7 blocks and 6 treatments. A first
# analysis with blocks as a factor shows block to be significant and treatments not. 

data(Hald17.4)
summary(aovp(Y~T+block,Hald17.4))

# Using the block to define a separate error strata tells a different story.

summary(aovp(Y~T+Error(block),Hald17.4))

# There appears to be a linear trend in the blocks. This may be investigated by
# extracting a linear component. The factor L was created by copying the block
# factor and assigning it a linear contrast, like this
# contrasts(L,1)<-contr.poly(7). The analysis then becomes.

summary(aovp(Y~T+L+Error(block),Hald17.4))

# The L factor is not significant under permutation. It is significant when aov()
# is used and the test is the usual F-test.


## From Venables and Ripley (2000)
# This is a 2^3 factorial in the variables N,P,K. It is fractioned by using the 
# three way interaction, NPK, into two fractions of 4. Each of these fractions is 
# allocated to 3 blocks, making 6 blocks in all. An analysis with block as a 
# variable is the following. As may be seen, aovp() discards the confounded NPK interaction.

data(NPK)
summary(aovp(yield ~ block + N*P*K, NPK))

# Since the NPK interaction was confounded with blocks, the experimenter no doubt judged 
# it of lesser interest. It may however be examined by including blocks as an additional
# error term as follows. The basic error level between blocks is of course larger than 
# that within blocks, so the NPK interaction would have to be substantially larger that
# it would have had to be were it tested within blocks.
summary(aovp(yield ~  N*P*K + Error(block), NPK))

# The SS calculated by aovp() are unique SS by default. That is,
# they are sums of squares for the difference of a model with and without the source. The
# resulting test is a test of the hypothesis that the source has no effect on the response.
# Sequential SS, which are those produced by aov() may be obtained by setting the
# parameter seqs=TRUE. simDesign is an unbalanced design created by the AlgDesign package.

data(simDesign)
summary(aovp(Y~.,simDesign))
summary(aovp(Y~.,simDesign,seqs=TRUE)) 

# Since there is only one stratum, these results are the same as would be obtained from
anova(lmp(Y~.,simDesign))

# ANOVA for numerical variables. First using contrasts, then numeric variables.
data(Federer276)
summary(aovp(Plants~Variety*Treatment+Error(Rep/Plot),Federer276))
data(Federer276Numeric)
summary(aovp(poly.formula(Plants~quad(Variety,Treatment)+Error(Rep/Plot)),Federer276Numeric))
# The coefficients and their p-values may be obtained by
summaryC(aovp(Plants~Variety*Treatment+Error(Rep/Plot),Federer276))


}


\references{
	\describe{
	\item{Hald, A. 1952}{Statistical theory with engineering applications. Wiley. Table 7.4}
	\item{Venables, W.N, and Ripley, B. D. (2000) }{Modern applied statistics with S-Plus, p176} 
	\item{Wald, A. (1947)}{Sequential analysis, Wiley, Sec. 5.3}
	}
}


\author{
	Bob Wheeler \email{bwheelerg@gmail.com}
}

\keyword{htest}

