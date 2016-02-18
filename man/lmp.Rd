\name{lmp}
\alias{lmp}
\title{Fitting and testing linear models with permutation tests.}
\description{
\code{lmp} is \code{\link{lm}} modified to use permutation tests instead of normal 
theory tests. Like \code{lm}, it can be used to carry out regression, 
single stratum analysis of variance and analysis of covariance . Timing 
differences between \code{lmp} and \code{lm} are negligible.
}
\usage{
lmp(formula, data,  perm="Exact", seqs=FALSE, center=TRUE, subset, weights, na.action, 
	method = "qr", model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
   		singular.ok = TRUE, contrasts = NULL, offset, ...)

}

\arguments{
	The arguments are mostly the same as for \code{\link{lm}}. 
Additional parameters may be included. They are described in the 
"Additional Parameters" section below. These additional parameters
are the same as for \code{\link{aovp}}.

  \item{formula}{a symbolic description of the model to be fit.}
  \item{data}{an optional data frame containing the variables
    in the model.  If not found in \code{data}, the variables are taken
    from \code{environment(formula)}, typically the environment from which
    \code{lmp} is called.}
  \item{perm}{"Exact", "Prob", "SPR" will produce permutation probabilities. 
	Anyting else, such as "", will produce F-test probabilites.}
  \item{seqs}{If TRUE, will calculate sequential SS. If FALSE, unique SS.}
  \item{center}{If TRUE will center numerical variables}
  \item{subset}{an optional vector specifying a subset of observations
    to be used in the fitting process.}
  \item{weights}{an optional vector of weights to be used
    in the fitting process. If specified, weighted least squares is used
    with weights \code{weights} (that is, minimizing \code{sum(w*e^2)});
    otherwise ordinary least squares is used.}
  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset.  The \dQuote{factory-fresh}
    default is \code{\link{na.omit}}.  Another possible value is
    \code{NULL}, no action.  Value \code{\link{na.exclude}} can be useful.}
  \item{method}{the method to be used; for fitting, currently only
    \code{method = "qr"} is supported; \code{method = "model.frame"} returns
    the model frame (the same as with \code{model = TRUE}, see below).}
  \item{model, x, y, qr}{logicals.  If \code{TRUE} the corresponding
    components of the fit (the model frame, the model matrix, the
    response, the QR decomposition) are returned.}
  \item{singular.ok}{logical. If \code{FALSE} (the default in S but
    not in \R) a singular fit is an error.}
  \item{contrasts}{an optional list. See the \code{contrasts.arg}
    of \code{\link{model.matrix.default}}.}
  \item{offset}{this can be used to specify an \emph{a priori}
    known component to be included in the linear predictor
    during fitting.  An \code{\link{offset}} term can be included in the
    formula instead or as well, and if both are specified their sum is used.}
  \item{\dots}{additional arguments to be passed.}
}

\section{Additional parameters}{
	These are the same as for \code{\link{aovp}}. 

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
The usual regression model EY=Xb is assumed. The vector b is divided into sources 
with dfi degrees of freedom for the ith source, and \code{anova(lmp())} will produce an 
ANOVA table for these sources. Either permutation test p-values or the usual F-test 
p-values will be output. Polynomial model terms are collected into 
sources, so that \code{Y~A+B+I(A^2)} will contain two sources, one for A with 2 df,
and one for B with 1 df. Sources for factors are treated as usual, and polynomial
terms and factors may be mixed in one model. The function \code{\link{poly.formula}} may 
be used to create polynomial models, and the function \code{\link{multResp}} may
be used to create a multi-response matrix for the lhs from variables in \code{data}.

One may also use \code{summary(lm())} to obtain coefficient estimates and 
estimates of the permutation test p-values. The \code{Exact} method will permute
the values exactly. The \code{Prob} and \code{SPR} methods will approximate
the permutation distribution by randomly exchanging pairs of Y elements. The \code{Exact}
method will be used by default when the number of observations is less than 
or equal to \code{maxExact}, otherwise \code{Prob} will be used.

Prob:  Iterations terminate when the estimated standard error of the estimated 
proportion p is less than p*Ca. The iteration continues until all sources and 
coefficients meet this criterion or until \code{maxIter} is reached. See Anscome(1953) 
for the origin of the criterion. 

SPR: This method uses sequential probability ratio tests to decide between 
the hypotheses \code{p0} and \code{p1} for a strength \code{(alpha, beta)} test. The test terminates 
upon the acceptance or rejection of \code{p0} or if \code{maxIter} is reached. See Wald (1947).
The power of the SPR is beta at \code{p0} and increases to 1-beta at \code{p1}. Placing \code{p0} and
\code{p1} close together makes the cut off sharp.  

Exact: This method generates all permutations of Y. It will generally be found 
too time consuming for more than 10 or 11 observations, but note that \code{\link{aovp}}
may be used to divide the data into small enough blocks for which exact 
permutation tests may be possible.

For Prob and SPR, one may set \code{nCycle} to unity to exchange all elements instead 
of just pairs at each iteration, but there seems to be no advantage to doing this
unless the number of iterations is small -- say less than 100.

The SS will be calculated \emph{sequentially}, just as \code{lm()} does; or they may be
calculated \emph{uniquely}, which means that the SS for each source is calculated 
conditionally on all other sources. This is SAS type III, which is also what \code{drop1()}
produces, except that \code{drop1()} will not drop main effects when interactions are present. 
The parameter \code{seqs} may be used to override the default unique calculation behavior.




}

\note{
There is a vignette with more details and an example. To access it, type 

vignette("lmPerm")

The default contrasts are set internally to \code{(contr.sum, contr.poly)}, which means
that factor coefficients are either pairwise contrasts with the last level or polynomial contrasts.

Numerical variables should be centered in order to make them orthogonal to the constant when
ANOVA is to be done. 

This function will behave identically to \code{lm()} if the following parameters are set:
\code{perm="", seq=TRUE, center=FALSE}. An exception for multiple responses is that an 
ANOVA table for each response is output instead of a call to \code{anova.mlm()}.

}

\value{
	The usual output from \code{\link{lm}}, with permutation p-values or F-test 
	p-values. The p-values for the coefficients are of necessity, two-sided.
	\item{Iter}{For Prob and SPR: The number of iterations until the criterion 
		is met.}
	\item{Accept}{For SPR: 1 the p0 hypothesis is accepted, and 0 for rejection or when no decision before \code{maxIter} occurs.} 
	}



\examples{
# 3x3 factorial with ordered factors, each is average of 12. 
# This is a saturated design with no df for error. The results tend to support 
# Cochran and Cox who used a guessed residual SS for their analysis. The design
# is balanced, so the sequential SS are the same as the unique SS. 
data(CC164)
summary(lmp(y ~ N * P, data = CC164, perm="")) # F-value output as if lm() was used.
summary(lmp(y ~ N * P, data = CC164,)) # Default, using "Exact" if possible.
summary(lmp(y ~ N * P, data = CC164, perm="SPR"))
anova(lmp(y ~ N * P, data = CC164))

# A two level factorial. The artificial data is N(0,1) with an effect of
# 1.5 added to factor X4. When the number of iterations are small, as in
# this case, using nCycle=1 is advantageous. 
X<-expand.grid(X1=1:2,X2=1:2,X3=1:2,X4=1:2)
X$Y<-c(0.99,1.34,0.88,1.94,0.63,0.29,-0.78,-0.89,0.43,-0.03,0.50,1.66,1.65,1.78,1.31,1.51)
summary(lmp(Y~(X1+X2+X3+X4)^2,X,"SP"))  # The prob method is used because "SP" is not recognized.
summary(lmp(Y~(X1+X2+X3+X4)^2,X,"SPR"))
summary(lmp(Y~(X1+X2+X3+X4)^2,X,"SPR",nCycle=1)) #An additional parameter being passed.

# A saturated design with 15 variables in 16 runs. The orginal analysis by Quinlan pooled the mean
# squares from the 7 smallest effcts and found many variables to be significant. Box, reanalyzed
# the data using half-normal plots and found only variables E and G to be important. The permutation
# analysis agrees with this conclusion.
data(Quinlan)
summary(lmp(SN~.,Quinlan))

# A design containing both a polynomial variable and a factor
data(simDesignPartNumeric)
anova(lmp(poly.formula(Y~quad(A,B)+C),simDesignPartNumeric))


}

\seealso{\code{\link{summary.lmp}}, \code{\link{aovp}}}

\references{
	\describe{
	\item{Chochran, W. and Cox, G. (1957) p164}{Experimental Design, 2nd Ed.}
		John Wiley & Sons, New York. 
	\item{Wald, A. (1947)}{Sequential analysis, Wiley, Sec. 5.3}
	\item{Quinlan, J. (1985)}{"Product improvement by application of Taguchi methods." in American Supplier Institute News
		(special symposium ed.) Dearborn, MI. American Supplier Institute. 11-16.}
	\item{Box, G. (1988)}{Signal-to-noise ratios, performance criteria, and transformations. Technometics. 30-1. 1-17.}
	}
}



\author{
	Bob Wheeler \email{rwheeler@echip.com}
}

\keyword{htest}