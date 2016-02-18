# single response
"summary.lmp" <- 
function (object, correlation = FALSE, symbolic.cor = FALSE, ...)
{
    z <- object
    p <- z$rank
    if (p == 0) {
        r <- z$residuals
        n <- length(r)
        w <- z$weights
        if (is.null(w)) {
            rss <- sum(r^2)
        } else {
            rss <- sum(w * r^2)
            r <- sqrt(w) * r
        }
        resvar <- rss/(n - p)
        ans <- z[c("call", "terms")]
        class(ans) <- "summary.lm"
        ans$aliased <- is.na(coef(object))  # used in print method
        ans$residuals <- r
        ans$df <- c(0, n, length(ans$aliased))
        ans$coefficients <- matrix(NA, 0, 4)
        dimnames(ans$coefficients)<-
            list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
        ans$sigma <- sqrt(resvar)
        ans$r.squared <- ans$adj.r.squared <- 0
        return(ans)
    }
    Qr <- object$qr
    if (is.null(z$terms) || is.null(Qr))
	stop("invalid \'lm\' object:  no 'terms' nor 'qr' component")
    n <- NROW(Qr$qr)
    rdf <- n - p
    if(is.na(z$df.residual) || rdf != z$df.residual)
        warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
    p1 <- 1:p
    ## do not want missing values substituted here
    r <- z$residuals
    f <- z$fitted
    w <- z$weights
    if (is.null(w)) {
        mss <- if (attr(z$terms, "intercept"))
            sum((f - mean(f))^2) else sum(f^2)
        rss <- sum(r^2)
    } else {
        mss <- if (attr(z$terms, "intercept")) {
            m <- sum(w * f /sum(w))
            sum(w * (f - m)^2)
        } else sum(w * f^2)
        rss <- sum(w * r^2)
        r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    se <- sqrt(diag(R) * resvar)
    est <- z$coefficients[ind<-Qr$pivot[p1]]         #REW
    tval <- est/se
    ans <- z[c("call", "terms")]
    ans$residuals <- r
	doPerm<-!is.null(object$perm)														#REW begin
	if (doPerm) perm<-object$perm$perm
	if (doPerm && (perm=="Prob" || perm=="SPR" || perm=="Exact")) {
		ans$perm<-perm
		if (perm=="Prob") {
			ans$coefficients<-cbind(est,object$perm$Mnt[ind], object$perm$Pt[ind])
			dimnames(ans$coefficients)<-list(names(z$coefficients)[ind],
				c("Estimate", "Iter", "Pr(Prob)"))
		}
		else if (perm=="SPR") {
			ans$coefficients<-cbind(est,object$perm$Mnt[ind], object$perm$Pt[ind],object$perm$acceptt[ind])
			dimnames(ans$coefficients)<-list(names(z$coefficients)[ind],
				c("Estimate", "Iter",  "Pr(SPR)", "Accept"))			
		}
		else {
			ans$coefficients<-cbind(est,object$perm$Pt[ind])
			dimnames(ans$coefficients)<-list(names(z$coefficients)[ind],
				c("Estimate","Pr(Exact)"))			
		}

	}
	else {
		ans$coefficients <-
		cbind(est, se, tval, 2*pt(abs(tval), rdf, lower.tail = FALSE))
		dimnames(ans$coefficients)<-
		list(names(z$coefficients)[Qr$pivot[p1]],
			c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
	}                                                                                  #REW end
    ans$aliased <- is.na(coef(object))  # used in print method
    ans$sigma <- sqrt(resvar)
    ans$df <- c(p, rdf, NCOL(Qr$qr))
    if (p != attr(z$terms, "intercept")) {
	df.int <- if (attr(z$terms, "intercept")) 1 else 0
	ans$r.squared <- mss/(mss + rss)
	ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int)/rdf)
	ans$fstatistic <- c(value = (mss/(p - df.int))/resvar,
			    numdf = p - df.int, dendf = rdf)
    } else ans$r.squared <- ans$adj.r.squared <- 0
    ans$cov.unscaled <- R
    dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1,1)]
    if (correlation) {
	ans$correlation <- (R * resvar)/outer(se, se)
	dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.lmp"
    ans
}

# Cycle through multiple responses
summary.mlmp<-
function(object, ...)
{
    coef <- coef(object)
    ny <- NCOL(coef)
## Not right:    if(is.null(ny)) return(NextMethod("summary"))
    effects <- object$effects
    resid <- residuals(object)
    fitted <- fitted(object)
    ynames <- colnames(coef)
									# REW begin
	if (is.null(object$perm))  
		doPerm<-FALSE
	else {
		perm=object$perm$perm
		doPerm<-(perm=="Prob" || perm=="SPR" || perm=="Exact")
	}
	if (doPerm){
		P<-object$perm$P
		Pt<-object$perm$Pt								
		if (perm!="Exact") {
			Mn<-object$perm$Mn
			Mnt<-object$perm$Mnt
		}
		if (perm=="SPR") {
			accept=object$perm$accept
			acceptt=object$perm$acceptt
		}  
	}                                            
										#REW end
    if(is.null(ynames)) {
	lhs <- object$terms[[2]]
	if(mode(lhs) == "call" && lhs[[1]] == "cbind")
	    ynames <- as.character(lhs)[-1]
	else ynames <- paste("Y", seq(ny), sep = "")
    }
    value <- vector("list", ny)
    names(value) <- paste("Response", ynames)
    cl <- oldClass(object)
    class(object) <- cl[match("mlmp", cl):length(cl)][-1]
    for(i in seq(ny)) {
		object$coefficients <- coef[, i]
		object$residuals <- resid[, i]
		object$fitted.values <- fitted[, i]
		object$effects <- effects[, i]
#		object$call$formula[[2]] <- object$terms[[2]] <- as.name(ynames[i]) # REW deleted
#					this is apparently only used for decoration.
		if (doPerm) {
			object$perm$P<-as.vector(P[,i])                     #REW begin
			object$perm$Pt<-as.vector(Pt[,i])
			if (perm!="Exact") {
				object$perm$Mn<-as.vector(Mn[,i])
				object$perm$Mnt<-as.vector(Mnt[,i])
			}
			if (perm=="SPR") {
				object$perm$accept<-as.vector(accept[,i])
				object$perm$acceptt<-as.vector(acceptt[,i])
			}  
		}                                      #REW end
		value[[i]] <- summary(object, ...)
    }
    class(value) <- "listof"
    value
}


"print.summary.lmp" <-
function (x, digits = max(3, getOption("digits") - 3),
              symbolic.cor = x$symbolic.cor,
	      signif.stars= getOption("show.signif.stars"),	...)
{
    cat("\nCall:\n")#S: ' ' instead of '\n'
    cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
    resid <- x$residuals
    df <- x$df
    rdf <- df[2]
    cat(if(!is.null(x$w) && diff(range(x$w))) "Weighted ",
        "Residuals:\n", sep="")
    if (rdf > 5) {
	nam <- c("Min", "1Q", "Median", "3Q", "Max")
	rq <- if (length(dim(resid)) == 2)
	    structure(apply(t(resid), 1, quantile),
		      dimnames = list(nam, dimnames(resid)[[2]]))
	else  structure(quantile(resid), names = nam)
	print(rq, digits = digits, ...)
    }
    else if (rdf > 0) {
	print(resid, digits = digits, ...)
    } else { # rdf == 0 : perfect fit!
	cat("ALL", df[1], "residuals are 0: no residual degrees of freedom!\n")
    }
    if (length(x$aliased) == 0) {
        cat("\nNo Coefficients\n")
    } else {
        if (nsingular <- df[3] - df[1])
            cat("\nCoefficients: (", nsingular,
                " not defined because of singularities)\n", sep = "")
        else cat("\nCoefficients:\n")
        coefs <- x$coefficients
        if(!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), length(colnames(coefs)), dimnames=list(cn, colnames(coefs))) #REW
            coefs[!aliased, ] <- x$coefficients
        }

		if (1==pmatch("(Intercept)",rownames(coefs),0)) {					#REW begin
			coefs<-coefs[-1,,drop=FALSE]						#Intercept must always be sig: keep matrix if only one row
		}
		
		if (!is.null(x$perm)) 
			isSPR<-x$perm=="SPR"
		else isSPR<-FALSE
		printCoefmat(coefs, digits=digits, signif.stars=signif.stars,      
				has.Pvalue = ifelse(isSPR,FALSE,TRUE),P.values=ifelse(isSPR,FALSE,TRUE), 
					na.print="NA", ...)								    #REW  end		
    }
    
    cat("\nResidual standard error:",
	format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")
    if (!is.null(x$fstatistic)) {
	cat("Multiple R-Squared:", formatC(x$r.squared, digits=digits))
	cat(",\tAdjusted R-squared:",formatC(x$adj.r.squared,digits=digits),
	    "\nF-statistic:", formatC(x$fstatistic[1], digits=digits),
	    "on", x$fstatistic[2], "and",
	    x$fstatistic[3], "DF,  p-value:",
	    format.pval(pf(x$fstatistic[1], x$fstatistic[2],
                           x$fstatistic[3], lower.tail = FALSE), digits=digits),
	    "\n")
    }
    correl <- x$correlation
    if (!is.null(correl)) {
	p <- NCOL(correl)
	if (p > 1) {
	    cat("\nCorrelation of Coefficients:\n")
	    if(is.logical(symbolic.cor) && symbolic.cor) {# NULL < 1.7.0 objects
		print(symnum(correl, abbr.colnames = NULL))
	    } else {
                correl <- format(round(correl, 2), nsmall = 2, digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop=FALSE], quote = FALSE)
            }
	}
    }
    cat("\n")#- not in S
    invisible(x)
}



"permute" <-
function (N=4,K=1,initialize=0) 
{
# First time, set initialize to 1. Subsequently, allow the default of zero to be used.
# result will be false for the last permutation set, and count will be the number of pairs in that set.

# Note, the pairs in vec are for indices 1:N, not 0:(N-1).
# N: number of items to permute
# K: number of pairs to return on each access
# count: the number of pairs acturally returned on an access
# result TRUE if not all permutations have been found. FALSE for the last set.
	result<-0
	vec<-rep(0,2*K)
	count<-0

	.C("permute",result=as.integer(result), N=as.integer(N), K=as.integer(K),vec=as.integer(vec),
			initialize=as.integer(initialize),count=as.integer(count),PACKAGE="lmPerm")
}






