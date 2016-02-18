
summary.aovp<-
function(object, intercept = FALSE, split,
                        expand.split = TRUE, keep.zero.df = TRUE, ...)
{

    splitInteractions <- function(split, factors, names, asgn, df.names)
    {
        ns <- names(split)
        for(i in unique(asgn)) {
            if(i == 0 || names[i+1] %in% ns) next
            f <- rownames(factors)[factors[, i] > 0]
            sp <- f %in% ns
            if(any(sp)) {              # some marginal terms are split
                if(sum(sp) > 1) {
                    old <- split[ f[sp] ]
                    nn <- f[sp]
                    names(nn) <- nn
                    marg <- lapply(nn, function(x)
                                   df.names[asgn == (match(x, names) - 1)])
                    term.coefs <- strsplit(df.names[asgn == i], ":", fixed=TRUE)
                    ttc <- sapply(term.coefs, function(x) x[sp])
                    rownames(ttc) <- nn
                    splitnames <- apply(expand.grid(lapply(old, names)), 1,
                                        function(x) paste(x, collapse="."))
                    names(splitnames) <- splitnames
                    tmp <- sapply(nn, function(i)
                                  names(old[[i]])[match(ttc[i, ], marg[[i]])] )
                    tmp <- apply(tmp, 1, function(x) paste(x, collapse="."))
                    new <- lapply(splitnames, function(x) match(x, tmp))
                    split[[ names[i+1] ]] <-
                        new[sapply(new, function(x) length(x) > 0)]
                } else {
                    old <- split[[ f[sp] ]]
                    marg.coefs <- df.names[asgn == (match(f[sp], names) - 1)]
                    term.coefs <- strsplit(df.names[asgn == i], ":", fixed=TRUE)
                    ttc <- sapply(term.coefs, function(x) x[sp])
                    new <- lapply(old, function(x)
                                  seq(along=ttc)[ttc %in% marg.coefs[x]])
                    split[[ names[i+1] ]] <- new
                }
            }
        }
        split
    }
    asgn <- object$assign[object$qr$pivot[1:object$rank]]
    uasgn <- unique(asgn)
    nterms <- length(uasgn)
    effects <- object$effects
    if(!is.null(effects))
        effects <- as.matrix(effects)[seq(along=asgn),,drop=FALSE]
    rdf <- object$df.resid
	nmeffect <- c("(Intercept)", attr(object$terms, "term.labels"))
    coef <- as.matrix(object$coef)
    resid <- as.matrix(object$residuals)
    wt <- object$weights
    if(!is.null(wt)) resid <- resid * wt^0.5
    nresp <- NCOL(resid)
    ans <- vector("list", nresp)
    if(nresp > 1) {
        names(ans) <- character(nresp)
        for (y in 1:nresp) {
            cn <- colnames(resid)[y]
            if(is.null(cn) || cn == "") cn <- y
            names(ans)[y] <- paste(" Response", cn)
        }
    }

    if(!is.null(effects) && !missing(split)) {
        ns <- names(split)
        if(!is.null(Terms <- object$terms)) {
            if(!is.list(split))
                stop("the 'split' argument must be a list")
            if(!all(ns %in% nmeffect))
                stop("unknown name(s) in the 'split' list")
        }
        if(expand.split) {
            df.names <- names(coef(object))
            split <- splitInteractions(split, attr(Terms, "factors"),
                                       nmeffect, asgn, df.names)
            ns <- names(split)
        }
    }
    for (y in 1:nresp) {
        if(is.null(effects)) {
            nterms <- 0
            df <- ss <- ms <- numeric(0)
            nmrows <- character(0)
        } else {
            df <- ss <- numeric(0)
            nmrows <- character(0)
            for(i in seq(nterms)) {
                ai <- (asgn == uasgn[i])
                df <- c(df, sum(ai))
                ss <- c(ss, sum(effects[ai, y]^2))
                nmi <- nmeffect[1 + uasgn[i]]
                nmrows <- c(nmrows, nmi)
                if(!missing(split) && !is.na(int <- match(nmi, ns))) {
                    df <- c(df, unlist(lapply(split[[int]], length)))
                    if(is.null(nms <- names(split[[int]])))
                        nms <- paste("C", seq(along = split[[int]]), sep = "")
                    ss <- c(ss, unlist(lapply(split[[int]],
                                              function(i, e)
                                              sum(e[i]^2), effects[ai, y])))
                    nmrows <- c(nmrows, paste("  ", nmi, ": ", nms, sep = ""))
                }
            }
        }
 		doPerm<-!is.null(object$perm)                          #REW begin
		if (doPerm) perm<-object$perm$perm
		if (doPerm && (perm=="Prob" || perm=="SPR" || perm=="Exact")) {
			ms<-ifelse(df>0,ss/df,NA)
			if (rdf>0) {
				df <- c(df, rdf)
				rss<-sum(resid[, y]^2)
				ss <- c(ss, rss)	
				ms<-c(ms,rss/rdf)			
				nmrows <- c(nmrows,  "Residuals")
			}
			x <- list(Df = df, "R Sum Sq" = ss, "R Mean Sq" = ms)				
			class(x) <- c("anova", "data.frame")
			row.names(x) <- format(nmrows)
			pm<-pmatch("(Intercept)",row.names(x),0)
			const<-(!intercept && pm>0)


			if (object$perm$perm=="Prob") {
				Iter<-c(as.vector((object$perm$Mn)[,y]),if(rdf>0)NA)
				Prob<-c(as.vector((object$perm$P)[,y]),if(rdf>0)NA)
				x$"Iter"=Iter
				x$"Pr(Prob)"=Prob
			}
			else if (object$perm$perm=="SPR") {
				Iter<-c(as.vector((object$perm$Mn)[,y]),if(rdf>0)NA)
				Accept<-c(as.vector((object$perm$accept)[,y]),if(rdf>0)NA)
				Prob<-c(as.vector((object$perm$P)[,y]),if(rdf>0)NA)
				x$"Iter"=Iter
				x$"Pr(SPR)"=Prob
				x$"Accept"=Accept

			}
			else {
				Prob<-c(as.vector((object$perm$P)[,y]),if(rdf>0)NA)
				x$"Pr(Exact)"=Prob
			}

			if(!keep.zero.df) x <- x[df > 0, ]
			if (const) x<-x[-pm,]

		} 
		else {	
		   if(rdf > 0) {
				df <- c(df, rdf)
				ss <- c(ss, sum(resid[, y]^2))
				nmrows <- c(nmrows,  "Residuals")
			}
			nt <- length(df)
			ms <- ifelse(df > 0, ss/df, NA)
			x <- list(Df = df, "R Sum Sq" = ss, "R Mean Sq" = ms)

            #Residuals have no P-value
			if(rdf > 0 && nt>1) { # REW added nt>1 to suppress F value text
				TT <- ms/ms[nt]
				TP <- pf(TT, df, rdf, lower.tail = FALSE)
				TT[nt] <- TP[nt] <- NA
				x$"F value" <- TT
				x$"Pr(>F)" <- TP
			}

			class(x) <- c("anova", "data.frame")
			row.names(x) <- format(nmrows)
			if(!keep.zero.df) x <- x[df > 0, ]
			pm <- pmatch("(Intercept)", row.names(x), 0)
			if(!intercept && pm > 0) x <- x[-pm ,]
		}
                                                          #REW end
        ans[[y]] <- x
    }
    class(ans) <- c("summary.aovp", "listof")
    ans
}
