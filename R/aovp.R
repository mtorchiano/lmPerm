aovp<-
function (formula, data = NULL, perm="Exact", seqs=FALSE, center=TRUE,projections = FALSE, qr = TRUE, 
    contrasts = NULL, ...)									
{

#REW begin
# get.factors ***********************************
	get.factors<-function(data) {
		sapply(data,is.factor)
	}

# set.contrasts ********************************
	# sets contrasts with column labels of the form "#d"
	# this enables them to be found and removed in making sourcelabels
	# the contrasts are also changed to contr.sum()
	set.contrasts<-function(data){
		vc<-get.factors(data)
		if (any(vc)) {
			fac<-(1:length(vc))[vc]
			for (i in fac) {
				cont<-contrasts(data[,i])
				if (is.null(colnames(cont))) { # contr.poly will have named columns
					n<-ncol(cont)
					cont<-contr.sum(n+1)
					nms<-paste("#",1:n,sep="")
					colnames(cont)<-nms
					contrasts(data[,i])<-cont
				}
			}	

		}
	data
	}

#sourceVec ******************************
		# Scans sourceNames and returns a unique index for the variables contained in each source
	sourceVec<-function (sourceNames,varnames) {
		ns<-length(sourceNames)
		nv<-length(varnames)

		# now make the index
		ind<-rep(0,ns)
		val<-2^(1:nv-1)
		for (i in 1:ns) {
			nm<-sourceNames[i]
			for (j in 1:nv) {
				vv<-paste("(^|[^a-zA-Z])",varnames[j],"($|[^a-zA-Z])",sep="")
				if (-1!=regexpr(vv,nm)) {
					ind[i]<-ind[i]+val[j]
				}
			}
		}

		# Fix so indexes form a one step increasing sequence
		ord<-sort(ind,index.return=TRUE)$ix
		
		ind<-ind[ord]
		ind<-c(ind,0)!=c(0,ind)
		ind<-ind[1:ns]
		ind[ind>1]<-1
		ind<-cumsum(ind)

		ord<-sort(ord,index.return=TRUE)$ix # reverses the sort
		ind<-ind[ord]

		ind
	}

# getSourceNames **************************************
	# Pick out the shortest name for the source name
	getSourceNames<-function(dataVars,ind) {
		un<-unique(ind)
		un<-un[un>0] # Skip the intercept
		sourcenames<-vector("character",length(un))
		m<-1
		for (i in un) {
			value<-dataVars[i==ind]
			nc<-nchar(value)
			value<-value[min(nc)==nc]
			sourcenames[m]<-value[1]
			m<-m+1
		}
		sourcenames
	}


# cleanSource *************************************
	cleanSource<-function(sourcenames) {
		sourcenames<-gsub("#[0-9]*","",sourcenames) # remove the #d from factor lables
		sourcenames<-gsub("(\\.L|\\.Q|\\.C|\\^[0-9]*)($|:)","\\2",sourcenames) # remove contr.poly names
		sourcenames
	}

# cleanCoef **************************************
	cleanCoef<-function(coef) {
		if (is.vector(coef))
			coefNames<-names(coef)
		else
			coefNames<-rownames(coef)
		coefNames<-gsub("#([0-9]*)","\\1",coefNames)
		if (is.vector(coef)) 
			names(coef)<-coefNames
		else
			rownames(coef)<-coefNames
		coef
	}	

# begin program ************************************


	settings<-TRUE
	useF<-TRUE
	dots<-list(...) 
	if (length(dots)>0){
		if (!is.null(dots$settings)) settings<-dots$settings
		if (!is.null(dots$useF)) useF<-dots$useF
	}

	doPerm<-(perm=="Prob" || perm=="SPR" || perm=="Exact") 
	singTol<-1e-8
		# the following checks for a balanced design unless seqs is set


	# Make sure that contr.sum is used (ANOVA requires at least that contrast columns sum to zero.)
	# This is not actually used, because set.contrasts() is used on each factor
	opcons <- options("contrasts")
	if (missing(contrasts))									
		options(contrasts = c("contr.sum", "contr.poly"))   
    on.exit(options(opcons))


#REW end
    Terms <- if (missing(data)) 
        terms(formula, "Error")
    else terms(formula, "Error", data = data)

    indError <- attr(Terms, "specials")$Error

    if (length(indError) > 1) 
        stop(sprintf(ngettext(length(indError), "there are %d Error terms: only 1 is allowed", 
            "there are %d Error terms: only 1 is allowed"), length(indError)), 
            domain = NA)
    lmcall <- Call <- match.call()
		lmcall[[1]] <- as.name("lmp")     # REW
	  lmcall$singular.ok <- TRUE
    if (projections) 
        qr <- lmcall$qr <- TRUE
    lmcall$projections <- NULL
    if (is.null(indError)) {  # if there is no error term, lmp() is called
		lmcall$qr<-TRUE														#REW Need the qr for later use
        fit <- eval(lmcall, parent.frame())   		           
        if (projections) 
            fit$projections <- proj(fit)

		if (doPerm) {
			class(fit) <- if (inherits(fit, "mlmp")) 
				c("maov", "aovp", "aov", oldClass(fit))
			else c("aovp", "aov", oldClass(fit))
		}
		else {

			class(fit) <- if (inherits(fit, "mlm")) 
				c("maov", "aov", oldClass(fit))
			else c("aov", oldClass(fit))
		}
        fit$call <- Call
        return(fit)
    }
    else {
		if (pmatch("weights",names(list(...)),0L))
			stop("weights are not supported in a multistratum aovp() fit")
	# REW deleted
#        opcons <- options("contrasts")
#        options(contrasts = c("contr.helmert", "contr.poly")) 
#       on.exit(options(opcons))
		if (!missing(data)){ # REW This makes data avail for lhs multiple columns when multResp() is used
			attach(data,warn.conflicts=FALSE,name = "aovp")    # REW
      on.exit(detach(name = "aovp"))
		}
      
        allTerms <- Terms
        errorterm <- attr(Terms, "variables")[[1 + indError]]
        eTerm <- deparse(errorterm[[2]], width.cutoff = 500, backtick = TRUE) # the name of the error term
        intercept <- attr(Terms, "intercept")
        ecall <- lmcall
        ecall$formula <- as.formula(paste(deparse(formula[[2L]], 
          width.cutoff = 500, backtick = TRUE), "~", eTerm, if (!intercept) 
            "- 1"), env = environment(formula))
        ecall$method <- "qr"
        ecall$qr <- TRUE
        ecall$contrasts <- NULL
		ecall[[1]] <- as.name("lm")     # Do not need to call lmp() for error strata REW
       er.fit <- suppressWarnings(eval(ecall, parent.frame())) # fit to error strata. ecall is (response ~ block)
											# warnings suppressed because of perm parameter REW
 #       options(opcons) REW deleted
       nmstrata <- attr(terms(er.fit), "term.labels") # error strata names
        nmstrata <- sub("^`(.*)`$", "\\1", nmstrata) # removes backticks
		erStrataNames<-nmstrata # REW
        nmstrata <- c("(Intercept)", nmstrata)
        qr.e <- er.fit$qr
        rank.e <- er.fit$rank
      if (rank.e < length(er.fit$coef)) 
            warning("Error() model is singular")
        qty <- er.fit$resid
        maov <- is.matrix(qty)
        asgn.e <- er.fit$assign[qr.e$piv[1:rank.e]]
        maxasgn <- length(nmstrata) - 1
        nobs <- NROW(qty)
        if (nobs > rank.e) {
            result <- vector("list", maxasgn + 2) # empty list, to be filled with strata
            asgn.e[(rank.e + 1):nobs] <- maxasgn + 1 # extends asgn.e with a value marking within sources
            nmstrata <- c(nmstrata, "Within")
        }
        else result <- vector("list", maxasgn + 1)
        names(result) <- nmstrata # this will now have elements like "(Intercept)", strata name and "Within"

		lmcall$formula<-form <- update(formula, paste(". ~ .-", 
            deparse(errorterm, width.cutoff = 500, backtick = TRUE))) # Remove error term to get the within formula.
 
        Terms <- terms(form)

        lmcall$method <- "model.frame"

		if (!missing(data)){ # REW This makes data avail for lhs multiple columns when multResp() is used
			attach(data,warn.conflicts=FALSE,name="aovp")    # REW
        on.exit(detach(name = "aovp"))
		}
        mf <- eval(lmcall, parent.frame())  # The within model frame -- the call is to lmp()
#REW begin

		respVar<-deparse(formula[[2L]],width.cutoff=500, backtick=TRUE) # the response variable

		if (attr(attr(mf,"terms"),"intercept") ==0) {
			center<-FALSE
			if (missing(seqs))
				seqs<-FALSE
		}

		# center the numeric columns except for the response
		didCenter<-FALSE
		if (center) {
			numericColumn <- sapply(mf, is.numeric)
			numericColumn[respVar]<-FALSE
			if (any(numericColumn)) {
				didCenter<-TRUE
				means <- apply(mf[, numericColumn,drop=FALSE], 2, mean)
				mf[, numericColumn] <- sweep(mf[, numericColumn, 
					drop = FALSE], 2, means)
			}
		}	
		
	
		mf<-set.contrasts(mf) # use labeled contrasts so the factor names may be distinguished

#REW end
		xvars<-as.character(attr(Terms,"variables"))[-1L]
        if ((yvar <- attr(Terms, "response")) > 0) # remove duplicate responses ???
            xvars <- xvars[-yvar]
        if (length(xvars) > 0) {
            xlev <- lapply(mf[xvars], levels)
            xlev <- xlev[!sapply(xlev, is.null)]
        }
        else xlev <- NULL
        resp <- model.response(mf) # the numerical values of the response
        qtx <- model.matrix(Terms, mf, contrasts)

# begin REW
		# Create sourcenames which allow for numeric vars 
		formulaVars<-all.vars(formula)
		if (!is.na(match(".",formulaVars)))
			formulaVars<-colnames(mf)
		eTerm<-unlist(strsplit(eTerm," *[\\+|\\/|\\*] *")) # in case of multiple error terms

		formulaSources<-colnames(qtx)

		ind<-sourceVec(formulaSources,formulaVars)
		inds<-sort(ind,index.return=TRUE)$ix

		qtx<-qtx[,inds, drop=FALSE]  # rearrange both the x matrix and the sources
		ind<-ind[inds]
		attr(qtx,"assign")<-ind # Replace the old source ordering attribute

		dataVars<-colnames(qtx)  # the column order has changed
		sourcenames<-getSourceNames(dataVars,ind)
		sourcenames<-cleanSource(sourcenames)
# end REW
        cons <- attr(qtx, "contrasts")
        dnx <- colnames(qtx)
        asgn.t <- attr(qtx, "assign")
        if (length(wts <- model.weights(mf))) {
            wts <- sqrt(wts)
            resp <- resp * wts
            qtx <- qtx * wts
        }
        qty <- as.matrix(qr.qty(qr.e, resp)) # mult y by t(Q) from error fit 
        if ((nc <- ncol(qty)) > 1) {
            dny <- colnames(resp)
            if (is.null(dny)) 
                dny <- paste("Y", 1:nc, sep = "")
            dimnames(qty) <- list(seq(nrow(qty)), dny)
        }
        else dimnames(qty) <- list(seq(nrow(qty)), NULL) #(row indices, response col names if more than 1 column)
        qtx <- qr.qty(qr.e, qtx)   # mult x by t(Q) from error fit. The rows of qtx are coefficients of
								   # linear functionals estimating source effects. The first sources
								   # are effects within block sources. Block sources are in the order given 
								   # in the Error() term. Each source is orthogonal to all preceeding sources. 
								   # Orthogonality will be changed within groups of effects if seqs is FALSE
        dimnames(qtx) <- list(seq(nrow(qtx)), dnx)
        for (i in seq(along = nmstrata)) {
            select <- asgn.e == (i - 1)
            ni <- sum(select)
            if (!ni) 
                next
            xi <- qtx[select, , drop = FALSE]
            cols <- colSums(xi^2) > 1e-05
            if (any(cols)) {
                xi <- xi[, cols, drop = FALSE]
#REW begin
				attr(xi, "assign") <- sources<-asgn.t[cols]			  
				qtyi<-qty[select, , drop = FALSE]					  
				fiti <- lm.fit(xi, qtyi)
				sources<-sources[fiti$qr$pivot[1:fiti$rank]] # in case of singularity
				fiti$terms <- Terms
				dof <- unlist(lapply(split(sources,sources), length)) 
				if (useF) {
					Qr<-t(qr.Q(fiti$qr,complete=TRUE))
					Qt<-Qr[1:fiti$rank,,drop=FALSE]

					rk<-fiti$rank
					if (rk==NROW(Qr))
						Qr=NULL
					else 
						Qr<-Qr[(rk+1):NROW(Qr),]
					dof<-c(dof,NCOL(Qt)-sum(dof)) # add the residual df
				} else {
					Qt<-t(qr.Q(fiti$qr)[,1:fiti$rank])  # rank should equal the number of non aliased coefs 
					dof<-c(dof,0)
				}
				if (!seqs) {
					mxS<-max(sources)
					if (mxS>0) {
						R<-qr.R(fiti$qr)
						R<-R[1:fiti$rank,1:fiti$rank,drop=FALSE]			
						V<-solve(t(R)%*%R)          # inv(X'X)
						w<-lm.fit(R,Qt,singular.ok=TRUE) 
						Lc<-w$coefficients			# coef(b)=inv(R'R)(QR)'y = inv(R)Q'y.
													# The rows of Lc will be replaced with the vectors of the liner estimates of the sources
													# which are inv(sqrt(Vb))Lc[b,]. These will be used to calculate effectsto replace fiti$effects.
						if (is.vector(Lc))
							Lc<-as.matrix(Lc)
						for (j in 0:mxS) {  
							ind<-i==sources
							if (!any(ind))
							 next
							if (j==mxS) {
								Lc[ind,]<-Qt[ind,]  # the last source needs no correction
								ind<-(1:fiti$rank)[ind]						
							} else {
								ind<-(1:fiti$rank)[ind]
								Vb<-V[ind,ind]
								if (length(ind)==1) {
									T<-if (Vb>singTol) 1/sqrt(Vb) else 0
								}
								else {
									T<-eigen(Vb)
									Tv<-T$values
									nZero<-Tv>singTol
									Tv[nZero]<-1/sqrt(Tv[nZero])
									T<-(T$vectors)%*%(diag(Tv))%*%t(T$vectors)
								}
								Lc[ind,]<-T%*%Lc[ind,,drop=FALSE]
							}
							effects<-Lc[ind,,drop=FALSE]%*%qtyi # non-sequential effects
							if (is.matrix(fiti$effects))
								fiti$effects[ind,]<-effects
							else
								fiti$effects[ind]<-as.vector(effects)
						}
						Qt<-Lc
					} 
				}
				if (useF)
					Qt<-rbind(Qt,Qr)

				if (doPerm){ 
					fiti$perm<-permtst(Qt,qtyi,dof,perm, ...) 
				}
           }
			else {
				y <- qty[select, , drop = FALSE]
				fiti <- list(coefficients = numeric(0L), residuals = y, 
				  fitted.values = 0 * y, weights = wts, rank = 0L, 
				  df.residual = NROW(y))
			}
			if (projections) 
				fiti$projections <- proj(fiti)

			if (doPerm) {													
				class(fiti) <- c(if (maov) "maov", "aovp", "aov", oldClass(er.fit))
			}
			else {
				# REW added "aovp" here
				class(fiti) <- c(if (maov) "maov", "aovp", "aov", oldClass(er.fit))
			}
			if (!is.null(fiti$terms))
				attr(fiti$terms,"term.labels")<-sourcenames	
			fiti$coefficients<-cleanCoef(fiti$coefficients)	
#REW end
			result[[i]] <- fiti
        }
       result <- result[!sapply(result, is.null)] # keep only non null fits
		class(result) <- c("aovlist", "listof")
        if (qr) 
            attr(result, "error.qr") <- qr.e
        attr(result, "call") <- Call
        if (length(wts)) 
            attr(result, "weights") <- wts
        attr(result, "terms") <- allTerms
        attr(result, "contrasts") <- cons
        attr(result, "xlevels") <- xlev
		#REW begin
		if (settings)
			print(paste("Settings: ",if (seqs) "sequential SS" else "unique SS", if (didCenter) ": numeric variables centered"))
		#REW end
		result
    }
}
