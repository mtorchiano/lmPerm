lmp<-
function (formula, data, perm="Exact", seqs=FALSE,center=TRUE,subset, weights, na.action,   
		method = "qr", model = TRUE, x = FALSE, y = FALSE,
		qr = TRUE, singular.ok = TRUE, contrasts = NULL,
		offset, ...)   # REW
{

#REW begin
# get.factors ***********************************
	get.factors<-function(data) {
		sapply(data,is.factor)
	}


# set.contrasts ********************************
	# sets contrasts with column labels of the form "#d"
	set.contrasts<-function(data){
		vc<-get.factors(data)
		if (any(vc)) {
			fac<-(1:length(vc))[vc]
			for (i in fac) {
				cont<-contrasts(data[,i]) # contr.poly will have column names
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
	sourceVec<-function(sourceNames,varnames) {
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
		un<-un[un>0] # Delete the intercept
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
	opcons <- options("contrasts")
	if (missing(contrasts))									
		options(contrasts = c("contr.sum", "contr.poly"))   
    on.exit(options(opcons))

	
#REW end
    
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
#    mf$singular.ok <- mf$model <- mf$method <- NULL
#    mf$x <- mf$y <- mf$qr <- mf$contrasts <- mf$... <- NULL
    m <- match(c("formula", "data", "subset", "weights", "na.action",
                 "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    
    # MTk: using environment instead of attach!
    # Is it useful in any way?
    extenv = new.env(parent = parent.frame())
    attr(extenv,"name") <- "lmp ext env"
    if(!missing(data)){ 
      for(coln in names(data)){
        assign(coln,data[[coln]],envir=extenv)
      }
    }
    ##
    #  was
    ##
    # 		if (!missing(data)){ # REW This makes data avail for lhs multiple columns when multResp() is used
    # 			attach(data,warn.conflicts=FALSE,name = "aovp")    # REW
    #       on.exit(detach(name = "aovp"))
    # 		}    
    #mf <- eval(mf, parent.frame())  #this is the only use of data. Hereafter,mf is data with attachments
    if(!missing(data)){
      mf <- eval(mf, data,parent.frame())  #this is the only use of data. Hereafter,mf is data with attachments
    }else{
      mf <- eval(mf, parent.frame())
    }
#REW begin

	respVar<-deparse(formula[[2L]],width.cutoff =500, backtick=TRUE) # the response variable
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
    if (method == "model.frame")
	return(mf)
    else if (method != "qr")
	warning(gettextf("method = '%s' is not supported. Using 'qr'", method),
                domain = NA)
    mt <- attr(mf, "terms") # allow model.frame to update it
	termLabels<-attr(mt,"term.labels")  #REW
    y <- model.response(mf, "numeric")
    w <- model.weights(mf)
    offset <- model.offset(mf)
    if(!is.null(offset) && length(offset) != NROW(y))
	stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                     length(offset), NROW(y)), domain = NA)

    if (is.empty.model(mt)) {
		x <- NULL
		z <- list(coefficients = if (is.matrix(y))
                    matrix(,0,3) else numeric(0), residuals = y,
		fitted.values = 0 * y, weights = w, rank = 0,
		df.residual = if (is.matrix(y)) nrow(y) else length(y))
		if(!is.null(offset)) 
		z$fitted.values <- offset
		class(z) <- c(if(is.matrix(y) && ncol(y)>1) "mlm", "lm")
    }
    else {
		x <- model.matrix(mt, mf, contrasts)
# begin REW
		# Create sourcenames which allow for numeric vars 
		formulaVars<-all.vars(formula)
		if (!is.na(match(".",formulaVars))) # in case . used
			formulaVars<-colnames(mf)
		dataVars<-colnames(x)
		ind<-sourceVec(dataVars,formulaVars)
		inds<-sort(ind,index.return=TRUE)$ix
		x<-x[,inds, drop=FALSE]  # rearrange both the x matrix and the sources
		ind<-ind[inds]
		attr(x,"assign")<-ind # Replace the old source ordering attribute

		dataVars<-colnames(x)  # the column order has changed
		sourcenames<-getSourceNames(dataVars,ind)
		sourcenames<-cleanSource(sourcenames)
		z <- if(is.null(w)) 
			suppressWarnings(lm.fit(x, y, offset = offset,singular.ok=singular.ok, ...))  # suppressWarnings because of paramaters in ...
		else 
			suppressWarnings(lm.wfit(x, y, w, offset = offset,singular.ok=singular.ok, ...)) 

		# Modifying z$sources for non-sequantial SS, which are coef(b)' V(b)^-1 coef(b) where coef(b) indicates 
		# the lm coefficients corresponding to a source for columns b: ie z$coefficients[b]
		sources<-attr(x,"assign")
		sources<-sources[z$qr$pivot[1:z$rank]]		# The pivot colummns form a linearally independent set with rank z$rank
		dof <- unlist(lapply(split(sources,sources), length))  # degrees of freedom for the sources
		if (useF) {
			Qr<-t(qr.Q(z$qr,complete=TRUE))
			Qt<-Qr[1:z$rank,]
			rk<-z$rank
			if (rk==NROW(Qr))
				Qr=NULL
			else 
				Qr<-Qr[(rk+1):NROW(Qr),]
			dof<-c(dof,NCOL(Qt)-sum(dof)) # add the residual df
		} else {
			Qt<-t(qr.Q(z$qr)[,1:z$rank]) # rank should be equal to the number of non-aliased coefs
			dof<-c(dof,0)
		}

		if (!(seqs)) {
			mxS<-max(sources)
			if (mxS>1) {			
				R<-qr.R(z$qr)							
				R<-R[1:z$rank,1:z$rank,drop=FALSE]	
					
				V<-solve(t(R)%*%R)  # i.e. inv(X'X)

				w<-suppressWarnings(lm.fit(R,Qt,singular.ok=singular.ok) ) # regress Q' on R to get inv(R)Q'=Lc.
				Lc<-w$coefficients # coef(b)=inv(R'R)(QR)'y = inv(R)Q'y.
								  # The rows of Lc will be replaced with the vectors of the linear estimators of the sources,
								  # which are inv(sqrt(Vb))Lc[b,]. These will be used to calculate effects to replace those in z$effects

					# calculate the SS for each source
				if (is.vector(Lc))
					Lc<-as.matrix(Lc)
				for (i in 0:mxS) {
					ind<-i==sources
					if (!any(ind))
						next
					if (i==mxS) {
						Lc[ind,]<-Qt[ind,] # the last effect needs no correction.
						ind<-(1:z$rank)[ind]
					} else {
						ind<-(1:z$rank)[ind]
						Vb<-V[ind,ind]
						if (length(ind)==1) {
							T<-if (Vb>singTol) 1/sqrt(Vb) else 0
						}
						else {  # Get T=sqrt(Vb)
							T<-eigen(Vb)
							Tv<-T$values
							nZero<-Tv>singTol
							Tv[nZero]<-1/sqrt(Tv[nZero])
							T<-(T$vectors)%*%(diag(Tv))%*%t(T$vectors)
						}

						Lc[ind,]<-T%*%Lc[ind,,drop=FALSE]
					}

					effects<-Lc[ind,,drop=FALSE]%*%y     # non-sequential effects

						# replace the z$effects from the lm.fit()
					if (is.matrix(z$effects))			
						z$effects[ind,]<-effects
					else {
						z$effects[ind]<-as.vector(effects)
					}
				}
				Qt<-Lc
			} 
		}
		if (useF)
			Qt<-rbind(Qt,Qr)
		if (doPerm) {
			z$perm<-permtst(Qt,y,dof,perm, ...)  
 			class(z) <- c(if(is.matrix(y) && ncol(y)>1) "mlmp", "lmp", "lm") 
		} 
		else	# note, REW included mlmp and lmp here to force anova.lmp to be called instead of anova.lm 
				# in order to deal with numerical variables, if any, and to output a table for each col of y.
				# instead of the multivariate response from mlm.		                              							    
			class(z) <- c(if(is.matrix(y) && ncol(y)>1) "mlmp", "lmp", "lm")			 
#REW end
   }
#REW begin
		# remove "#d" from coefficient names
	coefficients<-z$coefficients
	multCol<-is.matrix(coefficients)	
	coefNames<-if (multCol) rownames(coefficients) else names(coefficients)
	coefNames<-gsub("#([0-9]*)","\\1",coefNames) # remove the # from coefficient names
	if (multCol) 
		rownames(coefficients)<-coefNames
	else
		names(coefficients)<-coefNames
	z$coefficients<-coefficients
#REW end
    z$na.action <- attr(mf, "na.action")
    z$offset <- offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
	attr(z$terms,"term.labels")<-sourcenames #REW replace term labels with source names
    if (model)
	z$model <- mf
    if (ret.x)
	z$x <- x
    if (ret.y)
	z$y <- y
    if (!qr) z$qr <- NULL
	#REW begin
	if (settings)
		print(paste("Settings: ",if (seqs) "sequential SS" else "unique SS", if (didCenter) ": numeric variables centered"))
	#REW end
    z
}
