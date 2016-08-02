multResp<-
function (...) {
# Makes a multicolumn response from variables in data

	# after Venables Ripley, S Programming
    nameargs <- function(dots) {
        nm <- names(dots)
        fixup <- if (is.null(nm)) 
            seq(along = dots)
        else nm == ""
        dep <- sapply(dots[fixup], function(x) deparse(x, width.cutoff = 500)[1])
        if (is.null(nm)) 
            nm <- dep
        else {
            nm[fixup] <- dep
        }
        nm
    }

	dots<-as.list(substitute(list(...)))[-1]

	## When evaluated inside another eval, uses that eval environment
	n=1
	if(sys.nframe()>1)
	for(i in 1:(sys.nframe()-1)){
	  if(sys.call(-i)[[1]]=="eval"){
	    n <- i
	    break
	  }
	}

	mt<-NULL
	for (i in seq(along=dots)) {
		mt<-cbind(mt,as.matrix(eval.parent(dots[i][[1]],n)))
	}
	colnames(mt)<-nameargs(dots)
	mt
}
