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


	mt<-NULL
	for (i in seq(along=dots)) {
		mt<-cbind(mt,as.matrix(eval(dots[i][[1]])))
	}
	colnames(mt)<-nameargs(dots)
	mt
}
