summaryC<-
function (object, ...) 
{
    if (!is.null(attr(object, "weights"))) 
        cat("Note: The results below are on the weighted scale\n")
    dots <- list(...)
    strata <- names(object)
    if (strata[1L] == "(Intercept)") {
        strata <- strata[-1L]
        object <- object[-1L]
    }
    x <- vector(length = length(strata), mode = "list")
    names(x) <- paste("Error:", strata)
    for (i in seq_along(strata)) {
		if (NROW(object[[i]]$coefficients)==0) next
		if (NCOL(object[[i]]$coefficients)==1) {
			class(object[[i]])<-"lmp"
		} else {
			class(object[[i]])<-c("mlmp","lmp")
		}
		x[[i]] <- do.call("summary",c(list(object = object[[i]]), dots))
	}
    class(x) <- "listof"
    x
}
