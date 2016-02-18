poly.formula<-
function (frml) 
{
	# augments a formula with quad and cubic polynomials
    env <- environment(frml)
    nameargs <- function(...) {
        dots <- as.list(substitute(list(...)))[-1]
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


    quad <- function(...) {
		nms<-nameargs(...)
		nVars <- nargs()
        strg <- paste(paste("(", paste("X", 1:nVars, sep = "", 
            collapse = "+"), ")^2", sep = ""), "+", paste("I(X", 
            1:nVars, "^2)", sep = "", collapse = "+"))
        for (i in 1:nVars) {
            ag <- paste("X", i, sep = "")
            strg <- gsub(ag, replacement = nms[i], strg)
        }
        strg
    }
    cubic <- function(...) {
		nms<-nameargs(...)
		nVars <- nargs()
        strg <- paste(paste("(", paste("X", 1:nVars, sep = "", 
            collapse = "+"), ")^3", sep = ""), "+", paste("I(X", 
            1:nVars, "^2)", sep = "", collapse = "+"), "+", paste("I(X", 
            1:nVars, "^3)", sep = "", collapse = "+"))
        for (i in 1:nVars) {
            ag <- paste("X", i, sep = "")
            strg <- gsub(ag, replacement = nms[i], strg)
        }
        strg
    }
    cubicS <- function(...) {
		nms<-nameargs(...)
		nVars <- nargs()
        strg <- paste("(", paste("X", 1:nVars, sep = "", collapse = "+"), 
            ")^3", sep = "")
        for (i in 1:(nVars - 1)) {
            var <- paste("X", i, sep = "")
            strg <- paste(strg, "+", paste(paste("I(", var, paste("*X", 
                (i + 1):nVars, sep = ""), sep = ""), paste("*(", 
                var, paste("-X", (i + 1):nVars, "))", sep = "")), 
                collapse = "+"), sep = "", collapse = "+")
        }
        for (i in 1:nVars) {
            ag <- paste("X", i, sep = "")
            strg <- gsub(ag, replacement = nms[i], strg)
        }
        strg
    }
    findFunction <- function(name, string) {
        if (-1 == (strt <- regexpr(name, string))) 
            return(c(0, 0))
        head <- substr(string, 1, strt - 1)
        tail <- substr(string, strt, nchar(string))
        if (-1 == (fin <- regexpr(")", tail))) 
            return(c(0, 0))
        c(strt, strt + fin - 1)
    }
    frml <- deparse(frml, width.cutoff = 500)
    while ((0 != (pos <- findFunction("quad", frml))[1]) || (0 != 
        (pos <- findFunction("cubicS", frml))[1]) || (0 != (pos <- findFunction("cubic", 
        frml))[1])) {
        prog <- substr(frml, pos[1], pos[2])
        strHead <- substr(frml, 1, pos[1] - 1)
        strTail <- substr(frml, pos[2] + 1, nchar(frml))
        prog <- eval(parse(text = prog))
        frml <- paste(strHead, prog, strTail, sep = "")
    }
    frml <- as.formula(frml)
    environment(frml) <- env
    frml
}
