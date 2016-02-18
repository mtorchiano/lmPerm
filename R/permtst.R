permtst<-
function (Q, Y, dof, perm="Prob", ... )
{
	maxIter<-5000
	Ca<-0.1
	p0<-0.05
	p1<-0.06
	alpha<-beta<-0.01
	maxExact<-10
	nCycle<-1000

	dots<-list(...) 

	if (length(dots)>0){
		if (!is.null(dots$maxIter)) maxIter<-dots$maxIter
		if (!is.null(dots$Ca)) Ca<-dots$Ca
		if (!is.null(dots$p0)) p0<-dots$p0 
		if (!is.null(dots$p1)) p1<-dots$p1 
		if (!is.null(dots$alpha)) alpha<-dots$alpha 
		if (!is.null(dots$beta)) maxIter<-dots$beta 
		if (!is.null(dots$maxExact)) maxExact<-dots$maxExact 
		if (!is.null(dots$nCycle)) nCycle<-dots$nCycle
	}	

# Given QR=X, where Q is symmetrical and  orthogonal, the SS are given by Q'Y,
# where Q is the first ncol(X) columns of the symmetric Q.  The argument Q is t(Q).

	nc<-1
	if (!is.vector(Y)) {
		nc<-ncol(Y)
	}

	N<-NCOL(Q)

    nS <- length(dof)
 	K<-sum(dof)
	nSe<-nS-1
	Km<-K-dof[nS]
	Mn <- rep(0, nc*nSe)
	P <- rep(0, nc*nSe)

	Pt<-rep(0,nc*Km)
	Mnt<-rep(0,nc*Km)
	totPerm<-0
	if (N>maxExact && perm=="Exact")
		perm<-"Prob"
	if (perm=="Prob") {
		value <- .C("permuteProb", N = as.integer(N), nc=as.integer(nc), Y = as.double(Y), 
			Q = as.double(Q),nS = as.integer(nS), dof = as.integer(dof), 
			Ca = as.double(Ca), maxIter = as.integer(maxIter), Mn = as.integer(Mn), 
			P = as.double(P), Mnt=as.integer(Mnt), Pt=as.double(Pt), nCycle=as.integer(nCycle), PACKAGE = "lmPerm")
			result<-list(perm=perm,P=matrix(value$P,nSe,nc),Mn=matrix(value$Mn,nSe,nc),
				Pt=matrix(value$Pt,Km,nc),Mnt=matrix(value$Mnt,Km,nc))
	}
	else if (perm=="SPR") {
	    accept <- rep(0, nc*nS)
		acceptt<-rep(0,nc*Km)
		value <- .C("permuteSPR", N = as.integer(N), nc=as.integer(nc), Y = as.double(Y), 
			Q = as.double(Q), nS = as.integer(nS), dof = as.integer(dof), 
			maxIter = as.integer(maxIter), Mn = as.integer(Mn), Mnt=as.integer(Mnt),
			p0 = as.double(p0), p1 = as.double(p1), alpha = as.double(alpha), beta = as.double(beta), 
			accept = as.integer(accept), P = as.double(P), acceptt=as.integer(acceptt),
			Pt=as.double(Pt), nCycle=as.integer(nCycle), PACKAGE = "lmPerm")
			accept<-temp<-(matrix(value$accept,nS,nc))[1:nSe,drop=FALSE]
			accept[temp==-1]<-0
			acceptt<-temp<-matrix(value$acceptt,Km,nc)
			acceptt[temp==-1]<-0
			result<-list(perm=perm,P=matrix(value$P,nSe,nc),Mn=matrix(value$Mn,nSe,nc),accept=accept,
				Pt=matrix(value$Pt,Km,nc),Mnt=matrix(value$Mnt,Km,nc),acceptt=acceptt)
	}
	else if (perm=="Exact") {
	    value <- .C("permuteExact", N = as.integer(N), nc=as.integer(nc),Y = as.double(Y), 
			Q= as.double(Q), nS = as.integer(nS), dof = as.integer(dof), 
			P = as.double(P), Pt=as.double(Pt), totPerm=as.double(totPerm),PACKAGE = "lmPerm")
			result<-list(perm=perm,P=matrix(value$P,nSe,nc),Pt=matrix(value$Pt,Km,nc),totPerm=value$totPerm)
	}
	
	result

}
