#include "wheeler.h"
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>

void permute(int* result, int* N, int* K, int* vec, int* initialize, int* count);
void permuteSPR(int* Ni, int* nci, double* Y, double* Q, int* nSi, int* dof, int* Mi, 
			int* Mn, int* Mnt,double* p0,	double* p1, double* alpha, double* beta, int* accept, 
			double* P, int* acceptt, double* Pt, int* nCyclei);
void permuteExact(int* Ni, int* nci, double* Y, double* Q, int* nSi, int* dof, double* P,
				  double* Pt,double* totPermi);
void permuteProb(int* Ni, int* nci, double* Y, double* Q, int* nSi, int* dof, double* ci,int* Mi, 
			int* Mn, double* P,int *Mnt, double* Pt, int* nCyclei);
void	permuteRand(double** Y,int n);


#define tol 1e-8



#define N (*Ni)
#define nS (*nSi)
#define M (*Mi)
#define c (*ci)
#define nc (*nci)
#define nCycle (*nCyclei)
#define totPerm (*totPermi)

#define IndQ(A,B) B+A*K

#define iterMin 50   /* Minimum number of pairwise iterations before a decision may be made */

/********************
 Permutation generator with minimal changes

 p() the permutation array 1 to n+2
 pi() its inverse, identical to p() initially, 1 to n+2
 d() array of change directions, 1 to n
 changed() the pair changed, 1 to 2

 Initialize:
 p(i)=pi(i)=i
 d(i)=-1
 d(1)=0, p(1)=p(n+2)=m=n+1

 Repeat permute() while m != 1
 Initially, and at each repeat, p will contain a new
 permutation the differs by single pairwise interchange
 from the previous one.

 See Reingold, E.M., Nievergelt,J, and Deo, N. (1977)
 Combinatorial Algorithms Theory and Practice,
 Prentice Hall page 170.

 *******************/

bool DllMain(void)					 
{

    return true;
}



/***********
* Returns true if new permutations are generated, and false at the end.
* initialize should be 1 on the first access, and 0 susequently
***********/
void permute(
	int* result,  /* return value */
	int* Ni,        /* number of items to permute */
	int* K,
	int* vec,		/* dim 2*K, for K pairs at a time */
	int* initialize,    /* 1 initialize; 0 continue*/
	int* count     /* number of pairs actually returned*/
)
{
	int i;
	int temp;
	int a;
	int b;
	int n=(*Ni);
	int k=(*K);
	int lcount=0;
	static int p[100];
	static int pi[100];
	static int d[100];
	static int m;


	if (k<1 || n<2) {
		*result=false;
		return;
	}

	/* Initialize */
	/* Note: This original permutation (1 to n) is not returned.*/
	if (*initialize==1) {
		for (i=1;i<=n;i++) {
			p[i]=i;
			pi[i]=i;
			d[i]=-1;
		}
		p[0]=p[n+1]=(n+1);
		d[1]=0;
		m=n+1;

		*initialize=0;  /* Reset initialize to start making permutations */
	}


	do {
		if (m!=1) {	

			m=n;
   
			while  (p[pi[m]+d[m]]>m) { /* m is at a limit, try moving m-1 */
				d[m]=-d[m];
				m=m-1;
			}


    
			a=pi[m];
			b=pi[m] + d[m];
    
			temp=p[a];
			p[a]=p[b];
			p[b]=temp;

				/* p[b] is equal to m */
    
			pi[m]=pi[p[a]];
			pi[p[a]]=a;

		}
		else {  /* all permutations have been found */
			goto finishUp;
		}

		if (m!=1) {
			vec[lcount++]=a;
			vec[lcount++]=b;
		}

	}while(--k);

finishUp:
    *result=(m==1)?false:true;
	(*count)=lcount/2; /* Number of pairs */
	return;

	
}




/************ 
* Exact Permutation test using all permutations.
* Assuming the usual regression model of Y=Zt, where the K parameter vector t is divided 
* into nS sources. The statistic for usual test of the significance of a source may be
* shown to be SS=Y'QQ'Y, where the NxK matrix Q is a function of Z only. This function 
* finds the permutation distribution of SS by exchanging all pairs of elements of Y. 
***********/


void permuteExact(
	int* Ni,	   /* scalar number of observations in Y*/
	int* nci,     /* scalar number of columns in Y */
	double* Y,	   /* nc*N matrix of responses */
	double* Q,     /* N*K matrix of linear functionals in row major form */
	int* nSi,      /* scalar number of sources into which the functionals are divided */
	int* dof,     /* nS vector of df for the sources: K=sum(df) */
	double* P,      /* nc*nS vector of observed p values at termination*/
    double* Pt,     /* a nc*K vector of observed p values for the terms at termination*/
	double* totPermi   /* total number of permutatons */
)
{
	double* SS;  /* nc*nS matrix  of current sum of squares */
	double* oSS;  /* nc*nS matrix of Sum of squares of observations */
	double* b;   /* nc*K matrix of term estimates */
	double* obSq;  /* nc*K matrix of square of original estimates */
	int* pairs; /* vector of Kp permuted pairs, Note it is 2*Kp in length.*/
	int* dVec;   /* nc*nSe matrix of number of source exceedences */
	int* dVect;   /* nc*Kt matrix of number of term exceedences */
	double* rss; /* residual sum of squares */
	int K=0;       /* total number of terms */
	int Kt=0;	 /* number of non residual terms */
	int Kp=5000;
	int i;
	int j;
	int k;
	int u1;
	int u2;
	int iter;
	double temp;
	double y;
	double x;
	double y1;
	double y2;
	double x1;
	double x2;
	double delta;
	int count=0;
	int result=0;
	int initialize=1;
	bool useF=true;  /* Use F ratio instead of SS */
	int nSe=nS-1;   /* Number of effects */


	int l;
	double* Yp;
	double* bp;
	double* oSSp;
	double* obSqp;
	double* SSp;
	int* dVecp;
	int* dVectp;
	double* Pp;
	double* Ptp;

	totPerm=1; /* Allowing for initial sum of squares */

	for (i=0;i<nS;i++) /* find number of columns in Q */
		K+=dof[i];
	if (dof[nSe]==0) {
		useF=false;
	}

	Kt=K-dof[nSe];
	SS=(double *)S_alloc(nc*nS,sizeof(double)); /* if !useF, unused cell on end for do loops */
	oSS=(double *)S_alloc(nc*nS,sizeof(double)); /* if !useF, unused cell on end for do loops */
	dVec=(int *)S_alloc(nc*nSe,sizeof(int));
	dVect=(int *)S_alloc(nc*Kt,sizeof(int));
	b=(double *)S_alloc(nc*K,sizeof(double));
	obSq=(double *)S_alloc(nc*K,sizeof(double));
	pairs=(int *)S_alloc(2*Kp,sizeof(int));
	rss=(double *)S_alloc(nc,sizeof(double));



	/* Calculate b */
	for (l=0;l<nc;l++) {
		Yp=Y+N*l;
		bp=b+K*l;
		k=0;
		for (i=0;i<N;i++) {
			y=Yp[i];
			for (j=0;j<K;j++) {
				bp[j]+=y*Q[k++];  /* Note Q is row major */
			}	
		}
	}


	/* Get oSS and obSq */
	for (l=0;l<nc;l++) {
		k=0;
		oSSp=oSS+nS*l;
		obSqp=obSq+K*l;
		bp=b+K*l;
		for (i=0;i<nS;i++) {
			oSSp[i]=0;
			for (j=0;j<dof[i];j++) {   /* if useF=false, dof[nSe] is zero */
				x=bp[k];				
				oSSp[i]+=obSqp[k++]=x*x;
			}
		}
	}

	/* Scale oSSp and obSqp by rss */
	if (useF) {
		for (l=0;l<nc;l++) {
			oSSp=oSS+nS*l;
			obSqp=obSq+K*l;
			rss[l]=oSSp[nSe];
			if (rss[l]>tol) {
				for (i=0;i<nSe;i++) {
					oSSp[i]/=rss[l];
				}
				for (i=0;i<K;i++) {
					obSqp[i]/=rss[l];
				}
			} else {
				useF=false;
				rss[l]=1;
			}
		}
	} else {
		for (l=0;l<nc;l++)
			rss[l]=1;
	}


		/* initial SS */
	memcpy(SS,oSS,nc*nS*sizeof(double));


	initialize=1; /* insurance */

	do {
		permute(&result,Ni,&Kp,pairs,&initialize,&count);

		totPerm+=count;

		/* Iterate over all permutations */
		for (iter=0;iter<count;iter++) {
			/* Get two indices */
			u1=(int)pairs[2*iter]-1; 
			u2=(int)pairs[2*iter+1]-1;

			/* Update the SS */
			for (l=0;l<nc;l++) {
				Yp=Y+N*l;
				SSp=SS+nS*l;
				oSSp=oSS+nS*l;
				bp=b+K*l;
				obSqp=obSq+K*l;
				dVecp=dVec+nSe*l;
				dVectp=dVect+Kt*l;

				y1=Yp[u1];
				y2=Yp[u2];
				delta=y1-y2;
				k=0;
				for (i=0;i<nS;i++) {
					for (j=0;j<dof[i];j++) {  /* dof[nSe] is zero for !useF*/
						SSp[i]-=bp[k]*bp[k]; /* Remove k term */
						x1=Q[IndQ(u1,k)];
						x2=Q[IndQ(u2,k)];
						bp[k]+=delta*(x2-x1); /* Update k term with exchanged observations*/
						SSp[i]+=bp[k]*bp[k];
						k++;
					}
				}
				/* Update Y */
				temp=Yp[u1];
				Yp[u1]=Yp[u2];
				Yp[u2]=temp;

				if (useF)
					rss[l]=SSp[nSe];

				/* Tabulates dVec and dVect */
				k=0;
				for (i=0;i<nSe;i++) {
					dVecp[i]+=(SSp[i]+tol>oSSp[i]*rss[l]);
					for (j=0;j<dof[i];j++) {
						dVectp[k]+=(bp[k]*bp[k]+tol>obSqp[k]*rss[l]);
						k++;
					}
				}
			}
		}
	}while(result);


	for (l=0;l<nc;l++) {
		Pp=P+nSe*l;
		Ptp=Pt+Kt*l;
		dVecp=dVec+nSe*l;
		dVectp=dVect+Kt*l;

		for (i=0;i<nSe;i++) {
			Pp[i]=((double)dVecp[i]+1)/totPerm; /* Plus 1 for initial sum of squares */
		}
		for (i=0;i<Kt;i++)
			Ptp[i]=((double)dVectp[i]+1)/totPerm; /* Plus 1 for initial sum of squares */

	}

}





/* See Anscombe, F.J. (1953). Sequential estimation. J.R. Statist. Soc. B 15, 1-29. */

/************ 
* Permutation test using fraction of p to stop.
* Assuming the usual regression model of Y=Zt, where the K parameter vector t is divided 
* into nS sources. The statistic for usual test of the significance of a source may be
* shown to be SS=Y'QQ'Y, where the NxK matrix Q is a function of Z only. This function 
* approximates the permutation distribution of SS by randomly exchanging pairs of elements
* of Y. 
*
* The interation continues until the estimated standard error is within a fraction, c,
* of the estimated p; i.e., until (d/m)(1-d/m)/m<=(cd/m)^2, where d is the number of times
* the permutated statistic exceeds the original statistic in a sample of m successive permutations.
* Algebraic manipulation leads to the criterion: p>=1/(mc^2+1), where p=d/m.
*
* Iterations are continued until all sources and terms reach this goal or M is exceeded.
* Estimates of p for both sources and terms are returned. (A sources iterations are continued
* untill both the source and all its terms reach the goal.
***********/



void permuteProb(
	int* Ni,	   /* scalar number of observations in Y*/
	int* nci,     /* scalar Number of Y columns */
	double* Y,	   /* nc*N matrix of responses */
	double* Q,     /* NxK matrix of linear functionals in row major form */
	int* nSi,      /* scalar number of sources into which the functionals are divided */
	int* dof,     /* nS vector of df for the sources: K=sum(df) */
	double* ci,    /* scalar proportion */
	int* Mi,       /* scalar maximum limit on the number of iterations*/
	int* Mn,      /* nc*nSe matrix of sample sizes at termination */
	double* P,      /* nc*nS matrix of observed p values for the sources at termination*/
	int* Mnt,     /* nc*K matrix of sample sizes for the terms at termination */
    double* Pt,     /* nc*K vector of observed p values for the terms at termination*/
	int* nCyclei    /* complete permutation every nCycle*N iterations*/
)
{
	double* SS;  /* nc*nS matrix of current sum of squares */
	double* oSS; /* nc*nS matrix of Original SS */
	double* b;   /* nc*K matrix of term estimates */
	double* obSq;  /* nc*K matrix of square of original estimates */
	int* d;   /* nc*nS matrix of number of exceedences */
	int* dt;   /* nc*K matrix of number of exceedences for the terms */
	bool* accept; /* nc*nS matrix True when the corresponding source and all its terms have met goal*/
	int* acceptC; /* nc*Nx marrix Count of the number of terms in a source that have met goal */
	bool* acceptt; /* nc*K True when the corresponding term has met goal */
	double* rss;  /* Residual sum of squares */
	double crit;   /* 1/(iter*Cs+1) */
	int K=0;       /* total number of terms */
	int Kt=0;      /* total non residual terms */
	double cs=c*c;
	double p;
	int i;
	int j;
	int k;
	int u1;
	int u2;
	int iter;
	double temp;
	double y;
	double x;
	double y1;
	double y2;
	double x1;
	double x2;
	double delta;
	bool nToDo; /* number of sources not decided */
	bool useF=true;  /* Use F ratio instead of SS */
	int nSe=nS-1;   /* Number of effects */

	int l;
	double* Yp;
	double* bp;
	double* oSSp;
	double* obSqp;
	double* SSp;
	int* dp;
	int* dtp;
	double* Pp;
	double* Ptp;
	bool* acceptp;
	int* acceptCp;
	bool* accepttp;
	int* Mnp;
	int* Mntp;


	for (i=0;i<nS;i++) /* find number of columns in Q */
		K+=dof[i];

	if (dof[nSe]==0) { /* No residuals */
		useF=false;
	}
	
	Kt=K-dof[nSe];

	nToDo=nc*(nSe+Kt);

	SS=(double *)S_alloc(nc*nS,sizeof(double));  /* These should be zeroed as allocated */
	oSS=(double *)S_alloc(nc*nS,sizeof(double));
	d=(int *)S_alloc(nc*nSe,sizeof(int));
	accept=(bool *)S_alloc(nc*nS,sizeof(bool));
	b=(double *)S_alloc(nc*K,sizeof(double));
	obSq=(double *)S_alloc(nc*K,sizeof(double));
	dt=(int *)S_alloc(nc*Kt,sizeof(int));
	acceptC=(int *)S_alloc(nc*nS,sizeof(int));
	acceptt=(bool *)S_alloc(nc*Kt,sizeof(bool));
	rss=(double *)S_alloc(nc,sizeof(double));

	/* Calculate b */
	for (l=0;l<nc;l++) {
		Yp=Y+N*l;
		bp=b+K*l;
		k=0;
		for (i=0;i<N;i++) {
			y=Yp[i];
			for (j=0;j<K;j++) {
				bp[j]+=y*Q[k++];  /* Note Q is row major */
			}	
		}
	}

	/* Get oSS and obSq */
	for (l=0;l<nc;l++) {
		k=0;
		oSSp=oSS+nS*l;
		obSqp=obSq+K*l;
		bp=b+K*l;
		for (i=0;i<nS;i++) {
			oSSp[i]=0;
			for (j=0;j<dof[i];j++) { /* dof[nSe] may be zero */
				x=bp[k];
				oSSp[i]+=obSqp[k++]=x*x;
			}
		}
	}

	/* Scale oSSp and obSqp by rss */
	if (useF) {
		for (l=0;l<nc;l++) {
			oSSp=oSS+nS*l;
			obSqp=obSq+K*l;
			rss[l]=oSSp[nSe];
			if (rss[l]>tol) {
				for (i=0;i<nSe;i++) {
					oSSp[i]/=rss[l];
				}
				for (i=0;i<K;i++) {
					obSqp[i]/=rss[l];
				}
			} else {
				useF=false;
				rss[l]=1;
			}
		}
	} else {
		for (l=0;l<nc;l++)
			rss[l]=1;
	}

		/* initial SS */
	memcpy(SS,oSS,nc*nS*sizeof(double));


	GetRNGstate();

	/* Iterate until M or all sources decided */
	for (iter=1;iter<=M && nToDo>0;iter++) {
		/* Get two random indices */
		u1=(int)((double)N*unif_rand()); 
		u2=(int)((double)N*unif_rand());

		/* Update the SS */
		for (l=0;l<nc;l++) {

			Yp=Y+N*l;
			SSp=SS+nS*l;
			oSSp=oSS+nS*l;
			bp=b+K*l;
			obSqp=obSq+K*l;
			dp=d+nSe*l;
			dtp=dt+Kt*l;
			acceptp=accept+nS*l;
			acceptCp=acceptC+nS*l;
			accepttp=acceptt+Kt*l;
			Pp=P+nSe*l;
			Ptp=Pt+Kt*l;
			Mnp=Mn+nSe*l;
			Mntp=Mnt+Kt*l;

			if (iter==1 || 0==iter%(nCycle*N)) { /* permute everything at start and then once in nCycle*N cycles */
				permuteRand(&Yp,N);
				memset(bp,0,K*sizeof(double));
				k=0;
				for (i=0;i<N;i++) {
					y=Yp[i];
					for (j=0;j<K;j++) {
						bp[j]+=y*Q[k++];  
					}	
				}
				k=0;
				for (i=0;i<nS;i++) {
					SSp[i]=0;
					for (j=0;j<dof[i];j++) {
						x=bp[k++];
						SSp[i]+=x*x;
					}
				}
			}
			else {
				y1=Yp[u1];
				y2=Yp[u2];
				delta=y1-y2;
				k=0;
				for (i=0;i<nS;i++) {
					if (!acceptp[i]) {  /* acceptp[nSe] will always be 0 */
						for (j=0;j<dof[i];j++) {
							SSp[i]-=bp[k]*bp[k]; /* Remove k term */
							x1=Q[IndQ(u1,k)];
							x2=Q[IndQ(u2,k)];
							bp[k]+=delta*(x2-x1); /* Update k term with exchanged observations*/
							SSp[i]+=bp[k]*bp[k];
							k++;
						}
					}
					else
						k+=dof[i];
				}
				/* Update Y */
				temp=Yp[u1];
				Yp[u1]=Yp[u2];
				Yp[u2]=temp;
			}

			if (useF)
				rss[l]=SSp[nSe];

			/* accept or continue */
			crit=1.0/((double)iter*cs+1.0);
			k=0;
			for (i=0;i<nSe;i++) {
				for (j=0;j<dof[i];j++) {
					if (!accepttp[k] ) {
						dtp[k]+=(bp[k]*bp[k]+tol)>obSqp[k]*rss[l];
						p=(double)dtp[k]/(double)iter;
						if (iter>iterMin &&  (p>=crit)) {
							accepttp[k]=true;
							acceptCp[i]++;
							Ptp[k]=p;
							Mntp[k]=iter;
							nToDo--;
						}
					}
					k++;
				}
				if (!acceptp[i] ) { 
					dp[i]+=(SSp[i]+tol)>oSSp[i]*rss[l];
					p=(double)dp[i]/(double)iter;
					if (iter>iterMin &&  (p>=crit)) {
						if (acceptCp[i]>=dof[i]) {  /* acceptCp[nSe] is never checked */
							acceptp[i]=true;
							Pp[i]=p;
							Mnp[i]=iter;
							nToDo--;
						}
					}
				}
			}
		}
	}

	for (l=0;l<nc;l++){
		acceptp=accept+nS*l;
		Pp=P+nSe*l;
		Mnp=Mn+nSe*l;
		accepttp=acceptt+Kt*l;
		Ptp=Pt+Kt*l;
		Mntp=Mnt+Kt*l;
		dp=d+nSe*l;
		dtp=dt+Kt*l;


		k=0;
		for (i=0;i<nSe;i++) {
			if (!acceptp[i]) {
				Pp[i]=(double)dp[i]/(double)M;
				Mnp[i]=M;
			}
			for (j=0;j<dof[i];j++) {
				if (!accepttp[k]) {
					Ptp[k]=(double)dtp[k]/(double)M;
					Mntp[k]=M;
				}
				k++;
			}
		}
	}
	PutRNGstate();

}


/* See Wald, A. (1947). Sequential analysis, Wiley, Section 5.3 */

/************ 
* Permutation test using SPR to decide significance.
* Assuming the usual regression model of Y=Zt, where the K parameter vector t is divided 
* into nS sources. The statistic for usual test of the significance of a source may be
* shown to be SS=Y'QQ'Y, where the NxK matrix Q is a function of Z only. This function 
* approximates the permutation distribution of SS by randomly exchanging pairs of elements
* of Y. A SPR test of strngth (alpha, beta) is used to decide between two hypotheses, p0 and p1
* for the significance level of the test. The decision is output in the vector accept.
* This function terminates either when all sources have been decided or when the maximum 
* number of iterations, M, is reached. The estimate of the p value at termination is output 
* in P, and the number of actual iterations is output in Mn.
***********/


void permuteSPR(
	int* Ni,	   /* scalar number of observations in Y*/
	int* nci,     /* scalar number of of columns in Y */
	double* Y,	   /* nc*Nx matrix  of responses */
	double* Q,     /* NxK matrix of linear functionals in row major form */
	int* nSi,      /* scalar number of sources into which the functionals are divided */
	int* dof,     /* nS vector of df for the sources: K=sum(df) */
	int* Mi,       /* scalar maximum limit on the number of iterations*/
	int* Mn,      /* nc*nS matrix of sample sizes at termination */
	int* Mnt,     /* nc*K matrix of sample sizes for the terms at termination */
	double* p0,    /* scalar null hyp p value */
	double* p1,    /* scalar alternative hyp p value */
	double* alpha, /* scalar test size */
	double* beta,  /* scalar type II size */
	int* accept,  /* nc*nS matrix:  -1: rejected, 0: undecided, 1: accepted */
	double* P,      /* nc*nS matrix of observed p values at termination*/
 	int* acceptt,  /* nc*nS matrix:  -1: rejected, 0: undecided, 1: accepted */
	double* Pt,     /* nc*K matrix  of observed p values for the terms at termination*/
	int* nCyclei    /* complete permutation every nCycle*N iterations*/
)
{
	double an;  /* acceptance number an+m*rn */
	double rn;
	double bn;  /* rejection number bn+m*rn */
	double* SS;  /* nc*nS matrix of current sum of squares */
	double* oSS;  /* nc*nS matrix of Sum of squares of observations */
	double* b;   /* nc*K matrix of term estimates */
	double* obSq;  /* nc*K matrix of square of original estimates */
	int* d;   /* nc*nS matrix of number of exceedences */
	int* dt;   /* nc*K matrix of number of exceedences for the terms */
	int* acceptC; /* nc*nS matrix of the count of the number of terms in a source that have been decided */
	double* rss;  /* Residual sum of squares */
	int K=0;       /* total number of terms */
	int Kt=0;		/* total non residual terms */
	double pt;
	int i;
	int j;
	int k;
	int u1;
	int u2;
	int iter;
	double temp;
	int crita;
	int critr;
	double y;
	double x;
	double y1;
	double y2;
	double x1;
	double x2;
	double delta;
	bool nToDo; /* number of sources not decided */
	bool useF=true;  /* Use F ratio instead of SS */
	int nSe=nS-1;

	int l;
	double* Yp;
	double* bp;
	double* oSSp;
	double* obSqp;
	double* SSp;
	int* dp;
	int* dtp;
	double* Pp;
	double* Ptp;
	int* acceptp;
	int* acceptCp;
	int* accepttp;
	int* Mnp;
	int* Mntp;


	for (i=0;i<nS;i++) /* find number of columns in Q */
		K+=dof[i];

	if (dof[nSe]==0) { /* No residuals */
		useF=false;
	}
	
	Kt=K-dof[nSe];

	nToDo=nc*(nSe+Kt);


	SS=(double *)S_alloc(nc*nS,sizeof(double));  /* These should be zeroed as allocated */
	oSS=(double *)S_alloc(nc*nS,sizeof(double)); 
	d=(int *)S_alloc(nc*nSe,sizeof(int));
	b=(double *)S_alloc(nc*K,sizeof(double));
	obSq=(double *)S_alloc(nc*K,sizeof(double));
	dt=(int *)S_alloc(nc*Kt,sizeof(int));
	acceptC=(int *)S_alloc(nc*nS,sizeof(int));
	rss=(double *)S_alloc(nc,sizeof(double));

	pt=log(*p1/(*p0))-log((1-*p1)/(1-*p0));
	rn=log((1-*p0)/(1-*p1))/pt;

	an=log(*beta/(1-*alpha))/pt;
	bn=log((1-*beta)/(*alpha))/pt;


	/* Calculate b */
	for (l=0;l<nc;l++) {
		Yp=Y+N*l;
		bp=b+K*l;
		k=0;
		for (i=0;i<N;i++) {
			y=Yp[i];
			for (j=0;j<K;j++) {
				bp[j]+=y*Q[k++];  /* Note Q is row major */
			}	
		}
	}

	/* Get oSS and obSq */
	for (l=0;l<nc;l++) {
		k=0;
		oSSp=oSS+nS*l;
		obSqp=obSq+K*l;
		bp=b+K*l;
		for (i=0;i<nS;i++) {
			oSSp[i]=0;
			for (j=0;j<dof[i];j++) {
				x=bp[k];
				oSSp[i]+=obSqp[k++]=x*x;
			}
		}
	}

	/* Scale oSSp and obSqp by rss */
	if (useF) {
		for (l=0;l<nc;l++) {
			oSSp=oSS+nS*l;
			obSqp=obSq+K*l;
			rss[l]=oSSp[nSe];
			if (rss[l]>tol) {
				for (i=0;i<nSe;i++) {
					oSSp[i]/=rss[l];
				}
				for (i=0;i<K;i++) {
					obSqp[i]/=rss[l];
				}
			} else {
				useF=false;
				rss[l]=1;
			}
		}
	} else {
		for (l=0;l<nc;l++)
			rss[l]=1;
	}
	/* initial SS */
	memcpy(SS,oSS,nc*nS*sizeof(double));


	GetRNGstate();

	/* Iterate until M or all sources decided */
	for (iter=1;iter<=M && nToDo>0;iter++) {
		/* Get two random indices */
		u1=(int)((double)N*unif_rand()); 
		u2=(int)((double)N*unif_rand());

		crita=floor(an+(iter)*rn);
		critr=ceil(bn+(iter)*rn);


		/* Update the SS */
		for (l=0;l<nc;l++) {

			Yp=Y+N*l;
			SSp=SS+nS*l;
			oSSp=oSS+nS*l;
			bp=b+K*l;
			obSqp=obSq+K*l;
			dp=d+nS*l;
			dtp=dt+Kt*l;
			acceptp=accept+nS*l;
			acceptCp=acceptC+nS*l;
			accepttp=acceptt+Kt*l;
			Pp=P+nSe*l;
			Ptp=Pt+K*l;
			Mnp=Mn+nSe*l;
			Mntp=Mnt+K*l;

			if (iter==1 || 0==iter%(nCycle*N)) { /* permute everything at the start and then once in nCycle*N cycles */
				permuteRand(&Yp,N);
				memset(bp,0,K*sizeof(double));
				k=0;
				for (i=0;i<N;i++) {
					y=Yp[i];
					for (j=0;j<K;j++) {
						bp[j]+=y*Q[k++];  
					}	
				}
				k=0;
				for (i=0;i<nS;i++) {
					SSp[i]=0;
					for (j=0;j<dof[i];j++) {
						x=bp[k++];
						SSp[i]+=x*x;
					}
				}
			}
			else {

				y1=Yp[u1];
				y2=Yp[u2];
				delta=y1-y2;
				k=0;
				for (i=0;i<nS;i++) {
					if (acceptp[i]==0) {  /* acceptp[nSe] is always 0 */
						for (j=0;j<dof[i];j++) {
							SSp[i]-=bp[k]*bp[k]; /* Remove k term */
							x1=Q[IndQ(u1,k)];
							x2=Q[IndQ(u2,k)];
							bp[k]+=delta*(x2-x1); /* Update k term with exchanged observations*/
							SSp[i]+=bp[k]*bp[k];
							k++;
						}
					}
					else
						k+=dof[i];
				}
				/* Update Y */
				temp=Yp[u1];
				Yp[u1]=Yp[u2];
				Yp[u2]=temp;
			}

			if (useF)
				rss[l]=SSp[nSe];

			/* Compare SS with oSS and accept, reject or continue */
			k=0;
			for (i=0;i<nSe;i++) {
				for (j=0;j<dof[i];j++) {
					if (accepttp[k]==0 ) {
						dtp[k]+=(bp[k]*bp[k]+tol)>obSqp[k]*rss[l];
						if (dtp[k]<=crita /*&& iter>iterMin*/) {
							accepttp[k]=1;
							acceptCp[i]++;
							Ptp[k]=(double)dtp[k]/(double)iter;
							Mntp[k]=iter;
							nToDo--;
						}
						else if (dtp[k]>=critr /*&& iter>iterMin*/) {
							accepttp[k]=-1;
							acceptCp[i]++;
							Ptp[k]=(double)dtp[k]/(double)iter;
							Mntp[k]=iter;
							nToDo--;
						}
					}
					k++;
				}
				if (acceptp[i]==0 ) {
					dp[i]+=(SSp[i]+tol)>oSSp[i]*rss[l];
					if (dp[i]<=crita /*&& iter>iterMin*/) {
						if (acceptCp[i]>=dof[i]) {  /* acceptCp[nSe] is never checked */
							acceptp[i]=1;
							Pp[i]=(double)dp[i]/(double)(iter);
							Mnp[i]=iter;
							nToDo--;
						}
					}
					else 
						if (dp[i]>=critr /*&& iter>iterMin*/) {
							if (acceptCp[i]>=dof[i]) {
								acceptp[i]=-1;
								Pp[i]=(double)dp[i]/(double)(iter);
								Mnp[i]=iter;
								nToDo--;
							}
						}
				}
			}
		}
	}

	for (l=0;l<nc;l++){
		acceptp=accept+nS*l;
		Pp=P+nSe*l;
		Mnp=Mn+nSe*l;
		accepttp=acceptt+K*l;
		Ptp=Pt+K*l;
		Mntp=Mnt+K*l;
		dp=d+nSe*l;
		dtp=dt+K*l;


		k=0;
		for (i=0;i<nSe;i++) {
			if (acceptp[i]==0) {
				Pp[i]=(double)dp[i]/(double)M;
				Mnp[i]=M;
			}
			for (j=0;j<dof[i];j++) {
				if (accepttp[k]==0) {
					Ptp[k]=(double)dtp[k]/(double)M;
					Mntp[k]=M;
				}
				k++;
			}
		}
	}
	PutRNGstate();

}

/* permuteRand ***************************************************************
|  Randomly pemutes the n values in Y[] using the Fike
|  algorithm.  See Fike, "A permutation generation method"  The Computer
|  Journal, 18-1, Feb 75, 21-22.
*/

void	permuteRand(
	double** Y,
	int n
)
{
   int i, j;
   double temp;

   for (i=0;i<n;i++) {
      j=floor((1+i)*unif_rand());
      temp=(*Y)[j];
      (*Y)[j]=(*Y)[i];
      (*Y)[i]=temp;
   }


}

