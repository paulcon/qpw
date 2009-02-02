#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <f2c.h>
#include <CLAPACK/clapack.h>
#include "spectral_functions.h"

int recur
(
	int n,
	struct param p,
	double *ab
)
/*
 * Compute the recurrence coefficients associated with a Jacobi weight function
 * supported on the interval [l,r] with parameters a and b.
 */
{
	char jacobi_str[] = "jacobi";
	double l,r,a,b,a0,b2a2,s,o;
	int k;

	if(strcmp(jacobi_str,p.name)==0){
		l=p.l; r=p.r; a=p.parm1; b=p.parm2;
		a0=(b-a)/(a+b+2);
		b2a2=b*b-a*a;
		s=(r-l)/2;
		o=l+(r-l)/2;
		if(n>0){
			ab[0] = s*a0+o;
			ab[1] = 1.0;
		}
		for(k=1;k<n;k++){
			ab[2*k] = s*b2a2/((2*k+a+b)*(2*(k+1) + a+b))+o;
			ab[2*k+1] = ((r-l)*(r-l)*k*(k+a)*(k+b)*(k+a+b)) /
							((2*k+a+b)*(2*k+a+b)*(2*k+a+b+1)*(2*k+a+b-1));
		}
		return 0;
	} else {
		printf("ERROR: Only Jacobi family of parameters available. Your type was %s.",p.name);
		return 1;
	}
}

int jacobi_matrix
(
	int n,
	struct param p,
	double *JD,
	double *JE
)
/*
 * Construct the symmetric, tridiagonal Jacobi matrix associated with the parameter p.
 */
{
	int i;
	double *ab;

	// Compute the recurrence coefficients of the monic orthogonal polynomials.
	ab = (double *)calloc(2*n,sizeof(double));
	recur(n,p,ab);

	JD[0] = ab[0];
	if(n==1){ return 0; }
	JE[0]=sqrt(ab[3]);
	for(i=1;i<n-1;i++){
		JD[i]=ab[2*i];
		JE[i]=sqrt(ab[2*i+3]);
	}
	JD[n-1]=ab[2*(n-1)];
	return 0;
}

int jacobi_eigenvectors
(
	int n,
	struct param p,
	double *Lambda,
	double *Q
)
/*
 * Construct the eigenvectors of the nxn Jacobi matrix associated with parameter p.
 */
{
	char flag='I';
	double *JE, *WORK;
	int INFO=1,N,LDZ;

	JE=(double *)calloc(n-1,sizeof(double));
	WORK=(double *)calloc(2*n-2,sizeof(double));

	jacobi_matrix(n,p,Lambda,JE);

	N=n; LDZ=n;
	dsteqr_(&flag,&N,Lambda,JE,Q,&LDZ,WORK,&INFO);

	if(INFO==0){
		return 0;
	} else {
		printf("ERROR: There was an error computing the eigenvectors of the Jacobi matrix.\n");
		return 1;
	}
}

int evaluate_ops
(
	int n,
	struct param p,
	double point,
	double *P
)
/*
 * Evaluate the vector of orthonormal polynomials P associated with p evaluated at a point.
 */
{
	int res,i;
	double *ab;

	ab=(double *)calloc(2*n,sizeof(double));

	res=recur(n,p,ab);

	P[0]=1.0/sqrt(ab[1]);
	if(n==1){ return 0; }
	P[1]=(point-ab[0])/sqrt(ab[3]);
	if(n==2){ return 0; }
	for(i=2;i<n;i++){
	    P[i]=((point-ab[2*(i-1)])*P[i-1]-sqrt(ab[2*(i-1)+1])*P[i-2])/sqrt(ab[2*i+1]);
	}

	return 0;
}
