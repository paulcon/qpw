#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "quad_functions.h"
#include "spectral_functions.h"

#define EPS 3.0e-14 //EPS is the relative precision.

int test_jacobi_recur()
{
	int n=5,res,i;
	double s,t1,t2;
	double *ab;
	double ab_true[10];
	struct param p = {"jacobi",-1.0,1.0,0.0,0.0};

	ab_true[0]=0.0; ab_true[1]=1.0;
	ab_true[2]=0.0; ab_true[3]=3.3333333333333331e-01;
	ab_true[4]=0.0; ab_true[5]=2.6666666666666666e-01;
	ab_true[6]=0.0; ab_true[7]=2.5714285714285712e-01;
	ab_true[8]=0.0; ab_true[9]=2.5396825396825395e-01;

	ab = (double *)calloc(2*n,sizeof(double));
	res = recur(n,p,ab);

	s=0.0;
	for(i=0;i<n;i++){
		t1=ab[2*i]-ab_true[2*i];
		t2=ab[2*i+1]-ab_true[2*i+1];
		s=s+t1*t1+t2*t2;
	}

	if(res==0 && sqrt(s)<sqrt(EPS)){
		return 0;
	} else {
		return 1;
	}

}

int test_jacobi_matrix()
{
	int n=5,res,i;
	double s,t1,t2;
	double *JD, *JE;
	double JD_true[5], JE_true[4];
	struct param p = {"jacobi",-1.0,1.0,0.0,0.0};

	JD_true[0]=0.0;
	JD_true[1]=0.0;
	JD_true[2]=0.0;
	JD_true[3]=0.0;
	JD_true[4]=0.0;
	JE_true[0]=0.577350269189626;
	JE_true[1]=0.516397779494322;
	JE_true[2]=0.507092552837110;
	JE_true[3]=0.503952630678970;

	JD=(double *)calloc(n,sizeof(double));
	JE=(double *)calloc(n-1,sizeof(double));
	res=jacobi_matrix(n,p,JD,JE);

	s=0.0;
	for(i=0;i<n-1;i++){
		t1=JD[i]-JD_true[i];
		t2=JE[i]-JE_true[i];
		s=s+t1*t1+t2*t2;
	}
	s=s+(JD[n-1]-JD_true[n-1])*(JD[n-1]-JD_true[n-1]);

	if(res==0 && sqrt(s)<sqrt(EPS)){
		return 0;
	} else {
		return 1;
	}
}

int test_jacobi_eigenvecs(){
	int n=5,res,i,j;
	double *Q, *Lambda;
	double Qtrue[25];
	double s,t;
	struct param p = {"jacobi",-1.0,1.0,0.0,0.0};

	Q=(double *)calloc(n*n,sizeof(double));
	Lambda=(double *)calloc(n,sizeof(double));

	Qtrue[0]=0.344185186386768;
	Qtrue[1]=-0.540215698889530;
	Qtrue[2]=0.563165025741707;
	Qtrue[3]=-0.456253217712814;
	Qtrue[4]=0.253735514371261;
	Qtrue[5]=-0.489197644362361;
	Qtrue[6]=0.456253217712814;
	Qtrue[7]=0.071185504167388;
	Qtrue[8]=-0.540215698889530;
	Qtrue[9]=0.505587073358043;
	Qtrue[10]=0.533333333333333;
	Qtrue[11]=-0.000000000000000;
	Qtrue[12]=-0.596284793999944;
	Qtrue[13]=-0.000000000000000;
	Qtrue[14]=0.600000000000000;
	Qtrue[15]=-0.489197644362361;
	Qtrue[16]=-0.456253217712814;
	Qtrue[17]=0.071185504167388;
	Qtrue[18]=0.540215698889530;
	Qtrue[19]=0.505587073358043;
	Qtrue[20]=0.344185186386768;
	Qtrue[21]=0.540215698889530;
	Qtrue[22]=0.563165025741706;
	Qtrue[23]=0.456253217712814;
	Qtrue[24]=0.253735514371261;

	res=jacobi_eigenvectors(n,p,Lambda,Q);

	s=0.0;
	for(i=0;i<n-1;i++){
		for(j=0;j<n;j++){
			t=Qtrue[i+n*j]-Q[i+n*j];
			s=s+t*t;
		}
	}

	if(res==0 && sqrt(s)<sqrt(EPS)){
		return 0;
	} else {
		return 1;
	}
}

int test_eval_ops(){
	double s,t,point=0.53;
	int n=5,i,res;
	double *P;
	double Ptrue[5];
	struct param p = {"jacobi",-1.0,1.0,0.0,0.0};

	P=(double *)calloc(n,sizeof(double));

	Ptrue[0]=1.000000000000000;
	Ptrue[1]=0.917986928011505;
	Ptrue[2]=-0.175866746430358;
	Ptrue[3]=-1.118643497452942;
	Ptrue[4]=-0.999499368750000;

	res=evaluate_ops(n, p, point, P);

	s=0.0;
	for(i=0;i<n;i++){
		t=P[i]-Ptrue[i];
		s=s+t*t;
	}

	if(res==0 && sqrt(s)<sqrt(EPS)){
		return 0;
	} else {
		return 1;
	}
}
