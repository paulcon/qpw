#ifndef SPECTRAL_FUNCTIONS_H_
#define SPECTRAL_FUNCTIONS_H_

#endif /*SPECTRAL_FUNCTIONS_H_*/

struct param{
	char name[40];
	double l;
	double r;
	double parm1;
	double parm2;
};

int recur
(
	int n,
	struct param p,
	double *ab
);

int jacobi_matrix
(
	int n,
	struct param p,
	double *JD,
	double *JE
);

int jacobi_eigenvectors
(
	int n,
	struct param p,
	double *Lambda,
	double *Q
);

int evaluate_ops
(
	int n,
	struct param p,
	double point,
	double *P
);
