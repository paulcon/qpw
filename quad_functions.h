#ifndef QUAD_FUNCTIONS_H_
#define QUAD_FUNCTIONS_H_

#endif /*QUAD_FUNCTIONS_H_*/

int gauss_1d_pts_wts
(
	double *x, 
	double *w, 
	int n
);

int clenshaw_curtis_1d_pts_wts
(
	double *x, 
	double *w, 
	int n
);

long nchoosek(long n, long k);

void enumerate_compositions
(
	int d,
	int q,
	int *comps
);

void update_index
(
	int *index,
	int *multi_index,
	int dim
);

int write_1d_points_and_weights
(
	char *rule_name,
	int *level_ptr
);

int write_tensor_points_and_weights
(
	char *rule_name,
	int *dim_ptr,
	int *level_ptr
);

int write_sparse_points_and_weights
(
	char *rule_name,
	int *dim_ptr,
	int *level_ptr
);

