#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "quad_functions.h"
#include "spectral_functions.h"
#include "unit_tests.h"

int main(int argc, char **argv){

	if(test_jacobi_recur()){
		printf("Jacobi recurrence fail.\n");
	} else {
		printf("Jacobi recurrence pass.\n");
	}

	if(test_jacobi_matrix()){
		printf("Jacobi matrix fail.\n");
	} else {
		printf("Jacobi matrix pass.\n");
	}

	if(test_jacobi_eigenvecs()){
		printf("Jacobi eigs fail.\n");
	} else {
		printf("Jacobi eigs pass.\n");
	}

	if(test_eval_ops()){
		printf("Evaluate polys fail.\n");
	} else {
		printf("Evaluate polys pass.\n");
	}

	return 0;
}
