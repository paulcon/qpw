#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "quad_functions.h"

int main(int argc, char **argv){
	int res, level, dim;
	long mem_lim = 2e8;
	long mem;
	char ch;
	char rule_name[100];
	char grid_type[100];
	char tensor_str[] = "tensor";
	char sparse_str[] = "sparse";
	char gauss_str[] = "gauss";
	char cc_str[] = "cc";

	if(argc != 5){ 
		printf("Usage: qpw ptype gtype dim level\n");
		printf("ptype: cc, gauss\n");
		printf("gtype: tensor, sparse\n");
		printf("dim: integer>0\n");
		printf("level: number of points\n");
		printf("\nView README for more details.\n");
		return 1;
	}
		
	strcpy(rule_name,argv[1]);
	strcpy(grid_type,argv[2]);
	dim=atoi(argv[3]);
	level=atoi(argv[4]);
		
	// Some error checking.
	if(dim<1){
		printf("ERROR: dim should be >0, but dim=%d.\n",dim); 
		return 1;
	}
	if(strcmp(gauss_str,rule_name)!=0 && strcmp(cc_str,rule_name)!=0){
		printf("ERROR: Point type should be gauss or cc, but ptype=%s.\n",rule_name); 
		return 1;
	}
	// One dimension
	if(dim==1){ 
		mem = 2*level*sizeof(double);
		if(mem>mem_lim){
			printf("You are about to allocate %ld bytes for points and weights. Are you sure you want to continue? (Y/N)\n",mem);
			scanf(" %c", &ch);
			ch = toupper(ch);
			if(ch!='Y'){ return 0; }
		}
		
		res = write_1d_points_and_weights(rule_name,&level); 
		printf("Wrote %d %s points.\n",res,rule_name);
		return 0;
	}
	
	// More error checking
	if(level<1){
		printf("ERROR: Level should be >0, but level=%d.\n",level); 
		return 1;
	}
	if(strcmp(tensor_str,grid_type)!=0 && strcmp(sparse_str,grid_type)!=0){
		printf("ERROR: Grid type should be tensor or sparse, but gtype=%s.\n",grid_type); 
		return 1;
	}
	
	if(strcmp(tensor_str,grid_type)==0){
		mem = (pow(level,dim)+level)*sizeof(double);
		if(mem>mem_lim){
			printf("You are about to allocate %ld bytes for points and weights. Are you sure you want to continue? (Y/N)\n",mem);
			scanf(" %c", &ch);
			ch = toupper(ch);
			if(ch!='Y'){ return 0; }
		}
		res = write_tensor_points_and_weights(rule_name,&dim,&level); 
		printf("Wrote %d %d-dimensional %s points on a %s grid.\n",res,dim,rule_name,grid_type);
	} else if(strcmp(sparse_str,grid_type)==0){
		mem = (pow(pow(2,level)+1,dim)+pow(2,level)+1)*sizeof(double);
		if(mem>mem_lim){
			printf("You are about to allocate %ld bytes for points and weights. Are you sure you want to continue? (Y/N)\n",mem);
			scanf(" %c", &ch);
			ch = toupper(ch);
			if(ch!='Y'){ return 0; }
		}
		res = write_sparse_points_and_weights(rule_name,&dim,&level); 
		printf("Wrote %d %d-dimensional %s points on a %s grid.\n",res,dim,rule_name,grid_type);
	} else {
		printf("ERROR: Something went terribly wrong.\n");	
	}
		
	return 0;	
}
