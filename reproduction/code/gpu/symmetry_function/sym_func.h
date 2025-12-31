/*
 * Header file for sym_func.cu
 */

#ifndef SYM_H
#define SYM_H

#include <string>
#include <iostream>

// cuda error check
#define HANDLE_ERROR(err) (HandleError(err, __FILE__, __LINE__))
static void HandleError(cudaError_t err, const char *file, int line) {
	if (err != cudaSuccess) {
		printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
		exit(EXIT_FAILURE);
	}
}

// Sym_Param structure
struct Sym_Param {
	int atomtype_num;		// the number of distinct elements
	double rrc;				// cutoff radius for radial environment
	double arc;				// cutoff radius for angular environment
	
	int rad_eta_num;		// number of parameter etas for radial
	int rad_rs_num;		    // number of Rs for radial
	int ang_zeta_num;
	int ang_theta_num;
	int ang_eta_num;
	int ang_rs_num;

	double *rad_etas;
	double *rad_rss;
	double *ang_zetas;
	double *ang_thetas;
	double *ang_etas;
	double *ang_rss;

	// radial, angular and total dimension
	int rad_dim;
	int ang_dim;
	int tot_dim;
};

// helper functions
__host__ __device__ void printParam(Sym_Param *params);
__host__ __device__ void printIndexmap(int *indexmap);

// host functions
__host__ void createSymParam(std::string fname, Sym_Param *parms);
__host__ void createSymParamAuto(std::string fname, Sym_Param *parms);
__host__ void memcpyParam(Sym_Param *params, Sym_Param *d_params);
__host__ void freeParam(Sym_Param *params, Sym_Param *d_params);
__host__ void buildIndexmap(int *indexmap, Sym_Param *params);

// device functions
__device__ double cutoff(double r, double rs);
__device__ int hash(char c);
__device__ int getIndexR(char c);
__device__ int getIndexA(char c, char d);
__device__ void buildrad(double *dists, char *d_atoms, int num, double *d_aevs, Sym_Param *d_params, int *d_indexmap, int offset);
__device__ void buildang(double *dists, double *d_coords, char *d_atoms, int num, double *d_aevs, Sym_Param *d_params, int *d_indexmap, int offset);

// kernel function
__global__ void buildaev(double *d_coords, char *d_atoms, int *nums, double *d_aevs, Sym_Param *d_params, int *d_indexmap);

#endif
