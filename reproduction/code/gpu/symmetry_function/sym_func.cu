/*
 * CUDA implementation of the radial and angular symmetry functions
 */

 #include <string>
 #include <fstream>
 #include <cassert>
 #include <iostream>
 
 #include "sym_func.h"
 
 #define ATOMTYPE_NUM 4
 
 using namespace std;
 
 /************************** Program parameters ****************************/
 const double PI = 3.14159265358979323846264338328;	// it's a long double, double = 3.141592653589793
 
 /*************************** Helper functions ****************************/
 __host__ __device__ void printParam(Sym_Param *params) {
	 printf("\nParameters: \n");
	 printf("%d\n", params->atomtype_num);
	 printf("%f\n", params->rrc);
	 printf("%f\n", params->arc);
	 printf("%d\n", params->rad_eta_num);
	 printf("%d\n", params->rad_rs_num);
	 printf("%d\n", params->ang_zeta_num);
	 printf("%d\n", params->ang_theta_num);
	 printf("%d\n", params->ang_eta_num);
	 printf("%d\n", params->ang_rs_num);
	 for (int i=0; i<params->rad_eta_num; i++)
		 printf("%f ", params->rad_etas[i]);
	 printf("\n");
	 for (int i=0; i<params->rad_rs_num; i++)
		 printf("%f ", params->rad_rss[i]);
	 printf("\n");
	 for (int i=0; i<params->ang_zeta_num; i++)
		 printf("%f ", params->ang_zetas[i]);
	 printf("\n");
	 for (int i=0; i<params->ang_theta_num; i++)
		 printf("%f ", params->ang_thetas[i]);
	 printf("\n");
	 for (int i=0; i<params->ang_eta_num; i++)
		 printf("%f ", params->ang_etas[i]);
	 printf("\n");
	 for (int i=0; i<params->ang_rs_num; i++)
		 printf("%f ", params->ang_rss[i]);
	 printf("\n");
 
	 printf("rad_dim: %d\n", params->rad_dim);
	 printf("ang_dim: %d\n", params->ang_dim);
	 printf("tot_dim: %d\n\n", params->tot_dim);
 }
 
 __host__ __device__ void printIndexmap(int *indexmap) {
	 printf("\nIndexmap: \n");
	 for (int i=0; i<ATOMTYPE_NUM + ATOMTYPE_NUM * (ATOMTYPE_NUM + 1) / 2; i++)
		 printf("%d ", indexmap[i]);
	 printf("\n\n");
 }
 
 /**************************** Host functions ****************************/
 // create the sym_param structure from a file, it doesn't check the validity of the input parameter file TODO
 __host__ void createSymParam(string fname, Sym_Param *params) {
	 ifstream paramf (fname);
	 // open the parameter file and initialize sym_param
	 if (paramf.is_open()) {
		 string line;
		 paramf >> params->atomtype_num;     getline(paramf, line);
		 paramf >> params->rrc;              getline(paramf, line);
		 paramf >> params->arc;              getline(paramf, line);
		 paramf >> params->rad_eta_num;      getline(paramf, line);
		 paramf >> params->rad_rs_num;       getline(paramf, line);
		 paramf >> params->ang_zeta_num;     getline(paramf, line);
		 paramf >> params->ang_theta_num;    getline(paramf, line);
		 paramf >> params->ang_eta_num;      getline(paramf, line);
		 paramf >> params->ang_rs_num;       getline(paramf, line);
		 // allocate space
		 params->rad_etas   = new double[params->rad_eta_num];
		 params->rad_rss    = new double[params->rad_rs_num];
		 params->ang_zetas  = new double[params->ang_zeta_num];
		 params->ang_thetas = new double[params->ang_theta_num];
		 params->ang_etas   = new double[params->ang_eta_num];
		 params->ang_rss    = new double[params->ang_rs_num];
		 // read in more parameters
		 for (int i=0; i<params->rad_eta_num; i++)
			 paramf >> params->rad_etas[i];
		 getline(paramf, line);
		 for (int i=0; i<params->rad_rs_num; i++)
			 paramf >> params->rad_rss[i];
		 getline(paramf, line);
		 for (int i=0; i<params->ang_zeta_num; i++)
			 paramf >> params->ang_zetas[i];
		 getline(paramf, line);
		 for (int i=0; i<params->ang_theta_num; i++)
			 paramf >> params->ang_thetas[i];
		 getline(paramf, line);
		 for (int i=0; i<params->ang_eta_num; i++)
			 paramf >> params->ang_etas[i];
		 getline(paramf, line);
		 for (int i=0; i<params->ang_rs_num; i++)
			 paramf >> params->ang_rss[i];
		 // initialize radial, angular, and total dimension
		 params->rad_dim = params->rad_eta_num * params->rad_rs_num;
		 params->ang_dim = params->ang_zeta_num * params->ang_theta_num * params->ang_eta_num * params->ang_rs_num;
		 params->tot_dim = params->rad_dim * ATOMTYPE_NUM + params->ang_dim * (ATOMTYPE_NUM * (ATOMTYPE_NUM+1) / 2);
	 }
	 else {
		 cout << fname << ": can not open" << endl;
		 assert(false);
	 }
 }
 
 __host__ void createSymParamAuto(string fname, Sym_Param *params) {
	 ifstream paramf (fname);
	 // open the parameter file and initialize sym_param
	 if (paramf.is_open()) {
		 string line;
		 paramf >> params->atomtype_num;     getline(paramf, line);
		 paramf >> params->rrc;              getline(paramf, line);
		 paramf >> params->arc;              getline(paramf, line);
		 paramf >> params->rad_eta_num;      getline(paramf, line);
		 paramf >> params->rad_rs_num;       getline(paramf, line);
		 paramf >> params->ang_zeta_num;     getline(paramf, line);
		 paramf >> params->ang_theta_num;    getline(paramf, line);
		 paramf >> params->ang_eta_num;      getline(paramf, line);
		 paramf >> params->ang_rs_num;       getline(paramf, line);
		 // allocate space
		 params->rad_etas   = new double[params->rad_eta_num];
		 params->rad_rss    = new double[params->rad_rs_num];
		 params->ang_zetas  = new double[params->ang_zeta_num];
		 params->ang_thetas = new double[params->ang_theta_num];
		 params->ang_etas   = new double[params->ang_eta_num];
		 params->ang_rss    = new double[params->ang_rs_num];
		 // automatically compute the parameters
		 assert(params->rad_rs_num > 1 && params->ang_theta_num > 1 && params->ang_rs_num > 1);
		 // assume there is only one radial eta, and the cross point is 0.6 high
		 double interval = params->rrc / params->rad_rs_num;
		 params->rad_etas[0] = -log(0.6) * 4 / (interval * interval);
		 // equally spaced Rss
		 for (int i=0; i<params->rad_rs_num; i++)
			 params->rad_rss[i] = (i+1) * interval;
		 // assume there is only one angular zeta, the cross point for zeta is 0.6
		 interval = 2 * PI / params->ang_theta_num;
		 params->ang_zetas[0] = log(0.6 / 2) / log((1 + cos(interval / 2)) / 2);
		 // equally spaced Angular thetas
		 for (int i=0; i<params->ang_theta_num; i++)
			 params->ang_thetas[i] = i * interval;
		 // assume there only one angular eta, the cross point is 0.6
		 interval = params->arc / params->ang_rs_num;
		 params->ang_etas[0] = -log(0.6) * 4 / (interval * interval);
		 // equally spaced Angular Rss
		 for (int i=0; i<params->ang_rs_num; i++)
			 params->ang_rss[i] = (i+1) * interval;
		 // initialize radial, angular, and total dimension
		 params->rad_dim = params->rad_eta_num * params->rad_rs_num;
		 params->ang_dim = params->ang_zeta_num * params->ang_theta_num * params->ang_eta_num * params->ang_rs_num;
		 params->tot_dim = params->rad_dim * ATOMTYPE_NUM + params->ang_dim * (ATOMTYPE_NUM * (ATOMTYPE_NUM+1) / 2);
	 }
	 else {
		 cout << fname << ": can not open" << endl;
		 assert(false);
	 }
 }
 
 // build the index map, which maps the atom-atom to corresponding index in the atomic environment vector
 __host__ void buildIndexmap(int *indexmap, Sym_Param *params) {	// Need to be done differently TODO
	 int n = ATOMTYPE_NUM;
	 int N = ATOMTYPE_NUM + ATOMTYPE_NUM * (ATOMTYPE_NUM + 1) / 2;
 
	 for (int i=0; i<n; i++) {
		 indexmap[i] = i * params->rad_dim;
	 }
	 for (int i=n; i<N; i++) {
		 indexmap[i] = n * params->rad_dim + (i-n) * params->ang_dim;
	 }
 }
 
 // copy parameters structure to Device
 __host__ void memcpyParam(Sym_Param *params, Sym_Param *d_params) {
	 // allocate space on device to store parameters, need to be freed using cudaFree in freeParam
	 double *d_rad_etas, *d_rad_rss, *d_ang_zetas, *d_ang_thetas, *d_ang_etas, *d_ang_rss;
	 HANDLE_ERROR(cudaMalloc((void **)&d_rad_etas,   params->rad_eta_num * sizeof(double)));
	 HANDLE_ERROR(cudaMalloc((void **)&d_rad_rss,    params->rad_rs_num * sizeof(double)));
	 HANDLE_ERROR(cudaMalloc((void **)&d_ang_zetas,  params->ang_zeta_num * sizeof(double)));
	 HANDLE_ERROR(cudaMalloc((void **)&d_ang_thetas, params->ang_theta_num * sizeof(double)));
	 HANDLE_ERROR(cudaMalloc((void **)&d_ang_etas,   params->ang_eta_num * sizeof(double)));
	 HANDLE_ERROR(cudaMalloc((void **)&d_ang_rss,    params->ang_rs_num * sizeof(double)));
	 // copy the parameters from Host to Device
	 HANDLE_ERROR(cudaMemcpy(d_rad_etas,   params->rad_etas,   params->rad_eta_num * sizeof(double), cudaMemcpyHostToDevice));
	 HANDLE_ERROR(cudaMemcpy(d_rad_rss,    params->rad_rss,    params->rad_rs_num * sizeof(double), cudaMemcpyHostToDevice));
	 HANDLE_ERROR(cudaMemcpy(d_ang_zetas,  params->ang_zetas,  params->ang_zeta_num * sizeof(double), cudaMemcpyHostToDevice));
	 HANDLE_ERROR(cudaMemcpy(d_ang_thetas, params->ang_thetas, params->ang_theta_num * sizeof(double), cudaMemcpyHostToDevice));
	 HANDLE_ERROR(cudaMemcpy(d_ang_etas,   params->ang_etas,   params->ang_eta_num * sizeof(double), cudaMemcpyHostToDevice));
	 HANDLE_ERROR(cudaMemcpy(d_ang_rss,    params->ang_rss,    params->ang_rs_num * sizeof(double), cudaMemcpyHostToDevice));
	 // free the space in params
	 delete [] params->rad_etas;
	 delete [] params->rad_rss;
	 delete [] params->ang_zetas;
	 delete [] params->ang_thetas;
	 delete [] params->ang_etas;
	 delete [] params->ang_rss;
	 // set the array pointers to the address on device
	 params->rad_etas   = d_rad_etas;
	 params->rad_rss	   = d_rad_rss;
	 params->ang_zetas  = d_ang_zetas;
	 params->ang_thetas = d_ang_thetas;
	 params->ang_etas   = d_ang_etas;
	 params->ang_rss	   = d_ang_rss;
	 // copy the Sym_Param structure from Host to Device
	 HANDLE_ERROR(cudaMemcpy(d_params, params, sizeof(Sym_Param), cudaMemcpyHostToDevice));
 }
 
 // cudaFree the space of parameter structure on device
 __host__ void freeParam(Sym_Param *params, Sym_Param *d_params) {
	 // free the space of the parameter arrays
	 HANDLE_ERROR(cudaFree(params->rad_etas));
	 HANDLE_ERROR(cudaFree(params->rad_rss));
	 HANDLE_ERROR(cudaFree(params->ang_zetas));
	 HANDLE_ERROR(cudaFree(params->ang_thetas));
	 HANDLE_ERROR(cudaFree(params->ang_etas));
	 HANDLE_ERROR(cudaFree(params->ang_rss));
	 // cuda free the Sym_Param structure
	 HANDLE_ERROR(cudaFree(d_params));
 }
 
 /*************************** Device functions ***************************/
 // cutoff function
 __device__ double cutoff(double r, double rs) {
	 return (r >= rs) ? 0.0 : (0.5 * cos(PI * r / rs) + 0.5);
 }
 
 // get the distance between two atoms
 __device__ double getDist(double *coords, int i, int j) {
	 return sqrt((coords[3*i]-coords[3*j]) * (coords[3*i]-coords[3*j]) +
				 (coords[3*i+1]-coords[3*j+1]) * (coords[3*i+1]-coords[3*j+1]) +
				 (coords[3*i+2]-coords[3*j+2]) * (coords[3*i+2]-coords[3*j+2]));
 }
 
 // hash function to get a value given an atom, assuing atom is represented using char TODO
 __device__ int myhash(char c) {
	 switch(c) { // assuming c is upper case
		 case 'C':
			 return 1;
		 case 'H':
			 return 2;
		 case 'N':
			 return 3;
		 case 'O':
			 return 4;
	 }
	 return -1;
 }
 
 // get the index of radial environment vector
 __device__ int getIndexR(char c) {
	 return myhash(c) - 1;
 }
 
 // get the index of angular environment vector
 __device__ int getIndexA(char c, char d) {
	 int a = myhash(c);
	 int b = myhash(d);
	 return a * (ATOMTYPE_NUM) + b - 1 - (a-1)*a/2; // tricky..., hard to get it right
 }
 
 // build the radial environment vector, on one thread? TODO
 __device__ void buildrad(double *dists, char *d_atoms, int num, double *d_aevs, Sym_Param *d_params, int *d_indexmap, int offset) {
	 for (int i=threadIdx.x + 1; i<num; i++) {
		 if (dists[i] >= d_params->rrc)
			 continue;
		 // compute the index for these two atoms
		 int index1 = (offset + threadIdx.x) * d_params->tot_dim  + d_indexmap[getIndexR(d_atoms[offset+i])];
		 int index2 = (offset + i) * d_params->tot_dim + d_indexmap[getIndexR(d_atoms[offset+threadIdx.x])];
		 // compute the radial environment
		 for (int j=0; j<d_params->rad_eta_num; j++) {
			 for (int k=0; k<d_params->rad_rs_num; k++) {
				 double tmp = cutoff(dists[i], d_params->rrc) * exp(-d_params->rad_etas[j] * (dists[i] - d_params->rad_rss[k]) * (dists[i] - d_params->rad_rss[k]));
				 d_aevs[index1++] += tmp;
				 d_aevs[index2++] += tmp;
			 }
		 }
	 }
 }
 
 // build the angular environment vector, on another thread? TODO current implement runs in O(n^3), could be improved to O(k^3), where k is # of neighbors within cutoff
 __device__ void buildang(double *dists, double *d_coords, char *d_atoms, int num, double *d_aevs, Sym_Param *d_params, int *d_indexmap, int offset) {
	 // correct the offset
	 char *atoms = d_atoms + offset;
	 for (int i=0; i<num-1; i++) {
		 // if the atom is the itself or out of cutoff radius, then skip it
		 if (i == threadIdx.x || dists[i] >= d_params->arc)
			 continue;
		 for (int j=i+1; j<num; j++) {
			 if (j == threadIdx.x || dists[j] >= d_params->arc)
				 continue;
			 // get the angle
			 double theta = acos((dists[i] * dists[i] + dists[j] * dists[j] - getDist(d_coords,offset+i,offset+j) * getDist(d_coords,offset+i,offset+j))
								  / (2 * dists[i] * dists[j]));
			 // get the index
			 int key = atoms[i] > atoms[j] ? getIndexA(atoms[j], atoms[i]) : getIndexA(atoms[i], atoms[j]);
			 int index = (offset + threadIdx.x) * d_params->tot_dim + d_indexmap[key];
			 // compute the angular environment
			 for (int ii=0; ii<d_params->ang_zeta_num; ii++) {
				 for (int jj=0; jj<d_params->ang_theta_num; jj++) {
					 for (int kk=0; kk<d_params->ang_eta_num; kk++) {
						 for (int ll=0; ll<d_params->ang_rs_num; ll++) {
							 d_aevs[index++] += cutoff(dists[i], d_params->arc) * cutoff(dists[j], d_params->arc) *
												 pow(2, 1-d_params->ang_zetas[ii]) *
												 pow((1 + cos(theta - d_params->ang_thetas[jj])), d_params->ang_zetas[ii]) *	// cos or cosf TODO
												 exp(-d_params->ang_etas[kk] * pow(((dists[i] + dists[j]) / 2 - d_params->ang_rss[ll]), 2));
						 }
					 }
				 }
			 }
		 }
	 }
 }
 
 /*************************** Kernel function ****************************/
 __global__ void buildaev(double *d_coords, char *d_atoms, int *d_nums, double *d_aevs, Sym_Param *d_params, int *d_indexmap) {	// char * for atoms is not enough TODO; is it good to sharevarialbe d_para?
	 // prepare the distance vector for individual atom
	 int num = d_nums[blockIdx.x];
	 if (threadIdx.x >= num)
		 return;
	 int offset = 0;
	 for (int i=0; i<blockIdx.x; i++) {
		 offset += d_nums[i];
	 }
	 double *dists = new double[num];
	 for (int i=0; i<num; i++) {
		 dists[i] = getDist(d_coords, offset + i, offset + threadIdx.x);
	 }
	 // compute the radial aev
	 buildrad(dists, d_atoms, num, d_aevs, d_params, d_indexmap, offset);
	 // compute the angular aev
	 buildang(dists, d_coords, d_atoms, num, d_aevs, d_params, d_indexmap, offset);
	 // free the space of dists
	 delete [] dists;
 } 