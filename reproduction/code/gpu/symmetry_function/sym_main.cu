#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include "sym_func.h"
#include "../chemdata/chemdata.h"

using namespace std;

#define BLOCKS 10000 // how many molecules to process one time, max is 65535 ? TODO
#define THREADS 512  // max number 1024?

int determineThreadNum(int n) {
	vector<int> nums = {32, 64, 128, 256, 512};
	for (int i=0; i<nums.size(); i++)
		if (n <= nums[i])
			return nums[i];
	return nums.back();
}

/*************************** Device functions ***************************/
int main(int argc, char *argv[])
{
	// check the input argument
	if (argc != 3) {
		cout << "usage: <parameter-file> <input-xyz>" << endl;
		return -1;
	}

	string xyzfnameB, atomsfnameB, metasfnameB, paramfname, aevsfname;
	paramfname  = paramfname + argv[1] + ".parm";
	xyzfnameB   = xyzfnameB + argv[2] + ".xyz";
	atomsfnameB = atomsfnameB + argv[2] + ".atom";
	metasfnameB = metasfnameB + argv[2] + ".meta";
	aevsfname   = aevsfname + argv[2] + ".aev";
	// set input and output streams
	ifstream xyzI (xyzfnameB, ios::binary), atomsI (atomsfnameB, ios::binary), metasI (metasfnameB, ios::binary);
	ofstream aevsO (aevsfname, ios::binary);
	// check whether the streams are currently associated with the corresponding files
	assert(xyzI.is_open() && atomsI.is_open() && metasI.is_open() && aevsO.is_open());
	// set the precision for output TODO binary not needed
	// aevsO.precision(10);
	// set the parameters
	Sym_Param params;
	createSymParamAuto(paramfname, &params);
	printParam(&params);
	// capture the start time
	cudaEvent_t start, stop;
	HANDLE_ERROR(cudaEventCreate(&start));
	HANDLE_ERROR(cudaEventCreate(&stop));
	HANDLE_ERROR(cudaEventRecord(start, 0));
	// allocate and initialize d_params on device
	Sym_Param *d_params;
	HANDLE_ERROR(cudaMalloc((void **)&d_params, sizeof(Sym_Param)));
	memcpyParam(&params, d_params);
	// build the index map and copy it to device
	int indexmap_size = params.atomtype_num + params.atomtype_num * (params.atomtype_num + 1) / 2;
	int *indexmap = new int[indexmap_size];
	buildIndexmap(indexmap, &params);
	printIndexmap(indexmap);
	int *d_indexmap;
	HANDLE_ERROR(cudaMalloc((void **)&d_indexmap, indexmap_size * sizeof(int)));
	HANDLE_ERROR(cudaMemcpy(d_indexmap, indexmap, indexmap_size * sizeof(int), cudaMemcpyHostToDevice));

	// read the ChemData structure from meta data
	ChemData dsmeta, dsatom, dsxyz;
	metasI.read((char *)&dsmeta, sizeof(ChemData));
	atomsI.read((char *)&dsatom, sizeof(ChemData));
	xyzI.read((char *)&dsxyz, sizeof(ChemData));
	// write the ChemData for aevs
	ChemData dsaevs(dsatom.row_num, params.tot_dim, sizeof(double));
	aevsO.write((char *)&dsaevs, sizeof(ChemData));
	// process the input xyz file
	unsigned long mol_num = 0, tot_mols = 0;
	while (tot_mols < dsmeta.row_num) {
		// get the number of molecules need to be read
		tot_mols += BLOCKS;
		mol_num  = (tot_mols > dsmeta.row_num) ? (dsmeta.row_num + BLOCKS - tot_mols) : BLOCKS;
		// get the total number of atoms need to be read
		int *metas = (int *)malloc(mol_num * sizeof(int));
		metasI.read((char *)metas, mol_num * sizeof(int));
		unsigned long atom_num = 0;
		int threads = 0;	// number of threads for kernel launch
		for (unsigned long i=0; i<mol_num; i++) {
			atom_num += metas[i];
			threads = max(threads, metas[i]);
		}
		cout << "atom_num: " << atom_num << endl;
		// define the sizes
		unsigned long size_coords = atom_num * 3 * sizeof(double);
		unsigned long size_atoms  = atom_num * sizeof(char);
		unsigned long size_aevs   = atom_num * params.tot_dim * sizeof(double);
		unsigned long size_metas  = mol_num * sizeof(int);
		// allocate space
		char *atoms    = (char *)malloc(size_atoms);
		double *coords = (double *)malloc(size_coords);
		double *aevs   = (double *)malloc(size_aevs);
		if (!coords || !atoms || !aevs)
			printf("malloc fail at line: %d\n", __LINE__);
		// read the atoms and cartesian coordinates
		atomsI.read((char *)atoms, atom_num * sizeof(char));
		xyzI.read((char *)coords, atom_num * 3 * sizeof(double));
		
		printf("tot_atom: %d, params.tot_dim: %d, sizeof(double): %d\n", atom_num, params.tot_dim, sizeof(double));
		// initializations for aevs
		memset(aevs, 0, size_aevs);
		// allocate memory and cpy data on device
		int *d_metas;
		char *d_atoms;
		double *d_coords, *d_aevs;
		HANDLE_ERROR(cudaMalloc((void **)&d_coords, size_coords));
		HANDLE_ERROR(cudaMalloc((void **)&d_atoms, size_atoms));
		HANDLE_ERROR(cudaMalloc((void **)&d_aevs, size_aevs));
		HANDLE_ERROR(cudaMalloc((void **)&d_metas, size_metas));
		HANDLE_ERROR(cudaMemcpy(d_coords, coords, size_coords, cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(d_atoms, atoms, size_atoms, cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(d_aevs, aevs, size_aevs, cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(d_metas, metas, size_metas, cudaMemcpyHostToDevice));
		// number of threads
		threads = determineThreadNum(threads);
		cout << "mol_num: " << mol_num << " threads: " << threads << endl;
		printf("aevs: %p, d_aevs: %p, size_aevs: %ld\n", aevs, d_aevs, size_aevs);
		// kernel launch
		buildaev<<<mol_num, threads>>>(d_coords, d_atoms, d_metas, d_aevs, d_params, d_indexmap);
		// copy back the result
		HANDLE_ERROR(cudaMemcpy(aevs, d_aevs, size_aevs, cudaMemcpyDeviceToHost));
		//HANDLE_ERROR(cudaMemcpy(aevs, d_aevs, 112 * sizeof(double), cudaMemcpyDeviceToHost));
		// write the result
		aevsO.write((char *)aevs, size_aevs);
		/*for (long i=0; i<atom_num; i++) {
			aevsO << atoms[i] << "\t";
			for (long j=0; j<params.tot_dim; j++) {
				aevsO << aevs[params.tot_dim*i + j] << "\t";
			}
			aevsO << "\n";
		}*/
		// clean up on Host and Device
		free(coords); free(atoms); free(aevs); free(metas);
		HANDLE_ERROR(cudaFree(d_coords));
		HANDLE_ERROR(cudaFree(d_atoms));
		HANDLE_ERROR(cudaFree(d_aevs));
		HANDLE_ERROR(cudaFree(d_metas));
	}

	// final clean up
	delete [] indexmap;
	HANDLE_ERROR(cudaFree(d_indexmap));
	// freeParam(&params, d_params);
	// VfreeParam(&params);

	// get stop time, and display the timing results
	HANDLE_ERROR(cudaEventRecord(stop, 0));
	HANDLE_ERROR(cudaEventSynchronize(stop));

	float elapsedTime;
	HANDLE_ERROR(cudaEventElapsedTime(&elapsedTime, start, stop));
	printf("Time to generate: %3.1f ms\n", elapsedTime);
	HANDLE_ERROR(cudaEventDestroy(start));
	HANDLE_ERROR(cudaEventDestroy(stop));

	return 0;
}
