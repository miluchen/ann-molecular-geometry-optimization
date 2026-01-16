#include "dataset.h"
#include "cuchemnet.h"
#include <cublas_v2.h>
#include <fstream>
#include <vector>
#include <random>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <cstring>

using namespace std;

/************* Member Functions of DataSet ****************/
DataSet::DataSet(): net_count(0), inputs(0), outputs(0), metas(0), _meta(0), presum_metas(0), d_inputs(0), d_outputs(0), 
					dim_meta(0), count(0), batch_size(0), pos(0), dim_ins(0), dim_outs(0),
					input_grads(0), output_grads(0), d_input_grads(0), d_output_grads(0), _gmeta(0), isGrad(false) {}

DataSet::~DataSet() {
	printf("START DESTORY DataSet!\n");
	for (int i=0; i<net_count; i++) {
		free(inputs[i]);
	}
	free(inputs);	
	free(outputs);
	free(metas);
	free(_meta);
	free(presum_metas);
	for (int i=0; i<net_count; i++) {
		HANDLE_ERROR(cudaFree(d_inputs[i]));
	}
	free(d_inputs);
	HANDLE_ERROR(cudaFree(d_outputs));
	// free the gradients info
	if (isGrad) {
		for (int i=0; i<net_count; i++)
			free(input_grads[i]);
		free(input_grads);
		free(output_grads);
		free(_gmeta);
		for (int i=0; i<net_count; i++)
			HANDLE_ERROR(cudaFree(d_input_grads[i]));
		free(d_input_grads);
		HANDLE_ERROR(cudaFree(d_output_grads));
		isGrad = false;
	}
	printf("END DESTORY DataSet!\n");
}

void DataSet::create(int _net_count, vector<string> &_inputfiles, string &_outputfile, string &_metafile, long *_dim1s, long *_dim2s, long _count) {
	net_count  = _net_count; // number of inputs == number of nets == dim_meta
	dim_meta   = _net_count;
	count      = _count; // count should be >= batch_size
	dim_ins    = _dim1s; // data are stored in config TODO
	dim_outs   = _dim2s; // data are stored in config TODO output number for each individual net
	batch_size = 0;
	pos        = 0;
	// allocate space
	inputs       = (double **)malloc(net_count * sizeof(double *));
	presum_metas = (long *)malloc(count * dim_meta * sizeof(long));
	d_inputs     = (double **)malloc(net_count * sizeof(double *));
	memset(d_inputs, 0, net_count * sizeof(double *));
	// read in data and allocate space
	auto start = std::chrono::high_resolution_clock::now();
	for (int i=0; i<net_count; i++) {
		readData(_inputfiles[i], (char **)(inputs + i));
	}
	readData(_outputfile, (char **)(&outputs));
	readData(_metafile, (char **)(&metas));
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Energy File: " << _outputfile << std::endl;
	std::cout << "Elapsed time for reading in data: " << elapsed.count() << " s\n";
	computePresum(presum_metas, metas, count, dim_meta);
	std::cout << "Computing presum complete!" << std::endl;
}

void DataSet::computePresum(long *a, int *b, long num, int dim) {
	for (int i=0; i<dim; i++)
		a[i] = b[i];
	for (long i=1; i<num; i++) {
		for (int j=0; j<dim; j++) {
			a[i*dim + j] = a[(i-1)*dim + j] + b[i*dim + j];
		}
	}
}

// read in double numerics from a file
void DataSet::readData(std::string &_fname, char **_location) {
	std::cout << _fname << std::endl;
	std::ifstream is (_fname, ios::binary);
	assert(is.is_open());

	ChemData ds;
	is.read((char *)&ds, sizeof(ChemData));
	*_location = (char *)malloc(ds.getSize());
	is.read(*_location, ds.getSize());
	cout << "print ds for " << _fname << endl;
	cout << ds.row_num << " " << ds.col_num << " " << ds.element_size << endl;
	std::cout << "size of content: " << _fname << " is " << ds.getSize() << std::endl;
}

void DataSet::printInfo(void) {
	std::cout << "PRINT INFO OF DATASET: " << std::endl;
	std::cout << "count: " << count << std::endl;
	std::cout << "net_count: " << net_count << std::endl;
	std::cout << "dim_meta: " << dim_meta << std::endl;
	std::cout << "batch_size: " << batch_size << std::endl;
	std::cout << "pos: " << pos << std::endl;
	std::cout << "dim_ins: ";
	for (int i=0; i<net_count; i++) {
		std::cout << dim_ins[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "dim_outs: ";
	for (int i=0; i<net_count; i++)
		std::cout << dim_outs[i] << " ";
	std::cout << std::endl;
}

// print the data of a data set
void DataSet::printData(void) {
	// return;
	std::cout << "PRINTING META DATA:" << std::endl;
	for (long i=0; i<count; i++) {
		for (long j=0; j<dim_meta; j++) {
			std::cout << metas[i * dim_meta + j] << "\t" << std::flush;
		}
		std::cout << std::endl;
	}
	std::cout << "PRINTING PRESUM-META DATA:" << std::endl;
	for (long i=0; i<count; i++) {
		for (long j=0; j<dim_meta; j++) {
			std::cout << presum_metas[i * dim_meta + j] << "\t" << std::flush;
		}
		std::cout << std::endl;
	}
	std::cout << "PRINTING INPUT DATA:" << std::endl;
	for (long i=0; i<count; i++)
		printMolecule(i);
	/*long *_sizes = presum_metas + (count - 1) * dim_meta;
	for (int i=0; i<net_count; i++) {
		std::cout << "INPUT DATA " << i << ":" << std::flush;
		for (long j=0; j<_sizes[i] * dim_ins[i]; j++) {
			std::cout << inputs[i][j] << " " << std::flush;
		}
		printf("\n");
	}*/
	std::cout << "PRINTING OUTPUT DATA:" << std::endl;
	for (long i=0; i<count; i++) {
		std::cout << outputs[i] << "\t" << std::flush;
	}
	std::cout << std::endl;
}

// print the ind-th molecule, index starts from 0, useful for error checking
void DataSet::printMolecule(long ind) {
	std::cout << "PRINTING MOLECULE: " << ind << std::endl;
	for (int i=0; i<dim_meta; i++) {
		for (long j=presum_metas[ind*dim_meta + i] - metas[ind * dim_meta + i]; j<presum_metas[ind*dim_meta + i]; j++) {
			for (long k=0; k<dim_ins[i]; k++) {
				std::cout << inputs[i][j*dim_ins[i] + k] << "\t" << std::flush;
			}
			std::cout << std::endl;
		}
	}
}

// update the d_inputs and d_outputs to po longthe next batch
void DataSet::nextBatch(long _size) {
	// batch size must be smaller than total count of molecules
	assert(_size <= count);
	// new batch size is different
	if (_size != batch_size) {
		std::cout << "New batch size: " << _size << ", Old batch size: " << batch_size << std::endl;
		batch_size = _size;
		free(_meta);
		_meta = (long *)malloc(_size * dim_meta * sizeof(long));
		HANDLE_ERROR(cudaFree(d_outputs));
		HANDLE_ERROR(cudaMalloc((void **)&d_outputs, _size * sizeof(double)));
	}
	if (pos == count) {
		std::cout << "NOTE: reached the end of inputs, need to apply random shuffle" << std::endl;
		randomShuffle();
		pos = 0;
	}
	// if there is no enough data for next batch, roll back
	if (pos + batch_size > count) {
		pos = count - batch_size;
	}
	// compute sizes of new inputs -> _meta data (presum)
	computePresum(_meta, metas + pos * dim_meta, _size, dim_meta);
	/*std::cout << "PRINT _meta: " << std::endl;
	std::cout << "_size : " << _size << std::endl;
	for (long i=0; i<_size; i++) {
		for (int j=0; j<net_count; j++)
			std::cout << _meta[i * net_count + j] << " " << std::flush;
		std::cout << std::endl;
	}*/
	long *_sizes = _meta + (_size - 1) * dim_meta;
	// free and allocate space of d_inputs, and copy data TODO it may not be necessary some time, too many free and malloc cuda
	// std::cout << "BEFORE D_OUTPUTS" << std::endl;
	// printData();
	for (int i=0; i<net_count; i++) {
		HANDLE_ERROR(cudaFree(d_inputs[i]));
		HANDLE_ERROR(cudaMalloc((void **)(d_inputs+i), _sizes[i] * dim_ins[i] * sizeof(double)));
		// double *_hostinput = inputs[i] + ((pos == 0) ? presum_metas[i] : presum_metas[(pos-1)*dim_meta + i]) * dim_ins[i]; // TODO check again!!! TODO TODO TODO
		double *_hostinput = inputs[i] + ((pos == 0) ? 0 : presum_metas[(pos-1)*dim_meta + i]) * dim_ins[i]; // TODO check again!!!
		HANDLE_ERROR(cudaMemcpy(d_inputs[i], _hostinput, _sizes[i] * dim_ins[i] * sizeof(double), cudaMemcpyHostToDevice));
	}
	// copy energies
	HANDLE_ERROR(cudaMemcpy(d_outputs, outputs + pos, _size * sizeof(double), cudaMemcpyHostToDevice));
	// next batch gradients
	if (isGrad) {
		for (int i=0; i<net_count; i++) {
			HANDLE_ERROR(cudaFree(d_input_grads[i]));
			HANDLE_ERROR(cudaMalloc((void **)(d_input_grads + i), _sizes[i] * dim_ins[i] * sizeof(double)));
			double *_hostgrad = input_grads[i] + ((pos == 0) ? 0 : presum_metas[(pos-1)*dim_meta + i]) * dim_ins[i];
			HANDLE_ERROR(cudaMemcpy(d_input_grads[i], _hostgrad, _sizes[i] * dim_ins[i] * sizeof(double), cudaMemcpyHostToDevice));
		}
		HANDLE_ERROR(cudaMemcpy(d_output_grads, output_grads + pos, _size * sizeof(double), cudaMemcpyHostToDevice));
	}
	// update pos
	pos += batch_size;
}

// random shuffle the inputs and outputs, but the corresponding data should be aligned
void DataSet::randomShuffle(void) {
	// TODO TODO TODO
	// return;
	auto start = std::chrono::high_resolution_clock::now();
	std::vector<long> perm(count);
	// initializes permutation to 0, 1, 2, ...
	iota(perm.begin(), perm.end(), 0);
	// random shuffle it
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(perm.begin(), perm.end(), g);
	// build random permutation of inputs of different nets accoding to perm
	long *_sizes = presum_metas + (count - 1) * dim_meta;
	for (int i=0; i<net_count; i++) {
		long *atom_perm = (long *)malloc(_sizes[i] * sizeof(long));
		long index = 0;
		for (long j=0; j<count; j++) {
			long offset = perm[j] * dim_meta + i;
			for (long k=presum_metas[offset] - metas[offset]; k<presum_metas[offset]; k++)
				atom_perm[index++] = k;
		}
		// apply this permutation to inputs and outputs
		for (long j=0; j<_sizes[i]; j++) {
			// check if element at index i has not been moved by checking if perm[i] is nonnegative
			long next = j;
			while (atom_perm[next] >= 0 && atom_perm[atom_perm[next]] >= 0) {
				dataSwap(inputs[i], next, atom_perm[next], dim_ins[i]);
				if (isGrad) {
					// dataSwap(input_grads[i], next, atom_perm[next], dim_ins[i]); TODO
				}
				long tmp = atom_perm[next];
				atom_perm[next] -= _sizes[i];
				next = tmp;
			}
		}
		free(atom_perm);
	}
	// permutate outputs (energies) and meta data: metas
	for (long i=0; i<count; i++) {
		long next = i;
		while (perm[next] >= 0 && perm[perm[next]] >= 0) {
			dataSwap(outputs, next, perm[next], 1); // dimension of energies is 1
			dataSwap(metas, next, perm[next], dim_meta);
			if (isGrad) {
				dataSwap(output_grads, next, perm[next], 1);
				dataSwap(GGmetas, next, perm[next], 1);
			}
			long tmp = perm[next];
			perm[next] -= perm.size();
			next = tmp;
		}
		if (perm[next] >= 0) // to make it negative so that we could recover it easily later
			perm[next] -= perm.size();
	}
	// compute presum
	computePresum(presum_metas, metas, count, dim_meta);
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "RandomShuffle time: " << elapsed.count() << " s\n";
	// CHECK
	/*std::cout << "Permutation: " << std::endl;
	for (long i=0; i<count; i++)
		std::cout << perm[i] + perm.size() << " " << std::flush;
	std::cout << std::endl;
	printData();*/
}

void DataSet::dataSwap(double *v, long i, long j, long dim) {
	for (long k=0; k<dim; k++) {
		double tmp = v[i * dim + k];
		v[i * dim + k] = v[j * dim + k];
		v[j * dim + k] = tmp;
	}
}

void DataSet::dataSwap(int *v, long i, long j, int dim) {
	for (long k=0; k<dim; k++) {
		long tmp = v[i * dim + k];
		v[i * dim + k] = v[j * dim + k];
		v[j * dim + k] = tmp;
	}
}

// read in the gradient files
void DataSet::readGradients(vector<string> &_inputfiles, string &_outputfile) {
	// need: gbatch_size, gmetas, presum_gmetas, _gmeta
	isGrad = true; // set the bool, which will be used when free the space
	input_grads   = (double **)malloc(net_count * sizeof(double *));
	d_input_grads = (double **)malloc(net_count * sizeof(double *));
	GGmetas       = (long *)malloc(count * dim_meta * sizeof(long));
	memset(d_input_grads, 0, net_count * sizeof(double *));
	// read in data
	for (int i=0; i<net_count; i++) {
		readData(_inputfiles[i], (char **)(input_grads + i));
	}
	readData(_outputfile, (char **)(&output_grads));
	// we can build GGmetas and gmetas from metas TODO or read it from a file?
	buildGGmetas(GGmetas, metas, count, dim_meta); // useful for computing _gm
	// readData(_metafile, (char **)(&gmetas));
}

void DataSet::buildGGmetas(long *a, int *b, long num, int dim) {
	for (int i=0; i<num; i++) {
		int _size = 0; // compute the number of atoms in current molecule
		for (int j=0; j<dim; j++)
			_size += b[i * dim + j];
		_size *= 3; // 3 correponds to x,y,z directions
		for (int j=0; j<dim; j++)
			a[i * dim + j] = b[i * dim + j] * _size + (i == 0 ? 0 : a[(i-1) * dim + j]);
	}
}

// a is _gmeta, b is metas, num is number of molecules, dim is the dimension of metadata
void DataSet::build_gmeta(long *a, int *b, long num, int dim) {
	long index = 0; // keep track the row number of 'a'
	for (long i=0; i<num; i++) {
		int _size = 0;
		for (int j=0; j<dim; j++)
			_size += b[i*dim + j];
		_size *= 3;
		while (_size-- > 0) {
			for (int j=0; j<dim; j++) {
				a[index * dim + j] = b[i * dim + j] + (index == 0 ? 0 : a[(index-1) * dim + j]);
			}
			index++;
		}
	}
}
