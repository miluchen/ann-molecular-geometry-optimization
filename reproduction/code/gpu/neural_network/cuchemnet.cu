/*
 *	Source code for back propagation and Adam algorithm of a fully connected neural net.
 */

#include <math.h>
#include <random>
#include <vector>
#include <cstring>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <limits>
#include <fstream>
#include <cublas_v2.h>

#include "cuchemnet.h"

#define THREADS 512
#define MIN(a, b) (((a) > (b)) ? (b) : (a))

// optimization TODO we can get rid of d_adam
#define ADAM_BETA1 0.9
#define ADAM_BETA2 0.999
#define ADAM_EPSILON (1e-8)
// parameter for exponential cost
#define TAU 0.5

// choose the activation function, the output layer is linear function (identity activation function)
#define SIGMOID true
#define GAUSSIAN false
#define RELU false
#define TANH false
#define ELU false	// exponential linear unit
#define LRELU false	// Leaky rectified linear unit
// choose the cost function
#define SQUARECOST true
#define EXPCOST false
// whether to use dropout and the probability
#define DROPOUT false
#define DROPOUT_P 1.0
// whether to use max norm regularization
#define MAXNORM true
#define MAXLEN 3.0

using namespace std;

/************* Helper Functions *************/
void printhelper(double *w, long line) {
	if (!DEBUG) {
		return;
	}
	double *hw = (double *)malloc(sizeof(double));
	HANDLE_ERROR(cudaMemcpy(hw, w, sizeof(double), cudaMemcpyDeviceToHost));
	std::cout << "printhelper in line " << line << ": " << *hw << std::endl;
	free(hw);
}

void printarray(double *w, long size, long line) {
	if (!DEBUG) {
		return;
	}
	double *hw = (double *)malloc(size * sizeof(double));
	std::cout << "LINE: " << line << std::endl;
	HANDLE_ERROR(cudaMemcpy(hw, w, size * sizeof(double), cudaMemcpyDeviceToHost));
	for (long i=0; i<size; i++) {
		std::cout << hw[i] << "\t" << std::flush;
	}
	std::cout << std::endl;
	free(hw);
}

void printarray(double **w, long size, long line) {
	double **hw = (double **)malloc(size * sizeof(double *));
	std::cout << "LINE: " << line << std::endl;
	HANDLE_ERROR(cudaMemcpy(hw, w, size * sizeof(double), cudaMemcpyDeviceToHost));
	for (long i=0; i<size; i++) {
		printf("%p\t", hw[i]);
	}
	std::cout << std::endl;
	free(hw);
}

void printarray(long *w, long size, long line) {
    if (!DEBUG) {
        return;
    }
    long *hw = (long *)malloc(size * sizeof(long));
    std::cout << "LINE: " << line << std::endl;
    HANDLE_ERROR(cudaMemcpy(hw, w, size * sizeof(long), cudaMemcpyDeviceToHost));
    for (long i=0; i<size; i++) {
        std::cout << hw[i] << "\t" << std::flush;
    }
    std::cout << std::endl;
    free(hw);
}

/************* Device Functions *************/
__device__ double sigmoid(double x) {
	return 1./(1. + exp(-x));
}

__device__ double dsigmoid(double x) {
	return sigmoid(x) * (1 - sigmoid(x));
}

__device__ double gaussian(double x) {
	return exp(-x * x);
}

__device__ double dgaussian(double x) {
	return -2 * x * exp(-x * x);
}

/************* Kernel Functions *************/
__global__ void activation(double *outputs, long _neuroncount) {
	// sigmoid function
	long index  = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _neuroncount) {
		outputs[index] = 1./(1. + exp(-outputs[index]));
		index += blockDim.x * gridDim.x;
	}
}

// sigmoid activation function
__global__ void activationSigmoid(double *outputs, long _size) {
	long index  = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _size) {
		outputs[index] = 1./(1. + exp(-outputs[index]));
		index += blockDim.x * gridDim.x;
	}
}

__global__ void activationSigmoid(double *outputs, double *derivatives, long _size) {
	double tmp;
	long index  = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _size) {
		tmp = 1./(1. + exp(-outputs[index]));
		outputs[index]     = tmp;
		derivatives[index] = tmp * (1 - tmp);
		index += blockDim.x * gridDim.x;
	}
}

// gaussian activation function
__global__ void activationGaussian(double *outputs, long _size) {
	long index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _size) {
		// is it faster this way? with one less memory access? TODO testings are needed
		outputs[index] = exp(-(outputs[index] * outputs[index]));
		index += blockDim.x * gridDim.x;
	}
}

__global__ void activationGaussian(double *outputs, double *derivatives, long _size) {
	double tmp;
	long index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _size) {
		// is it faster this way? with one less memory access? TODO testings are needed
		tmp = outputs[index];
		outputs[index]     = exp(-(tmp * tmp));
		derivatives[index] = -2 * tmp * exp(-(tmp * tmp));
		index += blockDim.x * gridDim.x;
	}
}

__global__ void activationRelu(double *outputs, long _size) {
	long index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _size) {
		if (outputs[index] <= 0)
			outputs[index] = 0;
		index += blockDim.x * gridDim.x;
	}
}

__global__ void activationRelu(double *outputs, double *derivatives, long _size) {
	long index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _size) {
		if (outputs[index] <= 0) {
			outputs[index] = 0;
			derivatives[index] = 0;
		}
		else {
			derivatives[index] = 1;
		}
		index += blockDim.x * gridDim.x;
	}
}

__global__ void activationLrelu(double *outputs, long _size) {
	long index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _size) {
		if (outputs[index] <= 0)
			outputs[index] = 0.01 * outputs[index];
		index += blockDim.x * gridDim.x;
	}
}

__global__ void activationLrelu(double *outputs, double *derivatives, long _size) {
	long index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _size) {
		if (outputs[index] <= 0) {
			outputs[index] = 0.01 * outputs[index];
			derivatives[index] = 0.01;
		}
		else {
			derivatives[index] = 1;
		}
		index += blockDim.x * gridDim.x;
	}
}

__global__ void activationElu(double *outputs, long _size) {
	long index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _size) {
		double tmp = outputs[index];
		if (tmp < 0)
			outputs[index] = exp(tmp) - 1;
		index += blockDim.x * gridDim.x;
	}
}

__global__ void activationElu(double *outputs, double *derivatives, long _size) {
	long index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _size) {
		double tmp = outputs[index];
		if (tmp < 0) {
			tmp = exp(tmp) - 1;
			outputs[index] = tmp;
			derivatives[index] = tmp + 1;
		}
		else {
			derivatives[index] = 1;
		}
		index += blockDim.x * gridDim.x;
	}
}

__global__ void activationTanh(double *outputs,long _size) {
	long index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _size) {
		outputs[index] = 2. / (1. + exp(-2 * outputs[index])) - 1;
		index += blockDim.x * gridDim.x;
	}
}

__global__ void activationTanh(double *outputs, double *derivatives, long _size) {
	double tmp;
	long index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _size) {
		tmp = 2. / (1. + exp(-2 * outputs[index])) - 1;
		outputs[index] = tmp;
		derivatives[index] = 1 - tmp * tmp;
		index += blockDim.x * gridDim.x;
	}
}

__global__ void outputActivation(double *outputs, long _neuroncount) {
	 long index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _neuroncount) {
		outputs[index] = outputs[index]; // a linear activation function
		index += blockDim.x * gridDim.x;
	}
}

// an identity activation function, outputs[index] = outputs[index]
__global__ void activationIdentity(double *outputs, long _size) {
	return;
}

__global__ void activationIdentity(double *outputs, double *derivatives, long _size) {
	 long index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _size) {
		// outputs[index] = outputs[index]; // a linear activation function
		derivatives[index] = 1.0;
		index += blockDim.x * gridDim.x;
	}
}

__global__ void Grad_activation(double *outputs, double *derivatives, long _size) {
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _size) {
		outputs[index] *= derivatives[index];
		index          += blockDim.x * gridDim.x;
	}
}

/*__global__ void updateWeightsAdam(double *weightsm, double *weightsv, double *weights, double *weights_grad, Adam_Setting *d_adam, double alphat, long _size) {
	 long index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _size) {
		weightsm[index] = d_adam->beta1 * weightsm[index] + (1-d_adam->beta1) * weights_grad[index];
		weightsv[index] = d_adam->beta2 * weightsv[index] + (1-d_adam->beta2) * weights_grad[index] * weights_grad[index];
		weights[index] -= alphat * weightsm[index] / (sqrt(weightsv[index]) + d_adam->epsilon);
		index          += blockDim.x * gridDim.x;
	}
}*/
__global__ void updateWeightsAdam(double *weightsm, double *weightsv, double *weights, double *weights_grad, double alphat, long _size) {
	 long index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _size) {
		weightsm[index] = ADAM_BETA1 * weightsm[index] + (1-ADAM_BETA1) * weights_grad[index];
		weightsv[index] = ADAM_BETA2 * weightsv[index] + (1-ADAM_BETA2) * weights_grad[index] * weights_grad[index];
		weights[index] -= alphat * weightsm[index] / (sqrt(weightsv[index]) + ADAM_EPSILON);
		index          += blockDim.x * gridDim.x;
	}
}

__global__ void propagateError(double *errorvalues_t, double *derivatives, long n) {
	 long index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < n) {
		errorvalues_t[index] *= derivatives[index];
		index += blockDim.x * gridDim.x;
	}
}

// accumulates the outputs from different nets to give the total energies
__global__ void sumEnergy(double **dc_outputs, double *e_outputs, long _m, long *d_meta) {
	long index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _m) {
		double _tmp = 0;
		for (int i=0; i<ATOMTYPE_NUM; i++) {
			for (int j = (index == 0) ? 0 : d_meta[(index-1) * ATOMTYPE_NUM  + i]; j < d_meta[index * ATOMTYPE_NUM + i]; j++) {
				_tmp += dc_outputs[i][j];
				// e_outputs[index] += dc_outputs[i][j];
			}
		}
		e_outputs[index] = _tmp;
		index += blockDim.x * gridDim.x;
	}
}

// need optimization TODO
__global__ void cudaOutputError(double *e_outputs, long *d_meta, double **d_errorvalues) {
	// thread index is the atom index of a molecule, block index is the index of molecule
	/*int tot = 0; // now necessary for now, since the molecule size is small TODO
	for (int i=0; i<ATOMTYPE_NUM; i++) {
		tot += d_meta[blockIdx.x * ATOMTYPE_NUM + i];
	}
	while (index1 < tot) {*/
	if (threadIdx.x == 0) {
		for (int i=0; i<ATOMTYPE_NUM; i++) {
			long j = (blockIdx.x == 0) ? 0 : d_meta[(blockIdx.x - 1) * ATOMTYPE_NUM + i];
			long size = (blockIdx.x == 0) ? d_meta[blockIdx.x * ATOMTYPE_NUM + i] : (d_meta[blockIdx.x * ATOMTYPE_NUM + i] - d_meta[(blockIdx.x-1) * ATOMTYPE_NUM + i]);
			if (size == 0)	// it's possible that size is 0, for example: CH4, there is no N or O
				continue;
			for (; j<d_meta[blockIdx.x * ATOMTYPE_NUM + i]; j++) {
				// d_errorvalues[i][j] = e_outputs[blockIdx.x] / size;
				d_errorvalues[i][j] = e_outputs[blockIdx.x];
			}
		}
	}
}

// derivatives will have all ones
__global__ void activationGradIdentity(double *outputs, double *derivatives, long _size, int _dim, long *_meta, int net) {
	int index     = threadIdx.x + blockIdx.x * blockDim.x;
	int mol_index = index / _dim;
	int num_atom  = 0;
	for (int i=0; i<ATOMTYPE_NUM; i++)
		num_atom += _meta[mol_index * ATOMTYPE_NUM + i] - (mol_index == 0 ? 0 : _meta[(mol_index-1) * ATOMTYPE_NUM + i]);
	int step      = (_meta[mol_index * ATOMTYPE_NUM + net] - (mol_index == 0 ? 0 : _meta[(mol_index-1) * ATOMTYPE_NUM + net])) * _dim;
	num_atom     *= 3;
	while (index < _size) {
		for (int i=0; i<num_atom; i++) {
			derivatives[step*i + index] = 1;
		}
		index += blockDim.x * gridDim.x;
	}
}

__global__ void activationGradSigmoid(double *outputs, double *derivatives, long _size, int _dim, long *_meta, int net) {
	int index     = threadIdx.x + blockIdx.x * blockDim.x;
	int mol_index = index / _dim;
	int num_atom  = 0;
	for (int i=0; i<ATOMTYPE_NUM; i++)
		num_atom += _meta[mol_index * ATOMTYPE_NUM + i] - (mol_index == 0 ? 0 : _meta[(mol_index-1) * ATOMTYPE_NUM + i]);
	int step      = (_meta[mol_index * ATOMTYPE_NUM + net] - (mol_index == 0 ? 0 : _meta[(mol_index-1) * ATOMTYPE_NUM + net])) * _dim;
	num_atom     *= 3;
	while (index < _size) {
		double tmp     = dsigmoid(outputs[index]);
		outputs[index] = sigmoid(outputs[index]);
		for (int i=0; i<num_atom; i++) {
			derivatives[index + step * i] = tmp;
		}
		index += blockDim.x * gridDim.x;
	}
}

__global__ void activationGradGaussian(double *outputs, double *derivatives, long _size, int _dim, long *_meta, int net) {
	int index     = threadIdx.x + blockIdx.x * blockDim.x;
	int mol_index = index / _dim;
	int num_atom  = 0;
	for (int i=0; i<ATOMTYPE_NUM; i++)
		num_atom += _meta[mol_index * ATOMTYPE_NUM + i] - (mol_index == 0 ? 0 : _meta[(mol_index-1) * ATOMTYPE_NUM + i]);
	int step      = (_meta[mol_index * ATOMTYPE_NUM + net] - (mol_index == 0 ? 0 : _meta[(mol_index-1) * ATOMTYPE_NUM + net])) * _dim;
	num_atom     *= 3;
	while (index < _size) {
		double tmp     = dgaussian(outputs[index]);
		outputs[index] = gaussian(outputs[index]);
		for (int i=0; i<num_atom; i++) {
			derivatives[index + step * i] = tmp;
		}
		index += blockDim.x * gridDim.x;
	}
}

__global__ void activationGradRelu(double *outputs, double *derivatives, long _size, int _dim, long *_meta, int net) {

}

__global__ void activationGradElu(double *outputs, double *derivatives, long _size, int _dim, long *_meta, int net) {

}

__global__ void activationGradTanh(double *outputs, double *derivatives, long _size, int _dim, long *_meta, int net) {

}

/************* Member Functions of cuchemnet ****************/
cuchemnet::cuchemnet(): c_layers(0), c_layercount(0), c_weights(0), c_wbias(0), c_deltas(0), c_outputs(0), c_weights_grad(0), c_wbias_grad(0),
						c_adam_weightsm(0), c_adam_weightsv(0), c_adam_wbiasm(0), c_adam_wbiasv(0), tot_weights(0), tot_wbias(0), tot_size(0), batch_size(0), c_ones(0), 
						c_derivatives(0), tot_derivatives(0), c_structure(0) {}

cuchemnet::~cuchemnet() {
	printf("START DESTORY cuchemnet!\n");
	delete [] c_layers;
	// free the space on device
	HANDLE_ERROR(cudaFree(c_weights));
	HANDLE_ERROR(cudaFree(c_deltas));
	HANDLE_ERROR(cudaFree(c_outputs));
	HANDLE_ERROR(cudaFree(c_weights_grad));
	HANDLE_ERROR(cudaFree(c_adam_weightsm));
	HANDLE_ERROR(cudaFree(c_adam_weightsv));
	HANDLE_ERROR(cudaFree(c_ones));
	HANDLE_ERROR(cudaFree(c_derivatives));
	HANDLE_ERROR(cudaFree(c_structure));
	printf("END DESTORY cuchemnet!\n");
}

// create the network structure, _inputcount is number of input features
void cuchemnet::create(long _inputcount, long _outputcount, long _hiddenlayercount, long *_hiddenlayers) {
	assert(_inputcount >=1 && _outputcount >=1 && _hiddenlayercount >=0);
	// initialization for the architecture
	c_layercount = _hiddenlayercount + 1;
	c_layers     = new culayer*[c_layercount];
	// build a vector to store the structure of the NN
	std::vector<long> _structure;
	_structure.push_back(_inputcount);
	for (int i=0; i<_hiddenlayercount; i++) {
		_structure.push_back(_hiddenlayers[i]);
	}
	_structure.push_back(_outputcount);
	_structure.push_back(1);	// assuming the output layer has fan_out = 1
	// allocate space for data
	tot_weights = 0, tot_wbias = 0;
	for (int i=0; i<c_layercount; i++) {
		tot_weights += _structure[i] * _structure[i+1];
		tot_wbias   += _structure[i+1];
	}
	std::cout << "tot_weights: " << tot_weights << " tot_wbias: " << tot_wbias << std::endl;
	tot_size = tot_weights + tot_wbias;
	// allocate space for data on device
	HANDLE_ERROR(cudaMalloc((void **)&c_weights, tot_size * sizeof(double)));         // a big chunck containing both weights and wbias
	HANDLE_ERROR(cudaMalloc((void **)&c_weights_grad, tot_size * sizeof(double)));    // will be initialized to 0 in every bpGradient/training step
	HANDLE_ERROR(cudaMalloc((void **)&c_deltas, tot_weights * sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void **)&c_adam_weightsm, tot_size * sizeof(double)));
	HANDLE_ERROR(cudaMalloc((void **)&c_adam_weightsv, tot_size * sizeof(double)));
	c_wbias       = c_weights + tot_weights;
	c_wbias_grad  = c_weights_grad + tot_weights; // will be initialized to '0' during every bpGradient step
	c_adam_wbiasm = c_adam_weightsm + tot_weights;
	c_adam_wbiasv = c_adam_weightsv + tot_weights;
	// initializations on deltas, wbias, and values for adam, all to 0
	HANDLE_ERROR(cudaMemset(c_deltas, 0, tot_weights * sizeof(double)));
	HANDLE_ERROR(cudaMemset(c_wbias, 0, tot_wbias * sizeof(double))); // TODO
	HANDLE_ERROR(cudaMemset(c_adam_weightsm, 0, tot_size * sizeof(double)));
	HANDLE_ERROR(cudaMemset(c_adam_weightsv, 0, tot_size * sizeof(double)));
	HANDLE_ERROR(cudaMemset(c_weights_grad, 0, tot_size * sizeof(double)));
	// create the layers and initialize weights
	createLayersWeights(_structure);
	// build the structure
	buildStructure();
}

void cuchemnet::buildStructure(void) {
	// build the structures info on Host
	int *structure = (int *)malloc((c_layercount + 1) * sizeof(int));
	structure[0] = c_layers[0]->inputcount;
	for (int i=0; i<c_layercount; i++) {
		structure[i+1] = c_layers[i]->neuroncount;
	}
	if (DEBUG) {
		cout << "STRUCTURE OF NET: " << endl;
		for (int i=0; i<=c_layercount; i++)
			cout << structure[i] << ", ";
		cout << endl;
	}
	// copy structure info to Device
	HANDLE_ERROR(cudaMalloc((void **)&c_structure, (c_layercount+1)*sizeof(int)));
	HANDLE_ERROR(cudaMemcpy(c_structure, structure, (c_layercount+1)*sizeof(int), cudaMemcpyHostToDevice));
	// free the space
	free(structure);
}

void cuchemnet::buildOnes(void) {
	double *__tmp = (double *)malloc(batch_size * sizeof(double));
	for (long i=0; i<batch_size; i++)
		__tmp[i] = 1.0;
	HANDLE_ERROR(cudaMemcpy(c_ones, __tmp, batch_size * sizeof(double), cudaMemcpyHostToDevice));
	free(__tmp);
}

// create the layers and initialize weights
void cuchemnet::createLayersWeights(std::vector<long> &_structure) {
	long offset_weights = 0, offset_wbias = 0, offset_outputs = 0;
	for (long i=0; i<c_layercount; i++) {
		c_layers[i] = new culayer;
		long _inputcount = _structure[i], _neuroncount = _structure[i+1], fan_out = _structure[i+2];
		// initializations for layer structure
		c_layers[i]->inputcount   = _inputcount;
		c_layers[i]->neuroncount  = _neuroncount;
		c_layers[i]->weights_size = _inputcount * _neuroncount;
		// initialize pointers
		c_layers[i]->weights       = c_weights + offset_weights;
		c_layers[i]->weights_grad  = c_weights_grad + offset_weights;
		c_layers[i]->deltas        = c_deltas + offset_weights;
		c_layers[i]->wbias         = c_wbias + offset_wbias;
		c_layers[i]->wbias_grad    = c_wbias_grad + offset_wbias;
		c_layers[i]->outputs       = c_outputs + offset_outputs; // TODO necessary? we are not allocating space for outputs at first
		c_layers[i]->adam_weightsm = c_adam_weightsm + offset_weights;
		c_layers[i]->adam_weightsv = c_adam_weightsv + offset_weights;
		c_layers[i]->adam_wbiasm   = c_adam_wbiasm + offset_wbias;
		c_layers[i]->adam_wbiasv   = c_adam_wbiasv + offset_wbias;
		// initialization for weights on host
		double *hweights = (double *)malloc(c_layers[i]->weights_size * sizeof(double));
		// set the uniform distribution, assuming sigmoid function TODO
		// std::default_random_engine gen;
		std::random_device rd; // TODO turn off the randomness just for testing
		std::mt19937 gen(rd());
		double r = 4 * sqrt(6) / sqrt(_inputcount + fan_out + 1);
		std::uniform_real_distribution<double> distribution(-r, r);
		for (long j=0; j<c_layers[i]->weights_size; j++) {
			hweights[j] = distribution(gen);
		}
		// copy weights from host to device
		HANDLE_ERROR(cudaMemcpy(c_layers[i]->weights, hweights, c_layers[i]->weights_size * sizeof(double), cudaMemcpyHostToDevice));
		free(hweights);
		// update offsets
		offset_weights += c_layers[i]->weights_size;
		offset_wbias   += _neuroncount;
		offset_outputs += _neuroncount * batch_size; // TODO remove it?
	}
}

// check _m and allocate required space
void cuchemnet::prepropagate(long _m, PropagateMode _mode, long _gm) {
    // check whether m == batch_size? otherwise cudaFree and Malloc a new space for outputs, need to be redo, since it may happen every run TODO
	if (_m != batch_size) {
        batch_size = _m;
        HANDLE_ERROR(cudaFree(c_outputs));
        // allocate new space
        HANDLE_ERROR(cudaMalloc((void **)&c_outputs, tot_wbias * _m * sizeof(double)));
        // update the output pointers in the layers
        c_layers[0]->outputs = c_outputs;
        for (int i=1; i<c_layercount; i++) {
            c_layers[i]->outputs = c_layers[i-1]->outputs + c_layers[i-1]->neuroncount * _m;
        }
        // also update c_ones
        HANDLE_ERROR(cudaFree(c_ones));
        HANDLE_ERROR(cudaMalloc((void **)&c_ones, _m * sizeof(double)));
        buildOnes();
    }
    if (_mode == WITH_DERIVATIVE && tot_derivatives != tot_wbias * _m) {
		if (DEBUG) {
			std::cout << "WITH_DERIVATIVE at: " << __LINE__ << std::endl;
		}
		allocateDerivatives(_m);
     }
     else if (_mode == GRAD_DERIVATIVE && tot_derivatives != tot_wbias * _gm) {
		if (DEBUG) {
			std::cout << "GRAD_DERIVATIVE at: " << __LINE__ << std::endl;
		}
		allocateDerivatives(_gm);
    }
}

// allocate space for the derivatives if required
void cuchemnet::allocateDerivatives(long _size) {
	assert(_size > 0);
	HANDLE_ERROR(cudaFree(c_derivatives));
	HANDLE_ERROR(cudaMalloc((void **)&c_derivatives, tot_wbias * _size * sizeof(double)));
	c_layers[0]->derivatives = c_derivatives;
	for (int i=1; i<c_layercount; i++) {
		c_layers[i]->derivatives = c_layers[i-1]->derivatives + c_layers[i-1]->neuroncount * _size;
	}
	tot_derivatives = tot_wbias * _size;
}

// calculate and propagate the network, m is the number of inputs
void cuchemnet::propagate(double *d_input, long _m, cublasHandle_t &handle, PropagateMode _mode, int net, long *_meta, long _gm) { // default value for _meta is NULL TODO _gmeta ?
	if (_m == 0)
		return;
	// check _m and batch_size and allocate space if needed
	prepropagate(_m, _mode, _gm);
	// check whether it's gradient propagation
	if (_mode == GRAD_PROPAGATION) {
		Grad_propagate(d_input, _m, handle);
		return;
	}
    // start from the input layer and go to the output layer
    const double _alpha = 1.0, _beta = 1.0, _zerobeta = 0.0;
    for (int i=0; i<c_layercount; i++) {
        // use c_ones to realize a matrix-matrix operation -> copy wbias to outputs, TODO: combine these two operations TODO better to have a kernel to run this or use DeviceToDevice copy
        HANDLE_CUBLAS(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, c_layers[i]->neuroncount, _m, 1, &_alpha,
                                    c_layers[i]->wbias, c_layers[i]->neuroncount, c_ones, _m, &_zerobeta, c_layers[i]->outputs, c_layers[i]->neuroncount));
        // add W * X to ouputs
        HANDLE_CUBLAS(cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, c_layers[i]->neuroncount, _m, c_layers[i]->inputcount, &_alpha,
                                    c_layers[i]->weights, c_layers[i]->inputcount, d_input, c_layers[i]->inputcount, &_beta, c_layers[i]->outputs, c_layers[i]->neuroncount));
        // apply the activation function
		applyActivation(i, _m, _mode, net);
		// update the input
		d_input = c_layers[i]->outputs;
    }
}

void cuchemnet::applyActivation(int i, long _m, PropagateMode _mode, int net, long *_meta) { // 'i' is the layer number
	long _num      = c_layers[i]->neuroncount * _m;
	long threadnum = MIN(THREADS, _num);
	long blocknum  = (_num + threadnum - 1) / threadnum;
	if (i == c_layercount - 1) {
		if (_mode == WITHOUT_DERIVATIVE) {
			activationIdentity<<<blocknum, threadnum>>>(c_layers[i]->outputs, _num);
		}
		else if (_mode == WITH_DERIVATIVE) {
			activationIdentity<<<blocknum, threadnum>>>(c_layers[i]->outputs, c_layers[i]->derivatives, _num);
		}
		else if (_mode == GRAD_DERIVATIVE) {
			activationGradIdentity<<<blocknum, threadnum>>>(c_layers[i]->outputs, c_layers[i]->derivatives, _num, c_layers[i]->neuroncount, _meta, net);
		}
	}
	else {
		if (_mode == WITHOUT_DERIVATIVE) {
			if (SIGMOID)
				activationSigmoid<<<blocknum, threadnum>>>(c_layers[i]->outputs, _num);
			else if (GAUSSIAN)
				activationGaussian<<<blocknum, threadnum>>>(c_layers[i]->outputs, _num);
			else if (RELU)
				activationRelu<<<blocknum, threadnum>>>(c_layers[i]->outputs, _num);
			else if (TANH)
				activationTanh<<<blocknum, threadnum>>>(c_layers[i]->outputs, _num);
			else if (ELU)
				activationElu<<<blocknum, threadnum>>>(c_layers[i]->outputs, _num);
			else if (LRELU)
				activationLrelu<<<blocknum, threadnum>>>(c_layers[i]->outputs, _num);
		}
		else if (_mode == WITH_DERIVATIVE) {
			if (SIGMOID)
				activationSigmoid<<<blocknum, threadnum>>>(c_layers[i]->outputs, c_layers[i]->derivatives, _num);
			else if (GAUSSIAN)
				activationGaussian<<<blocknum, threadnum>>>(c_layers[i]->outputs, c_layers[i]->derivatives, _num);
			else if (RELU)
				activationRelu<<<blocknum, threadnum>>>(c_layers[i]->outputs, c_layers[i]->derivatives, _num);
			else if (TANH)
				activationTanh<<<blocknum, threadnum>>>(c_layers[i]->outputs, c_layers[i]->derivatives, _num);
			else if (ELU)
				activationElu<<<blocknum, threadnum>>>(c_layers[i]->outputs, c_layers[i]->derivatives, _num);
			else if (LRELU)
				activationLrelu<<<blocknum, threadnum>>>(c_layers[i]->outputs, c_layers[i]->derivatives, _num);
		}
		else if (_mode == GRAD_DERIVATIVE) {
			if (SIGMOID)
				activationGradSigmoid<<<blocknum, threadnum>>>(c_layers[i]->outputs, c_layers[i]->derivatives, _num, c_layers[i]->neuroncount, _meta, net);
			else if (GAUSSIAN)
				activationGradGaussian<<<blocknum, threadnum>>>(c_layers[i]->outputs, c_layers[i]->derivatives, _num, c_layers[i]->neuroncount, _meta, net);
		}
	}
}

// for gradient propagation, the bias do not depend on the gradient inputs, so they are not included in the computation
void cuchemnet::Grad_propagate(double *d_input, long _m, cublasHandle_t &handle) {
	const double _alpha = 1.0, _zerobeta = 0.0;
    for (long i=0; i<c_layercount; i++) {
        // add W * X to ouputs
		if (DEBUG) {
			cout << "D_INPUT" << endl;
			printarray(d_input, c_layers[i]->inputcount * _m, __LINE__);
		}
        HANDLE_CUBLAS(cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, c_layers[i]->neuroncount, _m, c_layers[i]->inputcount, &_alpha,
                                    c_layers[i]->weights, c_layers[i]->inputcount, d_input, c_layers[i]->inputcount, &_zerobeta, c_layers[i]->outputs, c_layers[i]->neuroncount));
        // apply the activation function
        long _num      = c_layers[i]->neuroncount * _m;
        long threadnum = MIN(THREADS, _num);
        long blocknum  = (_num + threadnum - 1) / threadnum;
		// the activation for gradient propagation is element-wise multiplications
		if (DEBUG) {
			cout << "OUTPUTS: " << endl;
			printarray(c_layers[i]->outputs, _num, __LINE__);
			cout << "DERIVATIVES: " << endl;
			printarray(c_layers[i]->derivatives, _num, __LINE__);
		}
		Grad_activation<<<blocknum, threadnum>>>(c_layers[i]->outputs, c_layers[i]->derivatives, _num);
		// update offset and d_input pointer
		d_input = c_layers[i]->outputs;
	}
}

// compute the gradiensts using back-propagation given the output layer errors
void cuchemnet::helpbpGradient(double *d_input, double *errorvalues, cublasHandle_t &handle) {
	// apply back propagation
	const double _alpha = 1.0, _beta = 0.0;
	for (long i=c_layercount-1; i>=0; i--) {
		// update layerinput
		double *layerinput = i == 0 ? d_input : c_layers[i-1]->outputs;
		// errorvalues * layerinput -> gradient
		HANDLE_CUBLAS(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_T, c_layers[i]->inputcount, c_layers[i]->neuroncount, batch_size, &_alpha, 
									layerinput, c_layers[i]->inputcount, errorvalues, c_layers[i]->neuroncount, &_beta, c_layers[i]->weights_grad, c_layers[i]->inputcount));
		// wbias_grad = sum(errorvalus) build a (1,1,...,1) vector to speed up the calculation
		HANDLE_CUBLAS(cublasDgemv(handle, CUBLAS_OP_N, c_layers[i]->neuroncount, batch_size, &_alpha, errorvalues, c_layers[i]->neuroncount,
									c_ones, 1, &_beta, c_layers[i]->wbias_grad, 1));
		// weights' * errorvalues * layerinput * (1 - layerinput); not necessary when i == 0
		if (i > 0) {
            // back-propagate the error values to the next layer
			double *errorvalues_t; // TODO make use of space of outputs, no need to allocate new space.. first run kernel to have y*(1-y), then cublasDgemm?
			HANDLE_ERROR(cudaMalloc((void **)&errorvalues_t, c_layers[i]->inputcount * batch_size * sizeof(double))); // don't need to be set 0
			HANDLE_CUBLAS(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, c_layers[i]->inputcount, batch_size, c_layers[i]->neuroncount, &_alpha, c_layers[i]->weights, 
										c_layers[i]->inputcount, errorvalues, c_layers[i]->neuroncount, &_beta, errorvalues_t, c_layers[i]->inputcount)); 
			// call a kernel assuming sigmoid TODO f'(z) * (1-f'(z))
			long threadnum = MIN(c_layers[i]->inputcount * batch_size, THREADS);
			long blocknum  = (c_layers[i]->inputcount * batch_size + threadnum - 1) / threadnum;
			propagateError<<<blocknum, threadnum>>>(errorvalues_t, c_layers[i-1]->derivatives, c_layers[i]->inputcount * batch_size);
			// propagateError<<<blocknum, threadnum>>>(errorvalues_t, layerinput, c_layers[i]->inputcount * batch_size);

			HANDLE_ERROR(cudaFree(errorvalues));
			errorvalues = errorvalues_t;
		}
		else
			HANDLE_ERROR(cudaFree(errorvalues));
	}
}

// compute the gradients based on back-propagation
void cuchemnet::bpGradient(double *d_input, double *errorvalues, long _m, double lambda, cublasHandle_t &handle) {
	if (_m == 0) { // Memset grads to be 0, so that weights can be updated easily TODO if we have a large weights array in the EnergyNet, which contains everything
		HANDLE_ERROR(cudaMemset(c_weights_grad, 0, tot_size * sizeof(double)));
		return;
	}
	// initialize the weights_grad and wbias_grad to be 0 this seems unnecessary, but useful to check the gradients
	// HANDLE_ERROR(cudaMemset(c_weights_grad, 0, tot_size * sizeof(double)));
	// HANDLE_ERROR(cudaMemset(c_wbias_grad, 0, tot_wbias * sizeof(double)));
	// accumulate the gradients without weight decay term, GPU parallel, this is mini-batch, not SGD
	if (DEBUG) {
		std::cout << "ERRORVALUES: " << std::endl;	
		printarray(errorvalues, _m * c_layers[c_layercount-1]->neuroncount, __LINE__);
	}
	helpbpGradient(d_input, errorvalues, handle);

	if (DEBUG) {
		std::cout << "C_WEIGHTS_GRAD AFTER HELP: " << std::endl;
		printarray(c_weights_grad, tot_size, __LINE__);
	}
	// add the weight decay term and average the gradients
	// const double _alpha  = 1.0/_m; TODO???TODO
	// const double _alpha = 1.0;
	const double _lambda = lambda/_m;
	// weights_grad / m
	// HANDLE_CUBLAS(cublasDscal(handle, tot_weights, &_alpha, c_weights_grad, 1));
	// weights_grad / m + lambda * weights
	HANDLE_CUBLAS(cublasDaxpy(handle, tot_weights, &_lambda, c_weights, 1, c_weights_grad, 1));
	// wbias_grad / m
	// HANDLE_CUBLAS(cublasDscal(handle, tot_wbias, &_alpha, c_wbias_grad, 1));
	if (DEBUG) {
		std::cout << "C_WEIGHTS_GRAD AFTER BP: " << std::endl;
		printarray(c_weights_grad, tot_size, __LINE__);
	}
}

// update weights of the current net
/*void cuchemnet::updateAdam(double alphat, Adam_Setting *d_adam) {
	int threadnum = MIN(tot_size, THREADS);
	int blocknum  = (tot_size + threadnum - 1) / threadnum;
	updateWeightsAdam<<<blocknum, threadnum>>>(c_adam_weightsm, c_adam_weightsv, c_weights, c_weights_grad, d_adam, alphat, tot_size);
}*/
void cuchemnet::updateAdam(double alphat) {
	int threadnum = MIN(tot_size, THREADS);
	int blocknum  = (tot_size + threadnum - 1) / threadnum;
	updateWeightsAdam<<<blocknum, threadnum>>>(c_adam_weightsm, c_adam_weightsv, c_weights, c_weights_grad, alphat, tot_size);
}

/*void cuchemnet::updateAdam(double alphat, Adam_Setting *d_adam, long _size) {
	int threadnum = MIN(_size, THREADS);
	int blocknum  = (_size + threadnum - 1) / threadnum;
	updateWeightsAdam<<<blocknum, threadnum>>>(c_adam_weightsm, c_adam_weightsv, c_weights, c_weights_grad, d_adam, alphat, _size);
}*/
void cuchemnet::updateAdam(double alphat, long _size) {
	int threadnum = MIN(_size, THREADS);
	int blocknum  = (_size + threadnum - 1) / threadnum;
	updateWeightsAdam<<<blocknum, threadnum>>>(c_adam_weightsm, c_adam_weightsv, c_weights, c_weights_grad, alphat, _size);
}

// load parameters: weights and wbias, assuming the sizes of _weights and wbias are correct, and caller frees the _weights and _wbias
void cuchemnet::loadParam(double *_weights, long _weights_num, double *_wbias, long _wbias_num) {
	// assuming space are already created for weights, wbias on device
	assert(_weights_num == tot_weights && _wbias_num == tot_wbias);
	HANDLE_ERROR(cudaMemcpy(c_weights, _weights, _weights_num * sizeof(double), cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(c_wbias, _wbias, _wbias_num * sizeof(double), cudaMemcpyHostToDevice));
}

void cuchemnet::loadParam(double *_c_weights, long _c_tot_size) {
	assert(_c_tot_size == tot_size);
	HANDLE_ERROR(cudaMemcpy(c_weights, _c_weights, _c_tot_size * sizeof(double), cudaMemcpyHostToDevice));
}

// store weights and wbias, assuming the caller allocate the space and free _weights later
void cuchemnet::storeParam(double *_weights) {
	HANDLE_ERROR(cudaMemcpy(_weights, c_weights, (tot_size) * sizeof(double), cudaMemcpyDeviceToHost));
}

void cuchemnet::printParam(void) {
    std::cout << "layers weights: " << std::endl;
	double *_weights = (double *)malloc(tot_weights * sizeof(double));
	double *_wbias   = (double *)malloc(tot_wbias * sizeof(double));
	HANDLE_ERROR(cudaMemcpy(_weights, c_weights, tot_weights * sizeof(double), cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(_wbias, c_wbias, tot_wbias * sizeof(double), cudaMemcpyDeviceToHost));

	long index1 = 0, index2 = 0;
	std::cout.precision(std::numeric_limits<float>::digits10);
    for (long i=0; i<c_layercount; i++) {
		std::cout << "Layer " << i << ": " << std::endl;
        for (long j=0; j<c_layers[i]->neuroncount; j++) {
            for (long k=0; k<c_layers[i]->inputcount; k++) {
				std::cout << _weights[index1++] << "\t";
            }
            std::cout << "wbias: " << _wbias[index2++] << std::endl;
        }
        std::cout << std::endl;
    }
    return;
}

void cuchemnet::reset(void) {
	// HANDLE_ERROR(cudaMemset(c_deltas, 0, tot_weights * sizeof(double))); // not necessary for now
	HANDLE_ERROR(cudaMemset(c_adam_weightsm, 0, tot_size * sizeof(double))); // TODO put them together weightsm and weightsv
	HANDLE_ERROR(cudaMemset(c_adam_weightsv, 0, tot_size * sizeof(double)));
}

/************* Member Functions of EnergyNet ****************/
EnergyNet::EnergyNet(): net_count(0), e_nets(0), hc_outputs(0), dc_outputs(0), batch_size(0), e_outputs(0), d_meta(0) {}

EnergyNet::~EnergyNet() {
	// free the space on Host
	for (int i=0; i<net_count; i++) {
		delete e_nets[i];
	}
	delete [] e_nets;
	free(hc_outputs);
	// free the space on Device
	HANDLE_ERROR(cudaFree(dc_outputs));
	HANDLE_ERROR(cudaFree(e_outputs));
	HANDLE_ERROR(cudaFree(d_meta));
	// release hardware resources used by the CUBLAS library
	HANDLE_CUBLAS(cublasDestroy(handle));
}

void EnergyNet::create(int atomic_num, long *_inputcounts, long *_outputcounts, long *_hiddenlayercounts, long **_hiddenlayers) {
	net_count = atomic_num;
	e_nets    = new cuchemnet*[net_count];
	adam	  = Adam_Setting(); // use the default settings
	// build nets for individual atomic number
	for (int i=0; i<net_count; i++) {
		e_nets[i] = new cuchemnet;
		e_nets[i]->create(_inputcounts[i], _outputcounts[i], _hiddenlayercounts[i], _hiddenlayers[i]);
	}
	// allocate space
	int _size = net_count * sizeof(double *);
	hc_outputs = (double **)malloc(_size);
	HANDLE_ERROR(cudaMalloc((void ***)&dc_outputs, _size));
	// initialization
	memset(hc_outputs, 0, _size);
	HANDLE_ERROR(cudaMemset(dc_outputs, 0, _size));
	// initializes the CUBLAS library and creates the handle (allocates hardware resources on the host and device)
	HANDLE_CUBLAS(cublasCreate(&handle));
}

// propagate the inputs through different nets and sum the energy of individual molecule
void EnergyNet::propagate(double **d_inputs, long _m, long *_meta, PropagateMode _mode, long _gm) {
	if (_m == 0) {
		std::cout << "EnergyNet: _m == 0, LINE: " << __LINE__ << std::endl;
		return;
	}
	// check whether m == batch_size? otherwise cudaFree and Malloc a new space for outputs
	if (_m != batch_size) {
		// std::cout << "EnergyNet: _m is different from batch_size: _m = " << _m << ", batch_size = " << batch_size << std::endl;
		batch_size = _m;
		HANDLE_ERROR(cudaFree(e_outputs));
		HANDLE_ERROR(cudaFree(d_meta));
		// allocate new space
		HANDLE_ERROR(cudaMalloc((void **)&e_outputs, _m * sizeof(double)));
		HANDLE_ERROR(cudaMalloc((void **)&d_meta, _m * net_count * sizeof(long)));
	}
	// copy meta from Host to Device
	HANDLE_ERROR(cudaMemcpy(d_meta, _meta, _m * net_count * sizeof(long), cudaMemcpyHostToDevice));
	// propagate the inputs
	long *_sizes = _meta + (_m - 1) * net_count;
	for (int i=0; i<net_count; i++) {
		if (_sizes[i] == 0) {
			// std::cout << "_sizes[" << i << "] == 0" << std::endl;
			hc_outputs[i] = 0;
			continue;
		}
		e_nets[i]->propagate(d_inputs[i], _sizes[i], handle, _mode, i, _meta, _gm);
		hc_outputs[i] = e_nets[i]->getOutputs();
	}
	// copy outputs pointers from host to device
	HANDLE_ERROR(cudaMemcpy(dc_outputs, hc_outputs, net_count * sizeof(double *), cudaMemcpyHostToDevice));
	// update the e_outputs, which are the energies obtained after the network TODO, could make threadnum much smaller
	int threadnum = MIN(THREADS, _m);
	int blocknum  = (_m + threadnum - 1) / threadnum;
	// set e_outputs to 0's before compute the predicted energies
	HANDLE_ERROR(cudaMemset(e_outputs, 0, _m * sizeof (double)));
	sumEnergy<<<blocknum, threadnum>>>(dc_outputs, e_outputs, _m, d_meta);
}

// train the energy net using Adam algorithm, _meta need to be processed in trainer TODO
double EnergyNet::trainAdam(double **d_inputs, double *d_targets, long _m, long *_meta, double lambda) {
	// update adam.t and compute alphat
	adam.t++;
	double alphat = adam.alpha * sqrt(1 - pow(adam.beta2, adam.t)) / (1 - pow(adam.beta1, adam.t));
	// allocate adam and copy settings to device TODO put it in Device when create the EnergyNet
	//Adam_Setting *d_adam;
	//HANDLE_ERROR(cudaMalloc((void **)&d_adam, sizeof(Adam_Setting)));
	//HANDLE_ERROR(cudaMemcpy(d_adam, &adam, sizeof(Adam_Setting), cudaMemcpyHostToDevice));
	// propagate the network, _meta will be stored in d_meta during propagation
	propagate(d_inputs, _m, _meta, WITH_DERIVATIVE);
	// compute the total error
	double errort = bpGradient(d_inputs, d_targets, _m, _meta, lambda); // e_outputs are modified! or shall we just modify d_targets? TODO
	// update weights and wbias
	long *_sizes = _meta + (_m-1) * net_count;
	for (int i=0; i<net_count; i++) {
		if (_sizes[i] > 0) {
			// e_nets[i]->updateAdam(alphat, d_adam);
			e_nets[i]->updateAdam(alphat);
			if (MAXNORM) {
				// e_nets[i]->updateMaxNormHost();
				e_nets[i]->updateMaxNormDevice();
			}
		}
	}
	// e_nets[0]->updateAdam(alphat, d_adam);
	// printParam();
	// back-propagate errors through different neural nets
	// accumulate the gradients without decay term, GPU parallel, mini-batch
	// add the weight decay term and average the lambda and gradients
	//HANDLE_ERROR(cudaFree(d_adam));
	return errort;
}

// this function implements the max norm regularization method directly on Device, use CUDA dynamic parallelism
// there are perhaps hundreds of neurons in this ANN, so it's possible to allocate a block for a neuron
void cuchemnet::updateMaxNormDevice(void) {
	// launch the kernel
	updateNetMaxNorm<<<1, 1>>>(c_weights, c_wbias, c_structure, c_layercount);
	// cudaFree implies a device synchronization
}

void cuchemnet::updateWeightsMaxNormDevice(void) {
	updateWeightsNetMaxNorm<<<1, 1>>>(c_weights, c_structure, c_layercount);
}

// kernel function (used for launching other kernels) for max norm regularization, only occupies one block and one thread
__global__ void updateNetMaxNorm(double *weights, double *wbias, int *structure, int layercount) {
	int offset1 = 0;
	int offset2 = 0;
	for (int i=0; i<layercount; i++) {
		int inputcount  = structure[i];
		int neuroncount = structure[i+1]; // assume neuroncount is smaller than the largest allowed thread counts (1024 for current GPU)
		updateLayerMaxNorm<<<1, neuroncount>>>(weights + offset1, wbias + offset2, inputcount);
		offset1 += inputcount * neuroncount;
		offset2 += neuroncount;
	}
	cudaDeviceSynchronize();
}

__global__ void updateWeightsNetMaxNorm(double *weights, int *structure, int layercount) {
	int offset = 0;
	for (int i=0; i<layercount; i++) {
		int inputcount  = structure[i];
		int neuroncount = structure[i+1];
		updateWeightsLayerMaxNorm<<<1, neuroncount>>>(weights + offset, inputcount);
		offset += inputcount * neuroncount;
	}
	cudaDeviceSynchronize();
}

__global__ void updateLayerMaxNorm(double *weights, double *wbias, int dim) {
	double sum = 0;
	int offset = threadIdx.x * dim;
	for (int i=0; i<dim; i++) {
		sum += weights[offset + i] * weights[offset + i];
	}
	sum += wbias[threadIdx.x] * wbias[threadIdx.x];
	sum  = sqrt(sum);
	if (sum > MAXLEN) {
		double ratio = MAXLEN / sum;
		for (int i=0; i<dim; i++) {
			weights[offset + i] *= ratio;
		}
		wbias[threadIdx.x] *= ratio;
	}
}

__global__ void updateWeightsLayerMaxNorm(double *weights, int dim) {
	double sum = 0;
	int offset = threadIdx.x * dim;
	for (int i=0; i<dim; i++)
		sum += weights[offset+i] * weights[offset+i];
	sum = sqrt(sum);
	if (sum > MAXLEN) {
		double ratio = MAXLEN / sum;
		for (int i=0; i<dim; i++) {
			weights[offset + i] *= ratio;
		}
	}
}

/*// _weights are incoming weights, _size the number of them, bias is the corresponding bias
__global__ void updateNodeMaxNorm(double *_weights, int _size, double bias) {
	extern __shared__ double cache[]; // dynamic shared memory allocation
	// copy and square the weights into cache
	unsigned int index = threadIdx.x;
	while (index < _size) {
		cache[index] = _weights[index] * _weights[index];
		index       += blockDim.x;
	}
	__syncthreads();
	// sum the results
	index = threadIdx.x;
	for (unsigned int stride = blockDim.x/2; stride >= 1; stride >>= 1) {
		_syncthreads();
		if (index < 
	double sum = 0;
	for (int i=0; i<_size; i++) {
		sum += _weights[i] * _weights[i];
	}
	sum += (*bias) * (*bias);
	if (sum > max_len) {
		double ratio = max_len / sum;
		for (int i=0; i<_size; i++) {
			_weights[i] *= ratio;
		}
		*bias *= ratio;
	}
}*/

void cuchemnet::updateMaxNormHost(void) {
	double max_len = 3.0; // TODO better make it a parameter that can be tuned
	int index1 = 0, index2 = 0;
	// copy the weights from Device to Host
	double *h_weights = (double *)malloc(tot_size * sizeof(double));
	HANDLE_ERROR(cudaMemcpy(h_weights, c_weights, tot_size * sizeof(double), cudaMemcpyDeviceToHost));
	double *h_wbias = h_weights + tot_weights;
	if (DEBUG) {
		cout << "Before updateMaxNorm" << endl;
		printParam();
	}
	// perform the max norm updates
	for (int i=0; i<c_layercount; i++) {
		for (int j=0; j<c_layers[i]->neuroncount; j++) {
			if (DEBUG) {
				std::cout << "Index1: " << index1 << ", Index2: " << index2 << std::endl;
				std::cout << "inputcount: " << c_layers[i]->inputcount << std::endl;
				printf("h_weights: %p, h_weights[0]: %lf, h_wibas[0]: ", h_weights, h_weights[0]); std::cout << std::endl;
			}
			double norm = computeL2norm(h_weights + index1, c_layers[i]->inputcount, h_wbias[index2]);
			if (DEBUG) {
				std::cout << "Layer: " << i << ", Neuron: " << j << ", Norm: " << norm << std::endl;
			}
			if (norm > max_len) {
				if (DEBUG) {
					std::cout << "Norm exceeds the max_len" << std::endl;
				}
				double ratio = max_len / norm;
				applyRatio(h_weights + index1, c_layers[i]->inputcount, ratio);
				h_wbias[index2] *= ratio;
			}
			index1 += c_layers[i]->inputcount;
			index2++;
		}
	}
	// check the results from Host and Device max norm implementations
	if (DEBUG) {
		updateMaxNormDevice();
		double *d_tmp;
		HANDLE_ERROR(cudaMalloc((void **)&d_tmp, tot_size * sizeof(double)));
		HANDLE_ERROR(cudaMemcpy(d_tmp, h_weights, tot_size * sizeof(double), cudaMemcpyHostToDevice));
		cout << "relative error of Host and Device max norm regularization results: " << relativeError(d_tmp, c_weights, tot_size, true) << endl;
		HANDLE_ERROR(cudaFree(d_tmp));
	}
	HANDLE_ERROR(cudaMemcpy(c_weights, h_weights, tot_size * sizeof(double), cudaMemcpyHostToDevice));
	free(h_weights);
	if (DEBUG) {
		cout << "After updateMaxNorm" << endl;
		printParam();
	}
}

void cuchemnet::updateWeightsMaxNormHost(void) {
	int index1 = 0;
	// copy the weights from Device to Host
	double *h_weights = (double *)malloc(tot_weights * sizeof(double));
	HANDLE_ERROR(cudaMemcpy(h_weights, c_weights, tot_weights * sizeof(double), cudaMemcpyDeviceToHost));
	if (DEBUG) {
		cout << "Before updateWeightsMaxNorm" << endl;
		printParam();
	}
	// perform the max norm updates
	for (int i=0; i<c_layercount; i++) {
		for (int j=0; j<c_layers[i]->neuroncount; j++) {
			if (DEBUG) {
				std::cout << "Index1: " << index1 << std::endl;
				std::cout << "inputcount: " << c_layers[i]->inputcount << std::endl;
				printf("h_weights: %p, h_weights[0]: %lf\n", h_weights, h_weights[0]);
			}
			double norm = computeL2norm(h_weights + index1, c_layers[i]->inputcount);
			if (DEBUG) {
				std::cout << "Layer: " << i << ", Neuron: " << j << ", Norm: " << norm << std::endl;
			}
			if (norm > MAXLEN) {
				if (DEBUG) {
					std::cout << "Norm exceeds the MAXLEN" << std::endl;
				}
				double ratio = MAXLEN / norm;
				applyRatio(h_weights + index1, c_layers[i]->inputcount, ratio);
			}
			index1 += c_layers[i]->inputcount;
		}
	}
	// check the results from Host and Device max norm implementations
	if (DEBUG) {
		updateWeightsMaxNormDevice();
		double *d_tmp;
		HANDLE_ERROR(cudaMalloc((void **)&d_tmp, tot_weights * sizeof(double)));
		HANDLE_ERROR(cudaMemcpy(d_tmp, h_weights, tot_weights * sizeof(double), cudaMemcpyHostToDevice));
		cout << "relative error of Host and Device max norm regularization results: " << relativeError(d_tmp, c_weights, tot_weights, true) << endl;
		HANDLE_ERROR(cudaFree(d_tmp));
	}
	HANDLE_ERROR(cudaMemcpy(c_weights, h_weights, tot_weights * sizeof(double), cudaMemcpyHostToDevice));
	free(h_weights);
	if (DEBUG) {
		cout << "After updateWeightsMaxNorm" << endl;
		printParam();
	}
}

// bool _device tells whether v1 and v2 reside on device or host
double relativeError(double *vv1, double *vv2, long _m, bool _device) {
	cublasHandle_t handle;
	HANDLE_CUBLAS(cublasCreate(&handle));
	double *v1, *v2;
	if (!_device) {
		HANDLE_ERROR(cudaMalloc((void **)&v1, sizeof(double) * _m));
		HANDLE_ERROR(cudaMalloc((void **)&v2, sizeof(double) * _m));
		HANDLE_ERROR(cudaMemcpy(v1, vv1, sizeof(double) * _m, cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(v2, vv2, sizeof(double) * _m, cudaMemcpyHostToDevice));
	}
	else {
		v1 = vv1;
		v2 = vv2;
	}
	const double _palpha = 1.0, _nalpha = -1.0;
	double *diff, *sum;
	HANDLE_ERROR(cudaMalloc((void **)&diff, sizeof(double) * _m));
	HANDLE_ERROR(cudaMalloc((void **)&sum, sizeof(double) * _m));
	// store v2 - v1 in diff, and v2 + v1 in sum
	HANDLE_CUBLAS(cublasDcopy(handle, _m, v2, 1, diff, 1));
	HANDLE_CUBLAS(cublasDcopy(handle, _m, v2, 1, sum, 1));
	HANDLE_CUBLAS(cublasDaxpy(handle, _m, &_nalpha, v1, 1, diff, 1));
	HANDLE_CUBLAS(cublasDaxpy(handle, _m, &_palpha, v1, 1, sum, 1));
	// compute the L2 norms
	double norm1, norm2;
	HANDLE_CUBLAS(cublasDdot(handle, _m, diff, 1, diff, 1, &norm1));
	HANDLE_CUBLAS(cublasDdot(handle, _m, sum, 1, sum, 1, &norm2));
	// free the space
	if (!_device) {
		HANDLE_ERROR(cudaFree(v1));
		HANDLE_ERROR(cudaFree(v2));
	}
	HANDLE_ERROR(cudaFree(diff));
	HANDLE_ERROR(cudaFree(sum));

	HANDLE_CUBLAS(cublasDestroy(handle));
	return sqrt(norm1) / sqrt(norm2);
}

double cuchemnet::computeL2norm(double *a, int _size, double b) {
	double norm = 0;
	for (int i=0; i<_size; i++) {
		norm += a[i] * a[i];
	}
	norm += b * b;
	return sqrt(norm);
}

double cuchemnet::computeL2norm(double *a, int _size) {
	double norm = 0;
	for (int i=0; i<_size; i++) {
		norm += a[i] * a[i];
	}
	return sqrt(norm);
}

void cuchemnet::applyRatio(double *a, int _size, double ratio) {
	for (int i=0; i<_size; i++)
		a[i] *= ratio;
}

// function that evaluates the network error, default value for _outputs is NULL/0
double EnergyNet::evaluate(double **d_inputs, double *d_targets, long _m, long *_meta, double *_outputs) {
	// assume data is not too large to fit in memory
	propagate(d_inputs, _m, _meta, WITHOUT_DERIVATIVE);
	// copy the result back to host
	if (_outputs) {
		HANDLE_ERROR(cudaMemcpy(_outputs, e_outputs, _m * sizeof(double), cudaMemcpyDeviceToHost));
	}
	// compute cost -> RMSE TODO?
	return rootMSE(d_targets, _m);
}

double EnergyNet::evaluateRelativeError(double **d_inputs, double *d_targets, long _m, long *_meta) {
	// assume data is not too large to fit in memory
	propagate(d_inputs, _m, _meta, WITHOUT_DERIVATIVE);
	// compute cost
	return relativeErrorGPU(d_targets, e_outputs, _m, true);
}

// bool _device tells whether v1 and v2 reside on device or host
double EnergyNet::relativeErrorGPU(double *vv1, double *vv2, long _m, bool _device) {
	double *v1, *v2;
	if (!_device) {
		HANDLE_ERROR(cudaMalloc((void **)&v1, sizeof(double) * _m));
		HANDLE_ERROR(cudaMalloc((void **)&v2, sizeof(double) * _m));
		HANDLE_ERROR(cudaMemcpy(v1, vv1, sizeof(double) * _m, cudaMemcpyHostToDevice));
		HANDLE_ERROR(cudaMemcpy(v2, vv2, sizeof(double) * _m, cudaMemcpyHostToDevice));
	}
	else {
		v1 = vv1;
		v2 = vv2;
	}
	const double _palpha = 1.0, _nalpha = -1.0;
	double *diff, *sum;
	HANDLE_ERROR(cudaMalloc((void **)&diff, sizeof(double) * _m));
	HANDLE_ERROR(cudaMalloc((void **)&sum, sizeof(double) * _m));
	// store v2 - v1 in diff, and v2 + v1 in sum
	HANDLE_CUBLAS(cublasDcopy(handle, _m, v2, 1, diff, 1));
	HANDLE_CUBLAS(cublasDcopy(handle, _m, v2, 1, sum, 1));
	HANDLE_CUBLAS(cublasDaxpy(handle, _m, &_nalpha, v1, 1, diff, 1));
	HANDLE_CUBLAS(cublasDaxpy(handle, _m, &_palpha, v1, 1, sum, 1));
	// compute the L2 norms
	double norm1, norm2;
	HANDLE_CUBLAS(cublasDdot(handle, _m, diff, 1, diff, 1, &norm1));
	HANDLE_CUBLAS(cublasDdot(handle, _m, sum, 1, sum, 1, &norm2));
	// free the space
	if (!_device) {
		HANDLE_ERROR(cudaFree(v1));
		HANDLE_ERROR(cudaFree(v2));
	}
	HANDLE_ERROR(cudaFree(diff));
	HANDLE_ERROR(cudaFree(sum));

	return sqrt(norm1) / sqrt(norm2);
}

// compute the root mean square error, don't call it twice or call it after expCost, since e_outputs is modified TODO
double EnergyNet::rootMSE(double *d_targets, long _m) {
	double rmse = 0;
	const double _alpha = -1.0;
	HANDLE_CUBLAS(cublasDaxpy(handle, _m, &_alpha, d_targets, 1, e_outputs, 1));
	HANDLE_CUBLAS(cublasDdot(handle, _m, e_outputs, 1, e_outputs, 1, &rmse));
	return sqrt(rmse / _m);
}

// mean square cost
double EnergyNet::squareCost(double *d_targets, long _m) {
	double errorsquare = 0;
	const double _alpha = -1.0;
	HANDLE_CUBLAS(cublasDaxpy(handle, _m, &_alpha, d_targets, 1, e_outputs, 1));
	HANDLE_CUBLAS(cublasDdot(handle, _m, e_outputs, 1, e_outputs, 1, &errorsquare));
	return errorsquare / (2 * _m);
}

/*double EnergyNet::equilibriumCost(double *d_targets, double *d_scales, long _m) {
	const double _alpha = -1.0;
	HANDLE_CUBLAS(cublasDaxpy(handle, _m, &_alpha, d_targets, 1, e_outputs, 1));
}

__global__ void equCost(double *d_outputs, double *d_targets, double *d_gradients, long _m) {
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	while (index < _m) {
		d_outputs = (d_outputs[i] - d_targets[i]) * (d_outputs[i] - d_targets[i])*/

// compute the exponential cost function, assuming the inputs have propagated through the network
double EnergyNet::expCost(double *d_targets, long _m) {
	double errorexp = 0;
	// predicted energies are stored in e_outputs
	const double _alpha = -1.0;
	// const double tau    = 0.5; // set tau be 0.5 TODO
	// const double tau = 1000;
	HANDLE_CUBLAS(cublasDaxpy(handle, _m, &_alpha, d_targets, 1, e_outputs, 1));
	if (DEBUG) {
		printhelper(d_targets+1, __LINE__);
		printhelper(e_outputs+1, __LINE__);
	}
	HANDLE_CUBLAS(cublasDdot(handle, _m, e_outputs, 1, e_outputs, 1, &errorexp));
	if (DEBUG) {
		std::cout << "LINE: " << __LINE__ << ", errorexp: " << errorexp << std::endl;
	}
	// errorexp = tau * exp(1.0 / tau * errorexp);
	// errorexp = tau * exp(errorexp / (tau * _m)); // TODO added _m, is it correct?
	errorexp = TAU * exp(errorexp / (TAU * _m)); // TODO added _m, is it correct?
	if (DEBUG) {
		std::cout << "LINE: " << __LINE__ << ", errorexp: " << errorexp << std::endl;
	}
	
	return errorexp;
}

// a new cost function, which also takes the gradients info, so that the structures near the equilibrium have more important influence
// the gradients info are regarded as the root-mean-square TODO
/*double EnergyNet::equilibriumCost(double *d_targets, double *d_gradients, long _m) {
	
}*/

// compute the weight decay term for the cost function
double EnergyNet::weightDecay(double lambda, long *_sizes) {
	double weight_decay_term = 0;
	// std::cout << "LINE: " << __LINE__ << ", weight_decay_term: " << weight_decay_term << std::endl;
	if (lambda > 0) {
		for (int i=0; i<net_count; i++) {
			double tmp = 0;
			HANDLE_CUBLAS(cublasDdot(handle, e_nets[i]->tot_weights, e_nets[i]->c_weights, 1, e_nets[i]->c_weights, 1, &tmp));
			weight_decay_term += tmp * lambda / (2 * _sizes[i]);
			//std::cout << "LINE: " << __LINE__ << ", weight_decay_term: " << weight_decay_term << std::endl;
		}
	}
	// std::cout << "LINE: " << __LINE__ << ", weight_decay_term: " << weight_decay_term << std::endl;
	return weight_decay_term;
}

// compute and store the output layer errors for individual net, compute the total error, assuming the inputs have propagated through the network
double EnergyNet::outputError(double *d_targets, long _m, double **errorvalues) {
	// double tau = 0.5; // TODO
	// double tau = 1000;
	double errorexp, _alpha;
	if (SQUARECOST) {
		errorexp = squareCost(d_targets, _m); // this will modify e_outputs
		_alpha   = 1.0 / _m;
	}
	else if (EXPCOST) {
		errorexp = expCost(d_targets, _m); // this will modifty e_outputs TODO
		//_alpha   = errorexp * 2 / TAU;
		_alpha   = errorexp * 2 / (_m * TAU);
	}
	// double _alpha = errorexp * 2 / tau;
	// double _alpha = 1.0 / _m; // TODO now use square error to test
	HANDLE_CUBLAS(cublasDscal(handle, _m, &_alpha, e_outputs, 1));
	//std::cout << "PRINTING E_OUTPUTS: " << std::endl;
	//printarray(e_outputs, _m, __LINE__);
	// prepare d_errorvalues
	double **d_errorvalues;
	HANDLE_ERROR(cudaMalloc((void ***)&d_errorvalues, ATOMTYPE_NUM * sizeof(double *)));
	HANDLE_ERROR(cudaMemcpy(d_errorvalues, errorvalues, ATOMTYPE_NUM * sizeof(double *), cudaMemcpyHostToDevice));
	// launch the kernel to compute and store the output layer for each net
	int threadnum = MIN(THREADS, 1); // TODO
	int blocknum  = (_m + threadnum - 1) / threadnum;
	if (DEBUG) {
		std::cout << "E_OUTPUTS: " << std::endl;
		printarray(e_outputs, _m, __LINE__);
		std::cout << "D_META: " << std::endl;
		printarray(d_meta, _m * ATOMTYPE_NUM, __LINE__);
	}
	cudaOutputError<<<blocknum, threadnum>>>(e_outputs, d_meta, d_errorvalues);

	HANDLE_ERROR(cudaFree(d_errorvalues));	
	// TODO return errorexp;
	return sqrt(2 * errorexp);
}

// _m indicates the number of molecules
double EnergyNet::bpGradient(double **d_inputs, double *d_targets, long _m, long *_meta, double lambda) {
	// allocate space for the output layer errors for individual net, the Device space are freed in cuchemnet::helpbpGradient
	double **errorvalues = (double **)malloc(sizeof(double *) * net_count);
	long *_sizes         = _meta + (_m-1) * net_count;
	for (int i=0; i<net_count; i++) {
		HANDLE_ERROR(cudaMalloc((void **)(errorvalues+i), _sizes[i] * sizeof(double)));
	}
	// compute the total error without decay term and store the output layer errors for each individual net
	if (DEBUG) {
		std::cout << "_m: " << _m << " LINE: " << __LINE__ << std::endl;
		for (int i=0; i<net_count; i++)
			std::cout<< i << ": " << _sizes[i] << std::endl;
	}
	double errort = outputError(d_targets, _m, errorvalues);
	// compute the total error with weight decay term
	errort += weightDecay(lambda, _sizes);
	// back propagation on individual net
	for (int i=0; i<net_count; i++) {
		e_nets[i]->bpGradient(d_inputs[i], errorvalues[i], _sizes[i], lambda, handle);
	}

	free(errorvalues);
	return errort;
}

double EnergyNet::helpNumericalGradient(double *w, double **d_inputs, double *d_targets, long _m, long *_meta, double lambda) {
	double EPSILON = 0.0001;
	double *host_w = (double *)malloc(sizeof(double));
	HANDLE_ERROR(cudaMemcpy(host_w, w, sizeof(double), cudaMemcpyDeviceToHost));
	long *_sizes = _meta + (_m-1) * net_count;
	// compute output with w + epsilon
	*host_w += EPSILON;
	HANDLE_ERROR(cudaMemcpy(w, host_w, sizeof(double), cudaMemcpyHostToDevice));
	// printParam();
	propagate(d_inputs, _m, _meta, WITHOUT_DERIVATIVE);
	double j1;
	if (SQUARECOST) {
		j1 = squareCost(d_targets, _m) + weightDecay(lambda, _sizes);
	}
	else if (EXPCOST) {
		j1 = expCost(d_targets, _m) + weightDecay(lambda, _sizes);
	}
	// compute output with w - epsilon
	*host_w -= 2*EPSILON;
	HANDLE_ERROR(cudaMemcpy(w, host_w, sizeof(double), cudaMemcpyHostToDevice));
	// printParam();
	propagate(d_inputs, _m, _meta, WITHOUT_DERIVATIVE);
	double j2;
	if (SQUARECOST) {
		j2 = squareCost(d_targets, _m) + weightDecay(lambda, _sizes);
	}
	else if (EXPCOST) {
		j2 = expCost(d_targets, _m) + weightDecay(lambda, _sizes);
	}
	// convert w back to original value
	*host_w += EPSILON;
	HANDLE_ERROR(cudaMemcpy(w, host_w, sizeof(double), cudaMemcpyHostToDevice));
	free(host_w);

	std::cout << "j1: " << j1 << ", j2: " << j2 << std::endl;
	return (j1-j2)/(2*EPSILON);
}

// assume ans has been allocated and will be freed by the caller later
void EnergyNet::computeNumericalGradient(double *ans, double **d_inputs, double *d_targets, long _m, long *_meta, double lambda) {
	// printParam();
	long index = 0;
	for (int i=0; i<net_count; i++) {
		long _size = e_nets[i]->tot_size;
		for (long j=0; j<_size; j++) {
			ans[index++] = helpNumericalGradient(e_nets[i]->c_weights + j, d_inputs, d_targets, _m, _meta, lambda);
		}
	}
}

void EnergyNet::checkGradient(double **d_inputs, double *d_targets, long _m, long *_meta, double lambda) {
	// get the size of weights and wbias
	long _size = 0;
	for (int i=0; i<net_count; i++) {
		_size += e_nets[i]->tot_size;
	}
	// allocate space
	double *gradient_numerical = (double *)malloc(_size * sizeof(double));
	double *gradient_bp		   = (double *)malloc(_size * sizeof(double));
	// get the gradients numerically
	computeNumericalGradient(gradient_numerical, d_inputs, d_targets, _m, _meta, lambda);
	if (DEBUG) {
		cout << "computed numerical gradients" << endl;
	}
	// propagate the network, _meta will be stored in d_meta during propagation
	// propagate(d_inputs, _m, _meta, WITH_DERIVATIVE);
	propagate(d_inputs, _m, _meta, WITH_DERIVATIVE);
	if (DEBUG) {
		cout << "propagated the inputs" << endl;
	}
	// get the gradients from back-propagation
	bpGradient(d_inputs, d_targets, _m, _meta, lambda);
	if (DEBUG) {
		cout << "computed the back-propagation gradients" << endl;
	}
	long offset = 0;
	for (int i=0; i<net_count; i++) {
		HANDLE_ERROR(cudaMemcpy(gradient_bp + offset, e_nets[i]->c_weights_grad, e_nets[i]->tot_size * sizeof(double), cudaMemcpyDeviceToHost));
		offset += e_nets[i]->tot_size;
	}
	// print the weights
	for (long i=0; i<_size; i++)
		std::cout << gradient_numerical[i] << ", ";
	std::cout << std::endl;
	for (long i=0; i<_size; i++)
		std::cout << gradient_bp[i] << ", ";
	std::cout << std::endl;
	// compute the norm1 / norm2
	double norm1 = 0, norm2 = 0;
	for (long i=0; i<_size; i++) {
		norm1 += (gradient_bp[i] - gradient_numerical[i]) * (gradient_bp[i] - gradient_numerical[i]);
		norm2 += (gradient_bp[i] + gradient_numerical[i]) * (gradient_bp[i] + gradient_numerical[i]);
	}
	std::cout << "norm: " << sqrt(norm1) / sqrt(norm2) << std::endl;
	// free the space
	free(gradient_numerical);
	free(gradient_bp);
}

void EnergyNet::printParam(void) {
	std::cout << "PRINTING WEIGHTS OF ENERGY NET: " << std::endl;
	for (int i=0; i<net_count; i++) {
		std::cout << "NET " << i << ":" << std::endl;
		e_nets[i]->printParam();
	}
}

long EnergyNet::getParamSize(void) {
	long ans = 0;
	for (int i=0; i<net_count; i++)
		ans += e_nets[i]->getParamSize();
	return ans;
}

void EnergyNet::storeParam(double *_weights) {
	long offset = 0;
	for (int i=0; i<net_count; i++) {
		e_nets[i]->storeParam(_weights + offset);
		offset += e_nets[i]->getParamSize();
	}
}

// assuming the _weights are valid and the number of it is correct
void EnergyNet::loadParam(double *_weights) {
	long offset = 0;
	for (int i=0; i<net_count; i++) {
		e_nets[i]->loadParam(_weights + offset, e_nets[i]->getParamSize());
		offset += e_nets[i]->getParamSize();
	}
}

void EnergyNet::loadParam(const string &fname, bool Binary) {
	cout << "Read weights from " << fname << endl;
	unsigned long _size = 0;
	for (int i=0; i<net_count; i++) {
		_size += e_nets[i]->tot_size;
	}
	double *_weights = (double *)malloc(_size * sizeof(double));
	ifstream is;
	if (Binary) {
		is.open(fname, ios::binary);
		assert(is.is_open());
		is.seekg(0, ios::end);
		streampos len = is.tellg();
		assert(!(len % sizeof(double)) && (len / sizeof(double) == _size));
		is.seekg(0, ios::beg);
		is.read((char *)_weights, len);
	}
	else {
		is.open(fname);
		assert(is.is_open());
		int i = 0;
		double val;
		while (is >> val) { // if the file contains more than _size weights, it'll generate large i?
			assert( i < _size);
			_weights[i++] = val;
		}
		assert(i == _size);
	}
	if (DEBUG) {
		cout << "Weights read from " << fname << ":" << endl;
		cout << "Number of weights: " << _size << endl;
		for (int i=0; i<_size; i++)
			cout << _weights[i] << "\t";
		cout << endl;
	}
	loadParam(_weights);
	free(_weights);
}

/******************************* Gradient Propagation *********************************/
// compute the gradients analytically given the input data, assume d_grads and d_target_grads also point to device memory
double EnergyNet::evaluateGradient(double **d_inputs, long _m, long *_meta, double **d_grads, double *d_target_grads, long _gm, long *_gmeta) {
	// propagate the input and store the derivatives
	propagate(d_inputs, _m, _meta, GRAD_DERIVATIVE, _gm);
	// propagate the gradients, compute the analytical gradients
	propagate(d_grads, _gm, _gmeta, GRAD_PROPAGATION);
	// extract the output gradients, the results are stored in e_outputs
	return relativeErrorGPU(d_target_grads, e_outputs, _gm, true);
}

double EnergyNet::evaluateGradient(double **d_inputs, long _m, long *_meta, double **d_grads, double *d_target_grads, double *ret) {
	// propagate the input and store the derivatives
	propagate(d_inputs, _m, _meta, WITH_DERIVATIVE);
	// propagate the gradients, compute the analytical gradients
	propagate(d_grads, _m, _meta, GRAD_PROPAGATION);
	if (ret) {
		HANDLE_ERROR(cudaMemcpy(ret, e_outputs, _m * sizeof(double), cudaMemcpyDeviceToHost));
	}
	// extract the output gradients, the results are stored in e_outputs
	return relativeErrorGPU(d_target_grads, e_outputs, _m, true);
}

// d_inputs: every atom has an atomic environment vector
// d_grads: atom (x, y, z) has three gradients, to predict one, we need to use derivatives of all aevs of the molecules
// so naively, the size of d_grads is 3 * n * size(d_inputs), where n is the averaged number of atoms of a molecule
// however, since we want to have high precision near the equlibrium, we only use gradients near the equlibrium to refine the model
// but first we need to make sure refine with gradients does help to improve the accuracy of gradients prediction
// TODO TODO TODO now we assuming d_inputs has been made to match the number of d_grads, inefficient but easy to implement for now
// we need to separate the feed-forward propagation and back propagation TODO
double EnergyNet::Grad_trainAdam(double **d_inputs, long _m, long *_meta, double **d_grads, double *d_target_grads, long _gm, long *_gmeta, double lambda) {
	// first need to propagate the inputs, assume the derivatives of "y=f(z)" are computed and stored in "c_derivatives"
	propagate(d_inputs, _m, _meta, GRAD_DERIVATIVE, _gm); // TODO
	// propagate the gradient
	propagate(d_grads, _gm, _gmeta, GRAD_PROPAGATION);
	// perform the back propagation
	double errort = bpGradient(d_grads, d_target_grads, _gm, _gmeta, lambda);
	// update weights
	Grad_updateAdam(_meta + (_m - 1) * net_count);

	return errort;
}

double EnergyNet::Grad_trainAdam(double **d_inputs, long _m, long *_meta, double **d_grads, double *d_target_grads, double lambda) {
	// propagate the input
	propagate(d_inputs, _m, _meta, WITH_DERIVATIVE);
	// propagate the gradients
	propagate(d_grads, _m, _meta, GRAD_PROPAGATION);
	// perform the back propagation
	double errort = bpGradient(d_grads, d_target_grads, _m, _meta, lambda);
	zeroWbiasGrad();
	// update weights
	Grad_updateAdam(_meta + (_m - 1) * net_count);

	return errort;
}

double EnergyNet::Grad_Energy_trainAdam(double **d_inputs, double *d_targets, long _m, long *_meta, double **d_grads, double *d_target_grads, double lambda) {
	double errort_energy = trainAdam(d_inputs, d_targets, _m, _meta, lambda);
	double errort_grad   = Grad_trainAdam(d_inputs, _m, _meta, d_grads, d_target_grads, lambda);
	return errort_energy;
}

void EnergyNet::Grad_updateAdam(long *_sizes) {
	// update adam.t and alphat
	adam.t++;
	double alphat = adam.alpha * sqrt(1 - pow(adam.beta2, adam.t)) / (1 - pow(adam.beta1, adam.t));
	//Adam_Setting *d_adam;
    //HANDLE_ERROR(cudaMalloc((void **)&d_adam, sizeof(Adam_Setting)));
    //HANDLE_ERROR(cudaMemcpy(d_adam, &adam, sizeof(Adam_Setting), cudaMemcpyHostToDevice));
	// update weights
	for (int i=0; i<net_count; i++) {
		if (_sizes[i] > 0) {
			// only need to update weights
			//e_nets[i]->updateAdam(alphat, d_adam, e_nets[i]->tot_weights);
			e_nets[i]->updateAdam(alphat, e_nets[i]->tot_weights);
			if (MAXNORM) {
				if (DEBUG) {
					e_nets[i]->updateWeightsMaxNormHost();
				}
				// e_nets[i]->updateMaxNorm(); // The updateMaxNorm also updates wbias, is it required property for gradient back-propagation? TODO we have set c_wbias_grad to be zero and updateWeightsMaxNormDevice() only updates weights!
				e_nets[i]->updateWeightsMaxNormDevice();
			}
		}
	}
	// free the space
	//HANDLE_ERROR(cudaFree(d_adam));
}

void EnergyNet::reset(double _alpha) {
	adam.t = 0;
	adam.alpha = _alpha;
	// zero out some spaces
	for (int i=0; i<net_count; i++) {
		e_nets[i]->reset();
	}
}

double EnergyNet::helpNumericalGradient(double *w, double **d_inputs, long _m, long *_meta, double **d_grads, double *d_target_grads, double lambda) {
	double EPSILON = 0.0001;
    double *host_w = (double *)malloc(sizeof(double));
    HANDLE_ERROR(cudaMemcpy(host_w, w, sizeof(double), cudaMemcpyDeviceToHost));
    long *_sizes = _meta + (_m-1) * net_count;
	propagate (d_inputs, _m, _meta, WITH_DERIVATIVE);
    // compute output with w + epsilon
    *host_w += EPSILON;
    HANDLE_ERROR(cudaMemcpy(w, host_w, sizeof(double), cudaMemcpyHostToDevice));
	propagate (d_grads, _m, _meta, GRAD_PROPAGATION);
	double j1;
	if (SQUARECOST) {
		j1 = squareCost(d_target_grads, _m) + weightDecay(lambda, _sizes);
	}
	else if (EXPCOST) {
		j1 = expCost(d_target_grads, _m) + weightDecay(lambda, _sizes);
	}
	// compute output with w - epsilon
    *host_w -= 2*EPSILON;
    HANDLE_ERROR(cudaMemcpy(w, host_w, sizeof(double), cudaMemcpyHostToDevice));
	// propagate (d_inputs, _m, _meta, WITH_DERIVATIVE);
	propagate (d_grads, _m, _meta, GRAD_PROPAGATION);
	double j2;
	if (SQUARECOST) {
		j2 = squareCost(d_target_grads, _m) + weightDecay(lambda, _sizes);
	}
	else if (EXPCOST) {
		j2 = expCost(d_target_grads, _m) + weightDecay(lambda, _sizes);
	}
	// convert w back to original value
    *host_w += EPSILON;
    HANDLE_ERROR(cudaMemcpy(w, host_w, sizeof(double), cudaMemcpyHostToDevice));
    free(host_w);

	std::cout << "j1: " << j1 << ", j2: " << j2 << std::endl;
	return (j1 - j2) / (2 * EPSILON);
}

void EnergyNet::computeNumericalGradient(double *ans, double **d_inputs, long _m, long *_meta, double **d_grads, double *d_target_grads, double lambda) {
	long index = 0;
	for (int i=0; i<net_count; i++) {
		long _size = e_nets[i]->tot_size;
		for (long j=0; j<_size; j++) {
			ans[index++] = helpNumericalGradient(e_nets[i]->c_weights + j, d_inputs, _m, _meta, d_grads, d_target_grads, lambda);
		}
	}
}

void EnergyNet::zeroWbiasGrad(void) {
	for (int i=0; i<net_count; i++)
		HANDLE_ERROR(cudaMemset(e_nets[i]->c_wbias_grad, 0, e_nets[i]->tot_wbias * sizeof(double)));
}

void EnergyNet::Grad_checkGradient(double **d_inputs, long _m, long *_meta, double **d_grads, double *d_target_grads, double lambda) {
	// get the total size of weights and bias
	long _size = 0;
	for (int i=0; i<net_count; i++)
		_size += e_nets[i]->tot_size;
	// allocate space
	double *gradient_numerical = (double *)malloc(_size * sizeof(double));
	double *gradient_bp        = (double *)malloc(_size * sizeof(double));
	// get the gradients numerically
	computeNumericalGradient(gradient_numerical, d_inputs, _m, _meta, d_grads, d_target_grads, lambda);
	// get the gradients from back-propagation
	propagate(d_inputs, _m, _meta, WITH_DERIVATIVE);
	propagate(d_grads, _m, _meta, GRAD_PROPAGATION);
	bpGradient(d_grads, d_target_grads, _m, _meta, lambda);
	zeroWbiasGrad();
	long offset = 0;
	for (int i=0; i<net_count; i++) {
		HANDLE_ERROR(cudaMemcpy(gradient_bp + offset, e_nets[i]->c_weights_grad, e_nets[i]->tot_size * sizeof(double), cudaMemcpyDeviceToHost));
		offset += e_nets[i]->tot_size;
	}
	// print the weights
	std::cout << "numerical gradients: " << std::endl;
    for (long i=0; i<_size; i++)
        std::cout << gradient_numerical[i] << ", ";
    std::cout << std::endl;
	std::cout << "back-propagation gradients: " << std::endl;
    for (long i=0; i<_size; i++)
        std::cout << gradient_bp[i] << ", ";
    std::cout << std::endl;
    // compute the norm1 / norm2
    double norm1 = 0, norm2 = 0;
    for (long i=0; i<_size; i++) {
        norm1 += (gradient_bp[i] - gradient_numerical[i]) * (gradient_bp[i] - gradient_numerical[i]);
        norm2 += (gradient_bp[i] + gradient_numerical[i]) * (gradient_bp[i] + gradient_numerical[i]);
    }
    std::cout << "norm: " << sqrt(norm1) / sqrt(norm2) << std::endl;
    // free the space
    free(gradient_numerical);
    free(gradient_bp);
}
