/*
 *	Header file for back propagation in a fully connected neural net with Adam algorithm.
 */

#ifndef CHEMNET_H
#define CHEMNET_H

#include <vector>
#include <iostream>
#include <cublas_v2.h>

#define DEBUG 0
#define ATOMTYPE_NUM 2

// Propagate mode, which controls how the input propagates
// WITHOUT_DERIVATIVE: just propagate the input, not computing the derivatives
// WITH_DERIVATIVE: propagate the input and store the derivative infos for later back-propagation
// GRAD_DERIVATIVE: propagate the input and store the derivative infos for later gradient-propagation
// GRAD_PROPAGATION: use the derivative infos as the new activation functions
enum PropagateMode { WITHOUT_DERIVATIVE, WITH_DERIVATIVE, GRAD_DERIVATIVE, GRAD_PROPAGATION };
// different types of activation functions
// enum ActivationMode { SIGMOID, GAUSSIAN, TANH, RELU, LRELU, ELU };

/************************** Helper Functions ***********************/
// cuda error check
#define HANDLE_ERROR(err) (HandleError(err, __FILE__, __LINE__))
static void HandleError(cudaError_t err, const char *file, long line) {
    if (err != cudaSuccess) {
        printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
        exit(EXIT_FAILURE);
    }
}

#define HANDLE_CUBLAS(stat) (HandleCublas(stat, __FILE__, __LINE__))
static void HandleCublas(cublasStatus_t stat, const char *file, long line) {
	if (stat != CUBLAS_STATUS_SUCCESS) {
		printf("error in %s at line %d\n", file, line);
		exit(EXIT_FAILURE);
	}
}

// print out a value given the address of it on Device
void printhelper(double *w, long line);
void printarray(double *w, long size, long line);
double relativeError(double *, double *, long, bool);

/*************************** Data structures ***************************/
// adam setting parameters
struct Adam_Setting {
	double alpha, beta1, beta2, epsilon, t;
	// default values for the Adam algorithm
	Adam_Setting(): alpha(0.001), beta1(0.9), beta2(0.999), epsilon(1e-8), t(0) {}
	Adam_Setting(double _a, double _b1, double _b2, double _eps, double _t): 
					alpha(_a), beta1(_b1), beta2(_b2), epsilon(_eps), t(_t) {}
};

// layer structure, it doesn't store the real data
struct culayer {
	long inputcount;	// total number of elements in layerinput, assuming the network is fully connected
	long neuroncount;	// number of neurons
	long weights_size;	// inputcount * neuroncount
	double *weights;
	double *wbias;
	double *deltas;		// deltas of weights for the last update (used for momentum)
	double *outputs;
	double *derivatives;
	double *weights_grad;
	double *wbias_grad;
	// data used for adam algorithm
	double *adam_weightsm;
	double *adam_weightsv;
	double *adam_wbiasm;
	double *adam_wbiasv;
};

// the neural network structure
class cuchemnet {
	private:
		long batch_size;	// the number of inputs in mini-batch
		long tot_wbias;
		long tot_weights;
		long tot_size;
		long c_layercount;  // hiddenlayer count + 1
		// real data, the members in culayer polong to these data
		double *c_weights, *c_wbias, *c_deltas, *c_outputs, *c_weights_grad, *c_wbias_grad;
		double *c_adam_weightsm, *c_adam_weightsv, *c_adam_wbiasm, *c_adam_wbiasv;
		double *c_ones;	// all '1's in double precision

		int *c_structure; // store the structure of the neural net

		// used to store derivatives of f(z), useful with Gaussian activation function and Gradient propagation
		long tot_derivatives; // not useful for now
		double *c_derivatives;

		culayer **c_layers;  // pointers to all layers

		// build a vector with all 1's
		void buildOnes(void);
		void buildStructure(void);
		void helpbpGradient(double *d_input, double *errorvalues, cublasHandle_t &handle);

		friend class EnergyNet;

	public:
		// constructor: all values to 0
		cuchemnet();
		// destructor: free memories
		~cuchemnet();
		// create the network structure
		void create(long _inputcount, long _outputcount, long _hiddenlayercount, long *_hiddenlayers);
		void createLayersWeights(std::vector<long> &_structure);
		// computes the network values given an input pattern
		void propagate(double *d_input, long _m, cublasHandle_t &handle);
		void propagate(double *, long, cublasHandle_t &, PropagateMode, int net, long *_meta = 0, long _gm = 0);
		void Grad_propagate(double *, long, cublasHandle_t &);
		// check sizes and allocate required space
		void prepropagate(long, PropagateMode, long _gm = 0);
		void allocateDerivatives(long);
		void applyActivation(int, long, PropagateMode, int net, long *_meta = 0);
		// compute gradients based on back-propagation
		void bpGradient(double *d_input, double *errorvalues, long _m, double lambda, cublasHandle_t &handle); // TODO it's possible that _m == 0
		// update wegihts
		// void updateAdam(double alphat, Adam_Setting *d_adam);
		// void updateAdam(double alphat, Adam_Setting *d_adam, long _size);
		void updateAdam(double alphat);
		void updateAdam(double alphat, long _size);
		void updateMaxNormHost(void);
		void updateMaxNormDevice(void);
		// only concern the weights for max norm
		void updateWeightsMaxNormDevice(void);
		void updateWeightsMaxNormHost(void);
		// helper for max norm
		double computeL2norm(double *, int, double);
		double computeL2norm(double *, int);
		void applyRatio(double *, int, double);
		// get the point to the output layer outputs
		inline double *getOutputs(void) { return c_layers[c_layercount-1]->outputs; }
		// load weights and wbias
		void loadParam(double *_weights, long _weights_num, double *_wbias, long _wbias_num);
		void loadParam(double *_c_weights, long _c_weights_num);
		// store weights and wbias
		void storeParam(double *_weights);
		// print weights and wbias
		void printParam(void);
		// return the total size of parameters
		inline long getParamSize(void) { return tot_size; }
		// reset some space
		void reset(void);
};

// a class specifically designed for predicting molecular energy
class EnergyNet {
	private:
		int net_count;			// number of nets, == ATOMTYPE_NUM
		cuchemnet **e_nets;		// different nets for individual atomic number
		double **hc_outputs;	// pointers to outputs of different nets, on host
		double **dc_outputs;	// pointers to outputs of different nets, on device
		long batch_size;		// batch size, updated during the propagation process
		double *e_outputs;		// computed energies, on device, updated during the propagation process
		long *d_meta;			// meta data, stored on device, updated during feed-forward propagation
		// cublas handle, which points to an opaque structure holding the cuBLAS library context
		cublasHandle_t handle;
		// adam setting
		Adam_Setting adam;
		void zeroWbiasGrad(void);

	public:
		EnergyNet();
		~EnergyNet();

		// now assuming the neural net structures are the same per each individual atomic number
		void create(int atomic_num, long *_inputcounts, long *_outputcounts, long *_hiddenlayercounts, long **_hiddenlayers);
		// void initWeights(void);
		// void initBias(void);
		void propagate(double **d_inputs, long _m, long *_meta, PropagateMode _mode, long _gm = 0);
		double trainAdam(double **d_inputs, double *d_targets, long _m, long *_meta, double lambda);
		double evaluate(double **d_inputs, double *d_targets, long _m, long *_meta, double *_outputs = 0);
		void checkGradient(double **d_inputs, double *d_targets, long _m, long *_meta, double lambda);
		void Grad_checkGradient(double **d_inputs, long _m, long *_meta, double **d_grads, double *d_target_grads, double lambda);

		double weightDecay(double lambda, long *_sizes);
		double outputError(double *d_targets, long _m, double **errorvalues);
		double bpGradient(double **d_inputs, double *d_targets, long _m, long *_meta, double lambda);
		double rootMSE(double *d_targets, long _m);
		double squareCost(double *d_targets, long _m);
		double expCost(double *d_targets, long _m);
		double helpNumericalGradient(double *w, double **d_inputs, double *d_targets, long _m, long *_meta, double lambda);
		double helpNumericalGradient(double *w, double **d_inputs, long _m, long *_meta, double **d_grads, double *d_target_grads, double lambda);
		void computeNumericalGradient(double *ans, double **d_inputs, double *d_targets, long _m, long *meta, double lambda);
		void computeNumericalGradient(double *ans, double **d_inputs, long _m, long *_meta, double **d_grads, double *d_target_grads, double lambda);
		void printParam(void);
		long getParamSize(void);
		void storeParam(double *_weights);
		void loadParam(double *_weights);
		void loadParam(const std::string &fname, bool Binary = false);
		// evalute the relative error: sqrt(norm1) / sqrt(norm2), where norm1 and norm2 are the L2 norms of the differences and sums between two vectors
		double evaluateRelativeError(double **, double *, long, long *);
		// compute the relative error of two arrays, assuming they reside on Device memory
		double relativeErrorGPU(double *, double *, long, bool);
		// train the model using input gradients
		double Grad_trainAdam(double **d_inputs, long _m, long *_meta, double **d_grads, double *d_target_grads, long _gm, long *_gmeta, double lambda); // TODO
		// a version with one gradient per molecule
		double Grad_trainAdam(double **d_inputs, long _m, long *_meta, double **d_grads, double *d_target_grads, double lambda);
		double Grad_Energy_trainAdam(double **d_inputs, double *d_target, long _m, long *_meta, double **d_grads, double *d_target_grads, double lambda);
		void Grad_updateAdam(long *_sizes);
		// compute the gradients analytically given the input data
		// double *analyticGradient(double **d_inputs, long _m, long *_meta, double **d_grads, long _gm, long *_gmeta); // d_grad is dx/dqi
		double evaluateGradient(double **d_inputs, long _m, long *_meta, double **d_grads, double *d_target_grads, long _gm, long *_gmeta);
		double evaluateGradient(double **d_inputs, long _m, long *_meta, double **d_grads, double *d_target_grads, double *ret = 0);
		// reset some space
		void reset(double _alpha = 0.001);
};

/************************** Kernel Functions ***********************/
// activation kernels
__global__ void activationSigmoid(double *, long);
__global__ void activationSigmoid(double *, double *, long);
__global__ void activationGaussian(double *, long);
__global__ void activationGaussian(double *, double *, long);
__global__ void activationIdentity(double *, long);
__global__ void activationIdentity(double *, double *, long);
__global__ void Grad_activation(double *, double *, long);
__global__ void activationGradSigmoid(double *, double *, long, int, long *, int);
__global__ void activationGradGaussian(double *, double *, long, int, long *, int);
__global__ void activationGradIdentity(double *, double *, long, int, long *, int);

__global__ void activation(double *outputs, long _neuroncount);
__global__ void outputActivation(double *outputs, long _neuroncount);
// __global__ void updateWeightsAdam(double *, double *, double *, double *, Adam_Setting *, double, long);
__global__ void updateWeightsAdam(double *, double *, double *, double *, double, long);
__global__ void propagateError(double *, double *, long);
__global__ void sumEnergy(double **, double *, int, int *);
__global__ void cudaOutputError(double *e_outputs, long *d_meta, double **d_errorvalues);
__global__ void computeSigmoidNodeGradient(double *outputs, long _size);
// for max norm regularization
__global__ void updateNetMaxNorm(double *, double *, int *, int);
__global__ void updateLayerMaxNorm(double *, double *, int);
__global__ void updateWeightsNetMaxNorm(double *, int *, int);
__global__ void updateWeightsLayerMaxNorm(double *, int);

#endif
