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

#include "chemnet.h"
#include "cblas.h"

/************* Helper Functions *************/
inline double activation(double x) { // TODO: a faster implementation f(x) = 0.5 * (x * alpha / (1 + abs(x*alpha)) + 0.5
	// sigmoid function
	return 1./(1. + exp(-x));
}

inline double outputActivation(double x) {	// activation function for the output (regression/classification) layer TODO
	return x;	// a linear activation function
}

/************* Member Functions of layer ****************/
layer::layer(): inputcount(0), neuroncount(0), weights_size(0), weights(0), wbias(0), deltas(0), outputs(0), weights_grad(0), wbias_grad(0), adam_weightsm(0), adam_weightsv(0), adam_wbiasm(0), adam_wbiasv(0) {}

layer::~layer() {
	delete [] weights;
	delete [] wbias;
	delete [] deltas;
	delete [] outputs;
	delete [] weights_grad;
	delete [] wbias_grad;
	delete [] adam_weightsm;
	delete [] adam_weightsv;
	delete [] adam_wbiasm;
	delete [] adam_wbiasv;
}

// create a layer
void layer::create(int _inputcount, int _neuroncount, int fan_out) {
	if (_inputcount < 1 || _neuroncount < 1 || fan_out < 0)
		assert(false);

	inputcount   = _inputcount;
	neuroncount  = _neuroncount;
	weights_size = _neuroncount * _inputcount;
	// allocate space
	weights       = new double[weights_size];
	weights_grad  = new double[weights_size]; // will be initialized to 0 in every bpGradient/training step
	deltas        = new double[weights_size];
	wbias         = new double[_neuroncount];
	wbias_grad    = new double[_neuroncount]; // will be initialized to 0 in every bpGradient/training step
	outputs       = new double[_neuroncount]; // doesn't need initialization
	adam_weightsm = new double[weights_size];
	adam_weightsv = new double[weights_size];
	adam_wbiasm   = new double[_neuroncount];
	adam_wbiasv   = new double[_neuroncount];
	// initialize weights, deltas, and wbias, and data for adam
	memset(deltas, 0, weights_size * sizeof(double));
	memset(wbias, 0, _neuroncount * sizeof(double));
	memset(adam_weightsm, 0, weights_size * sizeof(double));
	memset(adam_weightsv, 0, weights_size * sizeof(double));
	memset(adam_wbiasm, 0, _neuroncount * sizeof(double));
	memset(adam_wbiasv, 0, _neuroncount * sizeof(double));
	// set the uniform distribution, assuming sigmoid function TODO
	std::default_random_engine gen((std::random_device())());
	double r = 4 * sqrt(6) / sqrt(_inputcount + fan_out + 1);
	std::uniform_real_distribution<double> distribution(-r, r);
	for (int i=0; i<weights_size; i++) {
		weights[i] = distribution(gen);
	}
}

// calculate the output vector of a layer
void layer::calculate(double *layerinput, bool isOutput) {
	// copy wbias to outputs
	cblas_dcopy(neuroncount, wbias, 1, outputs, 1);
	// add W * X to outputs
	cblas_dgemv(CblasRowMajor, CblasNoTrans, neuroncount, inputcount, 1.0, weights, inputcount, layerinput, 1, 1, outputs, 1);
	// activation function: TODO how to make it fast
	if (isOutput) {
		for (int i=0; i<neuroncount; i++)
			outputs[i] = outputActivation(outputs[i]);
	}
	else {
		for (int i=0; i<neuroncount; i++)
			outputs[i] = activation(outputs[i]);
	}
}

/************* Member Functions of chemnet ****************/
chemnet::chemnet(): c_layers(0), c_layercount(0) {}

chemnet::~chemnet() {
	for (int i=0; i<c_layercount; i++) {
		delete c_layers[i];
	}
	delete [] c_layers;
}

// create the network structure
void chemnet::create(int _inputcount, int _outputcount, int *_hiddenlayers, int _hiddenlayercount) {
	if (_inputcount < 1 || _outputcount < 1 || _hiddenlayercount < 0)
		assert(false);

	// initialization for the architecture
	c_layercount = _hiddenlayercount + 1;
	c_layers     = new layer*[c_layercount];
	adam         = Adam_Setting();
	// build a vector to store the structure of the NN
	std::vector<int> _structure;
	_structure.push_back(_inputcount);
	for (int i=0; i<_hiddenlayercount; i++)
		_structure.push_back(_hiddenlayers[i]);
	_structure.push_back(_outputcount);
	_structure.push_back(1);	// assuming the output layer has fan_out = 1
	// create the layers
	for (int i=0; i<c_layercount; i++) {
		c_layers[i] = new layer;
		c_layers[i]->create(_structure[i], _structure[i+1], _structure[i+2]);
	}
	// CHECK
	std::cout << "size of structure: " << _structure.size() << std::endl;
}

// initialize the output bias differently: output biass to optimal value as if weights were 0
void chemnet::initOutputWbias(double **desiredoutput, int m) {
	// TODO
}

// load parameters: weights and wbias, assuming the sizes of _weights and wbias are correct, and caller frees the _weights and _wbias
void chemnet::loadParam(double *_weights, int _weights_num, double *_wbias, int _wbias_num) {
	int index1 = 0, index2 = 0;
	for (int i=0; i<c_layercount; i++) {
		for (int j=0; j<c_layers[i]->weights_size; j++) {
			c_layers[i]->weights[j] = _weights[index1++];
		}
		for (int j=0; j<c_layers[i]->neuroncount; j++) {
			c_layers[i]->wbias[j] = _wbias[index2++];
		}
	}
	// CHECK
	assert((index1 == _weights_num) && (index2 == _wbias_num));
}

// calculate and propagate the network
void chemnet::propagate(double *input) {
	// start from the input layer and go to the output layer
	for (int i=0; i<c_layercount-1; i++) {
		c_layers[i]->calculate(input, false);
		input = c_layers[i]->outputs;
	}
	c_layers[c_layercount-1]->calculate(input, true);
}

// train neural net, alpha is learning rate, momentum is learning momentum
double chemnet::train(double **input, double **desiredoutput, int m, double alpha, double momentum, double lambda) {
	// total quadratic error
	double errort = 0;
	// compute gradients for weights and wbias, and total error
	bpGradient(input, desiredoutput, m, lambda, errort);
	// update weights and wbias
	for (int i=0; i<c_layercount; i++) {
		// momentum * deltas
		cblas_dscal(c_layers[i]->weights_size, momentum, c_layers[i]->deltas, 1);
		// delta = alpha * weights_grad + momentum * delta
		cblas_daxpy(c_layers[i]->weights_size, alpha, c_layers[i]->weights_grad, 1, c_layers[i]->deltas, 1);
		// weights - delta
		cblas_daxpy(c_layers[i]->weights_size, -1.0, c_layers[i]->deltas, 1, c_layers[i]->weights, 1);
		// wbias - alpha * wbias_grad
		cblas_daxpy(c_layers[i]->neuroncount, -alpha, c_layers[i]->wbias_grad, 1, c_layers[i]->wbias, 1);
	}
	// return the training error: 1/2 * MSE
	return errort / 2;
}

// train the neural net using Adam algorithm
double chemnet::train_adam(double **input, double **desiredoutput, int m, double lambda) {
	// update adam.t and compute alphat
	adam.t++;
	double alphat = adam.alpha * sqrt(1 - pow(adam.beta2, adam.t)) / (1 - pow(adam.beta1, adam.t));
	// total quadratic error
	double errort = 0;
	// compute gradients for weights and wbias, and total error
	bpGradient(input, desiredoutput, m, lambda, errort);
	// update weights and wbias
	for (int i=0; i<c_layercount; i++) {
		// update weights TODO GPU
		for (int j=0; j<c_layers[i]->weights_size; j++) {
			c_layers[i]->adam_weightsm[j] = adam.beta1 * c_layers[i]->adam_weightsm[j] + (1-adam.beta1) * c_layers[i]->weights_grad[j];
			c_layers[i]->adam_weightsv[j] = adam.beta2 * c_layers[i]->adam_weightsv[j] + (1-adam.beta2) * c_layers[i]->weights_grad[j] * c_layers[i]->weights_grad[j];
			c_layers[i]->weights[j]      -= alphat * c_layers[i]->adam_weightsm[j] / (sqrt(c_layers[i]->adam_weightsv[j]) + adam.epsilon);
		}
		// update wbias
		for (int j=0; j<c_layers[i]->neuroncount; j++) {
			c_layers[i]->adam_wbiasm[j]   = adam.beta1 * c_layers[i]->adam_wbiasm[j] + (1-adam.beta1) * c_layers[i]->wbias_grad[j];
			c_layers[i]->adam_wbiasv[j]   = adam.beta2 * c_layers[i]->adam_wbiasv[j] + (1-adam.beta2) * c_layers[i]->wbias_grad[j] * c_layers[i]->wbias_grad[j];
			c_layers[i]->wbias[j]        -= alphat * c_layers[i]->adam_wbiasm[j] / (sqrt(c_layers[i]->adam_wbiasv[j]) + adam.epsilon);
		}
	}
	// return the training error: 1/2 * MSE
	return errort / 2;
}

// given testing input data, propagate the network and return the outputs
std::vector<double> chemnet::test(double **input, int m) {
	std::vector<double> ans;
	for (int i=0; i<m; i++) {
		propagate(input[i]);
		for (int j=0; j<c_layers[c_layercount-1]->neuroncount; j++)
			ans.push_back(c_layers[c_layercount-1]->outputs[j]);
	}
	return ans;
}

// check the gradients obtained numerically and from back-propagation
void chemnet::checkGradient(double **input, double **desiredoutput, int m, double lambda) {
	// get the gradients from back-propagation
	std::vector<double> gradient_bp;
	double __tmp;
	bpGradient(input, desiredoutput, m, lambda, __tmp);
	for (int i=0; i<c_layercount; i++) {
		gradient_bp.insert(gradient_bp.end(), c_layers[i]->weights_grad, c_layers[i]->weights_grad + c_layers[i]->weights_size);
		gradient_bp.insert(gradient_bp.end(), c_layers[i]->wbias_grad, c_layers[i]->wbias_grad + c_layers[i]->neuroncount);
	}
	// get the gradients numerically
	std::vector<double> gradient_numerical = computeNumericalGradient(input, desiredoutput, m, lambda);
	// print the weights
	for (int i=0; i<gradient_bp.size(); i++) 
		std::cout << gradient_bp[i] << ", ";
	std::cout << std::endl;
	for (int i=0; i<gradient_numerical.size(); i++) 
		std::cout << gradient_numerical[i] << ", ";
	std::cout << std::endl;
	// compute the norm1 / norm2
	double norm1 = 0, norm2 = 0;
	for (int i=0; i<gradient_bp.size(); i++) {
		norm1 += (gradient_bp[i] - gradient_numerical[i]) * (gradient_bp[i] - gradient_numerical[i]);
		norm2 += (gradient_bp[i] + gradient_numerical[i]) * (gradient_bp[i] + gradient_numerical[i]);
	}
	std::cout << "norm: " << sqrt(norm1) / sqrt(norm2) << std::endl;
	std::cout << "size: " << (gradient_bp.size() == gradient_numerical.size() ? gradient_bp.size() : 0) << std::endl;
}

// compute the cost function: square error + regulization
double chemnet::squareCost(double **input, double **desiredoutput, int m, double lambda) {
	double errort = 0;
	for (int i=0; i<m; i++) {
		propagate(input[i]);
		for (int j=0; j<c_layers[c_layercount-1]->neuroncount; j++) {
			double output = c_layers[c_layercount-1]->outputs[j];
			errort += (output - desiredoutput[i][j]) * (output - desiredoutput[i][j]);
		}
	}

	// compute the weight decay term
	double weight_decay_term = 0;
	for (int i=0; i<c_layercount; i++) {
		weight_decay_term += cblas_ddot(c_layers[i]->weights_size, c_layers[i]->weights, 1, c_layers[i]->weights, 1);
	}

	return errort / (2*m) + lambda * weight_decay_term / 2;
}

double chemnet::helpNumericalGradient(double *w, double **input, double **desiredoutput, int m, double lambda) {
	double EPSILON = 0.0001;
	*w += EPSILON;
	double j1 = squareCost(input, desiredoutput, m, lambda);
	*w -= 2*EPSILON;
	double j2 = squareCost(input, desiredoutput, m, lambda);
	*w += EPSILON;
	return (j1-j2)/(2*EPSILON);
}
	
// compute numerical gradient
std::vector<double> chemnet::computeNumericalGradient(double **input, double **desiredoutput, int m, double lambda) {
	std::vector<double> ans;

	for (int i=0; i<c_layercount; i++) {
		for (int j=0; j<c_layers[i]->weights_size; j++) {
			ans.push_back(helpNumericalGradient(c_layers[i]->weights + j, input, desiredoutput, m, lambda));
		}
		for (int j=0; j<c_layers[i]->neuroncount; j++) {
			ans.push_back(helpNumericalGradient(c_layers[i]->wbias + j, input, desiredoutput, m, lambda));
		}
	}

	return ans;
}

double *chemnet::outputError(double *input, double *desiredoutput, double &errort) {
	double *errorvalues = new double[c_layers[c_layercount-1]->neuroncount];
	// propagate the network
	propagate(input);
	// errorvalues = outputs - desiredoutput TODO assuming an identity output activation function
	cblas_dcopy(c_layers[c_layercount-1]->neuroncount, c_layers[c_layercount-1]->outputs, 1, errorvalues, 1);
	cblas_daxpy(c_layers[c_layercount-1]->neuroncount, -1.0, desiredoutput, 1, errorvalues, 1);
	// compute total error TODO assuming a square error and identity output activation function
	errort += cblas_ddot(c_layers[c_layercount-1]->neuroncount, errorvalues, 1, errorvalues, 1);

	return errorvalues;
}

// compute the gradiensts with the weight decay term
void chemnet::helpbpGradient(double *input, double *desiredoutput, double &errort) {
	double *errorvalues = outputError(input, desiredoutput, errort); // need to free the memory later
	// apply back propagation
	for (int i=c_layercount-1; i>=0; i--) {
		// update layerinput
		double *layerinput = i == 0 ? input : c_layers[i-1]->outputs;
		// errorvalues * layerinput -> gradient
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, c_layers[i]->neuroncount, c_layers[i]->inputcount, 1,
					1.0, errorvalues, 1, layerinput, c_layers[i]->inputcount, 1.0, c_layers[i]->weights_grad, c_layers[i]->inputcount);
		// wbias_grad + errorvalus
		cblas_daxpy(c_layers[i]->neuroncount, 1.0, errorvalues, 1, c_layers[i]->wbias_grad, 1);
		// weights' * errorvalues * layerinput * (1 - layerinput); not necessary when i == 0
		if (i > 0) {
            // back-propagate the error values to the next layer, assuming sigmoid function TODO
			double *errorvalues_t = new double[c_layers[i]->inputcount];
			cblas_dgemv(CblasRowMajor, CblasTrans, c_layers[i]->neuroncount, c_layers[i]->inputcount,
						1.0, c_layers[i]->weights, c_layers[i]->inputcount, errorvalues, 1, 0.0, errorvalues_t, 1); // beta == 0 so as to initilize errorvalues_t to zeros
			for (int j=0; j<c_layers[i]->inputcount; j++) {
				errorvalues_t[j] *= layerinput[j] * (1 - layerinput[j]);
			}
			delete [] errorvalues;
			errorvalues = errorvalues_t;
		}
		else
			delete [] errorvalues;
	}
}
	
// compute the gradients based on back-propagation
void chemnet::bpGradient(double **input, double **desiredoutput, int m, double lambda, double &errort) {
	// initialize the weights_grad and wbias_grad to be 0
	for (int i=0; i<c_layercount; i++) {
		memset(c_layers[i]->weights_grad, 0, c_layers[i]->weights_size * sizeof(double));
		memset(c_layers[i]->wbias_grad, 0, c_layers[i]->neuroncount * sizeof(double));
	}
	// accumulate the gradients without weight decay term, TODO GPU parallel
	for (int i=0; i<m; i++)
		helpbpGradient(input[i], desiredoutput[i], errort);

	// average the total error
	errort /= m;

	// add the weight decay term and average the gradients
	for (int i=0; i<c_layercount; i++) {
		// weights_grad / m
		cblas_dscal(c_layers[i]->weights_size, 1.0/m, c_layers[i]->weights_grad, 1);
		// weights_grad / m + lambda * weights
		cblas_daxpy(c_layers[i]->weights_size, lambda, c_layers[i]->weights, 1, c_layers[i]->weights_grad, 1);
		// wbias_grad / m
		cblas_dscal(c_layers[i]->neuroncount, 1.0/m, c_layers[i]->wbias_grad, 1);
	}
}

void chemnet::printParam(chemnet &net) {
    std::cout << "layers weights: " << std::endl;
    for (int i=0; i<c_layercount; i++) {
		std::cout << "Layer " << i << ": " << std::endl;
        for (int j=0; j<c_layers[i]->neuroncount; j++) {
            for (int k=0; k<c_layers[i]->inputcount; k++) {
               std::cout << c_layers[i]->weights[j * c_layers[i]->inputcount + k] << "\t";
            }
            std::cout << "wbias: " << c_layers[i]->wbias[j] << std::endl;
        }
        std::cout << std::endl;
    }
    return;
}
