/*
 *	Header file for back propagation in a fully connected neural net with Adam algorithm.
 */

#ifndef CHEMNET_H
#define CHEMNET_H

#include <vector>

struct layer {
	int inputcount;		// total number of elements in layerinput, assuming the network is fully connected for now
	int neuroncount;	// number of neurons
	int weights_size;	// inputcount * neuroncount
	double *weights;
	double *wbias;
	double *deltas;		// deltas of weights for the last update (used for momentum)
	double *outputs;
	double *weights_grad;
	double *wbias_grad;

	// data used for adam algorithm
	double *adam_weightsm;
	double *adam_weightsv;
	double *adam_wbiasm;
	double *adam_wbiasv;

	// constructor. all initialized to 0
	layer();
	// frees the memory used by the layer
	~layer();
	// create the layer and allocates memory
	void create(int _inputcount, int _neuroncount, int fan_out);
	// computes all neurons performing the network formula
	void calculate(double *layerinput, bool isOutput);
};

struct Adam_Setting {
	double alpha, beta1, beta2, epsilon, t;
	// default values for the Adam algorithm
	Adam_Setting(): alpha(0.001), beta1(0.9), beta2(0.999), epsilon(1e-8), t(0) {}
};

// the neural network structure
class chemnet {
	private:
		layer **c_layers;  // pointers to all layers
		int c_layercount;  // hiddenlayer count + 1
		Adam_Setting adam;

		// need an activation function to process the output layer (regression layer) TODO for now it's an identity function
		double outputActivation(double x);
		// computes the network values given an input pattern
		void propagate(double *input);	
		// compute gradients based on back-propagation
		void bpGradient(double **input, double **desiredoutput, int m, double lambda, double &errort);
		void helpbpGradient(double *input, double *desiredoutput, double &errort);
		double *outputError(double *input, double *desiredoutput, double &errort);
		// compute gradients numerically
		std::vector<double> computeNumericalGradient(double **input, double **desiredoutput, int m, double lambda);
		double helpNumericalGradient(double *w, double **input, double **desiredoutput, int m, double lambda);
		double squareCost(double **input, double **desiredoutput, int m, double lambda);

	public:
		// constructor: all values to 0
		chemnet();
		// destructor: free memories
		~chemnet();
		// create the network structure
		void create(int _inputcount, int _outputcount, int *_hiddenlayers, int _hiddenlayercount);
		// load weights and wbias
		void loadParam(double *_weights, int _weights_num, double *_wbias, int _wbias_num);
		// initialize the wbias for output layer differently
		void initOutputWbias(double **desiredoutput, int m);
		// Updates the weights given a desired output using backpropagation algorithm
		double train(double **input, double **desiredoutput, int m, double alpha, double momentum, double lambda);
		// train the model using Adam algorithm
		double train_adam(double **input, double **desiredoutput, int m, double lambda);
		// test the neural network
		std::vector<double> test(double **input, int m);
		// check the gradients
		void checkGradient(double **input, double **desiredoutput, int m, double lambda);
		// return the output layer
		inline layer *getOutput() { return c_layers[c_layercount-1]; }
		// print the parameters of the network
		void printParam(chemnet &net);
};

#endif
