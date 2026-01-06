#include <iostream>
#include <algorithm>
#include "chemnet.h"

#define PATTERN_COUNT 4
#define PATTERN_SIZE 2
#define NETWORK_OUTPUT 1
#define HIDDEN_LAYERS 2
#define EPOCHS 20000

double lambda = 0.0001;	// weight decay parameter
double alpha = 1;		// learning rate
double momentum = 0.5;	// momentum rate

using namespace std;

int main()
{
	// create some patterns w.r.t. XOR
	//double pattern[PATTERN_COUNT][PATTERN_SIZE] = { {0, 0}, {0, 1}, {1, 0}, {1, 1} };
	double **pattern = new double*[PATTERN_COUNT];
	double **desiredout = new double*[PATTERN_COUNT];
	for (int i=0; i<PATTERN_COUNT; i++) {
		pattern[i] = new double[PATTERN_SIZE];
		desiredout[i] = new double[NETWORK_OUTPUT];
	}
	pattern[0][0] = 0, pattern[0][1] = 0, pattern[1][0] = 0, pattern[1][1] = 1, pattern[2][0] = 1, pattern[2][1] = 0, pattern[3][0] = 1, pattern[3][1] = 1;
	desiredout[0][0] = 0, desiredout[1][0] = 1, desiredout[2][0] = 1, desiredout[3][0] = 0;
	// desired output values
	//double desiredout[PATTERN_COUNT][NETWORK_OUTPUT] = { {0}, {1}, {1}, {0} };
	// hidden layers
	int hiddenlayers[HIDDEN_LAYERS] = { 3, 4 };

	// create the neural network object
	chemnet net;
	// create the netork
	net.create(PATTERN_SIZE, NETWORK_OUTPUT, hiddenlayers, HIDDEN_LAYERS);
	//net.create(PATTERN_SIZE, NETWORK_INPUTNEURONS, NETWORK_OUTPUT, 0, 0); // TODO: need to provide the hiddenlayers *

	cout << "print initialized weights: " << endl;
	net.printParam(net);

	net.checkGradient(pattern, desiredout, 4, lambda);

	// start the training
	for (int i=0; i<EPOCHS; i++) {
		double error = 0;
		//error += net.train(pattern, desiredout, PATTERN_COUNT, alpha, momentum, lambda);
		error += net.train_adam(pattern, desiredout, PATTERN_COUNT, lambda);
		error /= PATTERN_COUNT;
		// display error
		cout << "ERROR: " << error << "\n";
		//if (error < 1e-10)
		//	break;
	}

	// test all patterns
	vector<double> ret = net.test(pattern, PATTERN_COUNT);
	for (int i=0; i<PATTERN_COUNT; i++) {
		cout << "TESTED PATTERN " << i << " DESIRED OUTPUT: " << *desiredout[i] << " NET RESULT: " << ret[i] << endl;
	}

	cout << "print trained weights: " << endl;
	net.printParam(net);

	return 0;
}
