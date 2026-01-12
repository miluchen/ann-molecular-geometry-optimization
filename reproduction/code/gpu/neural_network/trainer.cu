#include "trainer.h"
#include <limits>
#include <cassert>
#include <fstream>

/*********************** Member functions for Trainer ****************************/
Trainer::Trainer(): epochs(0), error_threshold(0), batch_size(0), lambda(0), alpha(0), momentum(0), weights(0) {}

Trainer::~Trainer() { 
	printf("START DESTORY Trainer!\n");
	free(weights);
	printf("END DESTORY Trainer!\n");
}

void Trainer::create(ConFig &configs) {
	epochs          = configs.epochs;
	error_threshold = configs.error_threshold;
	batch_size      = configs.batch_size;
	lambda          = configs.lambda;
	alpha           = configs.alpha;
	momentum        = configs.momentum;
}

// the main training function
void Trainer::train_main(cuchemnet &net, DataSet &dataset) {
	// to make sure the model is trained on all training data at least once
	assert(epochs * batch_size >= dataset.training_set->count);
	// space to hold the optimal weights of the net
	weights = (double *)malloc(net.getParamSize() * sizeof(double));
	double min_valid_error = std::numeric_limits<double>::max();
	// after N cycles of mini-batch gradient descent, compute validation error, training error. (they have the same size
	long N = dataset.validation_set->count / batch_size; 
	N = N > 1 ? N : 1;
	// set the validation inputs and outputs on Device
	dataset.validation_set->nextBatch(dataset.validation_set->count);
	// starting training
	double running_errors = 0;
	for (long i=1; i<=epochs; i++) {
		dataset.training_set->nextBatch(batch_size);
		double errort = net.train_adam(dataset.training_set->d_inputs, dataset.training_set->d_outputs, dataset.training_set->batch_m, lambda);
		if ((i % 100) == 0) {
			printf("TARINING ERROR: %.15f\n", errort);
		}
		running_errors += errort;
		// compute validation and training errors
		if ((i % N) == 0 || i == epochs) {
			validation_errors.push_back(net.evaluate(dataset.validation_set->d_inputs, dataset.validation_set->d_outputs, dataset.validation_set->count));
			if (validation_errors.back() < min_valid_error) {
				printf("Found better parameters at epochs: %d, validation error: %f\n", i, validation_errors.back());
				// store the new best weights
				min_valid_error = validation_errors.back();
				net.storeParam(weights);
				// write weights into file
				std::fstream os ("__weights__", std::ios::out);
				for (int w=0; w<net.getParamSize(); w++)
					os << weights[w] << " ";
				os.close();
			}
			if ((i % N) == 0)
				training_errors.push_back(running_errors / N);
			else
				training_errors.push_back(running_errors / (epochs % N));
			running_errors = 0;
			// check the threshold
			if (min_valid_error < error_threshold) {
				printf("Validation error smaller than threshold: validation error: %f\n", validation_errors.back());
				break;
			}
		}
	}
}
