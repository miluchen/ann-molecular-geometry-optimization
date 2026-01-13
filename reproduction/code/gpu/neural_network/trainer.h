#ifndef TRAINER_H
#define TRAINER_H

#include "config.h"
#include "dataset.h"
#include "cuchemnet.h"

class Trainer {
	public:
		// parameters
		long epochs;
		double error_threshold;
		long batch_size;
		double lambda; // weight decay parameter
		double alpha; // learning rate
		double momentum; // momentum rate (not necessary for adam
		// weights (initialized to be 0's): store the optimal weights along the training
		double *weights;
		// different kinds of errors along the training
		std::vector<double> training_errors, validation_errors;

		Trainer();
		~Trainer();

		void create(ConFig &);
		// main training function, it should load datasets first, partition the sets into training_set, validation_set, testing_set
		void train_main(EnergyNet &net, DataSet *training_set, DataSet *validation_set);
		// train the gradients
		void train_grad(EnergyNet &net, DataSet *training_set, DataSet *validation_set);
		// save the weights, training errors and validation errors
		void saveWeights(const std::string &, int);
		void saveErrors(const std::string&, const std::string&);
		// restart the training, "weightsfname" contains the weights?
		void restart(EnergyNet &net, DataSet *training_set, DataSet *validation_set, std::string &weightsfname);
};

#endif
