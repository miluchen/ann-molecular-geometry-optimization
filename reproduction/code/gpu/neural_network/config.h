#ifndef CONFIG_H
#define CONFIG_H

#include <string>

struct ConFig {
	// configs for dataset
	long input_dim;
	long output_dim;
	long training_count;
	long validation_count;
	long testing_count;
	std::string training_inputf, training_outputf, validation_inputf, validation_outputf, testing_inputf, testing_outputf;
	// for the neural net hidden layers
	long hiddenlayercount;
	long *hiddenlayers;
	// for the training algorithm
	long epochs;
	double error_threshold;
	long batch_size;
	double lambda; // weight decay parameter
	double alpha; // learning rate
	double momentum; // momentum rate (not necessary for adam
	// for adam algorithm TODO not contained in chemnet
	double adam_alpha, adam_beta1, adam_beta2, adam_epsilon;

	ConFig();
	~ConFig();

	void create(std::string &configfile);
};

#endif
