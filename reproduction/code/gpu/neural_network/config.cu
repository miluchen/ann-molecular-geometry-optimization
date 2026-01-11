#include <fstream>
#include "config.h"

/*********************** Member functions for ConFig ****************************/
ConFig::ConFig(): input_dim(0), output_dim(0), training_count(0), validation_count(0), testing_count(0),
					hiddenlayercount(0), hiddenlayers(0), 
					epochs(0), error_threshold(0), batch_size(0), lambda(0), alpha(0), momentum(0),
					adam_alpha(0), adam_beta1(0), adam_beta2(0), adam_epsilon(0) {}

ConFig::~ConFig() {
	printf("START DESTORY ConFig!\n");
	if (hiddenlayercount > 0)
		delete [] hiddenlayers;
	printf("END DESTORY ConFig!\n");
}

// read in all the configuratoin data from a file
void ConFig::create(std::string &_configfile) {
	// TODO omit adam parameters for now
	std::ifstream is (_configfile);
	if (is) {
		is >> input_dim >> output_dim >> training_count >> validation_count >> testing_count >> hiddenlayercount;
		hiddenlayers = new long[hiddenlayercount];
		for (long i=0; i<hiddenlayercount; i++)
			is >> hiddenlayers[i];
		is >> epochs >> error_threshold >> batch_size >> lambda >> alpha >> momentum;
		is >> training_inputf >> training_outputf;
		is >> validation_inputf >> validation_outputf;
		is >> testing_inputf >> testing_outputf;
	}
	is.close();
}
