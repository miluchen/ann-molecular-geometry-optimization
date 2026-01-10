/*
 * This program need input: configfile and optional file of weights
 * It starts training the ANN or restart the training
 */

#include <iostream>
#include <fstream>
#include <algorithm>
#include "cuchemnet.h"
#include "trainer.h"
#include "dataset.h"
#include "config.h"

using namespace std;

int main(int argc, char *argv[])
{
	if (argc != 2 && argc != 3) {
		printf("Input Error: ./a.out <configfile> OR ./a.out <configfile> <weights>\n");
		return -1;
	}
	/* read in the config parameters */
	ConFig configs;
	string configfile (argv[1]);
	configs.create(configfile);
	configs.printParam();

	/* prepare the dataset */
	printf("start creating dataset ... \n");
	DataSet *training_set   = new DataSet;
	training_set->create(configs.net_count, configs.training_inputfs, configs.training_outputf, configs.training_metaf, 
						 configs.input_dim, configs.output_dim, configs.training_count);
	training_set->printInfo();
	DataSet *validation_set = new DataSet;
	validation_set->create(configs.net_count, configs.validation_inputfs, configs.validation_outputf, configs.validation_metaf,
						   configs.input_dim, configs.output_dim, configs.validation_count);
	printf("training set and validation set created!\n");
	if (DEBUG) {
        training_set->printMolecule(2);
        for (int i=0; i<5; i++)
            training_set->nextBatch(2);
        printf("Done testing DataSet!\n");
    }

	/* create and build the neural net */
	std:: cout << "start creating neural net ... " << std::endl;
	EnergyNet net;
	net.create(configs.net_count, configs.input_dim, configs.output_dim, configs.hiddenlayercount, configs.hiddenlayers);
	if (argc == 3) {
		printf("Restart from weights in %s\n", argv[2]);
		string weightsfname (argv[2]);
		net.loadParam(weightsfname, false);
	}
	std::cout << "energy net created!" << std::endl;
	if (DEBUG) {
		net.printParam();
	}

	/* create the trainer */
	Trainer trainer;
	trainer.create(configs);
	printf("trainer created\n");
	/* capture the start time */
    cudaEvent_t start, stop;
    HANDLE_ERROR(cudaEventCreate(&start));
    HANDLE_ERROR(cudaEventCreate(&stop));
    HANDLE_ERROR(cudaEventRecord(start, 0));
	/* start training */
	trainer.train_main(net, training_set, validation_set);
	printf("training completed\n");
	/* get stop time, and display the timing results */
    HANDLE_ERROR(cudaEventRecord(stop, 0));
    HANDLE_ERROR(cudaEventSynchronize(stop));
    float elapsedTime;
    HANDLE_ERROR(cudaEventElapsedTime(&elapsedTime, start, stop));
    printf("Time for training: %3.1f ms\n", elapsedTime);
    HANDLE_ERROR(cudaEventDestroy(start));
    HANDLE_ERROR(cudaEventDestroy(stop));

	/* save the weights and training, validation errors */
	cout << "paramsize: " << net.getParamSize() << endl;
	trainer.saveWeights("__weights__", net.getParamSize());
	trainer.saveErrors("__training_errors__", "__validation_errors__");
	/* free space */
	delete training_set;
	delete validation_set;

	cout << "Reach the End of Main!\n" << endl;
	return 0;
}
