/* 
 * this program takes a configfile and weights as input
 * and do the testing on energy and gradient if possible
 */

#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>
#include "cuchemnet.h"
#include "trainer.h"
#include "dataset.h"
#include "config.h"

using namespace std;

int main(int argc, char *argv[])
{
	if (argc != 3) {
		printf("Error: ./neural <configfile> <weights>\n");
		return -1;
	}
	/* read in the config parameters */
	ConFig configs;
	string configfile (argv[1]);
	configs.create(configfile);
	configs.printParam();

	/* create and build the neural net */
	std:: cout << "start creating neural net ... " << std::endl;
	EnergyNet net;
	net.create(configs.net_count, configs.input_dim, configs.output_dim, configs.hiddenlayercount, configs.hiddenlayers);
	printf("Loading weights from %s\n", argv[2]);
	string weightsfname (argv[2]);
	net.loadParam(weightsfname, false);
	std::cout << "energy net created!" << std::endl;
	if (DEBUG) {
		net.printParam();
	}

	/* print and plot the errors */
	DataSet *testing_set = new DataSet;
	testing_set->create(configs.net_count, configs.testing_inputfs, configs.testing_outputf, configs.testing_metaf,
						configs.input_dim, configs.output_dim, configs.testing_count);
	testing_set->readGradients(configs.testing_grad_inputfs, configs.testing_grad_outputf);
	testing_set->nextBatch(testing_set->count);

	double *ret = (double *)malloc(configs.testing_count * sizeof(double));
	/* testing the energy error */
	double testing_error = net.evaluate(testing_set->d_inputs, testing_set->d_outputs, testing_set->count, testing_set->_meta, ret);
	printf("testing_error: %.15f\n", testing_error);
	/* compute the relative error for testing data */
    double relative_error = net.evaluateRelativeError(testing_set->d_inputs, testing_set->d_outputs, testing_set->count, testing_set->_meta);
    printf("relative testing_error: %.15f\n", relative_error);
	/* print out the results for energy predictions */
	cout.precision(numeric_limits<double>::digits10);
	for (long i=0; i<configs.testing_count; i++) {
		cout << testing_set->outputs[i] << "\t" << ret[i] << endl;
	}
	free(ret);

	/* testing the gradient error TODO */
	if (testing_set->isGrad) {
		double *retgrads = (double *)malloc(testing_set->batch_size * sizeof(double));
		double testing_grad_error = net.evaluateGradient(testing_set->d_inputs, testing_set->batch_size, testing_set->_meta, testing_set->d_input_grads, testing_set->d_output_grads, retgrads);
		printf("\n\nrelative GRADIENT testing_error: %.15f\n", testing_grad_error);
		for (long i=0; i<configs.testing_count; i++) {
			cout << testing_set->output_grads[i] << "\t" << retgrads[i] << endl;; 
		}
		free(retgrads);
	}


	cout << "\nReach the End of Main!\n" << endl;
	return 0;
}
