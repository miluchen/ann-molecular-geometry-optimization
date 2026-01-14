#ifndef DATASET_H
#define DATASET_H

#include "config.h"
#include <rawdata.h>

// data structure that holds the input data and target values
// for now assuming we could read in all the data TODO
struct DataSet {
	int net_count;
	double **inputs;   // pointer on host
	double *outputs;  // pointer on host, it holds all the energies
	int *metas;	  // meda data that store info molecules
	long *_meta;	  // what passed to the neural net
	long *presum_metas; // TODO ??
	double **d_inputs; // pointer to device
	double *d_outputs;
	int dim_meta;	// dimension of meta data
	long count;		// total count of molecules, should be >= batch_size
	long pos;		// keep track of the current position of the input
	long batch_size;	// mini-batch size TODO maybe we don't need it here if we already have it in Trainer
	long *dim_ins;	// dimension of inputs
	long *dim_outs;	// dimension of outputs
	// for gradient back-propagation
	bool isGrad;
	double **input_grads; // pointer on host to gradients inputs
	double *output_grads; // pointer on host, it holds all gradients
	long *_gmeta;
	double **d_input_grads;
	double *d_output_grads;

	// std::streampos pos_in;	// position of current character in the input stream
	// std::streampos pos_out;
	// std::string inputfile;
	// std::string outputfile;
	DataSet();
	~DataSet();

	void create(int _net_count, std::vector<std::string> &_inputfiles, std::string &_outputfile, std::string &_metafile, long *_dim1, long *_dim2, long _count);
	void nextBatch(long _size);
	void randomShuffle(void);
	void dataSwap(double *v, long i, long j, long dim);
	void dataSwap(int *v, long i, long j, int dim);
	void readData(std::string &fname, char **_loc);
	void computePresum(long *a, int *b, long num, int dim);
	void printInfo(void);
	void printData(void);
	void printMolecule(long ind);
	void readGradients(std::vector<std::string> &_inputfiles, std::string &_outputfile); // _dim1, _dim2 are assumed to be the same
};

#endif
