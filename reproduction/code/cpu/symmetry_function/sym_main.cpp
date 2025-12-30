#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include "sym_func.h"

using namespace std;

int main(int argc, char *argv[])
{
	// check the input argument
	if (argc != 3) {
		cout << "usage: <parameter-file> <input-xyz>" << endl;
		return -1;
	}

	string inputfile, outputfile, paramfile;
	paramfile	= paramfile  + argv[1];
	inputfile	= inputfile  + argv[2] + ".xyz";
	outputfile	= outputfile + argv[2] + ".aev";

	// set input and output streams
	ifstream xyzfile (inputfile);
	ofstream aevfile (outputfile);
	// check whether the streams are currently associated with the corresponding files
	assert(aevfile.is_open() && xyzfile.is_open());
	// set the precision for output
	aevfile.precision(10);

	int atom_num;
	string line;
	AEV maev;
	maev.create(paramfile);

	// read individual xyz file and process it, store the values inside of AEV or write it to a file
	while (xyzfile >> atom_num) {
		getline(xyzfile, line), getline(xyzfile, line);

		double *coords = new double[atom_num * 3 * sizeof(double)];
		char *atoms    = new char[atom_num * sizeof(char)];
		int i = 0;
		while (i < atom_num) {
			xyzfile >> atoms[i];
			xyzfile >> coords[3*i] >> coords[3*i+1] >> coords[3*i+2];
			i++;
		}
		double **aev = maev.buildaev(coords, atoms, atom_num);
		// write the atomic environment vector
		for (i=0; i<atom_num; i++) {
			aevfile << atoms[i] << "\t";
			for (int j=0; j<maev.tot_dim; j++) {
				aevfile << aev[i][j] << "\t";
			}
			aevfile << "\n";
		}
		delete [] coords;
		delete [] atoms;
		for (i=0; i<atom_num; i++)
			delete [] aev[i];
		delete [] aev;
	}

	return 0;
}
