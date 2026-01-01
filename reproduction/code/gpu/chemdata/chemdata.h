/* 
* This header contains the header definition for binary data structures used.
*/

#ifndef CHEMDATA_H
#define CHEMDATA_H

#define ENERGYSIZE sizeof(double)
#define XYZSIZE sizeof(double)
#define GRADSIZE sizeof(double)

struct ChemData {
	unsigned long row_num;	// total number of rows
	unsigned int col_num;	// total number of columns
	unsigned int element_size;	// size in byte of each element, assuming all element have the same size/type

	ChemData(): row_num(0), col_num(0), element_size(0) {}
	ChemData(unsigned long r, unsigned int c, unsigned int s): row_num(r), col_num(c), element_size(s) {}

	unsigned long getSize() { return row_num * col_num * element_size; }
};

#endif
