#include <H5Cpp.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <iomanip>
#include "../chemdata/chemdata.h"

using namespace std;

// ------------------------------------------------------------
// Read one molecule group and write to streams
// ------------------------------------------------------------
void readMolecule(const H5::Group& group, const string& groupName,
                  ofstream& coordStream, ofstream& atomStream, 
                  ofstream& energyStream, ofstream& metaStream)
{
    // ---------- coordinates ----------
    H5::DataSet coordDS = group.openDataSet("coordinates");
    H5::DataSpace coordSpace = coordDS.getSpace();

    if (coordSpace.getSimpleExtentNdims() != 3)
        throw runtime_error(groupName + ": coordinates must be 3D");

    hsize_t coordDims[3];
    coordSpace.getSimpleExtentDims(coordDims, nullptr);

    hsize_t nFrames = coordDims[0];
    hsize_t nAtoms  = coordDims[1];
    hsize_t xyzDim  = coordDims[2];

    if (xyzDim != 3)
        throw runtime_error(groupName + ": coordinates third dimension must be 3");

    float *coordinates = new float[nFrames * nAtoms * 3];
    coordDS.read(coordinates, H5::PredType::NATIVE_FLOAT);
    coordStream.write((char *)coordinates, nFrames * nAtoms * 3 * sizeof(float));
    delete [] coordinates;

    // ---------- energies ----------
    H5::DataSet energyDS = group.openDataSet("energies");
    H5::DataSpace energySpace = energyDS.getSpace();

    if (energySpace.getSimpleExtentNdims() != 1)
        throw runtime_error(groupName + ": energies must be 1D");

    hsize_t energyDims[1];
    energySpace.getSimpleExtentDims(energyDims, nullptr);

    if (energyDims[0] != nFrames)
        throw runtime_error(groupName + ": energies dimension (" + 
                           to_string(energyDims[0]) + ") must match coordinates first dimension (" + 
                           to_string(nFrames) + ")");

    double *energies = new double[nFrames];
    energyDS.read(energies, H5::PredType::NATIVE_DOUBLE);
    energyStream.write((char *)energies, nFrames * sizeof(double));
    delete [] energies;

    // ---------- species ----------
    H5::DataSet speciesDS = group.openDataSet("species");
    H5::DataSpace speciesSpace = speciesDS.getSpace();

    if (speciesSpace.getSimpleExtentNdims() != 1)
        throw runtime_error(groupName + ": species must be 1D");

    hsize_t speciesDims[1];
    speciesSpace.getSimpleExtentDims(speciesDims, nullptr);

    if (speciesDims[0] != nAtoms)
        throw runtime_error(groupName + ": species dimension (" + 
                           to_string(speciesDims[0]) + ") must match coordinates second dimension (" + 
                           to_string(nAtoms) + ")");

    // Read species as fixed-size strings (STRSIZE 1)
    H5::DataType dtype = speciesDS.getDataType();
    size_t str_len = dtype.getSize(); // Should be 1 for STRSIZE 1
    char *speciesBuf = new char[nAtoms * str_len];
    speciesDS.read(speciesBuf, dtype);
    for (size_t i = 0; i < nFrames; i++) {
        atomStream.write((char *)speciesBuf, nAtoms * str_len * sizeof(char));
    }
    delete [] speciesBuf;
    
    int *meta = new int[nFrames];
    for (size_t i = 0; i < nFrames; i++)
        meta[i] = nAtoms;
    metaStream.write((char *)meta, nFrames * sizeof(int));
    delete [] meta;

    cout << "Molecule group: " << groupName << endl;
    cout << "  Frames:  " << nFrames << endl;
    cout << "  Atoms:   " << nAtoms << endl;
}

// Process the H5 file, find the number of conformations and atoms of all molecules, and
// write the ChemData header for the output files
void writeChemDataHeader(const string& filename, ofstream& coordStream, ofstream& atomStream, 
    ofstream& energyStream, ofstream& metaStream) {
    int numberOfConformations = 0;
    int numberOfAtoms = 0;
    H5::H5File file(filename, H5F_ACC_RDONLY);
    H5::Group root = file.openGroup("/");
    hsize_t nObjs = root.getNumObjs();
    for (hsize_t i = 0; i < nObjs; ++i) {
        string topLevelName = root.getObjnameByIdx(i);
        H5G_obj_t type = root.getObjTypeByIdx(i);

        if (type == H5G_GROUP) {
            H5::Group topLevelGroup = root.openGroup(topLevelName);
            
            // Iterate over molecule groups (e.g., "gdb11_s01-0", "gdb11_s01-1", ...)
            hsize_t nMolObjs = topLevelGroup.getNumObjs();
            
            for (hsize_t j = 0; j < nMolObjs; ++j) {
                string molName = topLevelGroup.getObjnameByIdx(j);
                H5G_obj_t molType = topLevelGroup.getObjTypeByIdx(j);
                
                if (molType == H5G_GROUP) {
                    H5::Group molGroup = topLevelGroup.openGroup(molName);
                    // ---------- coordinates ----------
                    H5::DataSet coordDS = molGroup.openDataSet("coordinates");
                    H5::DataSpace coordSpace = coordDS.getSpace();
                    hsize_t coordDims[3];
                    coordSpace.getSimpleExtentDims(coordDims, nullptr);
                    numberOfConformations += coordDims[0];
                    numberOfAtoms += coordDims[1] * coordDims[0];
                }
            }
        }
    }
    ChemData header(numberOfConformations, 1, sizeof(int));
    metaStream.write((char *)&header, sizeof(ChemData));
    header = ChemData(numberOfConformations, 1, sizeof(double));
    energyStream.write((char *)&header, sizeof(ChemData));
    header = ChemData(numberOfAtoms, 1, sizeof(char));
    atomStream.write((char *)&header, sizeof(ChemData));
    header = ChemData(numberOfAtoms, 3, sizeof(float));
    coordStream.write((char *)&header, sizeof(ChemData));

    file.close();
}

void writeData(const string& filename, ofstream& coordStream, ofstream& atomStream, 
    ofstream& energyStream, ofstream& metaStream) {
    H5::H5File file(filename, H5F_ACC_RDONLY);
    // Open root group
    H5::Group root = file.openGroup("/");
    // Iterate over root-level groups (e.g., "gdb11_s01")
    hsize_t nObjs = root.getNumObjs();

    for (hsize_t i = 0; i < nObjs; ++i) {
        string topLevelName = root.getObjnameByIdx(i);
        H5G_obj_t type = root.getObjTypeByIdx(i);

        if (type == H5G_GROUP) {
            H5::Group topLevelGroup = root.openGroup(topLevelName);
            
            // Iterate over molecule groups (e.g., "gdb11_s01-0", "gdb11_s01-1", ...)
            hsize_t nMolObjs = topLevelGroup.getNumObjs();
            
            for (hsize_t j = 0; j < nMolObjs; ++j) {
                string molName = topLevelGroup.getObjnameByIdx(j);
                H5G_obj_t molType = topLevelGroup.getObjTypeByIdx(j);
                
                if (molType == H5G_GROUP) {
                    H5::Group molGroup = topLevelGroup.openGroup(molName);
                    string fullPath = "/" + topLevelName + "/" + molName;
                    readMolecule(molGroup, fullPath, coordStream, atomStream, 
                                energyStream, metaStream);
                }
            }
        }
    }
    file.close();
}

// ------------------------------------------------------------
// Main
// ------------------------------------------------------------
int main(int argc, char* argv[])
{
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <file.h5> <output_prefix>" << endl;
        return 1;
    }

    const string filename = argv[1];
    const string outputPrefix = argv[2];

    // Open output streams in binary mode
    ofstream coordStream(outputPrefix + ".xyz", ios::binary);
    ofstream atomStream(outputPrefix + ".atom", ios::binary);
    ofstream energyStream(outputPrefix + ".e", ios::binary);
    ofstream metaStream(outputPrefix + ".meta", ios::binary);

    if (!coordStream.is_open() || !atomStream.is_open() || 
        !energyStream.is_open() || !metaStream.is_open()) {
        cerr << "Error: Could not open output files" << endl;
        return 1;
    }

    writeChemDataHeader(filename, coordStream, atomStream, energyStream, metaStream);
    writeData(filename, coordStream, atomStream, energyStream, metaStream);

    coordStream.close();
    atomStream.close();
    energyStream.close();
    metaStream.close();

    return 0;
}

