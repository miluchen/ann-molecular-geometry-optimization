#include <H5Cpp.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <iomanip>

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

    vector<float> coordinates(nFrames * nAtoms * 3);
    coordDS.read(coordinates.data(), H5::PredType::NATIVE_FLOAT);

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

    vector<double> energies(nFrames);
    energyDS.read(energies.data(), H5::PredType::NATIVE_DOUBLE);

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
    vector<char> speciesBuf(nAtoms * str_len);
    speciesDS.read(speciesBuf.data(), dtype);
    
    // Convert to vector of strings
    vector<string> species;
    for (size_t i = 0; i < nAtoms; ++i) {
        string str(speciesBuf.data() + i * str_len, str_len);
        // Remove null terminators and trailing whitespace
        str.erase(find(str.begin(), str.end(), '\0'), str.end());
        while (!str.empty() && str.back() == ' ') str.pop_back();
        species.push_back(str);
    }

    // ---------- Write to streams ----------
    // Write coordinates: format: x y z (one line per atom per frame)
    coordStream << fixed << setprecision(6);
    for (hsize_t f = 0; f < nFrames; ++f) {
        for (hsize_t a = 0; a < nAtoms; ++a) {
            size_t idx = f * nAtoms * 3 + a * 3;
            coordStream << coordinates[idx] << " "
                       << coordinates[idx + 1] << " "
                       << coordinates[idx + 2] << "\n";
        }
    }

    // Write atoms/species: format: species (one per atom, repeated for each frame)
    for (hsize_t f = 0; f < nFrames; ++f) {
        for (hsize_t a = 0; a < nAtoms; ++a) {
            atomStream << species[a] << "\n";
        }
    }

    // Write energies: format: energy (one per frame)
    energyStream << fixed << setprecision(10);
    for (hsize_t f = 0; f < nFrames; ++f) {
        energyStream << energies[f] << "\n";
    }

    // Write metadata: format: groupName nAtoms
    metaStream << groupName << " " << nAtoms << "\n";
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

    // Open output streams
    ofstream coordStream(outputPrefix + ".xyz");
    ofstream atomStream(outputPrefix + ".atom");
    ofstream energyStream(outputPrefix + ".e");
    ofstream metaStream(outputPrefix + ".meta");

    if (!coordStream.is_open() || !atomStream.is_open() || 
        !energyStream.is_open() || !metaStream.is_open()) {
        cerr << "Error: Could not open output files" << endl;
        return 1;
    }

    try {
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

        coordStream.close();
        atomStream.close();
        energyStream.close();
        metaStream.close();
    }
    catch (H5::Exception& e) {
        e.printErrorStack();
        return 2;
    }
    catch (exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 3;
    }

    return 0;
}

