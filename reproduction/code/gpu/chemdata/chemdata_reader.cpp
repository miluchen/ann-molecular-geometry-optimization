#include "chemdata.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>

using namespace std;

// Read and print data from a binary file with ChemData format
void readAndPrintChemData(const string& filename, const string& type) {
    ifstream file(filename, ios::binary);
    
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }
    
    // Read ChemData header
    ChemData header;
    file.read((char*)&header, sizeof(ChemData));
    
    if (file.gcount() != sizeof(ChemData)) {
        cerr << "Error: Could not read ChemData header from " << filename << endl;
        file.close();
        return;
    }
    
    cout << "File: " << filename << endl;
    cout << "  Header:" << endl;
    cout << "    row_num:     " << header.row_num << endl;
    cout << "    col_num:     " << header.col_num << endl;
    cout << "    element_size: " << header.element_size << " bytes" << endl;
    cout << "    total size:   " << header.getSize() << " bytes" << endl;
    cout << endl;
    
    // Validate element size matches expected type
    unsigned long totalElements = header.row_num * header.col_num;
    unsigned long expectedSize = totalElements * header.element_size;
    
    if (expectedSize != header.getSize()) {
        cerr << "Warning: Size mismatch in header" << endl;
    }
    
    // Read and print data based on specified type
    cout << "  Data:" << endl;
    
    if (type == "char") {
        // Validate element size
        if (header.element_size != sizeof(char)) {
            cerr << "Error: Type 'char' specified but element_size is " << header.element_size 
                 << " (expected " << sizeof(char) << ")" << endl;
            file.close();
            return;
        }
        // Character data (atoms/species)
        char* data = new char[totalElements];
        file.read((char*)data, header.getSize());
        
        if (file.gcount() != (streamsize)header.getSize()) {
            cerr << "Error: Could not read all data from " << filename << endl;
            delete[] data;
            file.close();
            return;
        }
        
        cout << "    Type: char" << endl;
        cout << "    Values (showing first 20):" << endl;
        unsigned long toShow = (totalElements < 20) ? totalElements : 20;
        for (unsigned long i = 0; i < toShow; ++i) {
            if (data[i] == '\0' || data[i] == ' ') {
                cout << "'" << data[i] << "' ";
            } else {
                cout << data[i] << " ";
            }
            if ((i + 1) % 20 == 0) cout << endl << "      ";
        }
        if (totalElements > 20) {
            cout << "... (total " << totalElements << " elements)" << endl;
        }
        cout << endl;
        
        delete[] data;
    }
    else if (type == "int") {
        // Validate element size
        if (header.element_size != sizeof(int)) {
            cerr << "Error: Type 'int' specified but element_size is " << header.element_size 
                 << " (expected " << sizeof(int) << ")" << endl;
            file.close();
            return;
        }
        
        // Integer data (metadata)
        int* data = new int[totalElements];
        file.read((char*)data, header.getSize());
        
        if (file.gcount() != (streamsize)header.getSize()) {
            cerr << "Error: Could not read all data from " << filename << endl;
            delete[] data;
            file.close();
            return;
        }
        
        cout << "    Type: int" << endl;
        cout << "    Values (showing first 20):" << endl;
        unsigned long toShow = (totalElements < 20) ? totalElements : 20;
        for (unsigned long i = 0; i < toShow; ++i) {
            cout << data[i] << " ";
            if ((i + 1) % 10 == 0) cout << endl << "      ";
        }
        if (totalElements > 20) {
            cout << "... (total " << totalElements << " elements)" << endl;
        }
        cout << endl;
        
        delete[] data;
    }
    else if (type == "float") {
        // Validate element size
        if (header.element_size != sizeof(float)) {
            cerr << "Error: Type 'float' specified but element_size is " << header.element_size 
                 << " (expected " << sizeof(float) << ")" << endl;
            file.close();
            return;
        }
        
        // Float data (coordinates)
        float* data = new float[totalElements];
        file.read((char*)data, header.getSize());
        
        if (file.gcount() != (streamsize)header.getSize()) {
            cerr << "Error: Could not read all data from " << filename << endl;
            delete[] data;
            file.close();
            return;
        }
        
        cout << "    Type: float" << endl;
        cout << "    Values (showing first 20):" << endl;
        cout << fixed << setprecision(6);
        unsigned long toShow = (totalElements < 20) ? totalElements : 20;
        for (unsigned long i = 0; i < toShow; ++i) {
            cout << data[i] << " ";
            if ((i + 1) % 3 == 0) cout << endl << "      ";
        }
        if (totalElements > 20) {
            cout << "... (total " << totalElements << " elements)" << endl;
        }
        cout << endl;
        
        delete[] data;
    }
    else if (type == "double") {
        // Validate element size
        if (header.element_size != sizeof(double)) {
            cerr << "Error: Type 'double' specified but element_size is " << header.element_size 
                 << " (expected " << sizeof(double) << ")" << endl;
            file.close();
            return;
        }
        
        // Double data (energies, coordinates)
        double* data = new double[totalElements];
        file.read((char*)data, header.getSize());
        
        if (file.gcount() != (streamsize)header.getSize()) {
            cerr << "Error: Could not read all data from " << filename << endl;
            delete[] data;
            file.close();
            return;
        }
        
        cout << "    Type: double" << endl;
        cout << "    Values (showing first 20):" << endl;
        cout << fixed << setprecision(10);
        unsigned long toShow = (totalElements < 20) ? totalElements : 20;
        for (unsigned long i = 0; i < toShow; ++i) {
            cout << data[i] << " ";
            if (header.col_num == 3 && (i + 1) % 3 == 0) {
                cout << endl << "      ";
            } else if (header.col_num == 1 && (i + 1) % 5 == 0) {
                cout << endl << "      ";
            }
        }
        if (totalElements > 20) {
            cout << "... (total " << totalElements << " elements)" << endl;
        }
        cout << endl;
        
        delete[] data;
    }
    else {
        cerr << "Error: Unsupported type: " << type << endl;
        cerr << "  Supported types: char, int, float, double" << endl;
        file.close();
        return;
    }
    
    file.close();
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <binary_file> <type>" << endl;
        cerr << "  <type> must be one of: char, int, float, double" << endl;
        cerr << "  Example: " << argv[0] << " output.xyz float" << endl;
        cerr << "  Example: " << argv[0] << " output.atom char" << endl;
        cerr << "  Example: " << argv[0] << " output.e double" << endl;
        cerr << "  Example: " << argv[0] << " output.meta int" << endl;
        return 1;
    }
    
    string filename = argv[1];
    string type = argv[2];
    
    // Validate type
    if (type != "char" && type != "int" && type != "float" && type != "double") {
        cerr << "Error: Invalid type '" << type << "'" << endl;
        cerr << "  Supported types: char, int, float, double" << endl;
        return 1;
    }
    
    readAndPrintChemData(filename, type);
    
    return 0;
}

