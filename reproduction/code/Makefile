# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++14 -Wall

# Directories
BIN_DIR = bin
DATA_DIR = data

# Build executable for CPU version of symmetry function
aev_cpu:
	$(CXX) $(CXXFLAGS) cpu/symmetry_function/sym_main.cpp cpu/symmetry_function/sym_func.cpp -o $(BIN_DIR)/sym_cpu

# Build executable for GPU version of symmetry function
aev_gpu:
	$(CXX) $(CXXFLAGS) gpu/symmetry_function/sym_main.cu gpu/symmetry_function/sym_func.cu -o $(BIN_DIR)/sym_gpu

# Create bin directory if it doesn't exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Clean build artifacts
clean:
	rm -rf $(BIN_DIR)

.PHONY: clean aev_cpu