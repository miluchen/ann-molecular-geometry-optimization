# Reproduce Existing Researches
Here we focus on reproducing existing research studies, particularly those related to neural network potentials for organic molecules. Our goal is to implement and validate previously published methods and results, ensuring their reliability and correctness. This aids in advancing the field by providing transparent benchmarks and fostering further developments in data-driven molecular modeling.

## ANI-1
We aim to reproduce results leveraging the ANI-1 neural network potential, following the methodology described in Smith et al. [[1]](#1). The relevant dataset is available as described in [[2]](#2). Our first step will be to develop AEVLib, a GPU-accelerated library for computing atomic environment vectors. Next, we will build a C/C++ CUDA-based GPU-accelerated software for model training. Afterward, we plan to replicate additional related results.

### Code
The CPU implementation of the radial and angular symmetry functions can be found in [sym_func.cpp](./code/cpu/symmetry_function/sym_func.cpp). For an example of input parameters that define the behavior and dimensionality of these symmetry functions, refer to [param_sym.parm](./code/data/param_sym.parm).

## References

<a id="1">[1]</a>
Smith, J. S., Isayev, O., & Roitberg, A. E. (2017).
ANI-1: an extensible neural network potential with DFT accuracy at force field computational cost.
Chemical Science, 8, 3192–3203.
DOI: [10.1039/C6SC05720A](https://doi.org/10.1039/C6SC05720A)

<a id="2">[2]</a>
Smith, J. S., Nebgen, B., Lubbers, N., Isayev, O., & Roitberg, A. E. (2017).
Less is more: Sampling chemical space with active learning.
Scientific Data, 4, 170193.
DOI: [10.1038/sdata.2017.193](https://doi.org/10.1038/sdata.2017.193)