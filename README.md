# GNBG-Generator
**MATLAB Source Code for Generalized Numerical Benchmark Generator (GNBG)

This repository contains the MATLAB source code for the Generalized Numerical Benchmark Generator (GNBG), an advanced tool for creating diverse and challenging numerical optimization problems.

Full description of GNBG can be found in "D. Yazdani, M. N. Omidvar, D. Yazdani, K. Deb, and A. H. Gandomi, "GNBG: A Generalized and Configurable Benchmark Generator for Continuous Numerical Optimization," arXiv prepring arXiv:2312.07083, 2023."

Features:

GNBG Problem Instance Generator: A versatile module to generate a wide array of optimization problems with user-defined characteristics and difficulty levels.
Sample Optimizer - Differential Evolution (DE): Includes an implementation of the Differential Evolution algorithm, serving as an example optimizer to test the generated problem instances.

Usage:

The GNBG generator allows for the customization and generation of unique problem instances tailored to specific research needs.
The included DE optimizer can be used to benchmark these instances.

Getting Started:

Clone or download the repository.
Open the MATLAB source files.
Customize the GNBG settings in BenchmarkGenerator.m as per your requirements to generate the desired problem instances.
Use the included DE algorithm in main.m to solve these instances, or integrate your own optimization methods for benchmarking.

Additional Information:

The repository aims to facilitate research in the field of numerical optimization by providing a flexible and user-friendly platform for generating and testing optimization problems.
For detailed instructions and more information about the capabilities of GNBG, refer to the documentation provided within the code.
