# HF-SHELL
**by W. Ryssens & Y. Alhassid** 
   
[![DOI](https://zenodo.org/badge/288576572.svg)](https://zenodo.org/badge/latestdoi/288576572)
[![License: GPL v3](https://img.shields.io/github/license/wryssens/hf-shell)](https://www.gnu.org/licenses/gpl-3.0)

We present the code HF-SHELL for solving the self-consistent mean-field equations for configuration-interaction shell model Hamiltonians in the proton-neutron formalism. The code can calculate both ground-state and finite-temperature properties in the Hartree-Fock (HF), HF+Bardeen-Cooper-Schrieffer (HF+BCS), and the Hartree-Fock-Bogoliubov (HFB) mean-field approximations. Particle-number projection after variation is incorporated to reduce the grand-canonical ensemble to the  canonical ensemble, making the code particularly suitable for the calculation of nuclear state densities. The code does not impose axial symmetry and allows for  triaxial quadrupole deformations. The self-consistency cycle is particularly robust through the use of the heavy-ball optimization technique and the implementation of different options to constrain the quadrupole degrees of freedom.


### Compilation
---
The code can be compiled by executing 

`make` 

in the main directory. If you don't have mod/ and obj/ folders in your current directory, the compilation process will create them. 

The compilation process can be modified through the Makefile; relevant options are

* CXX     :  compiler
* CXXFLAGS:  flags to pass to the compiler (optimization level, ...)

Executing 

`make clean`

will clean all object and module files, forcing a complete new compilation next time around.


### Input 
-----

The code needs the following input

* A file detailing the single-particle model space. 
* A file that includes the (angular-momentum coupled) two-body matrix elements of the nuclear two-body interaction.
* (Optional) A file containing the reduced matrix elements of r^2 in the single-particle model space; these are used to calculate quadrupole matrix 
  elements.

Examples are included in the **examples/** folder. 

The details of the requested calculation should be provided to the code through STDIN. This input is built out of two different namelists, **&modelspace/** and **&config**. We refer the reader to the paper and the examples provided for the complete details on all input keywords.

### Execution 
-----

When the relevant input data is found in the **in.dat** file, execute the 
code with the following command

`./hf_shell.exe < in.dat ` 

### Output
-----

The code prints detailed output to STDOUT. We recommend the user pipes this output to a separate file using

`./hf_shell.exe < in.dat > output.out`

or

`./hf_shell.exe < in.dat | tee output.out`

for post-processing as desired. 

This output might be too detailed, and the code also outputs a table that summarizes the final results of all calculations. The table can be found in a file that is named following the **outfile** variable.

For detailed information on the format of the output, as well as formulas for the more complicated quantities calculated by the code, please refer to the paper.


### Auxiliary script
-----

We also provide an auxiliary python script **level_densities.py** (which relies on a working installation of Numpy) to perform post-processing of tabulated output files of HF-SHELL for the calculation of the nuclear state density.  This script is included in the **scripts/** folder.

### Examples
----

The examples folder includes ready-to-run example scripts for the following applications:

* Free energy surface calculation (including triaxial deformations) for Mg24
* Calculation of the pairing transition for 144Nd
* Calculation of the shape transition for 162Dy 
* Calculation of the state density for 162Dy

These are the same examples that are used to demonstrate the use of the code in the paper. To reproduce the results, run the corresponding example bash script. For more details, please see the comments in run.NUCXY.sh scripts.

All examples can be easily executed on a laptop computer: the Mg24 calculations take a few seconds, while the 144Nd and 162Dy calculations take a few minutes. The state density calculation for 162Dy takes somewhat longer (about ten minutes).

After running the calculations, the figures of the paper can be reproduced by running the plot.NUCXY.py scripts in the corresponding subfolders. These scripts are Python scripts and rely on working installations of the Numpy and the Matplotlib Python libraries. They produce the corresponding figures of the paper in eps format. Only for 24Mg do we not include a plotting script, as the post-processing of the code output is nontrivial; we have produced the plots in the paper with an in-house interpolation and plotting library.

### Citing
-----

If you use this code in research work, please consider citing us: 

>W. Ryssens & Y. Alhassid,  
>*Finite temperature mean-field approximations for shell model Hamiltonians; the code HF-SHELL*, submitted to EPJA.
