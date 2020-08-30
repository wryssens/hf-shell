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

The compilation process can easily be modified through the Makefile; relevant options are

* CXX     :  compiler
* CXXFLAGS:  flags to pass to the compiler (optimization level, ...)

Finally, executing 

`make clean`

will clean all object and module files, forcing a complete fresh compilation next time around.


### Input 
-----

The code needs the following pieces of input

* A file detailing the single-particle model space. 
* A file comprised of the (angular momentum coupled) two-body matrix elements of the nuclear two-body interaction.
* (Optional) A file containing the reduced matrix elements of r^2 in the single-particle model space; these are used to calculated quadrupole matrix 
  elements.

Examples of all of these are included in the **examples/** folder. 

Finally, the details of the requested calculation should be fed to the code through STDIN. This input should be built out of two different namelists, **&modelspace/** and **&config**. We refer the reader to the paper and the examples provided for the full details on all input keywords.

### Execution 
-----

If the relevant input data is found on the **in.dat** file, one can execute the 
code with the following command

`./hf_shell.exe < in.dat ` 

### Output
-----

While running, the code prints extremely detailed output to STDOUT. We recommend the user pipes this output to a separate file with either

`./hf_shell.exe < in.dat > output.out`

or

`./hf_shell.exe < in.dat | tee output.out`

for post-processing as desired. 

This output is in practice often too complete, and to aid the user the code also outputs an extra table that summarizes the end result of all calculations, in a file that is named following the **outfile** variable.

For detailed information on the format of the output, as well as formulas for the more complicated quantities calculated by the code, please refer to the paper.


### Auxiliary script
-----

We also provide an auxiliary python script **level_densities.py** (which relies on a working installation of Numpy) to perform postprocessing of tabulated output files of HF-SHELL for the calculation of the nuclear state density.  It can be found in the **scripts/** folder.

### Examples
----

The examples folder includes ready-to-run example scripts for the following applications:

* Triaxial deformation surface calculation for Mg24.
* Calculation of the pairing transition for 144Nd.
* Calculation of the shape transition for 162Dy 
* Calculation of the state density for 162Dy.

These are the same calculations that are showcased in the paper. Reproducing these calculations can be done by simply running the appropriate example bash script. For more details, please see the comments in run.NUCXY.sh scripts.

All examples can be easily attempted on a modest laptop: the Mg24 calculations take a few seconds, the 144Nd and 162Dy calculations take a few minutes. The state density calculation for 162Dy is somewhat longer, it executes on one of the authors laptop in roughly ten minutes.

After running the calculations, the figures of the paper can be recreated by running the plot.NUCXY.py scripts in the relevant subfolders. These scripts are Python scripts and rely on working installations of the Numpy and the Matplotlib Python libraries. They produce the corresponding figure of the paper in .eps format. Only for 24Mg do we not include a plotting script, as the postprocessing of the codes output is nontrivial; we have produced the plots in the paper with an in-house interpolation and plotting library.

### Citing
-----

If you use this code in research work, please consider citing us. 

>W. Ryssens & Y. Alhassid,  
>*Finite temperature mean-field approximations for shell model Hamiltonians; the code HF-SHELL*, submitted to EPJA.
