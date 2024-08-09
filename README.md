# HydBCsForOF

**HydBCsForOF** is a set of hydraulic engineering boundary conditions for the Volume-of-Fluid solver "interfoam" of OpenFOAM. This is neither a part of openfoam nor endorsed by the owners of OpenFOAM. The provided boundary conditions allow the specification of water flow rates for variable water levels and the prescribtions of water levels with variable flow rates.

## Publications

The numerical background is described in this paper, which we kindly ask you to cite if you are using **HydBCsForOF**:  
Thorenz, C. (2024): 'Boundary Conditions for Hydraulic Structures Modelling with OpenFOAM',
10th International Symposium on Hydraulic Structures, Zürich. ISSN 0374-0056 , DOI: [10.3929/ethz-b-000675949    ](https://doi.org/10.3929/ethz-b-000675949    )

## Installation

You have to install OpenFOAM in the ESI variant, version v2212: https://www.openfoam.com/news/main-news/openfoam-v2212. Further details are given in the accompanying paper (see above).

You can copy and paste the following code directly into your command line, after you activated your OpenFOAM environment:


```
git clone https://github.com/baw-de/HydBCsForOF.git
cd HydBCsForOF/boundaryConditions
. ./AllWmake

```

THERE WAS A TYPO IN THE PAPER! YOU HAVE TO EXECUTE `. ./AllWmake` !

The compilation should finish with a line which specifies the location and name of the compiled library.

## Example

A simple example is provided in the directory 'testcase'. Further details are given in the accompanying paper (see above).

For parallel execution, you can execute following commands, starting from the testcase folder:

```
cd inter 
cp 0/bak* 0 
setFields 
decomposePar 
mpirun -np 16 interFoam -parallel
```

This will change into the computation folder, initialize the variable fields, decompose the case into 16 subdomains and run the testcase on 16 CPU cores in parallel. The results can be checked by loading 'results.foam' in Paraview already while the job is running.

## Transfer to your own example

If you want to use this example as a blueprint for your own cases, it is recommendable to copy the whole "inter" folder as a starting point. In the example, "xmin" is defined as an inlet boundary for a given flowrate while "xmax" and "zmax" are defined as a outlet boundaries with a given water level. If your boundaries have different names, you have to adapt all files in the inter/0/bak folder. For each file, all boundary name entries (xmin,xmax, ... , column) must be changed to reflect your own setup. Take care: The entries for e.g. "xmin" in all variable files belong to each other and must be changed in all files, not only in one!

## Development

**HydBCsForOF** is developed at the [Federal Waterways Engineering and Research Institute](https://www.baw.de/). Updates will be incorporated here in the repository by BAW. 

## History of changes

### 2024-08-09
- Updated inter/system/fvSchemes to use linearUpwind. More stable for div(rhoPhi,U).
- Updated README.md with further explanations


## License 

**HydBCsForOF** is distributed by the [Federal Waterways Engineering and Research Institute](https://www.baw.de/) 
and is freely available and open source, licensed under the 
[GNU General Public License 3](https://www.gnu.org/licenses/gpl.html). No 
See [LICENSE.txt](LICENSE.txt) for details.



