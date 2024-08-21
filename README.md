# HydBCsForOF

**HydBCsForOF** is a set of hydraulic engineering boundary conditions for the Volume-of-Fluid solver "interfoam" of OpenFOAM. This is neither a part of OpenFOAM nor endorsed by the owners of OpenFOAM. The provided boundary conditions allow the specification of water flow rates (variable water levels) and the prescribtions of water levels (variable flow rates).

## Publications

The numerical background is described in this paper, which we kindly ask you to cite if you are using **HydBCsForOF**:  
Thorenz, C. (2024): 'Boundary Conditions for Hydraulic Structures Modelling with OpenFOAM',
10th International Symposium on Hydraulic Structures, Zürich. ISSN 0374-0056 , DOI: [10.3929/ethz-b-000675949           ](https://doi.org/10.3929/ethz-b-000675949           )

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

A simple example is provided in the directory 'testcase'. Further details on the example are given in the accompanying paper (see above).

For parallel execution, you can execute the following commands, starting from the testcase folder:

```
cd inter 
cp 0/bak* 0 
setFields 
decomposePar 
mpirun -np 16 interFoam -parallel
```

This will change into the computation folder, copy and initialize the variable field files, decompose the case into 16 subdomains and run the testcase on 16 CPU cores in parallel. The results can already be checked while the job is running by loading 'results.foam' into Paraview (remember to switch from "Reconstructed Case" to "Decomposed Case" as you are running in parallel!) 

## Transfer to your own example

Basic knowledge in OpenFOAM is required to setup examples based on these boundary conditions. It is necessary that you know  how to

- generate a mesh with blockMesh and snappyHexMesh or any other tool of your choice
- setup an interFoam case (i.e. how to apply boundary conditions for all boundaries for all field variables)
- set the initial conditions
- run the case
- post-process the results

before going on.

If you want to use this example as a blueprint for your own cases, it is recommendable to copy the whole "inter" folder as a starting point. In the example, "xmin" is defined as an inlet boundary for a given flowrate while "xmax" and "zmax" are defined as outlet boundaries with a given water level. If your boundaries have different names, you have to adapt all files in the inter/0/bak folder. For each file, all boundary name entries (xmin,xmax, ... , column) must be changed to reflect the names of your own setup. 

If in your own exmple a patch named "inlet" should be your inlet and you if you have no patch "xmin", you can simply change the string "xmin" to "inlet" in all files of inter/0/bak and change the values according to your needs. If your example has patch a "outlet" which you want to act as a fixed water level instead of "xmax", just exchange "xmax" with "outlet" in all files and adapt the values where needed. If your case has multiple walls, you can duplicate the entries for "column" in all files and change the name from "column" to your choice as needed.
  
Take care: The entries for e.g. "xmin" (or any other boundary) in one of the files must fit to the entries for "xmin" (or any other boundary) in the other files. Entries for a boundary must be changed in all files, not only in one!

You can use the entries in all files of inter/0/bak for 

- "xmin" as an example for a fixed flowrate
- "xmax" as an example for a fixed water level
- "column" as an example for a wall

In order to facilitate the setting of boundary conditions, we provide a simple Bash script. First copy your mesh (i.e. the polyMesh folder) to the folder inter/constant. On the command line, change to the inter folder. There you execute the script 

```
./setBoundaries.sh
```

This will guide you through the process of setting the boundary conditions for your case.

Take care: This is just a basic setup. Check i.e. in the files in 0/bak the provided values for the turbulence model at your inlet: Are they reasonable? Adjust as necessary and copy the files to 0 before going on! Maybe you also want to or have to adjust some other, more special, entries.

## Development

**HydBCsForOF** is developed at the [Federal Waterways Engineering and Research Institute](https://www.baw.de/). Updates will be incorporated here in the repository by BAW. 

## History of changes

### 2024-08-21
- Provide a script for a simple, basic setup of boundary conditions.

### 2024-08-20
- Updated inter/system/fvSchemes to use GAMG preconditioner. Faster in some cases. 
- Updated README.md with minor corrections.

### 2024-08-09
- Updated inter/system/fvSchemes to use linearUpwind. More stable for div(rhoPhi,U).
- Updated README.md with further explanations


## License 

**HydBCsForOF** is distributed by the [Federal Waterways Engineering and Research Institute](https://www.baw.de/) 
and is freely available and open source, licensed under the 
[GNU General Public License 3](https://www.gnu.org/licenses/gpl.html). No 
See [LICENSE.txt](LICENSE.txt) for details.



