# HydBCsForOF

**HydBCsForOF** is a set of hydraulic engineering boundary conditions for the Volume-of-Fluid solver "interfoam" of OpenFOAM. this is neither a part of openfoam nor endorsed by the owners of OpenFOAM. The provided boundary conditions allow the specification of water flow rates for variable water levels and the prescribtions of water levels with variable flow rates.

## Installation

You have to install OpenFOAM in the ESI variant, version v2212: https://www.openfoam.com/news/main-news/openfoam-v2212. Further details are given in the accompanying paper (see below).
An the command line navigate to the folder "boundaryConditions" and execute "./AllWmake". THERE WAS A TYPO IN THE PAPER! 

## Development

**HydBCsForOF** is developed at the [Federal Waterways Engineering and Research Institute](https://www.baw.de/).

## Publications

The numerical background is described in this paper, which we kindly ask you to cite if you are using **HydBCsForOF**:  
Thorenz, C. (2024): 'Boundary Conditions for Hydraulic Structures Modelling with OpenFOAM',
10th International Symposium on Hydraulic Structures, Zürich. ISSN 0374-0056 , DOI: 10.3929/ethz-b-000675949

## License 

**HydBCsForOF** is distributed by the [Federal Waterways Engineering and Research Institute](https://www.baw.de/) 
and is freely available and open source, licensed under the 
[GNU General Public License 3](https://www.gnu.org/licenses/gpl.html). 
See [LICENSE.txt](LICENSE.txt) for details.



