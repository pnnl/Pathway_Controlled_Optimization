# PCO_methods

 ## Pathway  Constrained Opitimization For Maximization of Growth and Entropy Production Rates in Metabolism

 This repo contains code to reproduce the results in the paper: An Approach to Learn Regulation to Maximize Growth and Entropy Production Rates in Metabolism

 (paper link)




## Setup

 We recommend running this repo in a conda environment detailed below.


### Important Prerequisite

 Note that the optimization is formulated in pyomo a python based optimization and modeling language. At this time Pyomo is only supported on Linux, Mac OS/X and other Unix variants


 ### create the conda environment

 ```console
 user@machine:~$ conda config --add channels conda-forge 
 user@machine:~$ conda config --set channel_priority strict
 user@machine:~$ conda create -n PCO_methods python=3.9
 user@machine:~$ conda activate PCO_methods 
 (PCO_methods) user@machine:~$ conda install pyomo pandas numpy scipy
 (PCO_methods) user@machine:~$ conda install jupyter pyutilib
 (PCO_methods) user@machine:~$ conda install cobra
 (PCO_methods) user@machine:~$ conda install ipopt
 (PCO_methods) user@machine:~$ conda install equilibrator-api
 ```

## Optimization Example

An example for running the PCO optimization on the Rubrum metabolic model from the paper is given in the jupyter notebook: **growth_optimization.ipynb**

The details for the optimization formulation are contained in **PCO_methods.py**

The model inputs are contained in the csv files: **Stoichiometric_matrix.csv**, **Equilibrium_constants.csv**, **Metabolites.csv** 

### **Important Note**

 The notebook growth_optimization.ipynb will not directly reproduce the results of the paper. By default ipopt comes installed with the linear solver Mumps. The optimization hyperparameters have been selected to achieve better convergence with the Mumps solver but will in general produce lower quality solutions than what can be achieved with other solver choices

 To reproduce the paper results the HSL linear solver MA57 must be used with ipopt and the optimization hyperparameters must be set to those given in the paper. A license for MA57 can be requested here: https://www.hsl.rl.ac.uk/licencing.html

## Rubrum Model Construction

Details for the construction of the Rubrum metabolic model can be found in the jupyter notebook **Rubrum-model-generation.ipynb**. 

The notebook steps through construction of the metabolic model from the metabolic pathways in the folder **model_files**. Note that this script does not need to be run in order to run the optimization example. 

Where possible standard free energies for reactions are computed with eQuilibrator : https://equilibrator.weizmann.ac.il However some reactions are not in eQuilibrator, details for the computation in those cases are given in the notebook.


---

### Disclaimer

This material was prepared as an account of work sponsored by an agency of the United States Government.  Neither the United States Government nor the United States Department of Energy, nor Battelle, nor any of their employees, norany jurisdiction or organization that has cooperated in the development of these materials, makes any warranty, express or implied, or assumes any legal liability or responsibility for the accuracy, completeness, or usefulness or any information, apparatus, product, software, or process disclosed, or represents that its use would not infringe privately owned rights.
Reference herein to any specific commercial product, process, or service by trade name, trademark, manufacturer, or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United StatesGovernment or any agency thereof, or Battelle Memorial Institute. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or any agency thereof.

PACIFIC NORTHWEST NATIONAL LABORATORY
operated by
BATTELLE
for the
UNITED STATES DEPARTMENT OF ENERGY
under Contract DE-AC05-76RL01830
