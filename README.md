# 3D Navier Stokes Solver
3D Navier Stokes Solver using 5th Order WENO and Upwind scheme written in Julia with RK3 time stepping. 

## How to use
### Input
Please use the input.yaml to change the parameters to run the solver

#### Params
* nitermax: Max number of iterations
* nx,ny,nz: Grid number in each direction
* cfl: CFL number for the simulation
* Re: Reynolds number for the simulation
* nplot: Plot/Data save frequency
* tend: Total run time for simulation

#### Flags
* ivis: viscous terms enable/disable 0 or 1
* iflx: pick flux scheme 1 or 2  
* ihpc: 0 or 1 (some fixes to run on Headless machines)

### Output
Output format used is .vtk files and .pvd file corresponding to the snapshots generated during the run. 



## Packages used
* Plots
* OffsetArrays
* WriteVTK
* NPZ
* YAML
* Printf

## Tested with 
* Julia version 1.6.7 
* Julia version 1.9.2
