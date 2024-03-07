# APPM
Asymptotic-Preserving Plasma Modelling (APPM) is a project for simulating plasma under both non-neutral and quasi-neutral regimes. 
Spcecifically, the (rescaled) Euler-Maxwell system parametrized by scaled Deybe length is solved for both regimes in a unified code.

For details, check 
[Master_thesis_TianweiYu.pdf](https://github.com/TianweiCSE/APPM/files/9423681/Master_thesis_TianweiYu.pdf)

A quick visulization of the result (E-field lambda = 0.5)

![E_lambda-0-5_withInsulator (1)](https://user-images.githubusercontent.com/59021698/186666005-cef9b4e5-308e-434e-8b60-ff1dd65b4464.gif)

## Dependence
* Eigen 3.3.7 
* hdf5 1.10.9

## Compilation
* On Windows (with MinGW)

        mkdir build
        cd build
        cmake .. -G "MinGW Makefiles"
        cmake --build .

* On Linux (with gcc)

        mkdir build
        cd build
        cmake ..
        make

## Usage

The ***simulation parameters*** are inputed through `AppmSolverParams.txt` where you can set the following:
- `maxIterations`: max iteration times
- `maxTime`: end time of simulation
- `lambda`: rescaled Debye length
- `itersPerWrite`: iteration interval to store snapshots

The ***mesh parameters*** are inputed through `primalMeshParams.txt` where you can set the following:
- `axialLayers`: cell number along z-axis (axial axis)
- `refinements`: refinement level of the cross section of the plasma region (mesh size is halved per +1)
- `outerLayers`: refinement level of the outer non-conducting region (mesh size is haved per +1)
- `zMax`: length of the cylindrical region (=5 by default)
- `electrodeRadius`: radius of the contact boundary, e.g. radius of the plasma region (=1 be default)

These two files locate at the source directory.

To run the code, 
- on windows, use command:

        appm

- on Linux, use command:

        ./appm 
     
## Data Visualization
The solutions are stored in ``solutions_<meshType>_<entityType>.xdmf`` which rely on ``snapshot-<iteration>.h5`` generated during the simulation. 

> For instance, ``solutions_dual_cell.xdmf`` represents the quantities associated on dual cells, i.e. fluid quantities (density, momentum, energy) and reconstructed E-field and B-field.

To visualize the solutions, use **Paraview** (5.10.0):
>1. open ``solutions_dual_cell.xdmf`` 
>2. select **Xdmf3ReaderT** 
>3. apply and add filter **Cell Data to Point Data** 
>4. apply and add filer **Glyph**   
>5. select the variable to visualize 
>6. adjust visualization parameters and apply
