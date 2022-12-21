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
        
## Data Visualization
We use Paraview 
