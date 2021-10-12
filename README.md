# Bridget2021
## Bridget's Research Project

## Prerequisites

1. Operating System.
    *  SUSE Linux Enterprise Server 15.
    *  Ubuntu 14.10.
    *  MacOS.
    
2. GCC/G++ version 8.2.0 or above.

3. CMake 3.11 or above.

You can load these depencencies on Cori typing `source cori-modules.sh`.

## Dependencies
    
1. CombBLAS.
  * Download or clone CombBLAS from `https://github.com/PASSIONLab/CombBLAS.git`.
  * Export the path to this directory as an environment variable `COMBBLAS_HOME`.
   ```
      git clone https://github.com/PASSIONLab/CombBLAS.git
      export COMBBLAS_HOME=$PWD
   ```
  * The following commands can be used to build and install CombBLAS:
  ```
    cd $COMBBLAS_HOME/CombBLAS
    mkdir build
    mkdir install
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ../
    make -j4
    make install         
  ```
1. Bridget's Project.
  * The following commands can be used to build and install Bridget's Prject:
  ```
    mkdir build
    cd build
    cmake ../
    make`       
  ```
  * You can run it using ``srun -n 1 -c 1 ./bridget'' using a single core on a single node since in the example the matrix is of dimension 10. You can increase the number of cores and processes when increasing the size of the matrix.
