# SLICER: Simulation Light Cone Builder

[![C++ Build CI](https://github.com/TiagoBsCastro/SLICER/actions/workflows/ci.yml/badge.svg)](https://github.com/TiagoBsCastro/SLICER/actions/workflows/ci.yml)
![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)

**SLICER** stands for **Simulation Light Cone Builder**, a specialized tool crafted for astrophysicists and cosmologists. Developed in C++ and parallelized using MPI, SLICER efficiently generates mass maps from particles and past light cones of structures and substructures from SUBFIND catalogs, all from GADGET simulation snapshots (format 2).

## Key Features
- **High-Performance Computing**: Utilizes MPI for efficient parallel computing, enabling the handling of extensive cosmological datasets.
- **Versatile Data Processing**: It is capable of creating both mass maps and past light cones of structures. The mass maps can be written for all particle species or separately.
- **Outputs in FITS Format**: Delivers results for mass maps in the Flexible Image Transport System (FITS) format, which is extensively used in the astronomical community for storing scientific data.
- **PINOCCHIO like PLC**: The past-light cone (PLC) of (sub)-structures is written with the same structure used by PINOCCHIO ([Pinocchio on GitHub](https://github.com/pigimonaco/Pinocchio)).

## Dependencies
Before installing SLICER, ensure that the following dependencies are met:
- **C++ Compiler**: A modern C++ compiler (e.g., GCC, Clang)
- **MPI**: For parallel computation (e.g., MPICH, OpenMPI)
- **CMake**: For building the application (version 3.10 or higher)
- **GNU Scientific Library (GSL)**: Provides numerous mathematical routines
- **CFITSIO**: Library for reading and writing FITS files
- **CCfits** (optional): A C++ wrapper for CFITSIO, required only if using C++ interfaces

### Installing Dependencies
On Ubuntu, you can install these dependencies using:
```bash
sudo apt-get update
sudo apt-get install libgsl-dev libcfitsio-dev cmake mpich libccfits-dev
```

## Installing 

1. Clone the Repository
```bash
git clone https://github.com/TiagoBsCastro/SLICER.git
cd SLICER
```
2. Configure with CMake
Create a build directory and run CMake to configure the project:
```bash
mkdir build
cd build
cmake .. -DCMAKE_PREFIX_PATH=$HOME
```

3. Build the Project
Compile the project:
```bash
make
```

## Running the SLICER

Modify the InputParams.ini file to adjust the simulation parameters according to your requirements.

## To-Do List

* Remove numberOfLensPerPlane directive and add as an entry in InputParams.ini.

## Contributing

Contributions to SLICER are welcome!

## License

SLICER is licensed under the GNU GPL3 license. See LICENSE for more details.
