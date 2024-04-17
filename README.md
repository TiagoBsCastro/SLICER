# SLICER
SLICER (Simulation LIght conE buildeR) is designed to create mass maps (FITS files) from GADGET numerical simulation snapshot formats. It is compatible with both Dark Matter (DM) and Hydro-dynamic simulation runs.

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
