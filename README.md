# SLICER: Simulation Light Cone Builder

[![C++ Build CI](https://github.com/TiagoBsCastro/SLICER/actions/workflows/ci.yml/badge.svg)](https://github.com/TiagoBsCastro/SLICER/actions/workflows/ci.yml)
![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)

**SLICER** stands for **Simulation Light Cone Builder**, a specialized tool crafted for astrophysicists and cosmologists. Developed in C++ and parallelized using MPI, SLICER efficiently generates mass maps from particles and past light cones of structures and substructures from SUBFIND catalogs, all from GADGET simulation snapshots (format 2).

## Key Features
- **High-Performance Computing**: Utilizes MPI for efficient parallel computing, enabling the handling of extensive cosmological datasets.
- **Versatile Data Processing**: It is capable of creating both mass maps and past light cones of structures. The mass maps can be written for all particle species or separately.
- **Outputs in FITS Format**: Delivers results for mass maps in the Flexible Image Transport System (FITS) format, which is extensively used in the astronomical community for storing scientific data.
- **PINOCCHIO like PLC**: The past-light cone (PLC) of (sub)-structures is written with the same structure used by PINOCCHIO ([Pinocchio on GitHub](https://github.com/pigimonaco/Pinocchio)).
- **Lensing post-processing**: SLICER is distributed with several Python post-processing tools that efficiently extract lensing maps from the mass maps. These tools can be found in the [`Lens/` directory](./Lens/).

## Dependencies
Before installing SLICER, ensure that the following dependencies are met:
- **C++ Compiler**: A modern C++ compiler (e.g., GCC, Clang)
- **MPI**: For parallel computation (e.g., MPICH, OpenMPI)
- **CMake**: For building the application (version 3.10 or higher)
- **GNU Scientific Library (GSL)**: Provides numerous mathematical routines
- **CFITSIO**: Library for reading and writing FITS files
- **CCfits** (optional): A C++ wrapper for CFITSIO

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

To run SLICER, follow these simple steps:

1. **Modify the `InputParams.ini` File**: Configure the simulation parameters by editing the [InputParams.ini](./examples/InputParams.ini) file. Here is a breakdown of what each parameter means:

   - **Number of Map Pixels**: Defines the square root of the number of pixels of the output map. Example: `256`.
   - **Source Redshift**: The redshift of the source plane. Example: `0.5`.
   - **Field of View**: Angular aperture of the field of view in degrees. Example: `2.0`.
   - **File with Snapshots**: Path to the file listing the snapshots to use. Example: [snapshot_list.txt](./examples/snapshot_list.txt).
   - **Snapshots Directory**: Directory containing the simulation snapshots. Example: `/home/tcastro/L128N256/`.
   - **PLC Sim. Name**: The naming convention used for your simulation outputs. Example: `gadget`.
   - **Seed for Pos. Center**: Seed used for random positioning of the center. Example: `-229`.
   - **Seed for Pos. Reflec.**: Seed for random reflections in positioning. Example: `-230`.
   - **Seed for Axis Sel.**: Seed used for random axis selection. Example: `-231`.
   - **Part. in Planes**: Controls whether particle types are written to separate lens planes or collectively to the same lens plane. Set to `0` to write collectively, or `1` to separate by particle type.
   - **PLC Directory**: Output directory for the past light cones. Example: `/home/tcastro/test_`.
   - **PLC Suffix**: Suffix for the output files. Example: `0`.
   - **Part. Degradation**: `log_2` of the factor by which the particle number is reduced for faster processing. Example: `0` (no reduction).
   - **DE-EOS w**: Equation of state parameter for dark energy. Example: `-1.0`.

2. **Run SLICER**:
   Navigate to the build directory and execute SLICER:
   ```bash
   ./SLICER

## To-Do List

* Remove numberOfLensPerPlane directive and add as an entry in InputParams.ini.

## Contributing

Contributions to SLICER are welcome!

## License

SLICER is licensed under the GNU GPL3 license. See LICENSE for more details.
