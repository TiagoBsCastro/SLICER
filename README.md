# SLICER: Simulation Light Cone Builder

<img src="./slicer.png" width="440" alt="Logo">

[![C++ Build CI](https://github.com/TiagoBsCastro/SLICER/actions/workflows/ci.yml/badge.svg)](https://github.com/TiagoBsCastro/SLICER/actions/workflows/ci.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
[![Documentation](https://img.shields.io/badge/docs-available-green.svg)](https://tiagobscastro.github.io/SLICER/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11048430.svg)](https://doi.org/10.5281/zenodo.11048430)
[![Last Updated](https://img.shields.io/badge/updated-May%2024-orange.svg)](https://tiagobscastro.github.io/SLICER/)

<div align="justify">

**SLICER** stands for **Simulation Light Cone Builder**, a specialized tool crafted for astrophysicists and cosmologists. Developed in C++ and parallelized using MPI, SLICER efficiently generates mass maps from particles and past light cones of structures and substructures from SUBFIND catalogs, all from GADGET simulation snapshots (format 2).

## Key Features
- **High-Performance Computing**: Utilizes MPI for efficient parallel computing, enabling the handling of extensive cosmological datasets.
- **Versatile Data Processing**: It is capable of creating both mass maps and past light cones of structures. The mass maps can be written for all particle species or separately.
- **Outputs in FITS Format**: Delivers results for mass maps in the Flexible Image Transport System (FITS) format, which is extensively used in the astronomical community for storing scientific data.
- **PINOCCHIO-like PLC**: The past-light cone (PLC) of (sub)-structures is written with the same structure used by PINOCCHIO ([Pinocchio on GitHub](https://github.com/pigimonaco/Pinocchio)).
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
cmake .. -DCMAKE_PREFIX_PATH=$(Your Library Path)
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
   - **Snapshots Directory**: Directory containing the snapshots. Example: `/home/tcastro/L128N256/`.
   - **PLC Sim. Name**: The naming convention used for your simulation outputs. Example: `gadget`.
   - **Seed for Pos. Center**: Seed used for random positioning of the center. Example: `-229`.
   - **Seed for Pos. Reflec.**: Seed for random reflections in positioning. Example: `-230`.
   - **Seed for Axis Sel.**: Seed used for random axis selection. Example: `-231`.
   - **Part. in Planes**: Controls whether particle types are written to separate lens planes or collectively to the same lens plane. Set to `0` to write collectively, or `1` to separate by particle type.
   - **PLC Directory**: Output directory for the past light cones. Example: `/home/tcastro/test_`.
   - **PLC Suffix**: Suffix for the output files. Example: `0`.
   - **Part. Degradation**: `log₂` of the factor by which the particle number is reduced for faster processing. Example: `0` (no reduction).
   - **DE-EOS w**: Equation of state parameter for dark energy. Example: `-1.0`.

2. **Run SLICER**:
   Navigate to the build directory and execute SLICER:
   ```bash
   ./SLICER InputParams.ini
   ```

### Configuration Definitions and Compile-Time Directives

SLICER includes several predefined constants and compile-time options that you can adjust according to your specific research needs.

#### Special Definitions in Code
The following definitions within the SLICER code help control various functional aspects:

- `POS_U 1.0`: Converts units from BoxSize unit length to kpc/h, defined in [`gadget2io.h`](./SLICER/gadget2io.h).
- `MAX_M 1e3`: Establishes a mass threshold; particles exceeding this threshold are assigned zero mass. This limit is particularly useful as the chemical enrichment model may sometimes produce stellar particles with unreasonably high masses, defined in [`densitymaps.h`](./SLICER/densitymaps.h).
- `DO_NGP false`: Determines whether to use the Nearest Grid Point (NGP) instead of the Triangular-Shaped Cloud (TSC) for the Mass Assignment Scheme (MAS), defined in [`densitymaps.h`](./SLICER/densitymaps.h).
- `numberOfLensPerSnap 4`: Indicates the number of lens planes to construct from a single snapshot, defined in [`densitymaps.h`](./SLICER/densitymaps.h).

To modify these definitions, locate the respective `.h` files in the `SLICER` directory and adjust the `#define` statements accordingly.

#### Compile-Time Directives
These options can be defined at compile time to modify SLICER's behavior:

- **FixedPLCVertex**: When activated, centers snapshots consistently at \(x, y, z = 0.5, 0.5, 0.5\), ensuring a uniform center across snapshots.
- **ReplicationOnPerpendicularPlane**: By default, SLICER restricts replications of the simulation box in the direction perpendicular to the lens plane, limiting the maximum possible aperture. Enabling this directive permits such replications, although it's important to note that repeated replication might lead to wide-angle effects under the flat sky approximation used by SLICER.

To enable these compile-time directives, add them when running the CMake configuration command:
```bash
cmake -DUSE_FIXED_PLC_VERTEX=ON -DUSE_REPLICATION=ON ..
```

These settings empower advanced users to tailor SLICER’s functionality to better suit their computational and analytical requirements.

## Documentation

Visit our [Documentation](https://tiagobscastro.github.io/SLICER/) for full API details.

## Contributing

Contributions to SLICER are welcome! Check our [Contributing](./.github/CONTRIBUTING.md) session.

## Featured Projects

SLICER has been employed in a variety of research projects. This section highlights some of the significant projects and published papers that have utilized SLICER.

### The [BEHOMO](https://valerio-marra.github.io/BEHOMO-project/) Project

The BEHOMO (Beyond Homogeneity) project leverages Λ Lemaître-Tolman-Bondi N-body simulations to explore large-scale structures and their evolution. SLICER has been instrumental in producing the lensing maps from the BEHOMO simulations.

* Marra, V., Castro, T., Camarena, D., Borgani, S., & Ragagnin, A. (2022). The BEHOMO project: Λ Lemaître-Tolman-Bondi N-body simulations. Astronomy & Astrophysics, 664, A179. https://doi.org/10.1051/0004-6361/202243539

### [Magneticum](https://magneticum.org) Collaboration

Magneticum is a series of hydrodynamical simulations focusing on the impact of baryons on cosmological structures. Several studies under this collaboration have used SLICER to analyze the effects of baryonic processes on cosmic structures.

* Castro, T., Quartin, M., Giocoli, C., Borgani, S., & Dolag, K. (2018). The effect of baryons in the cosmological lensing PDFs. Monthly Notices of the Royal Astronomical Society, 478(1), 1305-1325. https://doi.org/10.1093/mnras/sty1117

* Martinet, N., Castro, T., Harnois-Déraps, J., Jullo, E., Giocoli, C., et al. (2021). Impact of baryons in cosmic shear analyses with tomographic aperture mass statistics. Astronomy & Astrophysics, 648, A115. https://doi.org/10.1051/0004-6361/202040155


### Collaborative Research

SLICER has also been used in other papers to quantify the impact of baryons on the cosmological analysis of weak lensing data.

* Harnois-Déraps, J., Martinet, N., Castro, T., Dolag, K., Giblin, B., et al. (2021). Cosmic shear cosmology beyond two-point statistics: A combined peak count and correlation function analysis of DES-Y1. Monthly Notices of the Royal Astronomical Society, 506(2), 1623-1650. https://doi.org/10.1093/mnras/stab1623

* Heydenreich, S., Brück, B., Burger, P., Harnois-Déraps, J., Unruh, S., et al. (2022). Persistent homology in cosmic shear - II. A tomographic analysis of DES-Y1. Astronomy & Astrophysics, 667, A125. https://doi.org/10.1051/0004-6361/202243868

* Burger, P. A., Friedrich, O., Harnois-Déraps, J., Schneider, P., Asgari, M., et al. (2023). KiDS-1000 cosmology: Constraints from density split statistics. Astronomy & Astrophysics, 669, A69. https://doi.org/10.1051/0004-6361/202244673
  
* Burger, P. A., Porth, L., Heydenreich, S., Linke, L., Wielders, N., et al. (2024). KiDS-1000 cosmology: Combined second- and third-order shear statistics. Astronomy & Astrophysics, 683, A103. https://doi.org/10.1051/0004-6361/202347986

* Alfradique, V., Castro, T., Marra, V., Quartin, M., Giocoli, C., & Monaco, P. (2024). A deconstruction of methods to derive one-point lensing statistics. [arXiv:2405.00147](https://arxiv.org/abs/2405.00147).

* Harnois-Deraps, J., & others. (2024). KiDS-1000 and DES-Y1 combined: Cosmology from peak count statistics. [arXiv:2405.10312](https://arxiv.org/abs/2405.10312).

### Ongoing Research

SLICER continues to be a tool of choice for upcoming studies and publications! Stay tuned for more ;)

## License

SLICER is licensed under the GNU GPL3 license. See [LICENSE](./LICENSE) for more details.

</div>
