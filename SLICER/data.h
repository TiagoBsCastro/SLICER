/**
 * @file data.h
 * @brief Header with different data structures and how to read them.
 *
 * This class provides methods to read the different data structures involved.
 */
#ifndef DATA_H_
#define DATA_H_
#include <gsl/gsl_spline.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdint> // For compatibility with format2 structures
#include <iostream>
#include "utilities.h"

using namespace std;

/**
 * @struct InputParams
 * @brief Stores parameters read from an input file for simulation configuration.
 *
 * This structure contains settings and parameters that control the behavior of the simulation,
 * such as pixel count, redshift information, and paths for input and output files.
 */
struct InputParams
{
  int npix;                           // Number of Pixels
  double zs;                          // Source Redshift
  double Ds;                          // Comoving Distance at zs (Will not be read from InputFiles)
  double fov;                         // Field of View in Degrees (Will not be read from InputFiles)
  bool hydro;                         // Hydro or DM only sim (Will not be read from InputFiles)
  string simType;                     // Gadget or SubFind
  double rgrid;                       // Physical grid for matter density (Will not be read from InputFiles)
  string filredshiftlist;             // File with the redshift list it may contain three columns: snap 1/(1+z) z
  string pathsnap;                    // Path where the snaphosts are located
  string simulation;                  // Simulation name (prefix infront at the snap file)
  int seedcenter, seedface, seedsign; // Random Seeds
  bool partinplanes;                  // True: Each gadget particle type will have its own Map; False: One Map for All particle types
  string directory;                   // Directory to save FITS files
  string suffix;                      // Suffix to FITS files
  int snopt;                          // Shot-noise option: 0-No random Degradation; 1-Half particles degradation; 2- Three quarters particle degradation ...
  string snpix;                       // Label string for the output according to npix and physical options
  bool physical;                      // How the pixelization should be done: True on physical distances. False on angular positions
  double w;                           // Dark-Energy EOS w
};

/**
 * @struct Header
 * @brief Represents the header information in a Gadget-2 simulation file.
 *
 * This structure stores various metadata about a Gadget-2 simulation, including particle counts,
 * mass array, time, redshift, and other simulation parameters.
 */
const int dummy = 14;
struct Header
{
  int32_t npart[6];
  double massarr[6];
  double time;
  double redshift;
  int32_t flag_sfr;
  int32_t flag_feedback;
  uint32_t npartTotal[6];
  int32_t flag_cooling;
  int32_t numfiles;
  double boxsize;
  double om0;
  double oml;
  double h;
  int32_t flag_sage;
  int32_t flag_metals;
  int32_t nTotalHW[6];
  int32_t flag_entropy;
  int32_t la[dummy];
};

/**
 * @struct Block
 * @brief Stores block metadata in a Gadget-2 data file.
 *
 * This structure represents the metadata for a block in a Gadget-2 file, detailing its size,
 * alignment, and other necessary data for file parsing.
 */
struct Block
{
  int32_t blocksize1;
  int8_t alignment[4];
  char name[4];
  int8_t padding[8];
  int32_t blocksize2;
};

/**
 * @struct Lens
 * @brief Represents lensing configuration for the simulation.
 *
 * This structure stores information about the lens planes used in the simulation,
 * including their indices, associated snapshots, and z positions.
 */
struct Lens
{
  int nplanes;              // Number of lens planes
  vector<int> replication;  // Number of repetitions of the i-th snapshot box
  vector<int> pll;          // Lens Plane indice
  vector<string> fromsnap;  // From which snapshot the i-th lens plane were build
  vector<int> fromsnapi;    // The indice inside redshift list from which snapshot the i-th lens plane were build
  vector<double> zsimlens;  // z of the i-th lens (z correspodent to d=1/2(ld+ld2))
  vector<double> ld;        // Start position of the i-th lens
  vector<double> ld2;       // End position of the i-th lens
  vector<double> zfromsnap; // z from the snapshot selected to build the i-th lens
  vector<bool> randomize;   // Bool variable to whether the positions should be re-randomized or not
  vector<int> nrepperp;     // Number of repetitions on the perpendicular plane
};

/**
 * @struct Random
 * @brief Holds the parameters for randomizing aspects of the simulation.
 *
 * This structure contains vectors that store values used to randomize the simulation's
 * initial conditions, such as the simulation box's center and orientation.
 */
struct Random
{
  vector<double> x0, y0, z0;    // ramdomizing the center of the simulation [0,1]
  vector<int> face;             // face of the dice
  vector<int> sgnX, sgnY, sgnZ; // randomizing the box axis signs
};

/**
 * @struct Gadget
 * @brief Stores the positions of different types of particles in a GADGET simulation.
 *
 * GADGET simulations categorize particles into different types, each represented by a separate vector.
 * This structure holds the positions of these particles across different types.
 */
struct Gadget
{
  // GADGET has 6 different particle type
  vector<float> xx0, yy0, zz0;
  vector<float> xx1, yy1, zz1;
  vector<float> xx2, yy2, zz2;
  vector<float> xx3, yy3, zz3;
  vector<float> xx4, yy4, zz4;
  vector<float> xx5, yy5, zz5;
};

/**
 * @brief Constructor for the SubFind object, initializing its internal state and preparing it for use.
 *
 * This constructor sets up the initial configuration of the SubFind object, based on the provided
 * parameters. It ensures that all necessary properties are set to default values if not specified
 * and prepares any internal resources needed for the object's operations.
 *
 * @param param1 Description of what param1 represents and how it's used in initialization.
 * @param param2 Description of what param2 does and its significance to the SubFind object.
 */
struct SubFind
{
  // SubFind has 1 particle type (Pinocchio PLC like structure)
  vector<uint32_t> id;         //                 1) Group ID
  vector<uint32_t> fsub;       //               2) FOF group index
  vector<double> truez;        //                 3) True Redshift
  vector<float> xx0, yy0, zz0; //          4-6) comoving position (Mpc/h)
  vector<float> vx0, vy0, vz0; //          7-9) velocity (km/s)
  vector<float> m;             //                      10) Mass (m200 crit. for Halos MSUB for subhalos)
  vector<double> theta;        //                 11) Theta (degree)
  vector<double> phi;          //                   12) Phi (degree)
  vector<double> vel;          //                   13) Peculiar velocity along the line-of-sight (km/s)
  vector<double> obsz;         //                  14) Observed redshift
  vector<uint32_t> nsub;       //               15) Number of Subhalos

  SubFind(int, bool); // Constructor declaration
};

/**
 * @brief Stream extraction operator for reading Header structure data from an input stream.
 * @param input Reference to the input stream from which to read.
 * @param header Reference to the Header structure to populate with read data.
 * @return Reference to the updated input stream.
 *
 * Reads binary data corresponding to the Header structure directly from the input stream into the provided Header object.
 */
inline istream &operator>>(istream &input, Header &header)
{
  input.read((char *)&header, sizeof(header));
  return input;
};

/**
 * @brief Stream extraction operator for reading Block structure data from an input stream.
 * @param input Reference to the input stream from which to read.
 * @param block Reference to the Block structure to populate with read data.
 * @return Reference to the updated input stream.
 *
 * Reads binary data corresponding to the Block structure directly from the input stream into the provided Block object.
 */
inline istream &operator>>(istream &input, Block &block)
{
  input.read((char *)&block, sizeof(block));
  return input;
};

/**
 * @brief Reads the input parameters file and stores the values in the provided InputParams structure.
 *
 * This function opens and reads the contents of a named input parameters file. It parses the file
 * and populates the provided InputParams structure with the values read from the file. If the file
 * cannot be opened, the program will print an error message and terminate.
 *
 * @param p Reference to an InputParams structure where read parameters will be stored.
 * @param name Name of the file to read. The file should be in the current working directory.
 * @return int Returns 0 on success, or exits with 1 if the file cannot be opened.
 */
int readInput(struct InputParams &p, string name);

#endif
