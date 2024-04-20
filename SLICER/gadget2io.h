/**
 * @file gadget2io.h
 * @brief Interfaces for handling Gadget-2 simulation data, including reading snapshot headers,
 * particle data, and performing unit conversions and other necessary preprocessing steps.
 */

#ifndef GADGET2IO_H
#define GADGET2IO_H

#include <cstdint>
#include "utilities.h"
#include "data.h"

#define POS_U 1.0 // Unit conversion from BoxSize unit lengh to kpc/h

/**
 * @brief Reads the header from a snapshot file and stores it in the provided Header structure.
 *
 * This function opens the specified input file stream, reads the header information into the
 * provided Header instance, and optionally closes the file stream.
 *
 * @param file_in Path to the input file containing the snapshot header.
 * @param header Reference to the Header structure where the header data will be stored.
 * @param fin Reference to the ifstream instance used for reading the file.
 * @param close Boolean flag indicating whether to close the ifstream after reading.
 * @return int Returns 0 on successful read, otherwise returns an error code.
 */
int readHeader(string file_in, Header &header, ifstream &fin, bool close);

/**
 * @brief Tests if the snapshot contains hydrodynamic (Hydro) particles.
 *
 * This function checks the Header of the snapshot to determine if hydrodynamic particles
 * are present based on the simulation parameters.
 *
 * @param p Reference to the InputParams structure containing simulation settings.
 * @param data Reference to the Header structure containing the snapshot header data.
 */
void testHydro(InputParams &p, Header &data);

/**
 * @brief Prints the details of a snapshot header.
 *
 * This function outputs the contents of the provided Header structure to the standard output,
 * typically used for debugging or information purposes.
 *
 * @param header The Header structure whose contents are to be printed.
 */
void printHeader(Header header);

/**
 * @brief Advances the file stream by skipping a specified number of variables.
 *
 * This function moves the read position of the given file stream forward by skipping over
 * a specified number of variables, each of a given size.
 *
 * @param fin Reference to the ifstream from which data is being read.
 * @param size The size of each variable to skip.
 * @param n The number of variables to skip.
 */
void fastforwardNVars(ifstream &fin, size_t size, size_t n);

/**
 * @brief Advances the file stream to the start of a specified block.
 *
 * This function searches for and advances the file stream to the beginning of a block
 * specified by its label. It also prints the file stream's metadata if the processor ID
 * (myid) is 0, typically used for debugging or logging.
 *
 * @param fin Reference to the ifstream from which data is being read.
 * @param BLOCK The label of the block to which the file stream is to be advanced.
 * @param myid Processor ID used to determine if metadata should be printed.
 */
void fastforwardToBlock(ifstream &fin, string BLOCK, int myid);

/**
 * @brief Reads a block of data from a file stream into a provided array.
 *
 * This template function reads a specified number of elements from a file stream,
 * where each element is of type T, and stores them in the provided array. It first
 * advances to the specified block before reading the data. It also produces monitoring
 * outputs if the processor ID (myid) is 0.
 *
 * @tparam T The type of the elements to read from the file.
 * @param fin Reference to the ifstream from which data is being read.
 * @param num The number of elements to read.
 * @param block The label of the block from which to read the elements.
 * @param scalar Pointer to the array where read elements are stored.
 * @param myid Processor ID used for conditional monitoring outputs.
 */
template <typename T>
void readBlock(ifstream &fin, size_t num, string block, T *scalar, int myid)
{

  T dummy; // Dummy vars to read x,y, and z
  fastforwardToBlock(fin, block, myid);
  /* Loop on different types */
  for (size_t pp = 0; pp < num; pp++)
  {

    fin.read((char *)&dummy, sizeof(dummy));
    scalar[pp] = dummy;
  }
};

/**
 * @brief Reads particle positions from a snapshot and applies randomization and replication adjustments.
 *
 * This function reads the position block from a snapshot file stream and stores the positions
 * in a multi-dimensional array according to the randomization and box replication parameters.
 * It also handles conditional monitoring outputs based on the processor ID.
 *
 * @param fin Reference to the ifstream to read data from.
 * @param data Reference to the Header structure containing snapshot headers.
 * @param p Reference to the InputParams structure with simulation settings.
 * @param random Reference to the Random structure detailing the randomization plan.
 * @param isnap Snapshot index for logging and processing.
 * @param xx Array to store the read positions.
 * @param rcase Replication factor for box dimensions.
 * @param myid Processor ID for conditional monitoring.
 */
void readPos(ifstream &fin, Header &data, InputParams &p, Random &random,
             int isnap, float *xx[6][3], float rcase, int myid);

/**
 * @brief Reads particle velocities from a snapshot and applies randomization.
 *
 * This function reads the velocity block from a snapshot file stream and stores the velocities
 * in a multi-dimensional array according to the randomization settings. It outputs monitoring
 * messages based on the processor ID.
 *
 * @param fin Reference to the ifstream to read data from.
 * @param data Reference to the Header structure containing snapshot headers.
 * @param p Reference to the InputParams structure with simulation settings.
 * @param random Reference to the Random structure detailing the randomization plan.
 * @param isnap Snapshot index for logging and processing.
 * @param vv Array to store the read velocities.
 * @param myid Processor ID for conditional monitoring.
 */
void readVel(ifstream &fin, Header &data, InputParams &p, Random &random,
             int isnap, float *vv[6][3], int myid);

/**
 * @brief Retrieves group velocities from subhalos, applying randomization.
 *
 * This function calculates the group velocities of halos by accessing subhalos data.
 * It includes randomized components and requires reading several snapshot files, which
 * can be computationally intensive.
 *
 * @param halos Reference to SubFind structure containing halo data.
 * @param p Reference to InputParams structure with simulation settings.
 * @param random Reference to Random structure detailing the randomization plan.
 * @param FILE Base filename for snapshot data.
 * @param isnap Snapshot index for processing.
 * @return int Status code of the function (0 for success).
 */
int getGVel(SubFind &halos, InputParams &p, Random &random, string FILE, int isnap);

/**
 * @brief Retrieves the group ID based on halo positions within a specific snapshot.
 *
 * This function identifies and returns the group ID of halos based on their positions
 * within a specified snapshot file.
 *
 * @param halos Reference to SubFind structure containing halo data.
 * @param File Filename of the snapshot.
 * @param ffmin Minimum file frame index.
 * @param ffi File frame index currently processed.
 * @param nhalos Number of halos processed (output parameter).
 * @return int The group ID of the halos.
 */
int getGID(SubFind &halos, string File, int ffmin, int ffi, int &nhalos);

/**
 * @brief Computes the LineOfSight (LOS) velocities for halos.
 *
 * This function calculates and updates the LOS velocities for each halo in the SubFind structure.
 *
 * @param halos Reference to SubFind structure containing halo data.
 */
void getLOSVel(SubFind &halos);

/**
 * @brief Calculates the cosmological redshift for halos, excluding peculiar velocities.
 *
 * This function computes the true cosmological redshift of halos based on their comoving
 * distances, using GSL spline interpolation for precision.
 *
 * @param halos Reference to SubFind structure containing halo data.
 * @param data Reference to the Header containing the simulation specifics.
 * @param getZl GSL spline for redshift interpolation.
 * @param accGetZl GSL interpolation accelerator.
 * @param lens Lens structure detailing the lensing configuration.
 * @param isnap Snapshot index for processing.
 */
void getTrueZ(SubFind &halos, Header &data, gsl_spline *getZl,
              gsl_interp_accel *accGetZl, Lens lens, int isnap);

/**
 * @brief Computes the angular positions of halos, matching the coordinate system used in Pinocchio.
 *
 * This function calculates and updates the angular positions for halos within the SubFind structure.
 *
 * @param halos Reference to SubFind structure containing halo data.
 */
void getAngular(SubFind &halos);

/**
 * @brief Reads a redshift list file and extracts snapshot paths and redshifts until the source depth.
 *
 * This function reads from a file listing redshifts and associated snapshot paths, storing them
 * in provided vectors up to the first redshift deeper than the source defined in InputParams.
 *
 * @param filredshiftlist Filename of the redshift list file.
 * @param snapred Vector to store extracted redshifts.
 * @param snappath Vector to store corresponding snapshot paths.
 * @param snapbox Vector to store box sizes (not mentioned in parameters but assumed needed from context).
 * @param p Reference to InputParams for source depth and other settings.
 * @return int Status code of the function (0 for success).
 */
int readRedList(string filredshiftlist, vector<double> &snapred, vector<string> &snappath, vector<double> &snapbox, InputParams &p);

#endif