#include "data.h"

/**
 * @brief Writes processed halo data to a binary output file in a specific format.
 *
 * This function takes halo data, applies transformations to their coordinates and velocities,
 * and writes them to a binary file. The data is filtered based on certain conditions related to
 * their angular positions and physical properties before being written. The output is structured
 * for compatibility with specific post-processing or visualization tools.
 *
 * @param halos Reference to the SubFind structure containing the halo data.
 * @param data Reference to the Header structure containing metadata like the box size.
 * @param p Reference to the InputParams structure containing parameters like field of view.
 * @param snappl Snapshot label used in generating the output file name.
 * @param ff Frame file index used in generating the output file name.
 *
 * ## Details
 * - The function computes the field of view in radians and uses it to filter halos based on their
 *   angular positions (`theta`, `phi`).
 * - Only halos within the specified field of view and with positive mass and redshift are processed.
 * - Each halo's position and velocity data are transformed and written to the file in a structured
 *   format that includes ID, position, velocity, mass, observed and true redshift, and angular coordinates.
 * - The angular coordinates are adjusted to a specific format (degrees, with adjustments for range and orientation).
 *
 * ## Output File Format
 * The data for each halo is written in a sequence of binary fields that include:
 * - Halo ID (unsigned long long int)
 * - True redshift (double)
 * - Positions (x, y, z as doubles)
 * - Velocities (vx, vy, vz as doubles)
 * - Mass (double)
 * - Angular coordinates (theta, phi as doubles)
 * - Peculiar velocity (double)
 * - Observed redshift (double)
 *
 * ## Usage Example
 * To use this function, ensure that the halo, header, and input parameters are properly set up and that
 * the snapshot label and frame file index are correctly specified.
 */
void writePLC (SubFind &halos, Header &data, InputParams &p, string snappl, int ff);
