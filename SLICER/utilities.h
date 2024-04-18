#ifndef UTILITIES_H_
#define UTILITIES_H_
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <gsl/gsl_errno.h>
#include <valarray>

using namespace std;

const double speedcunit = 2.99792458e+3; // speed of light / H0/h

//const double tiny=1e-6;
//const double tiny4 = 1e-4;
//const double conv = 648000/M_PI;

// conversion: double or int -> string
static const char fINT[] = "%i";
static const char fLONG[] = "%lli";
static const char fDP0[] = "%1.0f";
static const char fDP1[] = "%2.1f";
static const char fDP2[] = "%3.2f";
static const char fDP3[] = "%4.3f";
static const char fDP4[] = "%5.4f";
static const char fDP5[] = "%6.5f";
static const char ee3[] = "%4.3e";

/**
 * @brief Converts a value of any type to a formatted string based on the provided format specifier.
 *
 * This template function takes a value of a generic type T and converts it into a string
 * according to the specified format. It utilizes sprintf to perform the formatting, which
 * allows for precise control over the output format, similar to printf syntax. This function
 * is useful for creating formatted strings from numerical or other data types dynamically.
 *
 * @tparam T The type of the value to be converted into a string. Can be any type that is 
 *           compatible with the sprintf format specifiers.
 * @param val Reference to the value to be converted into a string.
 * @param fact The format specifier as a C-style string, which determines the formatting of the output string.
 *             It should be compatible with the type T (e.g., "%d" for integers, "%f" for floating-point numbers).
 * @return string A string object containing the formatted representation of 'val'.
 *
 * Example usage:
 * \code{.cpp}
 * int number = 42;
 * std::string result = sconv(number, "%d");  // result will be "42"
 * double pi = 3.14159;
 * std::string pi_str = sconv(pi, "%.2f");    // pi_str will be "3.14"
 * \endcode
 *
 * ## Note:
 * Ensure that the format specifier matches the type of the value being formatted to avoid runtime errors.
 * Be aware of buffer overflow risks associated with sprintf and ensure that the format string and the
 * buffer size are properly managed.
 */
template <class T>
string sconv (T &val, const char *fact)
{
  char VAL[20]; sprintf (VAL, fact, val);
  return string(VAL);
}

/**
 * @brief Calculates the mass assignment weight using the Triangular Shaped Cloud (TSC) kernel.
 *
 * This function computes the weight for a particle relative to a grid point based on the distance
 * between them, using the TSC kernel. The TSC kernel is a piecewise quadratic function that 
 * effectively assigns mass to grid points based on particle proximity, within a specified range.
 *
 * @param ixx The actual position of the particle.
 * @param ixh The position of the grid point.
 * @param dx The grid spacing, which defines the scale of the kernel.
 * @return float The calculated weight according to the TSC kernel. If the distance exceeds the range 
 *               of influence defined by the kernel, the weight will be zero.
 *
 * The weight is determined based on the relative distance between the particle and the grid point,
 * scaled by the grid spacing:
 * - If the distance is less than or equal to 0.5 times the grid spacing, the weight is calculated as:
 *   \f$\frac{3}{4} - x^2\f$, where \f$x\f$ is the scaled distance.
 * - For distances greater than 0.5 times and up to 1.5 times the grid spacing, the weight follows
 *   a quadratic drop-off: \f$0.5 \times (\frac{3}{2} - x)^2\f$.
 * - Beyond 1.5 times the grid spacing, the weight is 0, indicating no contribution to the grid point.
 */
float weight (float ixx, float ixh, double dx);

/**
 * @brief Converts Cartesian coordinates to polar coordinates or celestial coordinates (RA, Dec).
 *
 * This function converts Cartesian coordinates (x, y, z) to either polar coordinates (theta, phi, d)
 * or celestial coordinates (Right Ascension, Declination, distance) based on the 'radec' flag.
 *
 * @param x The x-coordinate in Cartesian coordinates.
 * @param y The y-coordinate in Cartesian coordinates.
 * @param z The z-coordinate in Cartesian coordinates.
 * @param ang1 Reference to the output angle 1; theta or RA depending on 'radec'.
 * @param ang2 Reference to the output angle 2; phi or Dec depending on 'radec'.
 * @param d Reference to the output distance from the origin.
 * @param radec If true, calculates Right Ascension (RA) and Declination (Dec); if false, calculates theta and phi.
 *
 * The function computes the distance 'd' as the Euclidean distance from the origin to the point (x, y, z).
 * Depending on the 'radec' flag, it calculates either:
 * - Polar coordinates: 
 *   - Theta (ang1) as the angle from the positive z-axis.
 *   - Phi (ang2) as the angle from the positive x-axis in the x-y plane.
 * - Celestial coordinates:
 *   - Right Ascension (ang1) from the y-z plane.
 *   - Declination (ang2) as the angle from the equatorial plane.
 */
void getPolar(double x, double y, double z, double &ang1, double &ang2, double &d, bool radec);

/**
 * @brief Distributes mass to a grid based on particle positions and weights using different mass assignment schemes.
 *
 * This function computes the distribution of mass onto a grid with specified dimensions. Particles are assigned
 * to the grid based on their positions (`x`, `y`) and weights (`w`). The mass assignment can either be
 * direct to the nearest grid point (NGP) if `do_NGP` is true, or distributed to neighboring points using a
 * TSC weight function if `do_NGP` is false. The function ensures particles are properly assigned within the
 * bounds of the grid, and that position and weight arrays are correctly matched in length.
 *
 * @param x Vector of x positions for particles, normalized between 0 and 1.
 * @param y Vector of y positions for particles, corresponding to `x`.
 * @param w Vector of weights for each particle.
 * @param nn The number of grid points per dimension (grid is nn x nn).
 * @param do_NGP Boolean flag to use Nearest Grid Point (NGP) method; if false, uses a TSC weight function for distribution.
 * @return valarray<float> A flat array representing the 2D grid, where mass values are assigned according to the particle weights and positions.
 *
 * ## Grid Distribution Details
 * The mass is distributed on a 2D grid based on the specified positions and weights:
 * - **NGP Mode**: Directly assigns the full weight of each particle to the nearest grid point.
 * - **Weighted Mode**: Assigns mass to a 3x3 block of grid points centered on the nearest point, adjusted by the weight function which accounts for the distance from each grid point to the actual particle position.
 *
 * ## Error Handling
 * Throws runtime errors if the lengths of `x`, `y`, and `w` do not match, ensuring data consistency.
 */
valarray<float> gridist_w (vector<float>, vector<float> , vector<float>, int, bool);

#endif
