#include <valarray>
#include <stdexcept>
#include <CCfits/CCfits>
#include "data.h"
#include "utilities.h"

using namespace std;
using namespace CCfits;

#define MAX_M 1e3    // Threshold for mass; Particler heavier than MAX_M will be attached zero mass
#define POS_U 1.0    // Unit conversion from BoxSize unit lengh to kpc/h
#define DO_NGP false // Use NGP as the MAS
#define numberOfLensPerSnap 2 // Number of Lens to be builded from a snap

/*
 * Reads the vector with snapshots redshifts "zsnap" and returns the position
 * of the snapshot that the comoving distance computed on its redshift is the
 * closest to "dlens". Comoving distance is Computed with the GetDl interpolated
 * function.
*/
int getSnap (vector <double> & zsnap, gsl_spline *GetDl,
                    gsl_interp_accel *accGetDl, double dlens);

/*
 * Creates the lens Planes plan according to the InputParams "p", the sim.
 * specification contained in Header "header". The plan is stored in the
 * "lens" Structure
 *
 * - snapred is the vector with the snapshots redshifts
 * - snappath is the vector with the snapshots paths
 * - GetDl, accGetDl, GetZl, and accGetZl are auxiliary interp. func. to compute
 *   the comoving distance given z and vice-versa
 * - Each lens will be BoxSize/numberOfLensPerSnap thick
 * - if myid == 0: monitoring messages are produced
 */
void buildPlanes(InputParams &p, Lens &lens,
  vector <double> & snapred, vector <string> & snappath, vector <double> & snapbox,
  gsl_spline *GetDl,  gsl_interp_accel *accGetDl, gsl_spline *GetZl, gsl_interp_accel *accGetZl,
  int numOfLensPerSnap, int myid);

/*
 * Creates the Randomization plan according to the "lens" planes
 * and the params inside the InputParams "p". The plan is stored
 * it in the "random" structure.
 * - Each lens will be BoxSize/numberOfLensPerSnap thick
 * - if myid == 0: monitoring messages are produced
 *
 */
void randomizeBox (Random & random, Lens & lens, InputParams & p,
                                int numOfLensPerSnap, int myid);

/*
* Compute the number of replications on the perpendicular plane are necessary
* !!!! ONLY USED IF THE DIRECTIVE ReplicationOnPerpendicularPlane == 1 !!!!
*/
void computeReplications(double fov, double boxl, double Ds, int myid, double & fovradiants, int & nrepperp);

/*
 * Test if the chosen angular aperture is allowed.
 * !! SLICER do not allow repetitions of the Box in the plane parallel !!
 * !!                            the PLC axis                          !!
 */
int testFov(double fov, double boxl, double Ds, int myid, double & fovradiants);

/*
 * Maps the Particles on the grids "mapxyi[6]" using the TSC MAS
 * (unless USE_DGP is True) according to the InputParams "p", the
 * sim. specifications read from Header "data", the lens plan "lens"
 * the particle positions at "xx".
 *
 * - fovradiants is the PLC angular aperture in radians
 * - isnap is the number of boxes repetitions so far
 * - fin is the fstream snapshot file. Used to read masses for hydro part.
 * - if myid == 0: monitoring messages are produced
 *
 */
int mapParticles(ifstream & fin, Header &data, InputParams &p, Lens &lens,
    float* xx[6][3], double fovradiants, int isnap, valarray<float>(& mapxyi)[6],
    int(& ntotxyi)[6], int myid);

/*
 * Creates the density maps mapping the particles in the planes mapxy.
 * See the description of the mapParticles routine
 */
int createDensityMaps (InputParams &p, Lens &lens, Random &random, int isnap,
  int ffmin, int ffmax, string File, double fovradiants, double rcase,
  gsl_spline *GetDl, gsl_interp_accel *accGetDl, gsl_spline *GetZl,
  gsl_interp_accel *accGetZl, valarray<float> &mapxytot,
  valarray<float> (&mapxytoti)[6], int (&ntotxyi)[6], int myid);

/*
 * Writes down the density maps.
 * See the description of the mapParticles routine
 */
void writeMaps (InputParams &p, Header &data, Lens &lens, int isnap, double zsim,
                 string snappl, string snpix, valarray<float>& mapxytotrecv,
                 valarray<float>(& mapxytotirecv)[6], int(& ntotxyi)[6], int myid);

/* Outputfile name */
string fileOutput (InputParams p, string snappl, int label = 0);
