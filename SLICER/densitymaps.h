#include <valarray>
#include <stdexcept>
#include <CCfits/CCfits>
#include "data.h"
#include "utilities.h"
#include "gadget2io.h"

using namespace std;
using namespace CCfits;

#define MAX_M 1e3    // Threshold for mass; Particler heavier than MAX_M will be attached zero mass
#define DO_NGP false // Use NGP as the MAS
#define numberOfLensPerSnap 4 // Number of Lens to be builded from a snap

/**
 * @brief Finds the snapshot index where the computed comoving distance is closest to the specified value.
 * 
 * This function reads a vector of snapshot redshifts and identifies the index of the snapshot
 * that has the comoving distance closest to 'dlens'. The comoving distance is calculated using
 * the provided GSL spline interpolation function.
 * 
 * @param zsnap Vector of redshifts for each snapshot.
 * @param GetDl GSL spline used for interpolating the comoving distance.
 * @param accGetDl GSL interpolation accelerator for the comoving distance calculation.
 * @param dlens The target comoving distance to match.
 * @return int The index of the snapshot closest to the specified comoving distance.
 */
int getSnap (vector <double> & zsnap, gsl_spline *GetDl,
                    gsl_interp_accel *accGetDl, double dlens);

/**
 * @brief Constructs lens planes based on simulation specifications and input parameters.
 * 
 * This function creates lens planes according to the provided simulation specifications
 * and stores them in the 'lens' structure. It uses various interpolation functions to
 * compute comoving distances and requires paths and redshifts of snapshots.
 * 
 * @param p Reference to the InputParams structure containing simulation parameters.
 * @param lens Reference to the Lens structure where the lens planes will be stored.
 * @param snapred Vector containing the redshifts of the snapshots.
 * @param snappath Vector containing the paths to the snapshot files.
 * @param snapbox Vector containing the box sizes for each snapshot (not used in the provided snippet).
 * @param GetDl GSL spline for comoving distance interpolation.
 * @param accGetDl GSL interpolation accelerator for GetDl.
 * @param GetZl GSL spline for redshift interpolation.
 * @param accGetZl GSL interpolation accelerator for GetZl.
 * @param numOfLensPerSnap Number of lens planes per snapshot.
 * @param myid Processor ID for conditional logging.
 * @return int Status code indicating success or failure of the function.
 */
int buildPlanes(InputParams &p, Lens &lens,
  vector <double> & snapred, vector <string> & snappath, vector <double> & snapbox,
  gsl_spline *GetDl,  gsl_interp_accel *accGetDl, gsl_spline *GetZl, gsl_interp_accel *accGetZl,
  int numOfLensPerSnap, int myid);

/**
 * @brief Randomizes simulation box conditions based on lens planes and input parameters.
 * 
 * Creates a plan for randomizing aspects of the simulation box, such as position and orientation,
 * based on the specified lens planes and stores it in the 'random' structure.
 * 
 * @param random Reference to the Random structure to store the randomization plan.
 * @param lens Reference to the Lens structure containing the lens planes.
 * @param p Reference to the InputParams structure containing the simulation parameters.
 * @param numOfLensPerSnap Number of lens planes per snapshot used for thickness calculations.
 * @param myid Processor ID for conditional logging.
 */
void randomizeBox (Random & random, Lens & lens, InputParams & p,
                                int numOfLensPerSnap, int myid);

/**
 * @brief Computes the number of necessary replications on the perpendicular plane.
 * 
 * Calculates how many times the simulation box should be replicated on the plane perpendicular
 * to the lens axis to match the specified field of view, given the dark-energy equation of state parameter w.
 * This is only used if the directive ReplicationOnPerpendicularPlane is enabled.
 * 
 * @param fov Field of view in degrees.
 * @param boxl Box size of the simulation.
 * @param Ds Comoving distance at source redshift.
 * @param myid Processor ID for conditional logging.
 * @param fovradiants Field of view in radians (output parameter).
 * @param nrepperp Number of replications on the perpendicular plane (output parameter).
 */
void computeReplications(double fov, double boxl, double Ds, int myid, double & fovradiants, int & nrepperp);

/**
 * @brief Tests if the chosen angular aperture is within allowed limits.
 *
 * This function checks if the specified field of view can be accommodated without needing to
 * replicate the simulation box in the plane parallel to the PLC (Particle Light Cone) axis, which
 * is not allowed by SLICER.
 *
 * @param fov Field of view in degrees.
 * @param boxl Length of the simulation box.
 * @param Ds Comoving distance to the source.
 * @param myid Processor identifier for logging purposes.
 * @param fovradiants Output parameter to store the field of view in radians.
 * @return int Returns 0 if the field of view is allowed, otherwise returns 1.
 */
int testFov(double fov, double boxl, double Ds, int myid, double & fovradiants);

/**
 * @brief Maps particles onto grids using the Triangular Shaped Cloud (TSC) Mass Assignment Scheme (MAS).
 *
 * This function processes particle data to map them onto specified grids based on their positions
 * and other simulation parameters. If USE_DGP is true, an alternative mapping scheme is used.
 *
 * @param fin Input stream to read particle masses for hydro simulations.
 * @param data Reference to the Header containing simulation specifics.
 * @param p Reference to InputParams containing the simulation parameters.
 * @param lens Reference to the Lens structure detailing the lens plan.
 * @param xx Arrays containing particle positions.
 * @param fovradiants The PLC angular aperture in radians.
 * @param isnap Current snapshot index or number of box repetitions so far.
 * @param mapxyi Reference to an array to store mapped particle data.
 * @param ntotxyi Reference to an array counting total particles mapped.
 * @param myid Processor identifier for conditional logging.
 * @return int Status of the mapping operation; typically returns 0 if successful.
 */
int mapParticles(ifstream & fin, Header &data, InputParams &p, Lens &lens,
    float* xx[6][3], double fovradiants, int isnap, valarray<float>(& mapxyi)[6],
    int(& ntotxyi)[6], int myid);

/**
 * @brief Creates density maps by mapping particles onto planes.
 *
 * Utilizes the mapping methodology described in the mapParticles routine to create density maps.
 * This involves interpolating comoving distances and other factors necessary for accurately
 * placing particles onto the simulation grid.
 *
 * @param p Reference to InputParams containing the simulation parameters.
 * @param lens Reference to the Lens structure detailing the lens plan.
 * @param random Reference to Random structure for accessing randomization plans.
 * @param isnap Index of the current snapshot.
 * @param ffmin Minimum file frame index.
 * @param ffmax Maximum file frame index.
 * @param File Name of the output file.
 * @param fovradiants Field of view in radians.
 * @param rcase A case identifier for processing.
 * @param GetDl GSL spline for comoving distance interpolation.
 * @param accGetDl GSL interpolation accelerator.
 * @param GetZl GSL spline for redshift interpolation.
 * @param accGetZl GSL interpolation accelerator.
 * @param mapxytot General map for total particle density.
 * @param mapxytoti Specific maps for different particle densities.
 * @param ntotxyi Total number of particles mapped.
 * @param myid Processor ID for logging.
 * @return int Returns status of the density map creation; 0 on success.
 */
int createDensityMaps (InputParams &p, Lens &lens, Random &random, int isnap,
  int ffmin, int ffmax, string File, double fovradiants, double rcase,
  gsl_spline *GetDl, gsl_interp_accel *accGetDl, gsl_spline *GetZl,
  gsl_interp_accel *accGetZl, valarray<float> &mapxytot,
  valarray<float> (&mapxytoti)[6], int (&ntotxyi)[6], int myid);

/**
 * @brief Writes density maps to disk.
 *
 * This function handles the output of density maps created by the mapParticles function to disk.
 * It formats the data appropriately and ensures that all necessary metadata and particle information
 * is correctly encoded.
 *
 * @param p Input parameters used in the simulation.
 * @param data Header information about the simulation.
 * @param lens Details of the lensing configuration used.
 * @param isnap Snapshot index used for this operation.
 * @param zsim Redshift at the simulation time.
 * @param snappl Path to the snapshot.
 * @param snpix Pixel label for the output file.
 * @param mapxytotrecv Array receiving the total mapped densities.
 * @param mapxytotirecv Arrays receiving specific mapped densities for different types.
 * @param ntotxyi Total number of particles involved in mapping.
 * @param myid Processor ID for conditional logging.
 */
void writeMaps (InputParams &p, Header &data, Lens &lens, int isnap, double zsim,
                 string snappl, string snpix, valarray<float>& mapxytotrecv,
                 valarray<float>(& mapxytotirecv)[6], int(& ntotxyi)[6], int myid);

/**
 * @brief Generates an output file name based on the simulation parameters and snapshot details.
 *
 * This function constructs a filename for output based on various simulation parameters and the specific
 * snapshot being processed, incorporating the provided label if specified.
 *
 * @param p InputParams structure with current simulation settings.
 * @param snappl Snapshot path label used as part of the filename.
 * @param label Optional label to include in
 */
string fileOutput (InputParams p, string snappl, int label = 0);
