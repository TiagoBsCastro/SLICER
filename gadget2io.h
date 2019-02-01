/* Header with definitions to read Gadget2 snapshot format */
#include "utilities.h"
#include "data.h"

/* Reads the "file_in" snapshot and stores it on the "header" instance
   of Header. "fin" is a ifstream instance that is leaved open for further
   reading if "close" == false.
*/
int readHeader (string file_in, Header &header, ifstream & fin, bool close);

/* Tests if the snapshot has Hydro particles */
void testHydro(InputParams & p, Header & data);

/* Prints the Snapshot Header */
void printHeader (Header header);

/* Fastforward "n" variables of size "size" of the stream "fin". */
void fastforwardNVars (ifstream & fin, size_t size, size_t n);

/* Fastforward the fstream "fin" to variable labled "BLOCK".
   Prints "fin" metadata if myid == 0.
*/
void fastforwardToBlock (ifstream & fin, string BLOCK, int myid);

/* Reads the block "Block" of size "num x sizeof(T)" of snapshot "fin"
   and stores it on the "scalar" address.

   Produces monitoring outputs if myid == 0
*/
template<typename T>
void readBlock(ifstream & fin, size_t num, string block, T *scalar, int myid){

  T dummy; // Dummy vars to read x,y, and z
  fastforwardToBlock (fin, block, myid);
  /* Loop on different types */
  for (int pp=0; pp<num; pp++){

    fin.read((char *)&dummy, sizeof(dummy));
    scalar[pp] = dummy;

  }
};

/* Reads the position block of snapshot fin according to Header "data"
   and InputParams "p". The positions are stored in the address xx
   according to the randomization plane "random" and the box replication
   rcase.

   Produces monitoration messages if myid == 0.
*/
void readPos (ifstream & fin, Header &data, InputParams &p, Random &random,
                              int isnap, float* xx[6][3], float rcase,int myid);

/* Reads the velocity block of snapshot fin according to Header "data"
   and InputParams "p". The velocities are stored in the address vv
   according to the randomization plane "random" and the box replication
   rcase.

   Produces monitoration messages if myid == 0.
*/
void readVel (ifstream & fin, Header &data, InputParams &p, Random &random,
                              int isnap, float* vv[6][3], int myid);

/* Get the "halos" group velocity from the "subhalos" group. The velocity
   is randomized according to InputParams "p" and Randomizantion Plan
   "random".

   The code searches for the subhalos.id in snapshots "FILE".{0..nfiles}
   !! Be Aware that this evolves reading several files for each halo !!
   !! the complexity probably evolves with n^2 and should be avoided !!
   !! for very large boxes                                           !!
*/
int getGVel(SubFind &halos, InputParams &p, Random &random, string FILE, int isnap);

/* Get the "halos" group ID corresponding according to its position on the
   snapshot "FILE.ff"
*/
int getGID(SubFind &halos, string File, int ffmin, int ffi, int &nhalos);

/* Get the "halos" LineOfSight velocity
*/
void getLOSVel(SubFind &halos);

/* Get the "halos" cosmological redshift (not taking into account pec.vel)
*/
void getTrueZ(SubFind &halos, Header &data, gsl_spline *getZl,
  gsl_interp_accel *accGetZl);

/* Get the "halos" angular positions.
   Matches Pinocchio cordinate system.
*/
void getAngular(SubFind &halos);

/*
  Reads the snapshots available in "filredshiftlist" up to the first deeper
  than the source in InputParams "p". The path for the snapshots and its
  redshifts are stored in snappath and snapred respectively.
*/
int readRedList(string filredshiftlist, vector <double> & snapred, vector <string> & snappath, InputParams &p);
