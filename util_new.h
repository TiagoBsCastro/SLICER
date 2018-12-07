#include <vector>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include "utilities.h"
#include "cosmology.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;
const double speedcunit = 2.99792458e+3; // speed of light x H0/h

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
template <class T>
string sconv (T &val, const char *fact)
{
  char VAL[20]; sprintf (VAL, fact, val);
  return string(VAL);
}

// parameters in the input file
struct InputParams{
  int npix; // Number of Pixels
  double zs; // Source Redshift
  double Ds; // Comoving Distance at zs (Will not be read from InputFiles)
  double fov; // Field of View in Degrees
  string filredshiftlist; // File with the redshift list it may contain three columns: snap 1/(1+z) z
  string pathsnap; // Path where the snaphosts are located
  string simulation; // Simulation name (prefix infront at the snap file)
  int seedcenter, seedface, seedsign; // Random Seeds
  bool partinplanes; // True: Each gadget particle type will have its own Map; False: One Map for All particle types
  string directory; // Directory to save FITS files
  string suffix; // Suffix to FITS files
  int snopt; // Shot-noise option: 0-No random Degradation; 1-Half particles degradation; 2- Three quarters particle degradation ...
};

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

struct Block
{
  int32_t blocksize1;
  int8_t alignment[4];
  char name[4];
  int8_t padding[8];
  int32_t blocksize2;
};

struct Lens{
  int nsnaps;                   // Number of lens planes
  vector <int>    replication;  // Number of repetitions of the i-th snapshot box
  vector <int>    pll;          // Lens Plane indice
  vector <string> fromsnap;     // From which snapshot the i-th lens plane were build
  vector <double> zsimlens;     // z of the i-th lens (z correspodent to d=1/2(ld+ld2))
  vector <double> ld;           // Start position of the i-th lens
  vector <double> ld2;          // End position of the i-th lens
  vector <double> zfromsnap;    // z from the snapshot selected to build the i-th lens
  vector <bool>   randomize;    // Bool variable to whether the positions should be re-randomized or not

};
struct Random{
  vector<double> x0, y0, z0;   // ramdomizing the center of the simulation [0,1]
  vector<int> face;            // face of the dice
  vector<int> sgnX, sgnY,sgnZ; // randomizing the box axis signs
}
// Operators to Read Header and Block
inline istream & operator>>(istream &input, Header &header)
{
  input.read((char *)&header, sizeof(header));
  return input;
};
inline istream & operator>>(istream &input, Block &block)
{
  input.read((char *)&block, sizeof(block));
  return input;
};

int readInput(struct InputParams *, string, string &, bool *);
int read_header (string , Header *, ifstream &, bool close);
int read_redlist(string, vector <double> &, vector <string> &, InputParams *);
void build_planes(InputParams *, Header *, Lens &, vector <double> &, vector <string> &, gsl_spline *,
                              gsl_interp_accel *, gsl_spline *, gsl_interp_accel *, int, int);
void randomize_box (Random &, Lens *, InputParams *, int, int);
