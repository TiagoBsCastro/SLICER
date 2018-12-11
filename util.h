#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <valarray>
#include <CCfits/CCfits>
#include <string.h>
#include "cosmology.h"
#include "utilities.h"
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

using namespace std;
const double speedcunit = 2.99792458e+3; // speed of light / H0/h

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
  double fov; // Field of View in Degrees (Will not be read from InputFiles)
  bool hydro; // Hydro or DM only sim (Will not be read from InputFiles)
  string simType; // Gadget or SubFind
  double rgrid; // Physical grid for matter density (Will not be read from InputFiles)
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
  int nplanes;                   // Number of lens planes
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
};
struct Gadget{
  // GADGET has 6 different particle type
  vector<float> xx0, yy0, zz0;
  vector<float> xx1, yy1, zz1;
  vector<float> xx2, yy2, zz2;
  vector<float> xx3, yy3, zz3;
  vector<float> xx4, yy4, zz4;
  vector<float> xx5, yy5, zz5;
};
struct SubFind{
  // SubFind has 1 particle type (Pinocchio PLC like structure)
  vector<uint32_t> id; //                 1) Group ID
  vector<double> truez;//                 2) True Redshift
  vector<float> xx0, yy0, zz0;//          3-5) comoving position (Mpc/h)
  vector<float> vx0, vy0, vz0;//          6-8) velocity (km/s)
  vector<float> m;//                      9) Mass (m200 crit. for Halos MSUB for subhalos)
  vector<double> theta;//                 10) Theta (degree)
  vector<double> phi;//                   11) Phi (degree)
  vector<double> vel;//                   12) Peculiar velocity along the line-of-sight (km/s)
  vector<double> obsz;//                  13) Observed redshift
  vector<uint32_t> nsub; //               14) Number of Subhalos

  SubFind(int); //Constructor declaration

};
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
int test_fov(double , double , double , int , double *);
void test_hydro(InputParams *, Header *);
void fastforwardToPos (ifstream &, int, int,  bool);
void fastforwardNVars (ifstream &, size_t, size_t, int);
void print_header (Header);
void fastforwardToBlock (ifstream &, string, int);
void ReadPos (ifstream &, Header *, InputParams *, Random *, int, float* xx[6][3], float rcase, int);
void ReadVel (ifstream &, Header *, InputParams *, Random *, int , float* vv[6][3], int);
void getPolar(double, double, double, double *, double *, double *);
valarray<float> gridist_w (vector<float>, vector<float> , vector<float>, int, bool);
int MapParticles(ifstream &, Header *, InputParams *, Lens *,
    float* xx[6][3], double, int, valarray<float>(& mapxyi)[6],
    int(& ntotxyi)[6], int);
void write_maps (InputParams *, Header *, Lens *, int, double, string, string,
  valarray<float>&, valarray<float>(& mapxytotirecv)[6], int(& ntotxyi)[6], int);
template<typename T>
void ReadBlock(ifstream & fin, size_t num, string block, T *scalar, int myid){

  T dummy; // Dummy vars to read x,y, and z
  fastforwardToBlock (fin, block, myid);
  /* Loop on different types */
  for (int pp=0; pp<num; pp++){

    fin.read((char *)&dummy, sizeof(dummy));
    scalar[pp] = dummy;

  }

};

void GetGVel(SubFind & , SubFind *);
void GetGID(SubFind &, SubFind *);
void GetTrueZ(SubFind &, Header *, gsl_spline *, gsl_interp_accel *);
void GetLOSVel(SubFind &);
void GetAngular(SubFind &);
void CreatePLC (SubFind &, Header *, InputParams *, string, int);
