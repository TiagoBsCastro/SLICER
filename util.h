#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <valarray>
#include <CCfits/CCfits>
#include <string.h>
#include "cosmology.h"
#include "utilities.h"

using namespace std;
using namespace CCfits;

const double speedcunit = 2.99792458e+3; // speed of light x H0/h

// parameters in the input file
struct InputParams{
  int npix; // Number of Pixels
  double boxl; //Box Size
  double zs; // Source Redshift
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

const int dummy = 14;
struct Header
{
  // long substituted with int32_t
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

struct Randomization
{
  int face;
  int sgnX, sgnY, sgnZ;
  double xc, yc, zc;
  Randomization(){};

  Randomization(int face, int sgnX, int sgnY, int sgnZ, double xc, double yc, double zc)
  {
      this->face = face;
      this->sgnX = sgnX;
      this->sgnY = sgnY;
      this->sgnZ = sgnZ;
      this->xc = xc;
      this->yc = yc;
      this->zc = zc;
  }
};

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

template <class T>
string sconv (T &val, const char *fact)
{
  char VAL[20]; sprintf (VAL, fact, val);
  return string(VAL);
}

int getSnap (vector <double> &, vector <double> &, vector <double> &, double );

template <class T>
int locate (const std::vector<T> &v, const T);

//float weight (float, float, double);

valarray<float> gridist_w (vector<float>, vector<float> , vector<float>, int, bool);

void readInput(struct InputParams *p, std::string name);

void getPolar(double, double, double, double *ra, double *dec, double *d);

template <class ForwardIterator>
void min_element (ForwardIterator *first, ForwardIterator *last,
  ForwardIterator *min_x, ForwardIterator *min_y, ForwardIterator *min_z){
  *min_x=*first; *min_y=*(first+1); *min_z=*(first+2);
  if (first==last) {}
  else if ((last-first)%3){
    cout << "Bad pointers" << endl;
    exit(-1);
  }
  else{
    do{
      if (*first<*min_x)
          *min_x=*first;
        if (*(first+1)<*min_y)
          *min_y=*(first+1);
        if (*(first+2)<*min_z)
          *min_z=*(first+2);
        first+=3;
    }while(first!=last);
  }
}

template <class ForwardIterator>
void max_element (ForwardIterator *first, ForwardIterator *last,
  ForwardIterator *max_x, ForwardIterator *max_y, ForwardIterator *max_z){
  *max_x=*first; *max_y=*(first+1); *max_z=*(first+2);
  if (first==last) {}
  else if ((last-first)%3){
    cout << "Bad pointers" << endl;
    exit(-1);
  }
  else{
  do{
    if (*first>*max_x)
      *max_x=*first;
    if (*(first+1)>*max_y)
      *max_y=*(first+1);
    if (*(first+2)>*max_z)
      *max_z=*(first+2);
    first+=3;
  }while(first!=last);
  }
}

void readParameters(string file_name,int *npix, double *boxl,
            double *zs, double *fov, string *filredshiftlist,
            string *filsnaplist, string *pathsnap,string *idc,
            int *seedcenter, int *seedface, int *seedsign,
            string *simulation, int *nfiles,string *partinplanes,
            string *directory,string *suffix,int *sn_opt, bool *do_NGP);

void build_plans(double dlup, InputParams *p, int numberOfLensPerSnap, int nsnaps, vector <double> & lred, vector <double> & zl, vector <double> & dl, vector <string> & lsnap, vector <double> & ld, vector <double> & ld2, vector <int> & replication, vector <double> & zfromsnap, vector <string> & fromsnap, vector <double> & zsimlens, vector <bool> & randomize, vector <int> & pll,int myid);

void randomize_box (vector <double> & x0, vector <double> & y0, vector <double> & z0, vector<int> & face, vector<int> & sgnX, vector<int> & sgnY, vector<int> & sgnZ, vector <int> & replication, int nrandom, vector <bool> & randomize, InputParams *p, int numberOfLensPerSnap, int myid);

void box_randomize(vector <double> &, vector <double> &, vector <double> &,
            vector <int> &, vector <int> &, vector <int> &, vector <int> &,
            int , int , int , int );

void plans_builder (vector<double> &, vector<double> &, double , double ,
            vector<int> &, vector<double> &, vector<double> &, vector<double> &, vector<double> &,
            vector<int> &, vector<int> &);

void fastforwardToPos(ifstream &,int, int, bool);

void fastforwardToMASS (ifstream &, int, Header *, int);

void fastforwardToBHMASS (ifstream & , int , Header *, int);

void read_dl(string, vector <double> &, vector <double> &, double);

int read_redlist(string, vector <double> &, vector <string> &, double);

void read_particles (ifstream *, Header, int, int, float, float *,
                      vector<double> &, vector<double> &, vector<double> &, vector<int> &,
                      vector<int> &, vector<int> &, vector<int> &);

void map_particles (float *, double *, int, int, float, Header, vector <double> &,
                          vector <double> &, int , valarray <float> & , int * , bool , bool );

void write_map (string ,string , string ,
                      string , string , string , int , int * ,
                      valarray <float> & , valarray <float> *, float , float , vector <double> , vector <double> ,
                      float , float , float , double * );

int read_header (string , Header *, ifstream &, bool close);

void read_pos (ifstream *,Header *, float *, float *, float *, float *, float *, float *);

void rand_pos (float * , Header , int , int , float , Randomization);

void read_mass (ifstream *, Header * , double * , double * , double * , double * , double * , double * );

void print_header (Header);

template <class T>
void read_block (ifstream *fin, T *ptr, string block_str){

  Block block;
  int32_t blocksize;

  string str="NULL";

  while(strcmp(str.c_str(),block_str.c_str())){

    fin->read((char *)&block, sizeof(block));
    str={block.name[0],block.name[1],block.name[2],block.name[3]};
    fin->read((char *)&blocksize, sizeof(blocksize));
    cout << blocksize << endl;
    if (!strcmp(str.c_str(),block_str.c_str()))
      fin->read((char *) ptr, blocksize);
    else
      cout << "Skiping Block: " << str.c_str() << endl;
      fin->seekg(blocksize,ios_base::cur);
    fin->read((char *)&blocksize, sizeof(blocksize));

  }
}
