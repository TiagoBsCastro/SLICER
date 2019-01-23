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

int readInput(struct InputParams &, string, string &, bool &);
int read_header (string , Header &, ifstream &, bool);
int read_redlist(string, vector <double> &, vector <string> &, InputParams &);
void build_planes(InputParams &, Header &, Lens &, vector <double> &, vector <string> &,
   gsl_spline *, gsl_interp_accel *, gsl_spline *, gsl_interp_accel *, int, int);
void randomize_box (Random &, Lens &, InputParams &, int, int);
int test_fov(double , double , double , int , double &);
void test_hydro(InputParams &, Header &);
void fastforwardToPos (ifstream &, int, int,  bool);
void fastforwardNVars (ifstream &, size_t, size_t, int);
void print_header (Header);
void fastforwardToBlock (ifstream &, string, int);
void ReadPos (ifstream &, Header &, InputParams &, Random &, int, float* xx[6][3], float rcase, int);
void ReadVel (ifstream &, Header &, InputParams &, Random &, int , float* vv[6][3], int);
void getPolar(double, double, double, double &, double &, double &, bool);
valarray<float> gridist_w (vector<float>, vector<float> , vector<float>, int, bool);
int MapParticles(ifstream &, Header &, InputParams &, Lens &,
    float* xx[6][3], double, int, valarray<float>(& mapxyi)[6],
    int(& ntotxyi)[6], int);
void write_maps (InputParams &, Header &, Lens &, int, double, string, string,
  valarray<float>&, valarray<float>(& mapxytotirecv)[6], int(& ntotxyi)[6], int);


void GetGVel(SubFind &, SubFind &, InputParams &, Random &, string, int);
int GetGID(SubFind &, string, int);
void GetTrueZ(SubFind &, Header &, gsl_spline *, gsl_interp_accel *);
void GetLOSVel(SubFind &);
void GetAngular(SubFind &);
void CreatePLC (SubFind &, Header &, InputParams &, string, int);
