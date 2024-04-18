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
template <class T>
string sconv (T &val, const char *fact)
{
  char VAL[20]; sprintf (VAL, fact, val);
  return string(VAL);
}

float weight (float ixx, float ixh, double dx);

void getPolar(double x, double y, double z, double &ang1, double &ang2, double &d, bool radec);

valarray<float> gridist_w (vector<float>, vector<float> , vector<float>, int, bool);

#endif
