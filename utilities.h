#ifndef UTILITIES_H_
#define UTILITIES_H_
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <gsl/gsl_errno.h>
#include <valarray>

using namespace std;

/*
 * created by:  Matthias Bartelmann, MPA Garching, 2003; ZAH, U. Heidelberg, 2006 - (bartelmann@uni-heidelberg.de)
 * modified by: Carlo Giocoli, ZAH-ITA Heidelberg, 2010; INAF-OABO Bologna,  2012 - (carlo.giocoli@unibo.it)
 */

const double pi=3.141592653589793238462643383279502884197;
const double pio2=1.57079632679489661923132169163975144209858;
const double twopi=6.283185307179586476925286766559005768394;
const double sqrt2=1.41421356237309504880168872420969807856967;
const double euler=0.5772156649015328606065120900824024310422;
const double nXbin=128.;
const double speedcunit = 2.99792458e+3; // speed of light / H0/h

const double tiny=1e-6;
const double tiny4 = 1e-4;
const double conv = 648000/M_PI;

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

void error (const std::string s);
void error (const bool f, const std::string s);
void gsl_error_handler (const char * reason, const char * file, int line, int gsl_errno);
void warning (const std::string s);

template <class T>
void fill_linear ( std::vector<T> &v, size_t n, T min, T max ){
  v.resize ( n );
  for ( size_t i = 0; i < n; i++ )
    v[i] = min + (max - min) * T (i)/T (n-1);
}
template <class T>
void fill_logarithmic ( std::vector<T> &v, size_t n, T min, T max ){
  v.resize (n);
  for ( size_t i = 0; i < n; i++ )
    v[i] = exp ( log (min) + ( log (max) - log (min) ) * T (i)/T (n-1) );
}
template <class T>
int locate (const std::vector<T> &v, const T x){
  size_t n = v.size ();
  int jl = -1;
  int ju = n;
  bool as = (v[n-1] >= v[0]);
  while (ju-jl > 1){
    int jm = (ju+jl)/2;
    if ((x >= v[jm]) == as)
      jl=jm;
    else
      ju=jm;
  }
  if (x == v[0])
    return 0;
  else if (x == v[n-1])
    return n-2;
  else
    return jl;
}

float weight (float ixx, float ixh, double dx);

void getPolar(double x, double y, double z, double &ang1, double &ang2, double &d, bool radec);

valarray<float> gridist_w (vector<float>, vector<float> , vector<float>, int, bool);

#endif
