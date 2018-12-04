#ifndef UTILITIES_H_
#define UTILITIES_H_
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <valarray>
#include <assert.h>

/** 
 * created by:  Matthias Bartelmann, MPA Garching, 2003; ZAH, U. Heidelberg, 2006 - (bartelmann@uni-heidelberg.de)
 * modified by: Carlo Giocoli, ZAH-ITA Heidelberg, 2010; INAF-OABO Bologna,  2012 - (carlo.giocoli@unibo.it)
 */

const double pi=3.141592653589793238462643383279502884197;
const double pio2=1.57079632679489661923132169163975144209858;
const double twopi=6.283185307179586476925286766559005768394;
const double sqrt2=1.41421356237309504880168872420969807856967;
const double euler=0.5772156649015328606065120900824024310422;
const double nXbin=128.;

const double tiny=1e-6;
const double tiny4 = 1e-4;
const double conv = 648000/M_PI;


void error (const std::string s);
void error (const bool f, const std::string s); 
void gsl_error_handler (const char * reason, const char * file, int line, int gsl_errno);
void warning (const std::string s);
double sign (const double a, const double b);
int nint (const double x);

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

void swap (double &a, double &b);

void sort (double *a, const int n);

std:: vector<double> rotatexYz(double theta,std:: vector<double> v);
std:: vector<double> rotatexyZ(double theta,std:: vector<double> v);

double median(std:: vector<double> vec);

float weight (float ixx, float ixh, double dx);

void InertialT(std:: valarray<float> f,std:: vector<double> x,int n,double f1, double f2,double f3, double f4, double f5, 
	       double &e1,double &e2,double &e3,double &e4, double &e5);

double getY(std:: vector<double> x, std:: vector<double> y,double xi);

double cubicInterpolate (double p[4], double x);

double bicubicInterpolate (double p[4][4], double x, double y);

double tricubicInterpolate (double p[4][4][4], double x, double y, double z);

double nCubicInterpolate (int n, double* p, double coordinates[]);

double getZ(std::valarray<float> map,int npix,double x,double y, std:: vector<double> xi,std:: vector<double> yi);

/**
 * Returns the pointer of the radial profile of the maps of a given field.
 * Spherical simmetry is assumed 
 */
double *estprof(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, double dr0, 

		double xmax,std:: vector<int> &vi, std:: vector<int> &vj,int ngal);
/**
 * Returns the pointer of the variance of the radial profile of the maps of a given field 
 * Spherical simmetry is assumed 
 */
double *estsigmaprof(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, double dr0, 
		     double xmax, std:: vector<int> &vi, std:: vector<int> &vj, int ngal, double* qm);
/**
 * Returns the pointer of the cumulative radial profile of the maps of a given field
 * Spherical simmetry is assumed 
 */
double *estcprof(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, double dr0, 
		 double xmax, std:: vector<int> &vi, std:: vector<int> &vj, int ngal);

double *estsigmacprof2(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, std:: vector<double> rbin,
		       std:: vector<int> &vi, std:: vector<int> &vj, int ngal, double* qm);

double *estprof2(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, std:: vector<double> rbin,
		 std:: vector<int> &vi, std:: vector<int> &vj,int ngal);

double *estsigmaprof2(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, std:: vector<double> rbin,
		      std:: vector<int> &vi, std:: vector<int> &vj, int ngal, double* qm);
double *estcprof2(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, std:: vector<double> rbin,
		  std:: vector<int> &vi, std:: vector<int> &vj, int ngal);

double *estsigmacprof2(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, std:: vector<double> rbin,
		       std:: vector<int> &vi, std:: vector<int> &vj, int ngal, double* qm);

double fcos(double x);

double fsin(double x);

#endif
