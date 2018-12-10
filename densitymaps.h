#include "util.h"
int CreateDensityMaps (InputParams *p, Lens *lens, Random *random, int isnap, int ffmin, int ffmax, string File, double fovradiants, double rcase, gsl_spline *GetDl, gsl_interp_accel *accGetDl,
                        gsl_spline *GetZl, gsl_interp_accel *accGetZl, valarray<float> &mapxytot, valarray<float> (&mapxytoti)[6], int (&ntotxyi)[6], int myid);
