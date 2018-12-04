#ifndef   	COSMOLOGY_H_
# define   	COSMOLOGY_H_
#include "math.h"
#include <vector>
#include <gsl/gsl_spline.h>

/**
 * created  by: Matthias Bartelmann, MPA Garching, 2003; ZAH, U. Heidelberg, 2006 - (bartelmann@uni-heidelberg.de)
 * modified by: Carlo Giocoli, ZAH-ITA Heidelberg, 2010; INAF-OABO Bologna,  2012 - (carlo.giocoli@unibo.it)
 */

class cosmology{
  protected:
    std:: vector<double> vz;
    double om,la,h0,omb;
    double wq; // equation-of-state parameter
    bool cw;   // flag indicating whether wq is constant
    double ew; // evolution exponent for dark energy
    double ri; // internally used radiation-density parameter today
    double ng; // power of the scale factor in the early-time growth factor
    struct qaTable{
      std::vector<double> a;
      std::vector<double> q;
    };
    qaTable qa;
    double gslIntegrateQag (double (*)(double, void*), double, double);
    double timeEarly(double);
    gsl_interp_accel *acc;
    gsl_spline *splinezt;
    void initQTable ();
    void initialise();
  public:
    /**
     * Constructor of the class cosmology.
     * All attributes have default values. If they are omitted, they are set as follows:
     *
     * Omega  =  0.3
     * Lambda =  0.7
     * Hubble =  0.7
     * Omega_baryon = 0.04
     */
    cosmology(double omega=0.3,double lambda=0.7,double hubble=0.73,double wquint=-1,double baryon=0.0456);
    /**
     * Deconstructor of the class cosmology
     */
    ~cosmology();
    /**
     * Returns the redshift from the cosmic time (Gyr)
     */
    double zfromt(double t);
    /**
     * Returns the virial overdensity according to the spherical collapse model
     */
    double deltaVF (double zc=0.0, bool sf=false);
    /**
     * Returns the growth factor
     */
    double growthFactor (double a, bool ff=true);
    /**
     * Routine providing the coefficients of the second-order differential equation describing the
     * growth of the density contrast. This is used in the growth-factor calculation if no fitting
     * formula is or can be used.
     */
    void dDeltaDa (double x, const double *y, double *dy);
    /**
     * Returns the linear overdensity according to the spherical collapse model
     */
    double deltaC (double zc);
    /**
     * Returns the dynamical time scale of an halo at redshift z
     */
    double tdyn(double z);
    /**
     * Returns the factor by which the dark-energy density at the scale factor a differs from its
     * present value.
     */
    double q(double a);
    /**
     * Returns the expansion function, i.e. the right-hand side of Friedmann's equation
     * without H_0.
     */
    double e(double a);
    /**
     * Returns the matter-density parameter as a function of redshift z
     */
    double omegab(double z=0);
    /**
     * Returns the baryon-density parameter as a function of redshift z
     */
    double omega(double z=0);
    /**
     * Returns the dark-energy density parameter as a function of redshift z
     */
    double lambda(double z=0);
    /**
     * Returns the Hubble parameter at redshift z
     */
    double hubble(double z=0);
    /**
     * Returns the scale factor at matter-radiation equality
     */
    double aEqual();
    /**
     * Returns the derivative of the expansion function
     */
    double ePrime( double a );
    /**
     * Returns the cosmic time in unit of the Hubble time 1/H_0
     */
    double time(double z);
    /**
     * Returns the equation-of-state exponent ew=3(1+wq) of the cosmological model.
     */
    double quintExp (double a=1.0);
    /**
     * Returns the integral kernel for the cosmic time.
     */
    double timeKernel(double a);
    /**
     * Returns the integral kernel for the angular-diameter distance.
     */
    double angDistKernel(double x);
    /**
     * Returns the Comoving distance between redshifts 0 and z in units of c/H_0.
     */
    double comovDist( double z );
    /**
     * Returns the angular-diameter distance between redshifts z1 and z2 in units of c/H_0.
     */
    double angularDist(double z1, double z2);
    /**
     * Returns critical density in cgs units. If the flag uf is true, the critical density is returned
     * in units of solar masses per cubic Mpc.
     */
    double criticalDensity(double z=0, bool uf=false);
    /**
     * Returns a factor to convert the function time in Gyr
     */
    double CfactorT; // factor to convert time in Gyr
    struct constants
    {
      double kB;            // Boltzman constant erg/K;
      double tcmb;          // cmmb temperature
      double km,kpc,Mpc;
      double yr,Myr,Gyr;
      double Msol,lightspeed;
      double kt;
      double hbar,hc,er,H0,GNewton,GNewton2;
      double rc,nc,rd;
    };
    /**
     * Structure of useful constants for the cosmology class
     */
    constants cn;
};
#endif
