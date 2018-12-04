#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_poly.h>
#include "cosmology.h"
#include "utilities.h"

double gslTimeKernel (double, void*);
double gslQKernel (double, void*);
double gslAngDistKernel (double, void*);

const double del=1.0e-2;

const double tol=1.0e-4;
const double amin=1.0e-5;
const double amax=1.0;
const int nq=64;
double zi[nq],ti[nq];
double st,wt;

cosmology::cosmology(double omega,double lambda,double hubble,double wquint,double baryon)
  :om(omega),la(lambda),h0(hubble),wq(wquint),omb(baryon),cw(true){
  if(wq!=-1){
    std:: cout << " setting wq = " << wq << std:: endl;
    std:: cout << " you are running MOKA with w/=-1 !!! " << std:: endl;
    std:: cout << " if you are aware of eventual problems " << std:: endl;
    std:: cout << " you can comment this line in cosmology.cpp " << std:: endl;
    exit(1);
  }
  ew=3.0*(1.0+wq);
  ri=0;
  initialise();
  st=0;
  wt=0;
}

cosmology:: ~cosmology(){}

void cosmology:: initialise(){
  cn.kB=1.380622e-16;
  cn.tcmb=2.726;
  cn.km=1.e5;
  cn.kpc=3.0856775806e21;
  cn.Mpc=1.e3*cn.kpc;
  cn.yr=365.25*24.0*3600.0;
  cn.Myr=1.e6*cn.yr;
  cn.Gyr=1.e9*cn.yr;
  cn.Msol=1.98892e+33;
  cn.lightspeed=2.99792458e+10;
  cn.kt=cn.kB*cn.tcmb;
  cn.hbar=1.05457266e-27;;
  cn.hc=cn.hbar*cn.lightspeed;
  cn.er=gsl_pow_3( cn.kt/cn.hc );
  cn.H0=100.0*cn.km/cn.Mpc;
  cn.GNewton=6.6732e-8;
  cn.GNewton2=4.3011e-9;
  cn.rc=3.0*gsl_pow_2(cn.H0)/8.0/pi/cn.GNewton;
  cn.nc=21.0/8.0*pow(4.0/11.0, 4.0/3.0);
  cn.rd=pi*pi/15.0*cn.er*cn.kt/cn.rc/gsl_pow_2(cn.lightspeed);
  cn.rd*=(1.0+cn.nc);
  fill_logarithmic (qa.a, nq, amin, amax);
  initQTable();
  qa.q.resize (nq);
  fill_linear(vz,nq,0.,2.);
  double KminMpc=cn.Mpc/1.e5;
  CfactorT=1/(hubble()*100.*cn.Gyr/KminMpc);
  for( int i=0; i<nq; i++ ){
    qa.q[i]=-ew*log( qa.a[i] );
    zi[i]= -1 + pow(10,vz[i]);
    ti[i]=-time(zi[i])*CfactorT;
  }
  splinezt = gsl_spline_alloc (gsl_interp_cspline, nq);
  gsl_spline_init (splinezt, ti, zi, nq);
  /**
   * exponent of the scale factor in the early-time growth factor
   */
  double x0 = aEqual ();
  double e0 = e (x0);
  double s0;
  gsl_poly_solve_quadratic (1.0, 2.0+x0*ePrime (x0)/e0, -1.5*om/x0/x0/x0/e0/e0, &s0, &ng);
}

void cosmology::initQTable()
{
  fill_logarithmic (qa.a, nq, amin, amax);
  qa.q.resize (nq);
  if( cw ){
    for( int i=0; i<nq; i++ )
      qa.q[i]=-ew*log( qa.a[i] );
  }
  else{
    for( int i=0; i<nq; i++ )
      qa.q[i]=gslIntegrateQag( &gslQKernel, 1, qa.a[i] );
  }
}
double cosmology::deltaVF (double zc, bool sf)
{
  const double delv=18*M_PI*M_PI;
  if (cw)
  {
    if (sf){ // use fit formulae suggested by Felix Stoehr
      double af=0.7076;
      double bf=0.4403;
      double om=omega(zc);
      if (la<del){
        af=0.1210;
        bf=0.6756;
      }
      return 0.5*delv*(1+(om-1)*af+pow(om,bf));
    }
    else{
      if (la<del || fabs(om+la-1.0)<del){
        double dv;
        double wc=1.0/omega(zc)-1.0;
        if (la<del) dv=delv*(1.0+0.5584*pow(wc,0.9635));
        else{
          if (fabs(ew)<del) dv=delv*(1.0+0.4093*pow(wc,0.9052));
          else{
            double ac=0.399-1.309*(pow(fabs(wq),0.426)-1.0);
            double bc=0.941-0.205*(pow(fabs(wq),0.938)-1.0);
            dv=delv*(1.0+ac*pow(wc,bc));
          }
        }
        return dv/(wc+1.0);
      }
      else{
        warning ("in deltaV: unsupported cosmology");
        return delv;
      }
    }
  }
  else{
    if (st==0) return 0.0;
    else return 0;
  }
}

double cosmology::deltaC (double zc){
  const double delc=1.6864702;
  if (cw){
    if (la<del || fabs(om+la-1.0)<del){
      double oc=log(omega(zc))/log(10.0);
      if (la<del) return delc*(1.0+0.0406*oc);
      else{
        if (fabs(ew)<del) return delc*(1.0+0.0123*oc);
        else
          return delc*
	    (1.0+(0.353*gsl_pow_4(wq)+1.044*gsl_pow_3(wq)+
		  1.128*gsl_pow_2(wq)+0.555*wq+0.131)*oc);
      }
    }
    else{
      warning("in deltaC: unsupported cosmology");
      return delc;
    }
  }
  else{
    if (st==0) return 0.0;
    else return 0;
  }
}

double cosmology:: tdyn(double z){
  return 2*pow(deltaVF(z)/deltaVF(0.),-0.5)*hubble(0)/hubble(z);
}

double cosmology::zfromt (double t){
  return gsl_spline_eval (splinezt, -t, acc);
}

double cosmology::criticalDensity(double z, bool uf){
  double fh=hubble(z)/h0;
  double fu=1;
  if( uf ){
    fu=gsl_pow_3(cn.Mpc)/cn.Msol;
  }
  return cn.rc*fu*fh*fh;
}

double cosmology:: q(double a){
  if (a<amin || a>amax){
    return exp( gslIntegrateQag( &gslQKernel, 1, a ) );
  }
  else{
    int i=locate(qa.a,a);
    double f=(a-qa.a[i])/(qa.a[i+1]-qa.a[i]);
    return exp( f*qa.q[i+1]+(1-f)*qa.q[i] );
  }
}

double cosmology:: e( double a ){
  return sqrt(ri/gsl_pow_4(a)+om/gsl_pow_3(a)+la*q(a)+(1.0-om-la)/gsl_pow_2(a));
}

double cosmology:: omega( double z ){
  double a=1.0/(1.0+z);
  return om/( gsl_pow_3(a)*gsl_pow_2( e( a ) ) );
}

double cosmology:: omegab( double z ){
  // for now no redshift evolution is assumed!!!
  return omb;
}

double cosmology:: lambda( double z ){
  double a=1.0/(1.0+z);
  return la*q(a)/gsl_pow_2( e( a ) );
}

double cosmology:: hubble( double z){
  return h0*e( 1.0/(1.0+z));
}

double cosmology:: gslIntegrateQag( double (*fc)(double, void*), double a, double b ){
  gsl_function gf;
  gf.function=fc;
  gf.params=this;
  double e,y;
  const size_t n=64;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc (n);
  gsl_integration_qag (&gf, a, b, tiny, tol, n, 1, w, &y, &e);
  gsl_integration_workspace_free (w);
  return y;
}

double cosmology:: aEqual(){
  return cn.rd/om;
}

double cosmology:: angDistKernel( double x ){
  return 1.0/e(1.0/x);
}

double cosmology:: angularDist( double z1, double z2 ){
  double rk=om+la-1.0;
  int k=0;
  if (rk<0.0) k=-1;
  if (rk>0.0) k= 1;
  if (fabs(rk)<=del) k=0;
  rk=sqrt(fabs(rk));

  if (la<=del && 1.0/(1.0+GSL_MAX(z1,z2))>=5.0*aEqual()){
    if (fabs(om-1.0)<=del)
      return 2.0/(1.0+z2)*(1.0/sqrt(1.0+z1)-1.0/sqrt(1.0+z2));
    else
      return 2.0/gsl_pow_2(om)/(1.0+z1)/gsl_pow_2(1.0+z2)*
	(sqrt(1.0+om*z1)*(2.0-om+om*z2)-
	 sqrt(1.0+om*z2)*(2.0-om+om*z1));
  }
  else{
    double d=gslIntegrateQag(&gslAngDistKernel, 1.0+z1, 1.0+z2);
    if (rk*d>del && k==-1) d=sinh(rk*d)/rk;
    if (rk*d>del && k== 1) d=sin(rk*d)/rk;
    return d/(1.0+z2);
  }
}

double cosmology:: comovDist( double z ){

  if (la<=del && 1.0/(1.0+GSL_MAX(0,z))>=5.0*aEqual()){
    if (fabs(om-1.0)<=del)
      return 2.0/(1.0+z)*(1.0-1.0/sqrt(1.0+z));
    else
      return 2.0/gsl_pow_2(om)/gsl_pow_2(1.0+z)*((2.0-om+om*z)-sqrt(1.0+om*z)*(2.0-om));
  }
  else{
    double d=gslIntegrateQag(&gslAngDistKernel, 1.0, 1.0+z);
    return d;
  }
}

void cosmology::dDeltaDa( double x, const double *y, double *dy )
{
  double c[2]={3.0/x+ePrime( x )/e( x ),1.5*omega( 1.0/x-1.0 )/x/x};
  dy[0]=y[1];
  dy[1]=-y[1]*c[0]+y[0]*c[1];
}

int cosDyDx( double x, const double y[], double f[], void *p )
{
  cosmology *co=static_cast<cosmology *>(p);
  co->dDeltaDa( x, y, f );
  return GSL_SUCCESS;
}

double cosmology::growthFactor (double a, bool ff)
{
  double z=1.0/a-1.0;
  if (fabs(ew)<=del && cw && ff)
    {
      double oz=omega( z );
      double lz=lambda( z );
      return 2.5*oz/(pow( oz, 4.0/7.0 )-lz+(1.0+oz/2.0)*(1.0+lz/70.0));
    }
  else
    {
      double x0=aEqual();
      if(a<x0)
	return pow( x0, ng-1 );
      const gsl_odeiv_step_type *T=gsl_odeiv_step_rkck;
      gsl_odeiv_step *s=gsl_odeiv_step_alloc( T, 2 );
      gsl_odeiv_control *c=gsl_odeiv_control_y_new( tiny, 0.0 );
      gsl_odeiv_evolve *e=gsl_odeiv_evolve_alloc( 2 );
      gsl_odeiv_system sys={ cosDyDx, 0, 2, this };
      double h=tiny;
      double y[2]={ x0, ng*pow( x0, ng-1 ) };
      while (x0<a)
	gsl_odeiv_evolve_apply( e, c, s, &sys, &x0, a, &h, y );
      gsl_odeiv_evolve_free( e );
      gsl_odeiv_control_free( c );
      gsl_odeiv_step_free( s );
      return y[0]/a;
    }
}

double cosmology::timeEarly(double a){
  double r=aEqual()/a;
  return 2.0/3.0*a*sqrt(a)/sqrt(om)*((1-2*r)*sqrt(1+r)+2*r*sqrt(r));
}

double cosmology::timeKernel(double a){
  return 1.0/(a*e( a ));
}

double cosmology::ePrime (double a)
{
  // not working for w/=-1!!!
  return
    -0.5*(4.0*ri/gsl_pow_4 (a)+3.0*om/gsl_pow_3 (a)+la*0.+
	  2.0*(1.0-om-la)/gsl_pow_2 (a))/e (a)/a;
  //-0.5*(4.0*ri/gsl_pow_4 (a)+3.0*om/gsl_pow_3 (a)+la*qPrime (a)+
  //	  2.0*(1.0-om-la)/gsl_pow_2 (a))/e (a)/a;
}

double cosmology:: time(double z){
  double a=1.0/(1.0+z);
  double e=5.0*aEqual();
  if (a>=e){
    double t=gslIntegrateQag( &gslTimeKernel, e, a );
    return t+timeEarly( e );
  }
  else
    return timeEarly( a );
}

double cosmology::quintExp (const double a){
  if (cw) return ew;
  else return 0;
}

double gslTimeKernel(double x, void *p){
  cosmology *co=static_cast<cosmology *>(p);
  return co->timeKernel(x);
}

double gslQKernel(double a, void *p){
  cosmology *co=static_cast<cosmology *>(p);
  return -co->quintExp( a )/a;
}

double gslAngDistKernel (double x, void *p){
  cosmology *co=static_cast<cosmology *>(p);
  return co->angDistKernel (x);
}
