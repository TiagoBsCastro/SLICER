#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <CCfits/CCfits>
#include <readSUBFIND.h>
#include <cstring>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include "util.h"

/*****************************************************************************/
/*                                                                           */
/*    This code has been developed to create maps from files of              */
/*      numerical simuations. Up to now it has been optimized to run on      */
/*      "CoDECS-like" simulations and to read gadget1 format files           */
/*                          - it runs on multiple snapshots                  */
/*                          - it reads the SUBFIND and FOF catalogues        */
/*                                                                           */
/*     it returns a list of .fits file of the 2D mass map in each plane      */
/*                          subfindinfield                                   */
/*                          fofinfield                                       */
/*       giving the positions et al. of the fof and subs present in the cone */
/*                                                                           */
/*                                                                           */
/*                              dev. by Carlo Giocoli - cgiocoli@gmail.com   */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*    This code has been adapted to run on Magneticum simulations:           */
/*    other functionalities removed for simplicity                           */
/*                                                                           */
/*      - Proper mass assignment to Hydro particles (TSC-SPH convolution)    */
/*      - Proper mass assignment to  BH   particles                          */
/*	- Bug fixed: Proper construction of light-cones in case of few       */
/*                   snapshots                                               */
/*                                                                           */
/*                 adapted & co-dev. by Tiago Castro tiagobscastro@gmail.com */
/*****************************************************************************/

using namespace std;
using namespace CCfits;

const int bleft = 14;
const double speedcunit = 2.99792458e+3;

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

template <class generic>
void delete_duplicates (vector<generic> &x){

  sort(x.begin(), x.end());
  for(int i = 0; i < x.size() - 1; i++) {
      if (x[i] == x[i + 1]) {
          x.erase(x.begin() + i);
          i--;
      }
  }
  return;
}

double getY(vector<double> x, vector<double> y,double xi){  // Interpolated routine to calculate y(xi)
  int nn = x.size();
  if(x[0]<x[nn-1]){         
    if(xi>x[nn-1]) return y[nn-1];
    if(xi<x[0]) return y[0];
  }  
  else{
    if(xi<x[nn-1]) return y[nn-1];
    if(xi>x[0]) return y[0];  
  }
  int i = locate (x,xi);
  i = std::min (std::max (i,0), int (nn)-2);
  double f=(xi-x[i])/(x[i+1]-x[i]);
  if(i>1 && i<nn-2){
    double a0,a1,a2,a3,f2;                                     
    f2 = f*f;
    a0 = y[i+2] - y[i+1] - y[i-1] + y[i];
    a1 = y[i-1] - y[i] - a0;
    a2 = y[i+1] - y[i-1];
    a3 = y[i];
    return a0*f*f2+a1*f2+a2*f+a3;
  }                                                                      
  else return f*y[i+1]+(1-f)*y[i];           
} 

// weight functions
// --- TSC kernel
float weight (float ixx, float ixh, double dx) {
  float DD = ixx-ixh;
  float x=fabs(DD)/dx;
  float w;    
  
  if(fabs(DD)<=0.5*dx)
    w=3./4.-x*x;
  else if(fabs(DD)>0.5*dx && fabs(DD)<=0.5*3.0*dx)
    w=0.5*((3./2.-x)*(3./2.-x));
  else w=0.;

  return w;

}
// WC6 kernel projected on 2D
double WC6_projected_kernel (double r){

  double w;

  if(r==0)
    w=91./(8*M_PI);
  else if(r<=1 && r>0)
    w=(-91.*(sqrt(1 - pow(r,2))*(-128. + 1344.*pow(r,2) - 7264.*pow(r,4) + 33488.*pow(r,6) + 160290.*pow(r,8) + 37495.*pow(r,10)) + 
       3465.*pow(r,8)*(40. + 24.*pow(r,2) + pow(r,4))*log(r/(1. + sqrt(1. - pow(r,2))))))/(1024.*M_PI);
  else w=0.;

  return w;

}

struct params { double xi; double yi; double xj; double yj; double h_TSC; double h_SPH; };// Params for TSC-SPH conv. 
										          // _i.: Index of the Grid _j.: Index of the part.

double TSC_WC6_conv (double *r, size_t dim, void *p){

  struct params * fp = (struct params *) p;	

  (void)(dim); /* avoid unused parameter warnings */
  double ri[2],rj[2],h[2];

  ri[0]=fp->xi;
  ri[1]=fp->yi;
  rj[0]=fp->xj;
  rj[1]=fp->yj;
  h[0] = fp->h_SPH;
  h[1] = fp->h_TSC;

  double rrj=sqrt( pow(r[0]-rj[0],2)+pow(r[1]-rj[1],2) ); // Distance between the particle and the integrand position

  return WC6_projected_kernel(rrj/h[0])*weight(r[0],ri[0],h[1])*weight(r[1],ri[1],h[1]);

}

// grid points distribution function with != weights
valarray<float> gridist_w (vector<float> x, vector<float> y , vector<float> w, int nn){
  std::valarray<float> grxy( nn*nn );  
  int n0 = x.size();
  if(n0!=y.size()){cout << "x and y positions don't match!" << endl; exit(-1);}
  if(n0!=w.size()){cout << "positions and wheights don't match!" << endl; exit(-1);}

  double dl = 1./double(nn);    
  // --- - - - - - - - - - - - - - - - - - - - - - - ---
  //                               _ _ _ _ _ _
  //  The order of the points is: |_6_|_7_|_8_|
  //                              |_3_|_4_|_5_|
  //                              |_0_|_1_|_2_|
  // coordinate between 0 and 1 and mass particle = 1
  //
  // --- - - - - - - - - - - - - - - - - - - - - - - ---
  for (int i=0; i<n0; i++){
    int grx=int(x[i]/dl); // Grid position integer {0..nn-1}
    int gry=int(y[i]/dl);
  
    int   gridpointx[9], gridpointy[9];
    float posgridx[9], posgridy[9];
    float wfx[9], wfy[9];
    
    gridpointx[0] = grx-1; // Center index x
    gridpointy[0] = gry-1; // Center index y
    
    for (int j=0; j<9; j++){
      gridpointx[j]=gridpointx[0]+(j%3)-1;
      gridpointy[j]=gridpointy[0]+(j/3)-1;
            
      posgridx[j]=(gridpointx[j]+0.5)*dl; //Physical Position corresponding to the center of the bin
      posgridy[j]=(gridpointy[j]+0.5)*dl;

      wfx[j] = sqrt(w[i])*weight(x[i],posgridx[j],dl);
      wfy[j] = sqrt(w[i])*weight(y[i],posgridy[j],dl);

      int grxc = gridpointx[j];
      int gryc = gridpointy[j];
      
      if(grxc>=0 && grxc<nn && gryc>=0 && gryc<nn) 
	grxy[grxc+nn*gryc] = grxy[grxc+nn*gryc] + float(wfx[j])*float(wfy[j]); 
    }
  } 
  return grxy;
}

// grid points distribution function with != weights using hydro part. SPH-TSC kernel convolution
valarray<float> gridist_w (vector<float> x, vector<float> y , vector<float> w, vector<float> h, int nn){

  valarray<float> grxy( nn*nn );
  int n0 = x.size();
  if(n0!=y.size()){cout << "x and y positions don't match!" << endl; exit(-1);}
  if(n0!=w.size()){cout << "positions and weights don't match!" << endl; exit(-1);}
  if(n0!=h.size()){cout << "positions and smoothing lengths don't match!" << endl; exit(-1);}

  double xl[2] = {-1.5,-1.5}; // Limits of integration. TSC_kernel uses the last and the next bins {-1..1} + {-.5,-.5}
  double xu[2] = {1.5,1.5};  // due to the center of the bin 
  double res, err;
  const gsl_rng_type *T;
  gsl_rng *rand;
  struct params p;
  p.h_TSC=1.; // TSC kernel radius
  gsl_monte_function  W = {&TSC_WC6_conv, 2, &p};
  size_t calls = pow(10,3);
  gsl_rng_env_setup ();
  T = gsl_rng_default;
  rand= gsl_rng_alloc (T);
  gsl_monte_plain_state *s = gsl_monte_plain_alloc (2);

  double dl = 1./double(nn);

  for (int i=0; i<n0; i++){

    if(w[i]>1) continue;

    double rx=x[i]/dl; // double precision version of grx
    int grx=int(rx); // x pos on the grid int {0..nn-1}
    double ry=y[i]/dl; // double precision version of gry
    int gry=int(ry); // y pos on the grid int {0..nn-1}
    double r=h[i]/dl; // SPH radius
    int gr=int(r)+1; // Size of the SPH radius in bin's unit

    p.xj=rx;p.yj=ry; p.h_SPH=r;

    for(int ix=grx-gr-1; ix<=grx+gr+1; ix++) // All the possible bins in x-axis englobed by r (SPH kernel)

      for(int iy=gry-gr-1; iy<=gry+gr+1; iy++){ // All the possible bins in y-axis englobed by r (SPH kernel)

        if(ix<nn && iy<nn && ix>0 && iy>0){

         if(r<10 && r>.1){
            p.xi=ix+.5; p.yi=iy+.5;
            gsl_monte_plain_integrate (&W, xl, xu, 2, calls, rand, s,&res, &err);
            grxy[ix+nn*(iy)] += w[i]*res;
         }else if(r>10){

           for(int ii=-1;ii<=1;ii++)for(int jj=-1;jj<=1;jj++){ // Making a binary integration over the TSC neighbours

             float xi=ix+ii+.5; float yi=iy+jj+.5;
             double rij=sqrt( pow(rx-xi,2) + pow(ry-yi,2) );
             grxy[ix+nn*iy] += w[i]*WC6_projected_kernel( rij/r )*((4.-3.*fabs(ii))/6)*((4.-3.*fabs(jj))/6);

            }

         }else{

           for(int ii=-1;ii<=1;ii++)for(int jj=-1;jj<=1;jj++){ // Using SPH particles as point particles

             float xi=ix+ii+.5; float yi=iy+jj+.5;
             double rij=sqrt( pow(rx-xi,2) + pow(ry-yi,2) );
             grxy[ix+ii+nn*(iy+jj)] += w[i]*weight(rx,xi,p.h_TSC)*weight(ry,yi,p.h_TSC);     
       
           }
        }
      }
    }
  }

  gsl_monte_plain_free (s);
  gsl_rng_free (rand);
  return grxy;
}
// grid points distribution function with != weights using hydro part. SPH-TSC kernel convolution
valarray<float> gridist_w (vector<float> x, vector<float> y , vector<float> w, vector<float> h, int nn, interp_kernel *interp_ptr){

  valarray<float> grxy( nn*nn );
  int n0 = x.size();
  if(n0!=y.size()){cout << "x and y positions don't match!" << endl; exit(-1);}
  if(n0!=w.size()){cout << "positions and weights don't match!" << endl; exit(-1);}
  if(n0!=h.size()){cout << "positions and smoothing lengths don't match!" << endl; exit(-1);}

  double h_TSC=1.; // TSC kernel radius

  double dl = 1./double(nn);

  cout<< "h max: " <<*max_element(h.begin(),h.end())/dl << endl;
  cout<< "h min: " <<*min_element(h.begin(),h.end())/dl << endl;
  cout<< "x max: " <<*max_element(x.begin(),x.end()) << endl;
  cout<< "x min: " <<*min_element(x.begin(),x.end()) << endl;
  cout<< "y max: " <<*max_element(y.begin(),y.end()) << endl;
  cout<< "y min: " <<*min_element(y.begin(),y.end()) << endl;

  for (int i=0; i<n0; i++){

    if(w[i]>1) continue;

    double rx=x[i]/dl; // double precision version of grx
    int grx=int(rx); // x pos on the grid int {0..nn-1}
    double ry=y[i]/dl; // double precision version of gry
    int gry=int(ry); // y pos on the grid int {0..nn-1}
    double r=h[i]/dl; // SPH radius
    int gr=int(r)+1; // Size of the SPH radius in bin's unit

    for(int ix=grx-gr-1; ix<=grx+gr+1; ix++) // All the possible bins in x-axis englobed by r (SPH kernel)

      for(int iy=gry-gr-1; iy<=gry+gr+1; iy++){ // All the possible bins in y-axis englobed by r (SPH kernel)

        if(ix<nn && iy<nn && ix>0 && iy>0){

          if(r<10 && r>.1){

            double xi=ix+.5, yi=iy+.5;
            double rij=sqrt(pow(xi-rx,2)+pow(yi-ry,2));
            grxy[ix+nn*(iy)] += w[i]*gsl_spline2d_eval(interp_ptr->spline, rij, r, interp_ptr->racc, interp_ptr->hacc);

          }else if(r>10){

            for(int ii=-1;ii<=1;ii++)for(int jj=-1;jj<=1;jj++){ // Making a binary integration over the TSC neighbours

              double xi=ix+ii+.5; double yi=iy+jj+.5;
              double rij=sqrt( pow(rx-xi,2) + pow(ry-yi,2) );
              if(ix+ii<nn && iy+jj<nn) grxy[ix+ii+nn*(iy+jj)] += w[i]*WC6_projected_kernel( rij/r )*((4.-3.*fabs(ii))/6)*((4.-3.*fabs(jj))/6);

            }

          }else{

            for(int ii=-1;ii<=1;ii++)for(int jj=-1;jj<=1;jj++){ // Using SPH particles as point particles

              float xi=ix+ii+.5; float yi=iy+jj+.5;
              double rij=sqrt( pow(rx-xi,2) + pow(ry-yi,2) );
              if(ix+ii<nn && iy+jj<nn) grxy[ix+ii+nn*(iy+jj)] += w[i]*weight(rx,xi,h_TSC)*weight(ry,yi,h_TSC);     
       
            }
         }
      }
    }
  }

  return grxy;
}

struct DATA
{
  // long substituted with int32_t
  int32_t npart[6];	
  double massarr[6];
  double time;
  double redshift;
  int32_t flag_sfr;
  int32_t flag_feedback;
  int32_t npartTotal[6];	
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
  int32_t la[bleft];
};

struct BLOCK
{
  int32_t blocksize1;
  int8_t alignment[4];
  char name[4];
  int8_t padding[8];
  int32_t blocksize2;
};

istream & operator>>(istream &input, DATA &Data) 
{ 
  input.read((char *)&Data, sizeof(Data));
  return input;
};

istream & operator>>(istream &input, BLOCK &block) 
{ 
  input.read((char *)&block, sizeof(block));
  return input;
};

void readParameters(string file_name,int *npix, double *boxl,
		    double *zs, double *fov, string *filredshiftlist,
		    string *filsnaplist, string *pathsnap,
		    string *idc, 
		    int *seedcenter, int *seedface, int *seedsign,
		    string *subfiles, string *simulation, int *nfiles, 
		    string *partinplanes, int *noSNAP, double *bufferdeg,
                    string *directory,string *suffix,string *kernel_file){ 

  ifstream ifilin;
  ifilin.open(file_name.c_str());
  if(ifilin.is_open()){

    string butstr;
  
    ifilin >> butstr; // number of pixels
    ifilin >> *npix;
    ifilin >> butstr; // boxl
    ifilin >> *boxl;
    ifilin >> butstr; // source redshift
    ifilin >> *zs;
    ifilin >> butstr; // field of view in degrees
    ifilin >> *fov;
    ifilin >> butstr; // file with the redshift list it may contain three columns: snap 1/(1+z) z
    ifilin >> *filredshiftlist;
    ifilin >> butstr; // path where the snaphosts are located
    ifilin >> *pathsnap;
    ifilin >> butstr; // simulation name (prefix infront at the snap file)
    ifilin >> *simulation;
    ifilin >> butstr;  // number of files per snapshot
    ifilin >> *nfiles;
    ifilin >> butstr; // path and file name of the comoving distance file (if not available you may use CosmoLib)
    ifilin >> *idc;
    ifilin >> butstr; // path of SUBFIND if NO I will not look for them
    ifilin >> *subfiles;
    ifilin >> butstr; // seed for the random location of the center
    ifilin >> *seedcenter;
    ifilin >> butstr; // seed for the random selection of the dice face
    ifilin >> *seedface;
    ifilin >> butstr; // seed for the selection of the sign of the coordinates
    ifilin >> *seedsign;  
    ifilin >> butstr; // which particles in the planes (ALL: one file for all, or NO: one for each)
    ifilin >> *partinplanes;  
    ifilin >> butstr; // if 0 read also snaphosts if =/0 only subs and fof area read
    ifilin >> *noSNAP;  
    ifilin >> butstr; // size of the buffer region in degrees for subs and fof
    ifilin >> *bufferdeg;
    ifilin >> butstr;//directory to save FITS files  
    ifilin >> *directory;
    ifilin >> butstr;//Sufix to save FITS files  
    ifilin >> *suffix;
    ifilin >> butstr;//TSC-SPH interpolation data
    ifilin >> *kernel_file;
  }
  else exit(-1);
}; 

void getPolar(double x, double y, double z, double *ra, double *dec, double *d){
  *d = sqrt(x*x+y*y+z*z);
  *dec = asin(x/(*d));
  *ra = atan2(y,z);
}

int main(int argc, char** argv){

  if(argc!=2){cout << "No params!! Nothing to be done!" << endl; return -1;} // Should be run ./exe paramfile

  string inifile=argv[1]; // Paramfile name

  cout << "   ------------------------------------------------------ " << endl;
  cout << "   -                                                    - " << endl;
  cout << "   -           2D Mapping Simulation Snapshot           - " << endl;
  cout << "   -                                                    - " << endl;
  cout << "   -               collapsing one dimension             - " << endl;
  cout << "   ------------------------------------------------------ " << endl;

  // ... masses of the different type of particles
  double m0,m1,m2,m3,m4,m5;
  // ******************** to be read in the INPUT file ********************
  // ... project - number of pixels
  int npix;
  // ... set by hand the redshift of the source and than loop only up to 
  // ... the needed snaphost when creating the light cone!
  double zs,Ds;
  string filredshiftlist,filsnaplist;
  // ... loop on different snapshots
  string pathsnap; // = "/dati1/cgiocoli/CoDECS/";
  double boxl; // Mpc/h
  string idc; // comoving distance file
  int seedcenter, seedface, seedsign,nfiles;
  double fov;
  string subfiles;
  string simulation; // Simulation Name  
  string partinplanes; // ALL all part in one plane, NO each part type in different planes
  int noSNAP; 
  double bufferdeg;
  string directory, suffix,kernel_file;

  string FilePath;
   
  readParameters(inifile,&npix,&boxl,&zs,&fov,
		 &filredshiftlist,&filsnaplist,
		 &pathsnap,&idc,&seedcenter,
		 &seedface,&seedsign,
		 &subfiles,&simulation,&nfiles,
		 &partinplanes,
		 &noSNAP,&bufferdeg,&directory,&suffix,&kernel_file);  // Reading the parameters File

  int nplanes,ncubes;
  vector<int> isnap,ibox;
  vector<double> iDfinal;
  vector<string> Reflection1,Axes1;
  vector<string> Reflection2,Axes2;
  vector<string> Reflection3,Axes3;
  vector<double> xM0,yM0,zM0;
  interp_kernel *interp_ptr;

  if(kernel_file.c_str()!="None"){

    ifstream interp_tab(kernel_file.c_str(), ios::in|ios::binary|ios::ate);
    unsigned int size=sizeof(float),size_file;
    vector<double> r;
    vector<double> h;
    vector<double> w;

    size_file = interp_tab.tellg()/(3*sizeof(float));
    cout << "Total size of the Interpolation data: " <<size_file << endl;
    interp_tab.seekg (0, ios::beg);

    if (interp_tab.is_open()){

      for(int i=0;i<size_file;i++){

        float ri,hi,wi;
        interp_tab.read(reinterpret_cast<char*>(&ri), sizeof(float));r.push_back(ri);
        interp_tab.read(reinterpret_cast<char*>(&hi), sizeof(float));h.push_back(hi);
        interp_tab.read(reinterpret_cast<char*>(&wi), sizeof(float));w.push_back(wi);

      }

    interp_tab.close();

    }

    delete_duplicates(r);
    delete_duplicates(h);

    cout << "r array length: "  << r.size() << endl;
    cout << "h array length: "  << h.size() << endl;

    interp_ptr = new interp_kernel ( &r[0], &h[0], &w[0], r.size(), h.size());
  }else{

    delete interp_ptr;

  }
  
  double Omega_matter,Omega_lambda,Omega_baryon,hubble;
  string snpix = conv(npix,fINT);

  // ... read the redshift list and the snap_available
  // ... to build up the light-cone
  ifstream redlist;
  redlist.open(filredshiftlist.c_str());

  vector <double> dtsnaplist;
  vector <double> tredlist;
  vector <int> lsnap; //lsnap and lred are the selected snapshots selected from dtsnaplist and tredlist (z<zs)
  vector<double> lred;

  if(redlist.is_open()){
    int buta;
    double butb,butc;
    while(redlist >> buta >> butb >> butc){
      tredlist.push_back(butc);
      dtsnaplist.push_back(buta);
      if(butc<zs){
	lsnap.insert(lsnap.begin(),buta);
	lred.insert(lred.begin(),butc);
      }

    }

  }else{
    cout << " redshift list file redshift_list.txt does not " << endl;
    cout << " exist in the Code dir ... check this out      " << endl;
    cout << "    I will STOP here !!! " << endl;
    exit(1);
  }

  int nsnaps = lsnap.size();

  cout << "  " << endl;
  cout << " opening path for snapshots: " << endl;
  cout << pathsnap << endl;
  cout << " " << endl;
  cout << " I will look for comoving distance " << endl;
  cout << "      file = " << idc << endl;
  cout << " " << endl;
  
  ifstream infiledc;
  vector<double> zl, dl;
  infiledc.open(idc.c_str());
  if(infiledc.is_open()){
    double zi,dli;
    while(infiledc >> zi >> dli){
      zl.push_back(zi);
      dl.push_back(dli*speedcunit);
      cout << zi << "  " << dli*speedcunit << endl;
    }
    infiledc.close();
  }
  else{
    cout << "  " << endl;
    cout << " the comoving distance file: " << idc << endl;
    cout << " does not exists " << endl;
    cout << " I will STOP here!!! " << endl;
    exit(1);
  }

  if(zs>zl[zl.size()-1]){
    cout << " source redshift larger than the highest available redshift in the comoving distance file " << endl;
    cout << "  that is = " << zl[zl.size()-1] << endl;
    cout << " I will STOP here !!! " << endl;
    exit(1);
  }

  int nreplications;

  // cout << " Ds = " << Ds << endl;
  Ds = getY(zl,dl,zs);  // comoving distance of the last plane
  // cout << " Ds = " << Ds << endl;
  // exit(1);
  nreplications = int(Ds/boxl)+1;

  cout << "  nreplications " << nreplications << "  " << Ds << "  " << std:: endl;
  vector<int> replication;
  vector<int> fromsnap;
  vector<double> lD;
  vector<double> lD2;
  
  for(int j=0;j<nreplications;j++){
    for(int i=0;i<nsnaps;i++){
      double ldbut = getY(zl,dl,lred[i]); // Getting dl(lred[i])
      if(ldbut>=j*boxl && ldbut<=(j+1.)*boxl){
	std:: cout << " simulation snapshots = " << ldbut << "  " << lred[i] << "  " << j << " from snap " 
		     << lsnap[i] << "  " << boxl*j << "  " << boxl*(j+1) 
		     << std:: endl; // j[i]-j[i-1] is the number of repetitions needed for each snapshot!
	replication.push_back(j);
	fromsnap.push_back(lsnap[i]);
	lD.push_back(ldbut);
      }
    }
  }

  cout << "  " << endl;
  cout << "  reorganazing the planes " << endl;
  cout << " " << endl;

  for(int i=0;i<lD.size();i++){
    if(i<(lD.size()-1)) lD2.push_back(lD[i+1]); // lD2[i] is initialized as lD[i-1] if i<imax; if i=imax lD2 is the last plane distance
    else lD2.push_back(Ds);
  }

  for(int i=0;i<lD.size();i++){
    for(int k=1;k<=512;k++){
      if(lD[i]<double(k)*boxl && lD2[i]>double(k)*boxl){  // If lD[i] is less than k*box_size AND lD[i+1] is higher than k*box_size 
	  lD[i+1] = lD2[i]; // Redundant??
	  lD2[i]=double(k)*boxl; // Assigning lD2 (source plane??) position to the end of the box
      }
    }
    if(lD[i]<513*boxl && lD2[i]>513*boxl){
      cout << " exiting ... increase the number of replications by hand in the file it is now 512 !!!! " << endl;
      exit(1);
    }
  }

  std:: cout << " Comoving Distance of the last plane " << Ds << std:: endl;
  std:: cout << "  " << endl;
  vector<double> zsimlens(nsnaps);
  cout << " nsnaps = " << nsnaps << endl;
  
  for(int i=0;i<nsnaps;i++){
    if(i<nsnaps-1){
      if(lD[i+1]-lD2[i]>boxl*1e-9){
	  fromsnap[i] = -fromsnap[i];
      }
    }
    double dlbut = (lD[i] + lD2[i])*0.5;    // dlbut is the half distance between lD and lD2!

    zsimlens[i] = getY(dl,zl,dlbut); // Getting the lens plane redshift

    std:: cout << zsimlens[i] << " planes = " << lD[i] << "  " << lD2[i] << "  " << replication[i] 
		 << " from snap " << fromsnap[i] << std:: endl;
  } 

  vector<double> bfromsnap,blD,blD2,bzsimlens,blred;
  vector<int> breplication, blsnap;
  int pl=0;
  vector<int> pll;

  for(int i=0;i<nsnaps;i++){
    bfromsnap.push_back(fabs(fromsnap[i]));
    blD.push_back(lD[i]);
    blD2.push_back(lD2[i]);
    bzsimlens.push_back(zsimlens[i]);
    breplication.push_back(replication[i]);
    blsnap.push_back(lsnap[i]);
    blred.push_back(lred[i]);
    if(fromsnap[i]<=0){// Duplicate the number of planes of the snaps marked with negative signal
      bfromsnap.push_back(-fromsnap[i]);
      blD.push_back(lD2[i]);
      blD2.push_back(lD[i+1]);
      double dlbut = (lD[i+1] + lD2[i])*0.5;    
      // half distance between the two!
      bzsimlens.push_back(getY(dl,zl,dlbut));
      breplication.push_back(replication[i+1]);
      blsnap.push_back(lsnap[i]);
      blred.push_back(lred[i]);
    }   
  }

  nsnaps = bfromsnap.size(); // updating nsnaps. 

  for(int i=0;i<blD.size();i++){
    switch(i){
      case 0:{
        if(blD[i]>boxl){// Duplicate the first planes in case begins at dl bigger than Box_size
          nsnaps+=pow(2,int(log2((blD[i])/boxl))+1)-1; // Updating nsnaps acording with the number of replications needed     
        }
      }
      default:{
        if(blD2[i]-blD[i]>boxl){//Testing if there's no role in light-cone
          nsnaps+=pow(2,int(log2((blD2[i]-blD[i])/boxl))+1)-1; // Updating nsnaps acording with the number of replications needed 
        }    
      }
    }
  }

  for(int i=0;i<nsnaps;i++){

    auto pos_z=bzsimlens.begin();
    auto pos_rep=breplication.begin();
    auto pos_lsnap=blsnap.begin();
    auto pos_lred=blred.begin();
    auto pos_fromsnap=bfromsnap.begin();
    auto pos_lD=blD.begin();
    auto pos_lD2=blD2.begin();

    switch(i){
      case 0:{
        if(blD[i]>boxl){// Duplicate the first planes in case begins at dl bigger than Box_size

          blD.insert(pos_lD+i,blD[i]/2);
          blD2.insert(pos_lD2+i,blD[i+1]);

          double dlbut = (blD[i] + blD2[i])*0.5; // half distance between the two!
          bzsimlens.insert(pos_z+i,getY(dl,zl,dlbut));
          breplication.insert(pos_rep+i,replication[i]);
          breplication[i]=blD2[i]/boxl;
          bfromsnap.insert(pos_fromsnap+i,bfromsnap[i]);
          blsnap.insert(pos_lsnap+i,blsnap[i]);
          blred.insert(pos_lred+i,blred[i]);
          i--; // Run over the same plane again to test if it is requered another repetition
          break;
        } 
      }
      default:{
        if(blD2[i]-blD[i]>boxl){// Duplicate the number of planes with separation bigger than Box_size

          blD.insert(pos_lD+i,blD[i]);
          blD[i+1]=0.5*(blD[i]+blD2[i]); // Update  blD
          blD2.insert(pos_lD2+i,blD2[i]);
          blD2[i]=blD[i+1]; // Update  blD

          double dlbut = (blD[i] + blD2[i])*0.5; // half distance between the two!
          bzsimlens.insert(pos_z+i,getY(dl,zl,dlbut));
          breplication.insert(pos_rep+i,replication[i]);
          breplication[i]=blD2[i]/boxl;
          dlbut = (blD[i+1] + blD2[i+1])*0.5; // half distance between the two!
          bzsimlens[i+1]=getY(dl,zl,dlbut);
          bfromsnap.insert(pos_fromsnap+i,bfromsnap[i]);
          blsnap.insert(pos_lsnap+i,blsnap[i]);
          blred.insert(pos_lred+i,blred[i]);
          i--; // Run over the same plane again to test if it is requered another repetition
          break;
        } 
      }  
    }
  }

  cout << "  " << endl;
  cout << "  re-reorganazing the planes " << std:: endl;
  cout << " " << endl; 
  
  ofstream planelist;
  string planes_list;
  planes_list = directory+"planes_list_"+suffix+".txt";

  planelist.open(planes_list.c_str());

  for(int i=0;i<nsnaps;i++){
    std:: cout << bzsimlens[i] << " bplanes = " << blD[i] << "  " << blD2[i] << "  " << breplication[i] << " from snap " << bfromsnap[i] << std:: endl;
    pl++;
    planelist <<  pl << "   " << bzsimlens[i] << "   " << blD[i] << "   " << blD2[i] << "   " << breplication[i] << "   " << bfromsnap[i] << "   " << blred[i] << std:: endl;
    pll.push_back(pl);
  }
    
  if(blD2[nsnaps-1]<Ds){
    // we need to add one more snaphost
    double s;
    for(int i=dtsnaplist.size()-1;i>0;i--){
      s=dtsnaplist[i];
      if(s<bfromsnap[nsnaps-1]){
	bfromsnap.push_back(s);
	double dlbut = (Ds+blD2[nsnaps-1])*0.5;
	double zbut = getY(dl,zl,dlbut);
	blD.push_back(blD2[nsnaps-1]);
	blD2.push_back(Ds);
	bzsimlens.push_back(zbut);
	breplication.push_back(breplication[nsnaps-1]+1);
	double zn = getY(dtsnaplist,tredlist,s);
	blsnap.push_back(int(s));
	blred.push_back(zn);
	nsnaps++;
	break;
      }
    }
    pl++;
    std:: cout << bzsimlens[nsnaps-1] << " bplanes = " << blD[nsnaps-1] << "  " << blD2[nsnaps-1] << "  " << breplication[nsnaps-1] << " from snap " << bfromsnap[nsnaps-1] << std:: endl;	
    planelist << pl << "   " << bzsimlens[nsnaps-1] << "   " << blD[nsnaps-1] << "   " << blD2[nsnaps-1] << "   " << breplication[nsnaps-1] << "   " << bfromsnap[nsnaps-1] << "   " << blred[nsnaps-1] << std:: endl;
    pll.push_back(pl);	 
  } 

  planelist.close();

  std:: cout << "  " << endl; 
  // randomization of the box realizations :
  int nrandom = breplication[nsnaps-1]+1;
  vector<double> x0(nrandom), y0(nrandom), z0(nrandom); // ramdomizing the center of the simulation [0,1]
  vector<int> face(nrandom); // face of the dice
  vector<int> sgnX(nrandom), sgnY(nrandom),sgnZ(nrandom);
  
  for(int i=0;i<nrandom;i++){
    srand(seedcenter+i*13);
    x0[i] = rand() / float(RAND_MAX);
    y0[i] = rand() / float(RAND_MAX);
    z0[i] = rand() / float(RAND_MAX);
    cout << "  " << endl;
    cout << " random centers  for the box " << i << " = " << x0[i] << "  " << y0[i] << "  " << z0[i] << endl;
    face[i] = 7;
    srand(seedface+i*5);
    while(face[i]>6 || face[i]<1) face[i] = int(1+rand() / float(RAND_MAX)*5.+0.5);
    std:: cout << " face of the dice " << face[i] << std:: endl;
    sgnX[i] = 2;
    srand(seedsign+i*8);
    while(sgnX[i] > 1 || sgnX[i] < 0) sgnX[i] = int(rand() / float(RAND_MAX)+0.5);
    sgnY[i] = 2;
    while(sgnY[i] > 1 || sgnY[i] < 0) sgnY[i] = int(rand() / float(RAND_MAX)+0.5);
    sgnZ[i] = 2;
    while(sgnZ[i] > 1 || sgnZ[i] < 0) sgnZ[i] = int(rand() / float(RAND_MAX)+0.5);
    if(sgnX[i]==0) sgnX[i]=-1;
    if(sgnY[i]==0) sgnY[i]=-1;
    if(sgnZ[i]==0) sgnZ[i]=-1;
    std:: cout << " signs of the coordinates = " << sgnX[i] << "  " << sgnY[i] << " " << sgnZ[i] << endl;
  }

  cout << "  " << endl;
  cout << "  " << endl;
  cout << " set the field of view to be square in degrees " << endl;
  double h0,fovradiants,bufferrad;
  double om0, omL0;
  fovradiants = fov/180.*M_PI;
  bufferrad = bufferdeg/180.*M_PI;
  // check if the field of view is too large with respect to the box size
  if(fovradiants*Ds>boxl){
    std:: cout << " field view too large ... I will STOP here!!! " << std:: endl;
    std:: cout << " value set is = " << fov << std:: endl;
    std:: cout << " maximum value allowed " << boxl/Ds*180./M_PI << " in degrees " << std:: endl;
    exit(1);
  }
  // check if the fov + buffer region for subs and fof is too large than the field of view
  if((fovradiants+2.*bufferrad)*Ds>boxl){
    std:: cout << " field view + 2 x buffer region is too large ... I will STOP here!!! " << std:: endl;
    std:: cout << " value set is = " << fov + 2*bufferdeg << std:: endl;
    std:: cout << " maximum value allowed " << boxl/Ds*180./M_PI << " (fov+2 x buffer) in degrees " << std:: endl;
    exit(1);
  }
  // loop on snapshots ****************************************
  cout << " now loop on " << nsnaps << " snapshots " << endl;
  cout << "  " << endl;
  for(int nsnap=0;nsnap<nsnaps;nsnap++){
    if(blD2[nsnap]-blD[nsnap] < 0){
      cout << " comoving distance of the starting point " << blD[nsnap] << endl;
      cout << " comoving distance of the final    point " << blD2[nsnap] << endl;
      cout << " please check this out! I will STOP here!!! " << endl;
      exit(1);
    }
    int rcase = breplication[nsnap];
    std::valarray<float> mapxytot( npix*npix );
    // type 0
    int ntotxy0;
    ntotxy0 = 0;
    std::valarray<float> mapxytot0( npix*npix );
    // type 1
    int ntotxy1;
    ntotxy1 = 0;
    std::valarray<float> mapxytot1( npix*npix );
    // type 2
    int ntotxy2;
    ntotxy2 = 0;
    std::valarray<float> mapxytot2( npix*npix );
    // type 3
    int ntotxy3;
    ntotxy3 = 0;
    std::valarray<float> mapxytot3( npix*npix );
    // type 4
    int ntotxy4;
    ntotxy4 = 0;
    std::valarray<float> mapxytot4( npix*npix );
    // type 5
    int ntotxy5;
    ntotxy5 = 0;
    std::valarray<float> mapxytot5( npix*npix );
  
    string snapnum;
    if( blsnap[nsnap]<10) snapnum = "00"+conv(blsnap[nsnap],fINT);
    else if(blsnap[nsnap]>=10 && blsnap[nsnap]<100 ) snapnum = "0"+conv(blsnap[nsnap],fINT);
    else snapnum = conv(blsnap[nsnap],fINT); 
    string File = pathsnap+"/snapdir_"+snapnum+"/"+simulation+"_"+snapnum;
    string snappl;
    if( pll[nsnap]<10) snappl = "00"+conv(pll[nsnap],fINT);
    else if(pll[nsnap]>=10 && pll[nsnap]<100 ) snappl = "0"+conv(pll[nsnap],fINT);
    else snappl = conv(pll[nsnap],fINT);     
    string Filesub;

    if(ifstream(directory+simulation+"."+snappl+".plane_"+snpix+"_"+suffix+".fits") && partinplanes == "ALL"){

      cout << directory+simulation+"."+snappl+".plane_"+snpix+"_"+suffix+".fits" << " "<< "Already exists" <<endl;
      continue; // If this lens plane was already created go to the next one

    }
  
    float num_float1, num_float2, num_float3; 

    // redshift and dl of the simulation
    double zsim, dlsim;
    
    for (unsigned int ff=0; ff<nfiles; ff++){ 
      // map for each mass type
      std::valarray<float> mapxy0( npix*npix );
      std::valarray<float> mapxy1( npix*npix );
      std::valarray<float> mapxy2( npix*npix );
      std::valarray<float> mapxy3( npix*npix );
      std::valarray<float> mapxy4( npix*npix );
      std::valarray<float> mapxy5( npix*npix );

      // GADGET has 6 different particle type
      vector<float> xx0(0), yy0(0), zz0(0),smoothing_l(0);
      vector<float> xx1(0), yy1(0), zz1(0);
      vector<float> xx2(0), yy2(0), zz2(0);
      vector<float> xx3(0), yy3(0), zz3(0);
      vector<float> xx4(0), yy4(0), zz4(0);
      vector<float> xx5(0), yy5(0), zz5(0);
      
      string file_in = File+"."+conv(ff,fINT);
      
      ifstream fin(file_in.c_str());
      if (!fin) {cerr <<"Error in opening the file: "<<file_in<<"!\n\a"; exit(1);}
      	
      cout <<"    reading the input file: "<<file_in<<endl;
      
      int32_t blockheader[5];
      fin.read((char *)&blockheader, sizeof(blockheader));
      DATA data; fin >> data;

      if(ff==0){

      	cout << "Printing Header Data" << endl;

	cout << "N. part.: " << data.npart[0] << " "<< data.npart[1] << " "<< data.npart[2] << " "<< data.npart[3] << " "<< data.npart[4] << " "<< data.npart[5]  << endl;
        cout << "Mass Array: "<< data.massarr[0]<< " "<< data.massarr[1]<< " "<< data.massarr[2]<< " "<< data.massarr[3]<< " "<< data.massarr[4]<< " "<< data.massarr[5]<< " "<< endl;
        cout << "Time: "<< data.time<< endl;
        cout << "Z: "<< data.redshift<< endl;
        cout << "Flag SFR.: "<< data.flag_sfr<< endl;
        cout << "Flag Feedback: "<< data.flag_feedback<< endl;
        cout << "N. tot.: "<< data.npartTotal[0]<<" "<< data.npartTotal[1]<<" "<< data.npartTotal[2]<<" "<< data.npartTotal[3]<<" "<< data.npartTotal[4]<<" "<< data.npartTotal[5]<<" "<< endl;
        cout << "Flag cooling: "<< data.flag_cooling<< endl;
        cout << "N. files: "<< data.numfiles<< endl;
        cout << "Box size: "<< data.boxsize<< endl;
        cout << "Omega_matter: "<< data.om0<< endl;
        cout << "Omega_DE: "<< data.oml<< endl;
        cout << "h: "<< data.h<< endl;
        cout << "Flag sage: "<< data.flag_sage<< endl;
        cout << "Flag metals: "<< data.flag_metals<< endl;
        cout << "N. tot HW: "<< data.nTotalHW[0]<<" "<< data.nTotalHW[1]<<" "<< data.nTotalHW[2]<<" "<< data.nTotalHW[3]<<" "<< data.nTotalHW[4]<<" "<< data.nTotalHW[5]<<" "<< endl;
        cout << "Flag entropy: "<< data.flag_entropy<<endl;
      }

      if(data.nTotalHW[0]+data.nTotalHW[1]+data.nTotalHW[2]+data.nTotalHW[3]+data.nTotalHW[4]+data.nTotalHW[5]!=0){

	cout << "More than 2^32 particles!! The code has to be changed! Exiting now!" << endl;
	return -1;
      
      }

      BLOCK block; fin >> block;
      cout << "Size of Header is " << sizeof(data) << endl;
      cout << "Should be         " << block.blocksize1 << endl;

      cout << "Fast Fowarding next block. Name: ";
      cout << block.name[0];
      cout << block.name[1];
      cout << block.name[2];
      cout << block.name[3] << endl;

      int8_t blockheader2[block.blocksize2];  // Magneticum
      fin.read((char *)&blockheader2, sizeof(blockheader2));

      BLOCK block2; fin >> block2;
      cout << "reading next block. Name: ";
      cout << block2.name[0];
      cout << block2.name[1];
      cout << block2.name[2];
      cout << block2.name[3] << endl;
      cout << "Should be                 POS " << endl;
      
      // total number of particles as the sum of all of them
      int dim = data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]+data.npart[4]+data.npart[5];
      int dimmass0=0;
      int hydro=0;

      for(int i=0;i<=5;i++){
	
	if(data.massarr[i]==0){dimmass0+=data.npart[i];}

      }

      if(dimmass0==0){
	
	cout << "		@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
	cout << "		@			   @" << endl;	
	cout << "		@  !!DM only simulation!!  @" << endl;
	cout << "		@                          @" << endl;
	cout << "		@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

      }
      else{
	
	cout << "		@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
	cout << "		@                          @" << endl;	
	cout << "		@  !!Hydro   simulation!!  @" << endl;
	cout << "		@                          @" << endl;
	cout << "		@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
	hydro=1;
      }


      cout << " .......................................................... " << endl;
      cout << "   number of particles in this snapshot: " << endl;
      cout << data.npart[0] << " " << data.npart[1] << " " << data.npart[2] 
	   << " " << data.npart[3] << " " << data.npart[4] << " " << data.npart[5] << endl;
      
      if(ff==0){
	// compute the comoving angular diameter distance at simulation redshift
	zsim = data.redshift;
	dlsim = getY(zl,dl,zsim);
	h0 = data.h;
	cout << "  " << endl;
	cout << "      __________________ COSMOLOGY __________________  " << endl;
	cout << " " << endl;
	om0 = data.om0;
	omL0 = data.oml;
	cout << "      Omegam = " << data.om0 << " " << "Omegal = " << data.oml << endl;
	cout << "           h = " << data.h   << " " << "BoxSize = " << data.boxsize << endl;
	cout << "      redshift = " << zsim <<   " " << "Dl (comoving) = " << dlsim << endl;
	if(abs(boxl - data.boxsize/1.e+3)>1.e-2 ){
	  cout << " set boxl and data.size differ ... check it! " << std:: endl;
	  cout << "  boxl = " << boxl << "  " << " data.boxsize = " << data.boxsize/1.e+3 << endl;
	  exit(1);
	}
	
	cout << "      _______________________________________________  " << endl;
	cout << " " << endl;
	cout << "   total number of partiles in the simulation: " << endl; 
	cout << data.npartTotal[0] << " " << data.npartTotal[1] << " " << data.npartTotal[2] 
	     << " " << data.npartTotal[3] << " " << data.npartTotal[4] << " " << data.npartTotal[5] << endl;
	cout << " " << endl;
	cout << "   xparticle type mass array: " << endl; 
	cout << data.massarr[0] << " " << data.massarr[1] << " " << data.massarr[2] 
	     << " " << data.massarr[3] << " " <<  data.massarr[4] << " " <<  data.massarr[5] << endl; 
	m0 = data.massarr[0];
	m1 = data.massarr[1];
	m2 = data.massarr[2];
	m3 = data.massarr[3];
	m4 = data.massarr[4];
	m5 = data.massarr[5];
      }
      if(noSNAP==0){
	// type0
	for (int pp=0; pp<data.npart[0]; pp++){
	  fin.read((char *)&num_float1, sizeof(num_float1)); 
	  fin.read((char *)&num_float2, sizeof(num_float2)); 
	  fin.read((char *)&num_float3, sizeof(num_float3)); 
 
	  float x, y, z;

	  float xb, yb, zb;

	  xb = sgnX[rcase]*(((num_float1)/data.boxsize));
	  yb = sgnY[rcase]*(((num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*(((num_float3)/data.boxsize));

	  // wrapping periodic condition 
	  if(xb>1.) xb = xb - 1.;
	  if(yb>1.) yb = yb - 1.;
	  if(zb>1.) zb = zb - 1.;
	  if(xb<0.) xb = 1. + xb;
	  if(yb<0.) yb = 1. + yb;
	  if(zb<0.) zb = 1. + zb;
	  switch (face[rcase]){
	  case(1):
	    x = xb;
	    y = yb;
	    z = zb;
	    break;
	  case(2):
	    x = xb;
	    y = zb;
	    z = yb;
	    break;
	  case(3):
	    x = yb;
	    y = zb;
	    z = xb;
	    break;
	  case(4):
	    x = yb;
	    y = xb;
	    z = zb;
	    break;
	  case(5):
	    x = zb;
	    y = xb;
	    z = yb;
	    break;
	  case(6):
	    x = zb;
	    y = yb;
	    z = xb;
	    break;
	  }
	  // recenter
	  x = x - x0[rcase];
	  y = y - y0[rcase];
	  z = z - z0[rcase];	
	  // wrapping periodic condition again
	  if(x>1.) x = x - 1.;
	  if(y>1.) y = y - 1.;
	  if(z>1.) z = z - 1.;
	  if(x<0.) x = 1. + x;
	  if(y<0.) y = 1. + y;
	  if(z<0.) z = 1. + z;
	  z+=double(rcase); // pile the cones
          xx0.push_back(x);
	  yy0.push_back(y);
	  zz0.push_back(z);
	}
	
	// type1
	for (int pp=data.npart[0]; pp<data.npart[0]+data.npart[1]; pp++) {
	  fin.read((char *)&num_float1, sizeof(num_float1)); 
	  fin.read((char *)&num_float2, sizeof(num_float2)); 
	  fin.read((char *)&num_float3, sizeof(num_float3)); 
	  float x, y, z;
	  float xb, yb, zb;

	  xb = sgnX[rcase]*(((num_float1)/data.boxsize));
	  yb = sgnY[rcase]*(((num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*(((num_float3)/data.boxsize));

	  // wrapping periodic condition 
	  if(xb>1.) xb = xb - 1.;
	  if(yb>1.) yb = yb - 1.;
	  if(zb>1.) zb = zb - 1.;
	  if(xb<0.) xb = 1. + xb;
	  if(yb<0.) yb = 1. + yb;
	  if(zb<0.) zb = 1. + zb;
	  switch (face[rcase]){
	  case(1):
	    x = xb;
	    y = yb;
	    z = zb;
	    break;
	  case(2):
	    x = xb;
	    y = zb;
	    z = yb;
	    break;
	  case(3):
	    x = yb;
	    y = zb;
	    z = xb;
	    break;
	  case(4):
	    x = yb;
	    y = xb;
	    z = zb;
	    break;
	  case(5):
	    x = zb;
	    y = xb;
	    z = yb;
	    break;
	  case(6):
	    x = zb;
	    y = yb;
	    z = xb;
	    break;
	  }
	  // recenter
	  x = x - x0[rcase];
	  y = y - y0[rcase];
	  z = z - z0[rcase];	
	  // wrapping periodic condition again
	  if(x>1.) x = x - 1.;
	  if(y>1.) y = y - 1.;
	  if(z>1.) z = z - 1.;
	  if(x<0.) x = 1. + x;
	  if(y<0.) y = 1. + y;
	  if(z<0.) z = 1. + z;
	  z+=double(rcase); // pile the cones
	  xx1.push_back(x);
	  yy1.push_back(y);
	  zz1.push_back(z);
	}
	
	// type2
	for (int pp=data.npart[0]+data.npart[1]; pp<data.npart[0]+data.npart[1]+data.npart[2]; pp++) {
	  fin.read((char *)&num_float1, sizeof(num_float1)); 
	  fin.read((char *)&num_float2, sizeof(num_float2)); 
	  fin.read((char *)&num_float3, sizeof(num_float3)); 
	  float x, y, z;
	  float xb, yb, zb;

	  xb = sgnX[rcase]*(((num_float1)/data.boxsize));
	  yb = sgnY[rcase]*(((num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*(((num_float3)/data.boxsize));

	  // wrapping periodic condition 
	  if(xb>1.) xb = xb - 1.;
	  if(yb>1.) yb = yb - 1.;
	  if(zb>1.) zb = zb - 1.;
	  if(xb<0.) xb = 1. + xb;
	  if(yb<0.) yb = 1. + yb;
	  if(zb<0.) zb = 1. + zb;
	  switch (face[rcase]){
	  case(1):
	    x = xb;
	    y = yb;
	    z = zb;
	    break;
	  case(2):
	    x = xb;
	    y = zb;
	    z = yb;
	    break;
	  case(3):
	    x = yb;
	    y = zb;
	    z = xb;
	    break;
	  case(4):
	    x = yb;
	    y = xb;
	    z = zb;
	    break;
	  case(5):
	    x = zb;
	    y = xb;
	    z = yb;
	    break;
	  case(6):
	    x = zb;
	    y = yb;
	    z = xb;
	    break;
	  }
	  // recenter
	  x = x - x0[rcase];
	  y = y - y0[rcase];
	  z = z - z0[rcase];	
	  // wrapping periodic condition again
	  if(x>1.) x = x - 1.;
	  if(y>1.) y = y - 1.;
	  if(z>1.) z = z - 1.;
	  if(x<0.) x = 1. + x;
	  if(y<0.) y = 1. + y;
	  if(z<0.) z = 1. + z;
	  z+=double(rcase); // pile the cones
	  xx2.push_back(x);
	  yy2.push_back(y);
	  zz2.push_back(z);

	}
	
	// type3
	for (int pp=data.npart[0]+data.npart[1]+data.npart[2]; 
	     pp<data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]; pp++) {
	  fin.read((char *)&num_float1, sizeof(num_float1)); 
	  fin.read((char *)&num_float2, sizeof(num_float2)); 
	  fin.read((char *)&num_float3, sizeof(num_float3)); 
	  float x, y, z;
	  float xb, yb, zb;

	  xb = sgnX[rcase]*(((num_float1)/data.boxsize));
	  yb = sgnY[rcase]*(((num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*(((num_float3)/data.boxsize));

	  // wrapping periodic condition 
	  if(xb>1.) xb = xb - 1.;
	  if(yb>1.) yb = yb - 1.;
	  if(zb>1.) zb = zb - 1.;
	  if(xb<0.) xb = 1. + xb;
	  if(yb<0.) yb = 1. + yb;
	  if(zb<0.) zb = 1. + zb;
	  switch (face[rcase]){
	  case(1):
	    x = xb;
	    y = yb;
	    z = zb;
	    break;
	  case(2):
	    x = xb;
	    y = zb;
	    z = yb;
	    break;
	  case(3):
	    x = yb;
	    y = zb;
	    z = xb;
	    break;
	  case(4):
	    x = yb;
	    y = xb;
	    z = zb;
	    break;
	  case(5):
	    x = zb;
	    y = xb;
	    z = yb;
	    break;
	  case(6):
	    x = zb;
	    y = yb;
	    z = xb;
	    break;
	  }
	  // recenter
	  x = x - x0[rcase];
	  y = y - y0[rcase];
	  z = z - z0[rcase];	
	  // wrapping periodic condition again
	  if(x>1.) x = x - 1.;
	  if(y>1.) y = y - 1.;
	  if(z>1.) z = z - 1.;
	  if(x<0.) x = 1. + x;
	  if(y<0.) y = 1. + y;
	  if(z<0.) z = 1. + z;
	  z+=double(rcase); // pile the cones
	  xx3.push_back(x);
	  yy3.push_back(y);
	  zz3.push_back(z);
	}
	
	// type4
	for (int pp=data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]; 
	     pp<data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]+data.npart[4]; pp++) {
	  fin.read((char *)&num_float1, sizeof(num_float1)); 
	  fin.read((char *)&num_float2, sizeof(num_float2)); 
	  fin.read((char *)&num_float3, sizeof(num_float3)); 
	  float x, y, z;
	  float xb, yb, zb;

	  xb = sgnX[rcase]*(((num_float1)/data.boxsize));
	  yb = sgnY[rcase]*(((num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*(((num_float3)/data.boxsize));

	  // wrapping periodic condition 
	  if(xb>1.) xb = xb - 1.;
	  if(yb>1.) yb = yb - 1.;
	  if(zb>1.) zb = zb - 1.;
	  if(xb<0.) xb = 1. + xb;
	  if(yb<0.) yb = 1. + yb;
	  if(zb<0.) zb = 1. + zb;
	  switch (face[rcase]){
	  case(1):
	    x = xb;
	    y = yb;
	    z = zb;
	    break;
	  case(2):
	    x = xb;
	    y = zb;
	    z = yb;
	    break;
	  case(3):
	    x = yb;
	    y = zb;
	    z = xb;
	    break;
	  case(4):
	    x = yb;
	    y = xb;
	    z = zb;
	    break;
	  case(5):
	    x = zb;
	    y = xb;
	    z = yb;
	    break;
	  case(6):
	    x = zb;
	    y = yb;
	    z = xb;
	    break;
	  }
	  // recenter
	  x = x - x0[rcase];
	  y = y - y0[rcase];
	  z = z - z0[rcase];	
	  // wrapping periodic condition again
	  if(x>1.) x = x - 1.;
	  if(y>1.) y = y - 1.;
	  if(z>1.) z = z - 1.;
	  if(x<0.) x = 1. + x;
	  if(y<0.) y = 1. + y;
	  if(z<0.) z = 1. + z;
	  z+=double(rcase); // pile the cones
	  xx4.push_back(x);
	  yy4.push_back(y);
	  zz4.push_back(z);
	}
	
	// type5
	for (int pp=data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]+data.npart[4]; 
	     pp<data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]+data.npart[4]+data.npart[5]; pp++) {
	  fin.read((char *)&num_float1, sizeof(num_float1)); 
	  fin.read((char *)&num_float2, sizeof(num_float2)); 
	  fin.read((char *)&num_float3, sizeof(num_float3)); 
	  float x, y, z;
	  float xb, yb, zb;

	  xb = sgnX[rcase]*(((num_float1)/data.boxsize));
	  yb = sgnY[rcase]*(((num_float2)/data.boxsize));
	  zb = sgnZ[rcase]*(((num_float3)/data.boxsize));

	  // wrapping periodic condition 
	  if(xb>1.) xb = xb - 1.;
	  if(yb>1.) yb = yb - 1.;
	  if(zb>1.) zb = zb - 1.;
	  if(xb<0.) xb = 1. + xb;
	  if(yb<0.) yb = 1. + yb;
	  if(zb<0.) zb = 1. + zb;
	  switch (face[rcase]){
	  case(1):
	    x = xb;
	    y = yb;
	    z = zb;
	    break;
	  case(2):
	    x = xb;
	    y = zb;
	    z = yb;
	    break;
	  case(3):
	    x = yb;
	    y = zb;
	    z = xb;
	    break;
	  case(4):
	    x = yb;
	    y = xb;
	    z = zb;
	    break;
	  case(5):
	    x = zb;
	    y = xb;
	    z = yb;
	    break;
	  case(6):
	    x = zb;
	    y = yb;
	    z = xb;
	    break;
	  }
	  // recenter
	  x = x - x0[rcase];
	  y = y - y0[rcase];
	  z = z - z0[rcase];	
	  // wrapping periodic condition again
	  if(x>1.) x = x - 1.;
	  if(y>1.) y = y - 1.;
	  if(z>1.) z = z - 1.;
	  if(x<0.) x = 1. + x;
	  if(y<0.) y = 1. + y;
	  if(z<0.) z = 1. + z;
	  z+=double(rcase); // pile the cones
	  xx5.push_back(x);
	  yy5.push_back(y);
	  zz5.push_back(z);
	}

      }

      if(hydro){

	int tot = (data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]+data.npart[4]+data.npart[5]);

      	BLOCK block; fin >> block;
      	cout << "Size of Next Block is " << block.blocksize1 << endl;
	cout << "Should be             " << 3*sizeof(int32_t)*tot << endl;

        cout << "Fast Fowarding next block. Name: ";
        cout << block.name[0];
        cout << block.name[1];
        cout << block.name[2];
        cout << block.name[3] << endl;

        fin.seekg(block.blocksize2/sizeof(int8_t),fin.cur);

	fin >> block;
      	cout << "Size of Header is " << block.blocksize1 << endl;
	cout << "Should be         " << 12*tot*sizeof(int8_t) << endl;

        cout << "Fast Fowarding next block. Name: ";
        cout << block.name[0];
        cout << block.name[1];
        cout << block.name[2];
        cout << block.name[3] << endl;

        fin.seekg(block.blocksize2/sizeof(int8_t),fin.cur);

	int mass_pos = fin.tellg();

        fin >> block;
        cout << "Reading next block. Name: ";
        cout << block.name[0];
        cout << block.name[1];
        cout << block.name[2];
        cout << block.name[3] << endl;
        cout << "Should be                 MASS" << endl;
	cout << "Saving MASS position in the stream" << endl;

	// Reading smoothing lengths

	for(int i = 0; i<4;i++){
          fin.seekg(block.blocksize2/sizeof(int8_t),fin.cur);
          fin >> block;
   	  cout << "Fast Fowarding next block. Name: ";
          cout << block.name[0];
          cout << block.name[1];
          cout << block.name[2];
          cout << block.name[3] << endl;	  
        }

      fin.seekg(block.blocksize2/sizeof(int8_t),fin.cur);
      fin >> block;
      cout << "Reading next block. Name: ";
      cout << block.name[0];
      cout << block.name[1];
      cout << block.name[2];
      cout << block.name[3] << endl;
      cout << "Should be                 HSML" << endl;

      for(int i=0;i<xx0.size();i++){

	float hsml;
        fin.read((char *)&hsml, sizeof(hsml));
	smoothing_l.push_back(hsml/(data.boxsize));
 
      }

      fin >> block;
      cout << "Reading next block. Name: ";
      cout << block.name[0];
      cout << block.name[1];
      cout << block.name[2];
      cout << block.name[3] << endl;
      cout << "Should be                 SFR " << endl;

      fin.seekg (mass_pos, fin.beg);

      fin >> block;
      cout << "Returning to MASS block. Name: ";
      cout << block.name[0];
      cout << block.name[1];
      cout << block.name[2];
      cout << block.name[3] << endl;
      cout << "Should be                 MASS" << endl;

      }
      else{fin.clear(); fin.close();}

      int n0 = xx0.size();
      int n1 = xx1.size();
      int n2 = xx2.size();
      int n3 = xx3.size();
      int n4 = xx4.size();
      int n5 = xx5.size();
      
      cout << "  " << endl;
      cout << n0 <<"   type (0) particles selected until now"<<endl;
      cout << n1 <<"   type (1) particles selected until now"<<endl;
      cout << n2 <<"   type (2) particles selected until now"<<endl;
      cout << n3 <<"   type (3) particles selected until now"<<endl;
      cout << n4 <<"   type (4) particles selected until now"<<endl;
      cout << n5 <<"   type (5) particles selected until now"<<endl;
      cout << "  " << endl;
      
      int totPartxy0;
      int totPartxy1;
      int totPartxy2;
      int totPartxy3;
      int totPartxy4;
      int totPartxy5;
      
      if(n0>0){
	
	// quadrate box
	double xmin=double(*min_element(xx0.begin(), xx0.end()));
	double xmax=double(*max_element(xx0.begin(), xx0.end()));  
	double ymin=double(*min_element(yy0.begin(), yy0.end()));
	double ymax=double(*max_element(yy0.begin(), yy0.end()));  
	double zmin=double(*min_element(zz0.begin(), zz0.end()));
	double zmax=double(*max_element(zz0.begin(), zz0.end()));  
	cout << " " << endl;
	cout << " n0 particles " << endl;
	cout << "xmin = " << xmin << endl;
	cout << "xmax = " << xmax << endl;
	cout << "ymin = " << ymin << endl;
	cout << "ymax = " << ymax << endl;
	cout << "zmin = " << zmin << endl;
	cout << "zmax = " << zmax << endl;
	cout << "  " << endl;

	if(xmin<0 || ymin<0 || zmin< 0){
	  cout << "xmin = " << xmin << endl;
	  cout << "xmax = " << xmax << endl;
	  cout << "ymin = " << ymin << endl;
	  cout << "ymax = " << ymax << endl;
	  cout << "zmin = " << zmin << endl;
	  cout << "zmax = " << zmax << endl;
	  cout << "  0 type check this!!! I will STOP here!!! " << endl;
	  exit(1);
	}

	cout << " ... mapping type 0 particles on the grid with " << npix << " pixels" << endl;
	// 2Dgrid
	vector<float> xs(0),ys(0),ms(0),hs(0);
	// vector<double> ra(0),dec(0);
	for(int l=0;l<n0;l++){
	
		if(hydro && data.massarr[0]==0){

 			fin.read((char *)&num_float1, sizeof(num_float1)); 
             	}	
		else{num_float1=m0;}	

	  double di = sqrt(pow(xx0[l]-0.5,2) + pow(yy0[l]-0.5,2) + pow(zz0[l],2))*data.boxsize/1.e+3;
	  if(di>=blD[nsnap] && di<blD2[nsnap]){
	    double rai,deci,dd;
	    getPolar(xx0[l]-0.5,yy0[l]-0.5,zz0[l],&rai,&deci,&dd);
	    if(fabs(rai)<=fovradiants*0.5 && fabs(deci)<=fovradiants*0.5){	  
	      double fovinunitbox = fovradiants*di/(data.boxsize/1.e+3);
	      xs.push_back((xx0[l]-0.5)/fovinunitbox+0.5);
	      ys.push_back((yy0[l]-0.5)/fovinunitbox+0.5);
              ms.push_back(num_float1);
	      hs.push_back(smoothing_l[l]/fovinunitbox);
	    }	  
	  }
	}
	totPartxy0=xs.size();
	// cout << " n0: totPartxy0 " << totPartxy0 << endl;
	ntotxy0+=totPartxy0;
	
	if(totPartxy0>0){
        	if(kernel_file.c_str()=="None") mapxy0 = gridist_w(xs,ys,ms,hs,npix);
                else mapxy0 = gridist_w(xs,ys,ms,hs,npix,interp_ptr);
        }

	// re-normalize to the total mass!
	double mtot0=0.;
	//double mnorm=accumulate(ms.begin(), ms.end(), 0.);
        double mnorm=0;
        for(int i=0;i<totPartxy0;i++)if(ms[i]<1){mnorm+=ms[i];}

	if(totPartxy0>0){
	  for(int l=0;l<npix*npix;l++){
	    mtot0 += mapxy0[l];
	  }
	  for(int l=0;l<npix*npix;l++){
	    mapxy0[l]=mapxy0[l]/mtot0*mnorm;
	  }
	  //std:: cout << " total mass in the map " << mtot0*m0 << std:: endl; 
	  //std:: cout << " total mass in particles " << totPartxy0*m0 << std:: endl; 
	  //exit(1);
	}
      }
      
      if(n1>0){

	// quadrate box
	double xmin=double(*min_element(xx1.begin(), xx1.end()));
	double xmax=double(*max_element(xx1.begin(), xx1.end()));  
	double ymin=double(*min_element(yy1.begin(), yy1.end()));
	double ymax=double(*max_element(yy1.begin(), yy1.end()));  
	double zmin=double(*min_element(zz1.begin(), zz1.end()));
	double zmax=double(*max_element(zz1.begin(), zz1.end()));  
	cout << " " << endl;
	cout << " n1 particles " << endl;
	cout << "xmin = " << xmin << endl;
	cout << "xmax = " << xmax << endl;
	cout << "ymin = " << ymin << endl;
	cout << "ymax = " << ymax << endl;
	cout << "zmin = " << zmin << endl;
	cout << "zmax = " << zmax << endl;
	cout << "  " << endl;
	if(xmin<0 || ymin<0 || zmin< 0){
	  cout << "xmin = " << xmin << endl;
	  cout << "xmax = " << xmax << endl;
	  cout << "ymin = " << ymin << endl;
	  cout << "ymax = " << ymax << endl;
	  cout << "zmin = " << zmin << endl;
	  cout << "zmax = " << zmax << endl;
	  cout << "  1 type check this!!! I will STOP here!!! " << endl;
	  exit(1);
	}
	cout << " ... mapping type 1 particles on the grid with " << npix << " pixels" << endl;
	// 2Dgrid
	vector<float> xs(0),ys(0),ms(0);
	// vector<double> ra(0),dec(0);
        cout << "n1:" << n1 << endl;
	for(int l=0;l<n1;l++){
	  
		if(hydro && data.massarr[1]==0){

 			fin.read((char *)&num_float1, sizeof(num_float1)); 
             	}
		else{num_float1=m1;}

	  double di = sqrt(pow(xx1[l]-0.5,2) + pow(yy1[l]-0.5,2) + pow(zz1[l],2))*data.boxsize/1.e+3;
	  if(di>=blD[nsnap] && di<blD2[nsnap]){
	    double rai,deci,dd;
	    getPolar(xx1[l]-0.5,yy1[l]-0.5,zz1[l],&rai,&deci,&dd);
	    if(fabs(rai)<=fovradiants*0.5 && fabs(deci)<=fovradiants*0.5){	  
	      double fovinunitbox = fovradiants*di/(data.boxsize/1.e+3);
	      xs.push_back((xx1[l]-0.5)/fovinunitbox+0.5);
	      ys.push_back((yy1[l]-0.5)/fovinunitbox+0.5);
              ms.push_back(num_float1);
	    }	  
	  }
	}
	totPartxy1=xs.size();
	cout << " n1: totPartxy1 " << totPartxy1 << endl;
	ntotxy1+=totPartxy1;
	
	if(totPartxy1>0){
        	mapxy1 = gridist_w(xs,ys,ms,npix);
        }
	
	// re-normalize to the total mass!
	double mtot1=0;
	double mnorm=accumulate(ms.begin(), ms.end(), 0.);

	if(totPartxy1>0){
	  for(int l=0;l<npix*npix;l++){
	    mtot1 += mapxy1[l];
	  }
	  for(int l=0;l<npix*npix;l++){
	    mapxy1[l]=mapxy1[l]/mtot1*mnorm;
	  }
	  //std:: cout << " total mass in the map " << mtot1*m1 << std:: endl; 
	  //std:: cout << " total mass in particles " << totPartxy1*m1 << std:: endl; 
	  //exit(1);
	}
      }
      
      if(n2>0){

	// quadrate box
	double xmin=double(*min_element(xx2.begin(), xx2.end()));
	double xmax=double(*max_element(xx2.begin(), xx2.end()));  
	double ymin=double(*min_element(yy2.begin(), yy2.end()));
	double ymax=double(*max_element(yy2.begin(), yy2.end()));  
	double zmin=double(*min_element(zz2.begin(), zz2.end()));
	double zmax=double(*max_element(zz2.begin(), zz2.end()));  
	cout << " " << endl;
	cout << " n2 particles " << endl;
	cout << "xmin = " << xmin << endl;
	cout << "xmax = " << xmax << endl;
	cout << "ymin = " << ymin << endl;
	cout << "ymax = " << ymax << endl;
	cout << "zmin = " << zmin << endl;
	cout << "zmax = " << zmax << endl;
	cout << "  " << endl;
	if(xmin<0 || ymin<0 || zmin< 0){
	  cout << "xmin = " << xmin << endl;
	  cout << "xmax = " << xmax << endl;
	  cout << "ymin = " << ymin << endl;
	  cout << "ymax = " << ymax << endl;
	  cout << "zmin = " << zmin << endl;
	  cout << "zmax = " << zmax << endl;
	  cout << "  2 type check this!!! I will STOP here!!! " << endl;
	  exit(1);
	}
	cout << " ... mapping type 2 particles on the grid with " << npix << " pixels" << endl;
	// 2Dgrid
	vector<float> xs(0),ys(0),ms(0);
	// vector<double> ra(0),dec(0);
	for(int l=0;l<n2;l++){

		if(hydro && data.massarr[2]==0){

 			fin.read((char *)&num_float1, sizeof(num_float1)); 
             	}
		else{num_float1=m2;}

	  double di = sqrt(pow(xx2[l]-0.5,2) + pow(yy2[l]-0.5,2) + pow(zz2[l],2))*data.boxsize/1.e+3;
	  if(di>=blD[nsnap] && di<blD2[nsnap]){
	    double rai,deci,dd;
	    getPolar(xx2[l]-0.5,yy2[l]-0.5,zz2[l],&rai,&deci,&dd);
	    if(fabs(rai)<=fovradiants*0.5 && fabs(deci)<=fovradiants*0.5){	  
	      double fovinunitbox = fovradiants*di/(data.boxsize/1.e+3);
	      xs.push_back((xx2[l]-0.5)/fovinunitbox+0.5);
	      ys.push_back((yy2[l]-0.5)/fovinunitbox+0.5);
              ms.push_back(num_float1);
	    }	  
	  }
	}
	totPartxy2=xs.size();
	// cout << " n2: totPartxy2 " << totPartxy2 << endl;
	ntotxy2+=totPartxy2;

	if(totPartxy2>0){
        	mapxy2 = gridist_w(xs,ys,ms,npix);
	}

	// re-normalize to the total mass!
	double mtot2=0;
	double mnorm=accumulate(ms.begin(),ms.end(),0.);
	if(totPartxy2>0){
	  for(int l=0;l<npix*npix;l++){
	    mtot2 += mapxy2[l];
	  }
	  for(int l=0;l<npix*npix;l++){
	    mapxy2[l]=mapxy2[l]/mtot2*mnorm;
	  }
	  //std:: cout << " total mass in the map " << mtot2*m2 << std:: endl; 
	  //std:: cout << " total mass in particles " << totPartxy2*m2 << std:: endl; 
	  //exit(1);
	}
      }
      
      if(n3>0){

	// quadrate box
	double xmin=double(*min_element(xx3.begin(), xx3.end()));
	double xmax=double(*max_element(xx3.begin(), xx3.end()));  
	double ymin=double(*min_element(yy3.begin(), yy3.end()));
	double ymax=double(*max_element(yy3.begin(), yy3.end()));  
	double zmin=double(*min_element(zz3.begin(), zz3.end()));
	double zmax=double(*max_element(zz3.begin(), zz3.end()));  
	cout << " " << endl;
	cout << " n3 particles " << endl;
	cout << "xmin = " << xmin << endl;
	cout << "xmax = " << xmax << endl;
	cout << "ymin = " << ymin << endl;
	cout << "ymax = " << ymax << endl;
	cout << "zmin = " << zmin << endl;
	cout << "zmax = " << zmax << endl;
	cout << "  " << endl;
	if(xmin<0 || ymin<0 || zmin< 0){
	  cout << "xmin = " << xmin << endl;
	  cout << "xmax = " << xmax << endl;
	  cout << "ymin = " << ymin << endl;
	  cout << "ymax = " << ymax << endl;
	  cout << "zmin = " << zmin << endl;
	  cout << "zmax = " << zmax << endl;
	  cout << "  3 type check this!!! I will STOP here!!! " << endl;
	  exit(1);
	}
	cout << " ... mapping type 3 particles on the grid with " << npix << " pixels" << endl;
	// 2Dgrid
	vector<float> xs(0),ys(0),ms(0);
	// vector<double> ra(0),dec(0);
	for(int l=0;l<n3;l++){

		if(hydro && data.massarr[3]==0){

 			fin.read((char *)&num_float1, sizeof(num_float1)); 
             	}
		else{num_float1=m3;}

	  double di = sqrt(pow(xx3[l]-0.5,2) + pow(yy3[l]-0.5,2) + pow(zz3[l],2))*data.boxsize/1.e+3;
	  if(di>=blD[nsnap] && di<blD2[nsnap]){	   vector<float> mm5(0);
	    double rai,deci,dd;
	    getPolar(xx3[l]-0.5,yy3[l]-0.5,zz3[l],&rai,&deci,&dd);
	    if(fabs(rai)<=fovradiants*0.5 && fabs(deci)<=fovradiants*0.5){	  
	      double fovinunitbox = fovradiants*di/(data.boxsize/1.e+3);
	      xs.push_back((xx3[l]-0.5)/fovinunitbox+0.5);
	      ys.push_back((yy3[l]-0.5)/fovinunitbox+0.5);
              ms.push_back(num_float1);
	    }	  
	  }
	}
	totPartxy3=xs.size();
	// cout << " n3: totPartxy3 " << totPartxy3 << endl;
	ntotxy3+=totPartxy3;
	
	if(totPartxy3>0){
		mapxy3 = gridist_w(xs,ys,ms,npix);
	}
	
	// re-normalize to the total mass!
	double mtot3=0;
	double mnorm=accumulate(ms.begin(),ms.end(),0.);
	if(totPartxy3>0){
	  for(int l=0;l<npix*npix;l++){
	    mtot3 += mapxy3[l];
	  }
	  for(int l=0;l<npix*npix;l++){
	    mapxy3[l]=mapxy3[l]/mtot3*mnorm;
	  }
	  //std:: cout << " total mass in the map " << mtot3*m3 << std:: endl; 
	  //std:: cout << " total mass in particles " << totPartxy3*m3 << std:: endl; 
	  //exit(1);
	}
      }
      
      if(n4>0){

	// quadrate box
	double xmin=double(*min_element(xx4.begin(), xx4.end()));
	double xmax=double(*max_element(xx4.begin(), xx4.end()));  
	double ymin=double(*min_element(yy4.begin(), yy4.end()));
	double ymax=double(*max_element(yy4.begin(), yy4.end()));  
	double zmin=double(*min_element(zz4.begin(), zz4.end()));
	double zmax=double(*max_element(zz4.begin(), zz4.end()));  
	cout << " " << endl;
	cout << " n4 particles " << endl;
	cout << "xmin = " << xmin << endl;
	cout << "xmax = " << xmax << endl;
	cout << "ymin = " << ymin << endl;
	cout << "ymax = " << ymax << endl;
	cout << "zmin = " << zmin << endl;
	cout << "zmax = " << zmax << endl;
	cout << "  " << endl;
	if(xmin<0 || ymin<0 || zmin< 0){
	  cout << "xmin = " << xmin << endl;
	  cout << "xmax = " << xmax << endl;
	  cout << "ymin = " << ymin << endl;
	  cout << "ymax = " << ymax << endl;
	  cout << "zmin = " << zmin << endl;
	  cout << "zmax = " << zmax << endl;
	  cout << "  4 type check this!!! I will STOP here!!! " << endl;
	  exit(1);
	}
	cout << " ... mapping type 4 particles on the grid with " << npix << " pixels" << endl;
	// 2Dgrid
	vector<float> xs(0),ys(0),ms(0);
	// vector<double> ra(0),dec(0);
	for(int l=0;l<n4;l++){

		if(hydro && data.massarr[4]==0){

 			fin.read((char *)&num_float1, sizeof(num_float1)); 
             	}
		else{num_float1=m4;}

	  double di = sqrt(pow(xx4[l]-0.5,2) + pow(yy4[l]-0.5,2) + pow(zz4[l],2))*data.boxsize/1.e+3;
	  if(di>=blD[nsnap] && di<blD2[nsnap]){
	    double rai,deci,dd;
	    getPolar(xx4[l]-0.5,yy4[l]-0.5,zz4[l],&rai,&deci,&dd);
	    if(fabs(rai)<=fovradiants*0.5 && fabs(deci)<=fovradiants*0.5){	  
	      double fovinunitbox = fovradiants*di/(data.boxsize/1.e+3);
	      xs.push_back((xx4[l]-0.5)/fovinunitbox+0.5);
	      ys.push_back((yy4[l]-0.5)/fovinunitbox+0.5);
              ms.push_back(num_float1);
	    }	  
	  }
	}
	totPartxy4=xs.size();
	// cout << " n4: totPartxy4 " << totPartxy4 << endl;
	ntotxy4+=totPartxy4;
	
	if(totPartxy4>0){
		mapxy4 = gridist_w(xs,ys,ms,npix);
	}
	
	// re-normalize to the total mass!
	double mtot4=0;
        double mnorm = accumulate(ms.begin(),ms.end(),0.); 
	if(totPartxy4>0){
	  for(int l=0;l<npix*npix;l++){
	    mtot4 += mapxy4[l];
	  }
	  for(int l=0;l<npix*npix;l++){
	    mapxy4[l]=mapxy4[l]/mtot4*mnorm;
	  }
	  //std:: cout << " total mass in the map " << mtot4*m4 << std:: endl; 
	  //std:: cout << " total mass in particles " << totPartxy4*m4 << std:: endl; 
	  //exit(1);
	}
      }
      
      if(n5>0){

	// quadrate box
	double xmin=double(*min_element(xx5.begin(), xx5.end()));
	double xmax=double(*max_element(xx5.begin(), xx5.end()));  
	double ymin=double(*min_element(yy5.begin(), yy5.end()));
	double ymax=double(*max_element(yy5.begin(), yy5.end()));  
	double zmin=double(*min_element(zz5.begin(), zz5.end()));
	double zmax=double(*max_element(zz5.begin(), zz5.end()));  
	cout << " " << endl;
	cout << " n5 particles " << endl;
	cout << "xmin = " << xmin << endl;
	cout << "xmax = " << xmax << endl;
	cout << "ymin = " << ymin << endl;
	cout << "ymax = " << ymax << endl;
	cout << "zmin = " << zmin << endl;
	cout << "zmax = " << zmax << endl;
	cout << "  " << endl;
	if(xmin<0 || ymin<0 || zmin< 0){
	  cout << "xmin = " << xmin << endl;
	  cout << "xmax = " << xmax << endl;
	  cout << "ymin = " << ymin << endl;
	  cout << "ymax = " << ymax << endl;
	  cout << "zmin = " << zmin << endl;
	  cout << "zmax = " << zmax << endl;
	  cout << "  5 type check this!!! I will STOP here!!! " << endl;
	  exit(1);
	}
	cout << " ... mapping type 5 particles on the grid with " << npix << " pixels" << endl;
	// 2Dgrid
	vector<float> xs(0),ys(0),ms(0);

	for(int l=0;l<n5;l++){

		if(hydro && data.massarr[5]==0){
	          
                  if(l==0){

                    fin.seekg(n5*sizeof(int32_t)/sizeof(int8_t),fin.cur);
                    cout << "\n" <<endl;
                    cout << "Fast Fowarding blocks: ";
		    for(int i=0;i<9;i++){

			fin >> block;
		        cout << block.name[0];
		        cout << block.name[1];
		        cout << block.name[2];
		        cout << block.name[3] << ", ";

		        fin.seekg(block.blocksize2/sizeof(int8_t),fin.cur);
                    }

                    cout<< "\n" << endl;
                    fin >> block;
		    cout << "Reading BH masses from. Block: ";
		    cout << block.name[0];
		    cout << block.name[1];
		    cout << block.name[2];
		    cout << block.name[3] << endl;
                  }

 		  fin.read((char *)&num_float1, sizeof(num_float1)); 
             	}
		else{num_float1=m5;}

	  double di = sqrt(pow(xx5[l]-0.5,2) + pow(yy5[l]-0.5,2) + pow(zz5[l],2))*data.boxsize/1.e+3;
	  if(di>=blD[nsnap] && di<blD2[nsnap]){
	    double rai,deci,dd;
	    getPolar(xx5[l]-0.5,yy5[l]-0.5,zz5[l],&rai,&deci,&dd);
	    if(fabs(rai)<=fovradiants*0.5 && fabs(deci)<=fovradiants*0.5){	  
	      double fovinunitbox = fovradiants*di/(data.boxsize/1.e+3);
	      xs.push_back((xx5[l]-0.5)/fovinunitbox+0.5);
	      ys.push_back((yy5[l]-0.5)/fovinunitbox+0.5);
              ms.push_back(num_float1);
	    }	  
	  }
	}
	totPartxy5=xs.size();
	// cout << " n5: totPartxy5 " << totPartxy5 << endl;
	ntotxy5+=totPartxy5;

	if(totPartxy5>0){	
	        mapxy5 = gridist_w(xs,ys,ms,npix);
	}

	if(hydro){

		fin >> block;

        	cout << endl << "If masses were read correctly next box should be BHMD " << endl;
	        cout << "Name:                                            ";
	        cout << block.name[0];
	        cout << block.name[1];
	        cout << block.name[2];
	        cout << block.name[3] << endl;

	}

	// re-normalize to the total mass!
	double mtot5=0;
	double mnorm = accumulate(ms.begin(),ms.end(),0.);
	if(totPartxy5>0){
	  for(int l=0;l<npix*npix;l++){
	    mtot5 += mapxy5[l];
	  }
	  for(int l=0;l<npix*npix;l++){
	    mapxy5[l]=mapxy5[l]/mtot5*mnorm;
	  }
	  //std:: cout << " total mass in the map " << mtot5*m5 << std:: endl; 
	  //std:: cout << " total mass in particles " << totPartxy5*m5 << std:: endl; 
	  //exit(1);
	}
      }
      std:: cout << " maps done! " << std:: endl;

      if(hydro){
	fin.clear();
	fin.close(); // Close file here
      }

      // sum up all the maps
      mapxytot+=mapxy0+mapxy1+mapxy2+mapxy3+mapxy4+mapxy5;
      mapxytot0+=mapxy0;
      mapxytot1+=mapxy1;
      mapxytot2+=mapxy2;
      mapxytot3+=mapxy3;
      mapxytot4+=mapxy4;
      mapxytot5+=mapxy5;


      std:: cout << " done map*tot " << std:: endl;

    }
    if(noSNAP==0){
      if(partinplanes=="ALL"){
	/**
	 * write image array(s) to FITS files all particles in a FITS file!
	 */
	long naxis = 2;
	long naxes[2]={ npix,npix };
	string fileoutput;
	// fileoutput = simulation+"."+snapnum+".plane_"+snpix+".fits";
	fileoutput = directory+simulation+"."+snappl+".plane_"+snpix+"_"+suffix+".fits";
        cout << "Saving the maps on: " << fileoutput << endl;
	std::unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	std::vector<long> naxex( 2 );
	naxex[0]=npix;
	naxex[1]=npix;
	PHDU *phxy=&ffxy->pHDU();
	phxy->write( 1, npix*npix, mapxytot );
	// phxy->addKey ("x0",x0," unit of the boxsize");
	// phxy->addKey ("y0",y0," unit of the boxsize");  
	phxy->addKey ("REDSHIFT",zsim," "); 
	phxy->addKey ("PHYSICALSIZE",fov," "); 
	phxy->addKey ("PIXELUNIT",1.e+10/h0," "); 
	phxy->addKey ("DlLOW",blD[nsnap]/h0,"comoving distance in Mpc"); 
	phxy->addKey ("DlUP",blD2[nsnap]/h0,"comoving distance in Mpc"); 
	phxy->addKey ("nparttype0",ntotxy0," "); 
	phxy->addKey ("nparttype1",ntotxy1," "); 
	phxy->addKey ("nparttype2",ntotxy2," "); 
	phxy->addKey ("nparttype3",ntotxy3," "); 
	phxy->addKey ("nparttype4",ntotxy4," "); 
	phxy->addKey ("nparttype5",ntotxy5," "); 
	phxy->addKey ("HUBBLE",h0," "); 
	phxy->addKey ("OMEGAMATTER",om0," "); 
	phxy->addKey ("OMEGALAMBDA",omL0," "); 
	phxy->addKey ("m0",m0," "); 
	phxy->addKey ("m1",m1," "); 
	phxy->addKey ("m2",m2," "); 
	phxy->addKey ("m3",m3," "); 
	phxy->addKey ("m4",m4," "); 
	phxy->addKey ("m5",m5," "); 
      }else{
	/**
	 * write image array(s) to FITS files each particle type in different planes
	 */
	// type0
	if(ntotxy0>0){
	  long naxis = 2;
	  long naxes[2]={ npix,npix };
	  string fileoutput;
	  // fileoutput = simulation+"."+snapnum+".plane_"+snpix+".fits";
          fileoutput = directory+simulation+"."+snappl+".ptype0_plane_"+snpix+"_"+suffix+".fits";
	  std::unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	  std::vector<long> naxex( 2 );
	  naxex[0]=npix;
	  naxex[1]=npix;
	  PHDU *phxy=&ffxy->pHDU();
	  phxy->write( 1, npix*npix, mapxytot0 );
	  // phxy->addKey ("x0",x0," unit of the boxsize");
	  // phxy->addKey ("y0",y0," unit of the boxsize");  
	  phxy->addKey ("REDSHIFT",zsim," "); 
	  phxy->addKey ("PHYSICALSIZE",fov," "); 
	  phxy->addKey ("PIXELUNIT",1.e+10/h0," "); 
	  phxy->addKey ("DlLOW",blD[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("DlUP",blD2[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("nparttype0",ntotxy0," "); 
	  phxy->addKey ("HUBBLE",h0," "); 
	  phxy->addKey ("OMEGAMATTER",om0," "); 
	  phxy->addKey ("OMEGALAMBDA",omL0," "); 
	  phxy->addKey ("m0",m0," "); 
	}
	// type1
	if(ntotxy1>0){
	  long naxis = 2;
	  long naxes[2]={ npix,npix };
	  string fileoutput;
	  // fileoutput = simulation+"."+snapnum+".plane_"+snpix+".fits";
          fileoutput = directory+simulation+"."+snappl+".ptype1_plane_"+snpix+"_"+suffix+".fits";
	  std::unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	  std::vector<long> naxex( 2 );
	  naxex[0]=npix;
	  naxex[1]=npix;
	  PHDU *phxy=&ffxy->pHDU();
	  phxy->write( 1, npix*npix, mapxytot1 );
	  // phxy->addKey ("x0",x0," unit of the boxsize");
	  // phxy->addKey ("y0",y0," unit of the boxsize");  
	  phxy->addKey ("REDSHIFT",zsim," "); 
	  phxy->addKey ("PHYSICALSIZE",fov," "); 
	  phxy->addKey ("PIXELUNIT",1.e+10/h0," "); 
	  phxy->addKey ("DlLOW",blD[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("DlUP",blD2[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("nparttype1",ntotxy1," "); 
	  phxy->addKey ("HUBBLE",h0," "); 
	  phxy->addKey ("OMEGAMATTER",om0," "); 
	  phxy->addKey ("OMEGALAMBDA",omL0," "); 
	  phxy->addKey ("m1",m1," "); 
	}
	// type2
	if(ntotxy2>0){
	  long naxis = 2;
	  long naxes[2]={ npix,npix };
	  string fileoutput;
	  // fileoutput = simulation+"."+snapnum+".plane_"+snpix+".fits";
          fileoutput = directory+simulation+"."+snappl+".ptype2_plane_"+snpix+"_"+suffix+".fits";
	  std::unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	  std::vector<long> naxex( 2 );
	  naxex[0]=npix;
	  naxex[1]=npix;
	  PHDU *phxy=&ffxy->pHDU();
	  phxy->write( 1, npix*npix, mapxytot2 );
	  // phxy->addKey ("x0",x0," unit of the boxsize");
	  // phxy->addKey ("y0",y0," unit of the boxsize");  
	  phxy->addKey ("REDSHIFT",zsim," "); 
	  phxy->addKey ("PHYSICALSIZE",fov," "); 
	  phxy->addKey ("PIXELUNIT",1.e+10/h0," "); 
	  phxy->addKey ("DlLOW",blD[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("DlUP",blD2[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("nparttype2",ntotxy2," "); 
	  phxy->addKey ("HUBBLE",h0," "); 
	  phxy->addKey ("OMEGAMATTER",om0," "); 
	  phxy->addKey ("OMEGALAMBDA",omL0," "); 
	  phxy->addKey ("m2",m2," "); 
	}
	// type3
	if(ntotxy3>0){
	  long naxis = 2;
	  long naxes[2]={ npix,npix };
	  string fileoutput;
	  // fileoutput = simulation+"."+snapnum+".plane_"+snpix+".fits";
          fileoutput = directory+simulation+"."+snappl+".ptype3_plane_"+snpix+"_"+suffix+".fits";
	  std::unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	  std::vector<long> naxex( 2 );
	  naxex[0]=npix;
	  naxex[1]=npix;
	  PHDU *phxy=&ffxy->pHDU();
	  phxy->write( 1, npix*npix, mapxytot3 );
	  // phxy->addKey ("x0",x0," unit of the boxsize");
	  // phxy->addKey ("y0",y0," unit of the boxsize");  
	  phxy->addKey ("REDSHIFT",zsim," "); 
	  phxy->addKey ("PHYSICALSIZE",fov," "); 
	  phxy->addKey ("PIXELUNIT",1.e+10/h0," "); 
	  phxy->addKey ("DlLOW",blD[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("DlUP",blD2[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("nparttype3",ntotxy3," "); 
	  phxy->addKey ("HUBBLE",h0," "); 
	  phxy->addKey ("OMEGAMATTER",om0," "); 
	  phxy->addKey ("OMEGALAMBDA",omL0," "); 
	  phxy->addKey ("m3",m3," "); 
	}
	// type4
	if(ntotxy4>0){
	  long naxis = 2;
	  long naxes[2]={ npix,npix };
	  string fileoutput;
	  // fileoutput = simulation+"."+snapnum+".plane_"+snpix+".fits";
          fileoutput = directory+simulation+"."+snappl+".ptype4_plane_"+snpix+"_"+suffix+".fits";
	  std::unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	  std::vector<long> naxex( 2 );
	  naxex[0]=npix;
	  naxex[1]=npix;
	  PHDU *phxy=&ffxy->pHDU();
	  phxy->write( 1, npix*npix, mapxytot4 );
	  // phxy->addKey ("x0",x0," unit of the boxsize");
	  // phxy->addKey ("y0",y0," unit of the boxsize");  
	  phxy->addKey ("REDSHIFT",zsim," "); 
	  phxy->addKey ("PHYSICALSIZE",fov," "); 
	  phxy->addKey ("PIXELUNIT",1.e+10/h0," "); 
	  phxy->addKey ("DlLOW",blD[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("DlUP",blD2[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("nparttype4",ntotxy4," "); 
	  phxy->addKey ("HUBBLE",h0," "); 
	  phxy->addKey ("OMEGAMATTER",om0," "); 
	  phxy->addKey ("OMEGALAMBDA",omL0," "); 
	  phxy->addKey ("m4",m4," "); 
	}
	// type5
	if(ntotxy5>0){
	  long naxis = 2;
	  long naxes[2]={ npix,npix };
	  string fileoutput;
	  // fileoutput = simulation+"."+snapnum+".plane_"+snpix+".fits";
	  fileoutput = directory+simulation+"."+snappl+".ptype5_plane_"+snpix+"_"+suffix+".fits";
	  std::unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	  std::vector<long> naxex( 2 );
	  naxex[0]=npix;
	  naxex[1]=npix;
	  PHDU *phxy=&ffxy->pHDU();
	  phxy->write( 1, npix*npix, mapxytot5 );
	  // phxy->addKey ("x0",x0," unit of the boxsize");
	  // phxy->addKey ("y0",y0," unit of the boxsize");  
	  phxy->addKey ("REDSHIFT",zsim," "); 
	  phxy->addKey ("PHYSICALSIZE",fov," "); 
	  phxy->addKey ("PIXELUNIT",1.e+10/h0," "); 
	  phxy->addKey ("DlLOW",blD[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("DlUP",blD2[nsnap]/h0,"comoving distance in Mpc"); 
	  phxy->addKey ("nparttype5",ntotxy5," "); 
	  phxy->addKey ("HUBBLE",h0," "); 
	  phxy->addKey ("OMEGAMATTER",om0," "); 
	  phxy->addKey ("OMEGALAMBDA",omL0," "); 
	  phxy->addKey ("m5",m5," "); 
	}     
      }
    }         
  }
  cout << " end of work ... ;-)  " << endl;

  delete interp_ptr;
  return 0;
}
