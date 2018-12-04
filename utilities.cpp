#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_linalg.h>
#include "utilities.h"

void error (const std::string s){                                                            
  std::cerr << "error: " << s << std::endl;                         
  exit (-1);                                                                           
}

void error (const bool f, const std::string s){ 
  if (f){    
    std::cerr << "error: " << s << std::endl;        
    exit (-1);                        
  }                                                                      
}   

void gsl_error_handler (const char * reason, const char * file, int line, int gsl_errno){
  std::cerr
    << "GSL error handler: in "
    << file << ", line " << line << ", error code " << gsl_errno << ": "
    << reason << std::endl;
}

void warning (const std::string s){ 
  std::cout << "warning: " << s << std::endl; 
}

double sign (const double a, const double b){ 
  return b>0?fabs(a):-fabs(a); 
}

int nint (const double x){ 
  return (x<0)?int(x-0.5):int(x+0.5); 
}

void swap (double &a, double &b){
  double t=a;
  a=b;
  b=t;
}

void sort (double *arr, const int n){
  const int m=7,ns=50;
  int l=0,ir=n-1;
  int js=-1;
  int i,j,k;
  double a;
  int *is=new int[ns];
  for (;;){
    if (ir-l<m){
	for (j=l+1;j<=ir;j++){
	  a=arr[j];
	  for (i=j-1;i>=l;i--){
	    if (arr[i] <= a) break;
	    arr[i+1]=arr[i];
	  }
	  arr[i+1]=a;
	}
	if (js == -1) break;
	ir=is[js--];
	l=is[js--];
    }
    else{
      k=(l+ir)/2;
      swap(arr[k],arr[l+1]);
      if (arr[l] > arr[ir])   swap(arr[l],arr[ir]);
      if (arr[l+1] > arr[ir]) swap(arr[l+1],arr[ir]);
      if (arr[l] > arr[l+1])  swap(arr[l],arr[l+1]);
      i=l+1;
      j=ir;
      a=arr[l+1];
      for (;;){
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	swap(arr[i],arr[j]);
      }
      arr[l+1]=arr[j];
      arr[j]=a;
      js+=2;
      if (js > ns) error("stack too small in sort.");
      if (ir-i+1 >= j-l){
	is[js]=ir;
	is[js-1]=i;
	ir=j-1;
      }
      else{
	is[js]=j-1;
	is[js-1]=l;
	l=i;
      }
    }
  }
  delete[] is;
}

std:: vector<double> rotatexYz(double theta,std:: vector<double> v){
  std:: vector<double> u(3);
  u[0] = cos(theta)*v[0]+sin(theta)*v[2];
  u[1] = v[1];
  u[2] = -sin(theta)*v[0]+cos(theta)*v[2];
  return u;
}

std:: vector<double> rotatexyZ(double theta,std:: vector<double> v){
  std:: vector<double> u(3);
  u[0] = cos(theta)*v[0]-sin(theta)*v[1];
  u[1] = sin(theta)*v[0]+cos(theta)*v[1];
  u[2] = v[2];
  return u;
}

double median(std:: vector<double> vec){
  typedef std:: vector<double>::size_type vec_sz;
  vec_sz size = vec.size();
  if (size == 0){
    std:: cout << "median of an empty vector" << std:: endl;
    return 0;
  }
  sort(vec.begin(), vec.end());
  vec_sz mid = size/2;
  return size % 2 == 0 ? (vec[mid] + vec[mid-1]) / 2 : vec[mid];
}

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

double get_ellipticity(double e11,double e12,double e21,double e22){
  double bb = -e11-e22;
  double cc = -e12*e21 + e11*e22;
  double DD = bb*bb-4.*cc;
  double l1,l2;
  if(DD>0){
    l1 = (- bb + sqrt(DD))*0.5;
    l2 = (- bb - sqrt(DD))*0.5;
    if(l1<0){
      std:: cout << " lambda1 < 0  (" << l1 << ")" << std:: endl;
      l1 = 0;
    }
    else l1 = sqrt(l1);
    if(l2<0){
      std:: cout << " lambda2 < 0  (" << l2 << ")" << std:: endl;
      l2 = 0;
    }
    else l2 = sqrt(l2);

    if(l1>l2){
      return (1. - l2/l1);
    }
    else{
      return (1. - l1/l2);
    }
  }
  if(DD<0){
    std:: cout << " Inertial tensor with imaginary eigenvalues " << std:: endl;
    l1 = 0.;
    l2 = 0.;
    return 0.;
  }
  if(DD==0){
    l1 = -bb*0.5;
    l2 = -bb*0.5;
    return 0.;
  }
  return 0.;
}


void InertialT(std:: valarray<float> f,std:: vector<double> x,int n,double f1, double f2, double f3, double f4, double f5,
	       double &e1,double &e2,double &e3,double &e4, double &e5){
  int pcase = 0;
  if(n<0){
    n = -n;
    pcase = 1;
  }
  double e111 = 0;
  double e112 = 0;
  double e121 = 0;
  double e122 = 0;

  double e211 = 0;
  double e212 = 0;
  double e221 = 0;
  double e222 = 0;

  double e311 = 0;
  double e312 = 0;
  double e321 = 0;
  double e322 = 0;

  double e411 = 0;
  double e412 = 0;
  double e421 = 0;
  double e422 = 0;

  double e511 = 0;
  double e512 = 0;
  double e521 = 0;
  double e522 = 0;

  for(int i=0;i<n;i++) for(int j=0;j<n;j++){
      if(pcase==0){
	if(f[i+n*j] >= f1){
	  e111 += x[j]*x[j]; 
	  e112 += x[i]*x[j]; 
	  e121 += x[j]*x[i];
	  e122 += x[i]*x[i];  
	}
	
	if(f[i+n*j] >= f2){
	  e211 += x[j]*x[j]; 
	  e212 += x[i]*x[j]; 
	  e221 += x[j]*x[i];
	  e222 += x[i]*x[i];  
	} 
	
	if(f[i+n*j] >= f3){
	  e311 += x[j]*x[j]; 
	  e312 += x[i]*x[j]; 
	  e321 += x[j]*x[i];
	  e322 += x[i]*x[i];  
	} 
	
	if(f[i+n*j] >= f4){
	  e411 += x[j]*x[j]; 
	  e412 += x[i]*x[j]; 
	  e421 += x[j]*x[i];
	  e422 += x[i]*x[i];  
	} 

	if(f[i+n*j] >= f5){
	  e511 += x[j]*x[j]; 
	  e512 += x[i]*x[j]; 
	  e521 += x[j]*x[i];
	  e522 += x[i]*x[i];  
	} 
      }
      if(pcase==1){
	if(f[i+n*j] <= f1){ 
	  e111 += x[j]*x[j]; 
	  e112 += x[i]*x[j]; 
	  e121 += x[j]*x[i];
	  e122 += x[i]*x[i];  
	}
	
	if(f[i+n*j] <= f2){ 
	  e211 += x[j]*x[j]; 
	  e212 += x[i]*x[j]; 
	  e221 += x[j]*x[i];
	  e222 += x[i]*x[i];  
	} 
	
	if(f[i+n*j] <= f3){ 
	  e311 += x[j]*x[j]; 
	  e312 += x[i]*x[j]; 
	  e321 += x[j]*x[i];
	  e322 += x[i]*x[i];  
	} 
	
	if(f[i+n*j] <= f4){ 
	  e411 += x[j]*x[j]; 
	  e412 += x[i]*x[j]; 
	  e421 += x[j]*x[i];
	  e422 += x[i]*x[i];  
	} 

	if(f[i+n*j] <= f5){ 
	  e511 += x[j]*x[j]; 
	  e512 += x[i]*x[j]; 
	  e521 += x[j]*x[i];
	  e522 += x[i]*x[i];  
	} 
      }
    }
  e1 = get_ellipticity(e111,e112,e121,e122);
  e2 = get_ellipticity(e211,e212,e221,e222);
  e3 = get_ellipticity(e311,e312,e321,e322);
  e4 = get_ellipticity(e411,e412,e421,e422);
  e5 = get_ellipticity(e511,e512,e521,e522);
}

double getY(std:: vector<double> x, std:: vector<double> y,double xi){
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

double getZ(std::valarray<float> map,int npix,double x,double y, std:: vector<double> xi,std:: vector<double> yi){  
  size_t i=locate (xi,x);
  size_t j=locate (yi,y);
  size_t tnpix = npix;
  double fQ11 = map[i+tnpix*j];
  double fQ21 = map[(i+1)+tnpix*j];
  double fQ12 = map[i+tnpix*(j+1)];
  double fQ22 = map[(i+1)+tnpix*(j+1)];
  //double fR1 = (xi[i+1] - x)/(xi[i+1]-xi[i])*fQ11+(x-xi[i])/(xi[i+1]-xi[i])*fQ21;
  //double fR2 = (xi[i+1] - x)/(xi[i+1]-xi[i])*fQ12+(x-xi[i])/(xi[i+1]-xi[i])*fQ22;
  // double fP = (yi[j+1] - y)/(yi[j+1]-yi[j])*fR1 + (y - yi[j])/(yi[j+1]-yi[j])*fR2; 
  double fxy = fQ11*(xi[i+1]-x)*(yi[j+1]-y)/((xi[i+1]-xi[i])*(yi[j+1]-yi[j])) + 
    fQ21*(x-xi[i])*(yi[j+1]-y)/((xi[i+1]-xi[i])*(yi[j+1]-yi[j])) + 
    fQ12*(xi[i+1]-x)*(y-yi[j])/((xi[i+1]-xi[i])*(yi[j+1]-yi[j])) + 
    fQ22*(x-xi[i])*(y-yi[j])/((xi[i+1]-xi[i])*(yi[j+1]-yi[j]));
  return fxy;
}

double cubicInterpolate (double p[4], double x) {
  return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}

double bicubicInterpolate (double p[4][4], double x, double y) {
  double arr[4];
  arr[0] = cubicInterpolate(p[0], y);
  arr[1] = cubicInterpolate(p[1], y);
  arr[2] = cubicInterpolate(p[2], y);
  arr[3] = cubicInterpolate(p[3], y);
  return cubicInterpolate(arr, x);
}

double tricubicInterpolate (double p[4][4][4], double x, double y, double z) {
  double arr[4];
  arr[0] = bicubicInterpolate(p[0], y, z);
  arr[1] = bicubicInterpolate(p[1], y, z);
  arr[2] = bicubicInterpolate(p[2], y, z);
  arr[3] = bicubicInterpolate(p[3], y, z);
  return cubicInterpolate(arr, x);
}

double nCubicInterpolate (int n, double* p, double coordinates[]) {
  assert(n > 0);
  if (n == 1) {
    return cubicInterpolate(p, *coordinates);
  }
  else {
    double arr[4];
    int skip = 1 << (n - 1) * 2;
    arr[0] = nCubicInterpolate(n - 1, p, coordinates + 1);
    arr[1] = nCubicInterpolate(n - 1, p + skip, coordinates + 1);
    arr[2] = nCubicInterpolate(n - 1, p + 2*skip, coordinates + 1);
    arr[3] = nCubicInterpolate(n - 1, p + 3*skip, coordinates + 1);
    return cubicInterpolate(arr, *coordinates);
  }
}

double * estprof(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, 
		 double dr0, double xmax, std:: vector<int> &vi, std:: vector<int> &vj, int ngal){     
  int nbin = int(xmax/dr0); 
  std:: cout << " nbins (in estprof) = " << nbin << std:: endl;                                    
  double *kr = new double[nbin];                                                                   
  for (int k=0;k<nbin;k++){                                                                        
    int contapx=0;                                                     
    kr[k] = 0;                            
    if(ngal>0){
      int bvi,bvj;
      for( int i=0; i<ngal; i++ ){
	bvi=vi[i];
	bvj=vj[i];
	if(r[bvi+ny*bvj]>dr0*double(k) && r[bvi+ny*bvj]<=dr0*double(k+1)){
	  contapx = contapx + 1;                  
	  kr[k] = kr[k] + q[bvi+ny*bvj];
	}                                                                                          
      }                                                                                            
    }
    else{
      // for each bin in r estimate the mean value                                                   
      for( int i=0; i<nx; i++ )                                                                      
	for( int j=0; j<ny; j++ ){                                                                   
	  if(r[i+ny*j]>dr0*double(k) && r[i+ny*j]<=dr0*double(k+1)){                                 
	    contapx = contapx + 1;                                                                   
	    kr[k] = kr[k] + q[i+ny*j];                                                               
	  }                                                                                          
	}                                                                                            
    }
    kr[k] = kr[k]/double(contapx);                                                                 
    if(contapx==0) kr[k]=0.;                                                                       
  }                                                                                                
  return kr; // return the pointer                                                                 
}

double * estprof2(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, 
		  std:: vector<double> rbin, std:: vector<int> &vi, std:: vector<int> &vj, int ngal){     
  int nbin = rbin.size()-1;
  double *kr = new double[nbin];                                                                   
  for (int k=0;k<nbin;k++){                                                                        
    int contapx=0;                                                     
    kr[k] = 0;                            
    if(ngal>0){
      int bvi,bvj;
      for( int i=0; i<ngal; i++ ){
	bvi=vi[i];
	bvj=vj[i];
	if(r[bvi+ny*bvj]>rbin[k] && r[bvi+ny*bvj]<=rbin[k+1]){
	  contapx = contapx + 1;                  
	  kr[k] = kr[k] + q[bvi+ny*bvj];
	}                                                                                          
      }                                                                                            
    }
    else{
      // for each bin in r estimate the mean value                                                   
      for( int i=0; i<nx; i++ )                                                                      
	for( int j=0; j<ny; j++ ){                                                                   
	  if(r[i+ny*j]>rbin[k] && r[i+ny*j]<=rbin[k+1]){
	    contapx = contapx + 1;                                                                   
	    kr[k] = kr[k] + q[i+ny*j];                                                               
	  }                                                                                          
	}                                                                                            
    }
    kr[k] = kr[k]/double(contapx);                                                                 
    if(contapx==0) kr[k]=0.;                                                                       
  }                                                                                                
  return kr; // return the pointer                                                                 
}

// variance of the profile                                                                         
double * estsigmaprof(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, double dr0, 
		      double xmax, std:: vector<int> &vi, std:: vector<int> &vj, int ngal, double* qm){                                                                                      
  int nbin = int(xmax/dr0);                                                                        
  double *kr = new double[nbin];                                                                   
  for (int k=0;k<nbin;k++){                                                                        
    int contapx=0;                                                                                 
    kr[k] = 0;                      
    if(ngal>0){
      int bvi,bvj;
      for( int i=0; i<ngal; i++ ){
	bvi=vi[i];
	bvj=vj[i];
	if(r[bvi+ny*bvj]>dr0*double(k) && r[bvi+ny*bvj]<=dr0*double(k+1)){                                 
	  contapx = contapx + 1;                                                                   
	  kr[k] = gsl_pow_2(q[bvi+ny*bvj]-qm[k]) + kr[k];                                              
	  if(kr[k]<0) {                                                                            
	    std:: cout << "negative " << kr[k] << std:: endl;                                      
	    exit(1);                                                                               
	  }                                                                                        
	}                                                                                          
      }                   
    }
    else{
      // for each bin in r estimate the mean value                                                   
      for( int i=0; i<nx; i++ )                                                                  
	for( int j=0; j<ny; j++ ){                                                                   
	  if(r[i+ny*j]>dr0*double(k) && r[i+ny*j]<=dr0*double(k+1)){                                 
	    contapx = contapx + 1;                                                                   
	    kr[k] = gsl_pow_2(q[i+ny*j]-qm[k]) + kr[k];                                              
	    if(kr[k]<0) {                                                                            
	      std:: cout << "negative " << kr[k] << std:: endl;                                      
	      exit(1);                                                                               
	    }                                                                                        
	  }                                                                                          
	}                   
    }                                                                         
    kr[k] = sqrt(kr[k]/double(contapx));                                                           
    if(contapx==0) kr[k]=0.;                                                                       
  }                                                                                                
  return kr; // return the pointer                                                                 
}

// variance of the profile                                                                         
double * estsigmaprof2(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, std:: vector<double> rbin, 
		      std:: vector<int> &vi, std:: vector<int> &vj, int ngal, double* qm){                                                                                      
  int nbin = rbin.size()-1;
  double *kr = new double[nbin];                                                                   
  for (int k=0;k<nbin;k++){                                                                        
    int contapx=0;                                                                                 
    kr[k] = 0;                      
    if(ngal>0){
      int bvi,bvj;
      for( int i=0; i<ngal; i++ ){
	bvi=vi[i];
	bvj=vj[i];
	if(r[bvi+ny*bvj]>rbin[k] && r[bvi+ny*bvj]<=rbin[k+1]){
	  contapx = contapx + 1;                                                                   
	  kr[k] = gsl_pow_2(q[bvi+ny*bvj]-qm[k]) + kr[k];                                              
	  if(kr[k]<0) {                                                                            
	    std:: cout << "negative " << kr[k] << std:: endl;                                      
	    exit(1);                                                                               
	  }                                                                                        
	}                                                                                          
      }                   
    }
    else{
      // for each bin in r estimate the mean value                                                   
      for( int i=0; i<nx; i++ )                                                                  
	for( int j=0; j<ny; j++ ){                                                                   
	if(r[i+ny*j]>rbin[k] && r[i+ny*j]<=rbin[k+1]){
	    contapx = contapx + 1;                                                                   
	    kr[k] = gsl_pow_2(q[i+ny*j]-qm[k]) + kr[k];                                              
	    if(kr[k]<0) {                                                                            
	      std:: cout << "negative " << kr[k] << std:: endl;                                      
	      exit(1);                                                                               
	    }                                                                                        
	  }                                                                                          
	}                   
    }                                                                         
    kr[k] = sqrt(kr[k]/double(contapx));                                                           
    if(contapx==0) kr[k]=0.;                                                                       
  }                                                                                                
  return kr; // return the pointer                                                                 
}

// create cumulative profile of the maps for each lensing component - spherical simmetry is assumed          
double * estcprof(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, 
		  double dr0, double xmax, std:: vector<int> &vi, std:: vector<int> &vj, int ngal){     
  int nbin = int(xmax/dr0);                                                                        
  std:: cout << " nbins (in estprof) = " << nbin << std:: endl;                                    
  double *kr = new double[nbin];                                                                   
  for (int k=0;k<nbin;k++){                                                                        
    double ranulus = dr0*(double(k)+double(k+1))*0.5;
    int contapx=0;                                                                                 
    kr[k] = 0;                                                                                     
    if(ngal>0){
      int bvi,bvj;
      for( int i=0; i<ngal; i++ ){
	bvi=vi[i];
	bvj=vj[i];
	if(r[bvi+ny*bvj]<=ranulus){                                 
	  contapx = contapx + 1;                                                                   
	  kr[k] = kr[k] + q[bvi+ny*bvj];                                                               
	}                                                                                          
      }                                                                                            
    }
    else{
      // for each bin in r estimate the mean value                                                   
      for( int i=0; i<nx; i++ )                                                                      
	for( int j=0; j<ny; j++ ){                                                         
	  if(r[i+ny*j]<=ranulus){                                 
	    contapx = contapx + 1;                                                                   
	    kr[k] = kr[k] + q[i+ny*j];                                                               
	  }                                                                                          
	}  
    }                                                                                          
    kr[k] = kr[k]/double(contapx);                                                                 
    if(contapx==0) kr[k]=0.;                                                                       
  }                                                                                                
  return kr; // return the pointer                                                                 
}

// variance of the profile                                                                         
double * estsigmacprof(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, 
		       double dr0, double xmax, std:: vector<int> &vi, std:: vector<int> &vj, int ngal, double* qm){               
  int nbin = int(xmax/dr0);                                                                        
  double *kr = new double[nbin];                                                                   
  for (int k=0;k<nbin;k++){    
    double ranulus = dr0*(double(k)+double(k+1))*0.5;                                                                    
    int contapx=0;                                                                                 
    kr[k] = 0;                                                                                     
    if(ngal>0){
      int bvi,bvj;
      for( int i=0; i<ngal; i++ ){
	bvi=vi[i];
	bvj=vj[i];
	if(r[bvi+ny*bvj]<=ranulus){                                                                                           
	  contapx = contapx + 1;                                                                   
	  kr[k] = gsl_pow_2(q[bvi+ny*bvj]-qm[k]) + kr[k];                                              
	  if(kr[k]<0) {                                                                            
	    std:: cout << "negative " << kr[k] << std:: endl;                                      
	    exit(1);                                                                               
	  }                                                                                        
	}                                                                                          
      }                      
    }
    else{
      // for each bin in r estimate the mean value                                                   
      for( int i=0; i<nx; i++ )                                                                      
	for( int j=0; j<ny; j++ ){         
	  if(r[i+ny*j]<=ranulus){                                                                                           
	    contapx = contapx + 1;                                                                   
	    kr[k] = gsl_pow_2(q[i+ny*j]-qm[k]) + kr[k];                                              
	    if(kr[k]<0) {                                                                            
	      std:: cout << "negative " << kr[k] << std:: endl;                                      
	      exit(1);                                                                               
	    }                                                                                        
	  }                                                                                          
	}                      
    }                                                                      
    kr[k] = sqrt(kr[k]/double(contapx));                                                           
    if(contapx==0) kr[k]=0.;                                                                       
  }                                                                                                
  return kr; // return the pointer                                                                 
}

// create cumulative profile of the maps for each lensing component - spherical simmetry is assumed          
double * estcprof2(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, 
		   std:: vector<double> rbin, std:: vector<int> &vi, std:: vector<int> &vj, int ngal){     
  int nbin = rbin.size()-1;                                                            
  double *kr = new double[nbin];                                                                   
  for (int k=1;k<nbin+1;k++){                                                                        
    double ranulus = rbin[k];
    int contapx=0;                                                                                 
    kr[k-1] = 0;                                                                                     
    if(ngal>0){
      int bvi,bvj;
      for( int i=0; i<ngal; i++ ){
	bvi=vi[i];
	bvj=vj[i];
	if(r[bvi+ny*bvj]<=ranulus){                                 
	  contapx = contapx + 1;                                                                   
	  kr[k-1] = kr[k-1] + q[bvi+ny*bvj];                                                               
	}                                                                                          
      }                                                                                            
    }
    else{
      // for each bin in r estimate the mean value                                                   
      for( int i=0; i<nx; i++ )                                                                      
	for( int j=0; j<ny; j++ ){                                                         
	  if(r[i+ny*j]<=ranulus){                                 
	    contapx = contapx + 1;                                                                   
	    kr[k-1] = kr[k-1] + q[i+ny*j];                                                               
	  }                                                                                          
	}  
    }                                                                                          
    kr[k-1] = kr[k-1]/double(contapx);                                                                 
    if(contapx==0) kr[k-1]=0.;                                                                       
  }                                                                                                
  return kr; // return the pointer                                                                 
}

// variance of the profile                                                                         
double * estsigmacprof2(std:: valarray<float> q,int nx,int ny, std:: valarray<float> r, 
		       std:: vector<double> rbin, std:: vector<int> &vi, std:: vector<int> &vj, int ngal, double* qm){               
  int nbin = rbin.size()-1;                                                            
  double *kr = new double[nbin];                                                                   
  for (int k=1;k<nbin+1;k++){    
    double ranulus = rbin[k];
    int contapx=0;                                                                                 
    kr[k-1] = 0;                                                                                     
    if(ngal>0){
      int bvi,bvj;
      for( int i=0; i<ngal; i++ ){
	bvi=vi[i];
	bvj=vj[i];
	if(r[bvi+ny*bvj]<=ranulus){                                                                                           
	  contapx = contapx + 1;                                                                   
	  kr[k-1] = gsl_pow_2(q[bvi+ny*bvj]-qm[k-1]) + kr[k-1];                                              
	  if(kr[k-1]<0) {                                                                            
	    std:: cout << "negative " << kr[k-1] << std:: endl;                                      
	    exit(1);                                                                               
	  }                                                                                        
	}                                                                                          
      }                      
    }
    else{
      // for each bin in r estimate the mean value                                                   
      for( int i=0; i<nx; i++ )                                                                      
	for( int j=0; j<ny; j++ ){         
	  if(r[i+ny*j]<=ranulus){                                                                                           
	    contapx = contapx + 1;                                                                   
	    kr[k-1] = gsl_pow_2(q[i+ny*j]-qm[k-1]) + kr[k-1];                                              
	    if(kr[k-1]<0) {                                                                            
	      std:: cout << "negative " << kr[k-1] << std:: endl;                                      
	      exit(1);                                                                               
	    }                                                                                        
	  }                                                                                          
	}                      
    }                                                                      
    kr[k-1] = sqrt(kr[k-1]/double(contapx));                                                           
    if(contapx==0) kr[k-1]=0.;                                                                       
  }                                                                                                
  return kr; // return the pointer                                                                 
}

double fcos(double x){
  double fastcos;
  //always wrap input angle to -PI..PI
  if (x < -pi) x += twopi;
  else if (x > pi){x -= twopi;}

  //compute cosine: sin(x + PI/2) = cos(x)
  x += 1.57079632;
  if (x >  pi) x -= twopi;
  
  if (x < 0){
    fastcos = 1.27323954 * x + 0.405284735 * x * x;
    if (fastcos < 0)
      fastcos = .225 * (fastcos *-fastcos - fastcos) + fastcos;
    else
      fastcos = .225 * (fastcos * fastcos - fastcos) + fastcos;
  }
  else{
    fastcos = 1.27323954 * x - 0.405284735 * x * x;
    if (fastcos < 0)
      fastcos = .225 * (fastcos *-fastcos - fastcos) + fastcos;
    else
      fastcos = .225 * (fastcos * fastcos - fastcos) + fastcos;
  }
  return fastcos;
}

double fsin(double x){
  double fastsin;
  //always wrap input angle to -PI..PI
  if (x < -pi) x += twopi;
  else if (x > pi){x -= twopi;}
  
  //compute sine
  if (x < 0){
    fastsin = 1.27323954 * x + .405284735 * x * x;
    if (fastsin < 0)
      fastsin = .225 * (fastsin *- fastsin - fastsin) + fastsin;
    else
      fastsin = .225 * (fastsin * fastsin - fastsin) + fastsin;
  }
  else{
    fastsin = 1.27323954 * x - 0.405284735 * x * x;
    if (fastsin < 0)
      fastsin = .225 * (fastsin *-fastsin - fastsin) + fastsin;
    else
      fastsin = .225 * (fastsin * fastsin - fastsin) + fastsin;
  }
  return fastsin;
}
