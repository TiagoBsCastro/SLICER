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

void getPolar(double x, double y, double z, double &ang1, double &ang2, double &d, bool radec){
  if(radec){
    d = sqrt(x*x+y*y+z*z);
    ang2 = asin(x/(d)); // dec
    ang1 = atan2(y,z); // ra
  }else{
    d = sqrt(x*x+y*y+z*z);
    ang1 = acos(z/(d)); // theta
    ang2 = atan2(y,x); // phi
  }
}

// grid points distribution function with != wheights
valarray<float> gridist_w (vector<float> x, vector<float> y , vector<float> w, int nn, bool do_NGP){
  valarray<float> grxy( nn*nn );
  int n0 = x.size();
  if(n0!=y.size()){cout << "x and y positions don't match!" << endl; exit(-1);}
  if(n0!=w.size()){cout << "positions and wheights don't match!" << endl; exit(-1);}

  double dl = 1./double(nn);
  // --- - - - - - - - - - - - - - - - - - - - - - - ---
  //                               _ _ _ _ _ _
  //  The order of the points is: |_7_|_8_|_9_|
  //                              |_4_|_5_|_6_|
  //                              |_1_|_2_|_3_|
  // coordinate between 0 and 1 and mass particle = 1
  //
  // --- - - - - - - - - - - - - - - - - - - - - - - ---
  for (int i=0; i<n0; i++){

    int grx=floor(x[i]/dl)+1;
    int gry=floor(y[i]/dl)+1;

    int   gridpointx[9], gridpointy[9];
    float posgridx[9], posgridy[9];
    float wfx[9], wfy[9];

    gridpointx[0] = grx;
    gridpointy[0] = gry;

    if(do_NGP){
      if(grx>=0 && grx<nn && gry>=0 && gry<nn) grxy[grx+nn*gry] = grxy[grx+nn*gry] + w[i];
    }
    else{
      for (int j=0; j<9; j++){
        gridpointx[j]=gridpointx[0]+(j%3)-1;
        gridpointy[j]=gridpointy[0]+(j/3)-1;

        posgridx[j]=(gridpointx[j]+0.5)*dl;
        posgridy[j]=(gridpointy[j]+0.5)*dl;

        wfx[j] = sqrt(w[i])*weight(x[i],posgridx[j],dl);
        wfy[j] = sqrt(w[i])*weight(y[i],posgridy[j],dl);

        int grxc = gridpointx[j];
        int gryc = gridpointy[j];

        if(grxc>=0 && grxc<nn && gryc>=0 && gryc<nn) grxy[grxc+nn*gryc] = grxy[grxc+nn*gryc] + float(wfx[j])*float(wfy[j]);
      }
    }
  }
  return grxy;
}
