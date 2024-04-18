#include "utilities.h"

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

  // --- - - - - - - - - - - - - - - - - - - - - - - ---
  //                               _ _ _ _ _ _
  //  The order of the points is: |_6_|_7_|_8_|
  //                              |_3_|_4_|_5_|
  //                              |_0_|_1_|_2_|
  // coordinate between 0 and 1 and mass particle = 1
  //
  // --- - - - - - - - - - - - - - - - - - - - - - - ---

  valarray<float> grxy( nn*nn );
  int n0 = x.size();
  double dl = 1./double(nn);
  int   gridpointx[9], gridpointy[9];
  float posgridx, posgridy;
  float wfx, wfy;

  if(n0!=y.size()){cout << "x and y positions don't match!" << endl; exit(-1);}
  if(n0!=w.size()){cout << "positions and wheights don't match!" << endl; exit(-1);}

  for (int i=0; i<n0; i++){

    gridpointx[4] = floor(x[i]/dl);
    gridpointy[4] = floor(y[i]/dl);

    if(do_NGP){
      if(gridpointx[4]>=0 && gridpointx[4]<nn && gridpointy[4]>=0 && gridpointy[4]<nn) 
        grxy[gridpointx[4]+nn*gridpointy[4]] = grxy[gridpointx[4]+nn*gridpointy[4]] + w[i];
    }
    else{
      for (int j=0; j<9; j++){

        gridpointx[j]=gridpointx[4]+(j%3)-1;
        gridpointy[j]=gridpointy[4]+(j/3)-1;

        posgridx=(gridpointx[j]+0.5)*dl;
        posgridy=(gridpointy[j]+0.5)*dl;

        wfx = sqrt(w[i])*weight(x[i],posgridx,dl);
        wfy = sqrt(w[i])*weight(y[i],posgridy,dl);

        if(gridpointx[j]>=0 && gridpointx[j]<nn && gridpointy[j]>=0 && gridpointy[j]<nn) 
          grxy[gridpointx[j]+nn*gridpointy[j]] = grxy[gridpointx[j]+nn*gridpointy[j]] + wfx*wfy;
      }
    }
  }
  return grxy;
}