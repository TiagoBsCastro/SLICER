#include "util.h"

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

double getY(vector<double> x,vector<double> y,double xi){  // Interpolated routine to calculate y(xi)
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
  i = min (max (i,0), int (nn)-2);
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

int getSnap (vector <double> & zsnap, vector <double> & dl, vector <double> & zl, double dlens){

  int pos;
  double aux=99999;

  for(int i=0;i<zsnap.size();i++){
      float test= abs( getY(zl,dl,zsnap[i])-dlens );
    if( test < aux ){
      aux = test;
      pos = i;
    }
  }

  return pos;
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
// grid points distribution function with != wheights
valarray<float> gridist_w (vector<float> x, vector<float> y , vector<float> w, int nn){
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
  return grxy;
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

void getPolar(double x, double y, double z, double *ra, double *dec, double *d){
  *d = sqrt(x*x+y*y+z*z);
  *dec = asin(x/(*d));
  *ra = atan2(y,z);
}

void readInput(struct InputParams *p, std::string name){
  // read fileinput
  std:: ifstream fin (name.c_str());
  if(fin.is_open());
  else{
    std:: cout << " Params.ini file does not exist where you are running the code " << std:: endl;
    std:: cout << " I will STOP here!!! " << std:: endl;
    exit(1);
  }
  std:: string str;
  std::getline(fin, str);
  std::getline(fin, str);
  p->npix = std::stoi(str);//         1. Number of Pixels
  std::getline(fin, str);
  std::getline(fin, str);
  p->boxl = std::stof(str);//         2. Box Size
  std::getline(fin, str);
  std::getline(fin, str);
  p->zs = std::stof(str);//           3. Redshift Source
  std::getline(fin, str);
  std::getline(fin, str);
  p->fov = std::stof(str);//          4. Field of View
  std::getline(fin, str);
  std::getline(fin,
      p->filredshiftlist);//          5. File with the redshift list it may contain three columns: snap 1/(1+z) z
  std::getline(fin, str);
  std::getline(fin,
      p->pathsnap);//                 6. Path where the snaphosts are located
  std::getline(fin, str);
  std::getline(fin,
      p->simulation);//               7. Simulation name
  std::getline(fin, str);
  std::getline(fin, str);
  p->nfiles = std::stoi(str);//       8. Number of files on Gadget Snapshot
  std::getline(fin, str);
  std::getline(fin, str);
  p->seedcenter = std::stoi(str);//   9. Seed center
  std::getline(fin, str);
  std::getline(fin, str);
  p->seedface = std::stoi(str);//     10. Seed Face
  std::getline(fin, str);
  std::getline(fin, str);
  p->seedsign = std::stoi(str);//     11. Seed sign
  std::getline(fin, str);
  std::getline(fin, str);
  p->partinplanes = std::stoi(str);// 12. True: Each gadget particle type will have its own Map; False: One Map for All
  std::getline(fin, str);
  std::getline(fin,
       p->directory);//               13. Directory to save FITS files
  std::getline(fin, str);
  std::getline(fin, str);
  p->suffix;//                        14. Suffix to FITS files
  std::getline(fin, str);
  std::getline(fin, str);
  p->snopt;//                         15. Shot-noise option: 0-No random Degradation; 1-Half particles degradation ...

}

void readParameters(string file_name,int *npix, double *boxl,
                    double *zs, double *fov, string *filredshiftlist,
                    string *filsnaplist, string *pathsnap,string *idc,
                    int *seedcenter, int *seedface, int *seedsign,
                    string *simulation, int *nfiles,string *partinplanes,
                    string *directory,string *suffix,int *sn_opt,bool *do_NGP){

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
    ifilin >> butstr; // seed for the random location of the center
    ifilin >> *seedcenter;
    ifilin >> butstr; // seed for the random selection of the dice face
    ifilin >> *seedface;
    ifilin >> butstr; // seed for the selection of the sign of the coordinates
    ifilin >> *seedsign;
    ifilin >> butstr; // which particles in the planes (ALL: one file for all, or NO: one for each)
    ifilin >> *partinplanes;
    ifilin >> butstr;//directory to save FITS files
    ifilin >> *directory;
    ifilin >> butstr;//Sufix to save FITS files
    ifilin >> *suffix;
    ifilin >> butstr;//Shot-noise opt: 0-No random degradation; 1-One half degradation; 2- three quarters degradation...
    ifilin >> *sn_opt;
    ifilin >> butstr;//Mass Assignement Scheme - If 0 use TSC (Recomended) if 1 use NGP
    ifilin >> *do_NGP;
  }
  else exit(-1);
};

void read_redlist(string filredshiftlist, vector <double> & lred, vector <int> & lsnap, double zs){
  ifstream redlist;
  redlist.open(filredshiftlist.c_str());
  if(redlist.is_open()){
    int buta; //snapshot number
    double butb,butc; // snapshots scale factor and redshift
    do{
      redlist >> buta >> butb >> butc;
      lsnap.push_back(buta); // Selected snapshots number and redshift
      lred.push_back(butc);
    }while(butc<zs);
  }else{
    cout << " redshift list file redshift_list.txt does not " << endl;
    cout << " exist in the Code dir ... check this out      " << endl;
    cout << "    I will STOP here !!! " << endl;
    exit(1);
  }
}

void read_dl(string idc, vector <double> & zl, vector <double> & dl, double zs){
  ifstream infiledc;
  infiledc.open(idc.c_str());
  if(infiledc.is_open()){
    double zi,dli;
    while(infiledc >> zi >> dli){
      zl.push_back(zi);
      dl.push_back(dli*speedcunit); // on Mpc/h
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
}

void plans_builder (vector<double> & zl, vector<double> & dl, double boxl, double zs,
  vector<int> & lsnap, vector<double> & lred, vector<double> & ld, vector<double> & ld2, vector<double> & zsimlens,
  vector<int> & fromsnap, vector<int> & replication){

  double dllow=0;
  double dlup=getY(zl,dl,zs);

  double ldbut;
  int pos, nrepi=0;
  int nrep=0;
  do{
    nrep++;
    nrepi++;
    ldbut=dllow+nrep*boxl;
    double dlens=ldbut-0.5*boxl;
    int pos_temp = getSnap(lred, dl, zl, dlens);
    cout << " simulation snapshots = " << ldbut << "  " << getY(dl,zl,ldbut) << "  " << nrep << " from snap "
      << lsnap[pos_temp] << "  " << getY(dl,zl,ldbut-0.5*boxl)
      << endl;
      ld.push_back(ldbut-boxl);
      ld2.push_back(ldbut);
    if ( nrep != 1 && pos_temp != pos){
      for ( int i=0; i<nrepi-1; i++ ) replication.push_back(nrep-1);
      nrepi=1;
    }
    pos=pos_temp;
    zsimlens.push_back(getY(dl,zl,dlens));
    fromsnap.push_back(lsnap[pos]);
  }while(ldbut<dlup);
  for ( int i=0; i<nrepi+1; i++ ) replication.push_back(nrep); // Last plane replications
}

void box_randomize(vector <double> & x0, vector <double> & y0, vector <double> & z0,
  vector <int> & sgnX, vector <int> & sgnY, vector <int> & sgnZ, vector <int> & face,
  int nrandom, int seedsign, int seedface, int seedcenter){

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
    cout << " signs of the coordinates = " << sgnX[i] << "  " << sgnY[i] << " " << sgnZ[i] << endl;
  }
}

void rand_pos (float * xx, HEADER header, int ptype, int nsnap, float rcase, RANDOMIZATION randomization_plan){

  int face=randomization_plan.face;
  int sgnX=randomization_plan.sgnX, sgnY=randomization_plan.sgnY, sgnZ=randomization_plan.sgnZ;
  double xc=randomization_plan.xc, yc=randomization_plan.yc, zc=randomization_plan.zc;

  for (int pp=0; pp<header.npart[ptype]; pp++){

    float x, y, z;
    float xb, yb, zb;

    xb=xx[3*pp]; yb=xx[3*pp+1]; zb=xx[3*pp+2];

    xb = sgnX*(((xb)/header.boxsize));
    yb = sgnY*(((yb)/header.boxsize));
    zb = sgnZ*(((zb)/header.boxsize));

    // wrapping periodic condition
    if(xb>1.) xb = xb - 1.;
    if(yb>1.) yb = yb - 1.;
    if(zb>1.) zb = zb - 1.;
    if(xb<0.) xb = 1. + xb;
    if(yb<0.) yb = 1. + yb;
    if(zb<0.) zb = 1. + zb;
    switch (face){
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
    x = x - xc;
    y = y - yc;
    z = z - zc;
    // wrapping periodic condition again
    if(x>1.) x = x - 1.;
    if(y>1.) y = y - 1.;
    if(z>1.) z = z - 1.;
    if(x<0.) x = 1. + x;
    if(y<0.) y = 1. + y;
    if(z<0.) z = 1. + z;
    z+=double(rcase); // pile the cones
    xx[3*pp]=x;
    xx[3*pp+1]=y;
    xx[3*pp+2]=z;
  }
}

void map_particles (float * xx, double * m, int ptype, int sn_opt, float fovradiants, HEADER header, vector <double> &ld,
    vector <double> &ld2, int nsnap, valarray <float> & mapxy, int * ntotxyptype, bool do_NGP, bool hydro){

  float xmin, xmax, ymin, ymax, zmin, zmax;
  int npix = sqrt(mapxy.size());
  int npart = header.npart[ptype];

  min_element(xx,xx+3*npart,&xmin,&ymin,&zmin);
  max_element(xx,xx+3*npart,&xmax,&ymax,&zmax);

  cout << "xmin = " <<  xmin << endl;
  cout << "xmax = " <<  xmax << endl;
  cout << "ymin = " <<  ymin << endl;
  cout << "ymax = " <<  ymax << endl;
  cout << "zmin = " <<  zmin << endl;
  cout << "zmax = " <<  zmax << endl;
  cout << "  " << endl;

  if( xmin<0 ||  ymin<0 ||  zmin< 0){
    cout << "xmin = " <<  xmin << endl;
    cout << "xmax = " <<  xmax << endl;
    cout << "ymin = " <<  ymin << endl;
    cout << "ymax = " <<  ymax << endl;
    cout << "zmin = " <<  zmin << endl;
    cout << "zmax = " <<  zmax << endl;
    cout << "check this!!! I will STOP here!!! " << endl;
    exit(1);
  }

  cout << " ... mapping type "<< ptype <<" particles on the grid with " << npix << " pixels" << endl;
  // 2Dgrid
  vector<float> xs(0),ys(0),ms(0);

  cout << "Min distance: "<< ld[nsnap]/header.boxsize*1.e+3<< " "<<ld2[nsnap]/header.boxsize*1.e+3 << endl;

  for(int l=0;l<npart;l++){

    float x=xx[3*l],y=xx[3*l+1],z=xx[3*l+2];

    double di = sqrt(pow(x-0.5,2) + pow(y-0.5,2) + pow(z,2))*header.boxsize/1.e+3;
    if(di>=ld[nsnap] && di<ld2[nsnap]){
      double rai,deci,dd;
      getPolar(x-0.5,y-0.5,z,&rai,&deci,&dd);
      if(fabs(rai)<=fovradiants*0.5 && fabs(deci)<=fovradiants*0.5){
        double fovinunitbox = fovradiants*di/(header.boxsize/1.e+3);
        xs.push_back((x-0.5)/fovinunitbox+0.5);
        ys.push_back((y-0.5)/fovinunitbox+0.5);
        if(sn_opt==0){
          if(hydro)
            ms.push_back(m[l]);
          else
            ms.push_back(m[0]);
        }
        else{
          if(rand()/ float(RAND_MAX) < 1./pow(2,sn_opt))
            if(hydro)
              ms.push_back(pow(2,sn_opt)*m[l]);
            else
              ms.push_back(pow(2,sn_opt)*m[0]);
          else
            ms.push_back(0.);
        }
      }
    }
  }
  int totPartxy=xs.size();
  *ntotxyptype+=totPartxy;

  if(totPartxy>0){
    mapxy = gridist_w(xs,ys,ms,npix,do_NGP);
    // re-normalize to the total mass!
    double mtot=0.;
    double mnorm=0.;

    for(int i=0;i<ms.size();i++)
      if(ms[i]<1) // To avoid Box4 bugged mass
        mnorm+=ms[i];

    mtot = mapxy.sum();
    if(mtot==0.)
      mtot=1.; //To avoid NaN
    mapxy*=mnorm/mtot;
  }
}

void write_map (string partinplanes,string directory, string simulation,
      string snappl, string snpix, string suffix, int nsnap, int * ntotxyptype,
      valarray <float> & mapxytot, valarray <float> *mapxytotptype, float fov, float zsim, vector <double> ld, vector <double> ld2,
      float h0, float om0, float omL0, double * m){

  int npix = sqrt(mapxytot.size());
  if(partinplanes=="ALL"){
    /*
    * write image array(s) to FITS files all particles in a FITS file!
    */
    long naxis = 2;
    long naxes[2]={ npix,npix };
    string fileoutput;
    fileoutput = directory+simulation+"."+snappl+".plane_"+snpix+"_"+suffix+".fits";
    cout << "Saving the maps on: " << fileoutput << endl;
    unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
    vector<long> naxex( 2 );
    naxex[0]=npix;
    naxex[1]=npix;
    PHDU *phxy=&ffxy->pHDU();
    phxy->write( 1, npix*npix, mapxytot );
    phxy->addKey ("REDSHIFT",zsim," ");
    phxy->addKey ("PHYSICALSIZE",fov," ");
    phxy->addKey ("PIXELUNIT",1.e+10/h0,"Mass unit in M_Sun");
    phxy->addKey ("DlLOW",ld[nsnap]/h0,"comoving distance in Mpc");
    phxy->addKey ("DlUP",ld2[nsnap]/h0,"comoving distance in Mpc");
    phxy->addKey ("nparttype0",ntotxyptype[0]," ");
    phxy->addKey ("nparttype1",ntotxyptype[1]," ");
    phxy->addKey ("nparttype2",ntotxyptype[2]," ");
    phxy->addKey ("nparttype3",ntotxyptype[3]," ");
    phxy->addKey ("nparttype4",ntotxyptype[4]," ");
    phxy->addKey ("nparttype5",ntotxyptype[5]," ");
    phxy->addKey ("HUBBLE",h0," ");
    phxy->addKey ("OMEGAMATTER",om0," ");
    phxy->addKey ("OMEGALAMBDA",omL0," ");
    phxy->addKey ("m0",m[0]," ");
    phxy->addKey ("m1",m[1]," ");
    phxy->addKey ("m2",m[2]," ");
    phxy->addKey ("m3",m[3]," ");
    phxy->addKey ("m4",m[4]," ");
    phxy->addKey ("m5",m[5]," ");
  }else{
  /**
  * write image array(s) to FITS files each particle type in different planes
  */
    for(int i=0;i<6;i++){
      if(ntotxyptype[i]>0){
        long naxis = 2;
        long naxes[2]={ npix,npix };
        string fileoutput;
        fileoutput = directory+simulation+"."+snappl+".ptype0_plane_"+snpix+"_"+suffix+".fits";
        unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
        vector<long> naxex( 2 );
        naxex[0]=npix;
        naxex[1]=npix;
        PHDU *phxy=&ffxy->pHDU();
        phxy->write( 1, npix*npix, mapxytotptype[i] );
        phxy->addKey ("REDSHIFT",zsim," ");
        phxy->addKey ("PHYSICALSIZE",fov," ");
        phxy->addKey ("PIXELUNIT",1.e+10/h0,"Mass unit in M_Sun");
        phxy->addKey ("DlLOW",ld[nsnap]/h0,"comoving distance in Mpc");
        phxy->addKey ("DlUP",ld2[nsnap]/h0,"comoving distance in Mpc");
        phxy->addKey ("nparttype0",ntotxyptype[i]," ");
        phxy->addKey ("HUBBLE",h0," ");
        phxy->addKey ("OMEGAMATTER",om0," ");
        phxy->addKey ("OMEGALAMBDA",omL0," ");
        phxy->addKey ("m0",m[i]," ");
      }
    }
  }
}

void read_header (ifstream *fin, HEADER &header){

  BLOCK block;
  int32_t blocksize;

  fin->read((char *)&block, sizeof(block));
  fin->read((char *)&blocksize, sizeof(blocksize));
  fin->read((char *)&header, sizeof(header));
  fin->read((char *)&blocksize, sizeof(blocksize));
}
void print_header (HEADER header){
  cout << "Printing Header Data" << endl;

  cout << "N. part.: " << header.npart[0] << " "<< header.npart[1] << " "<< header.npart[2] << " "<< header.npart[3] << " "<< header.npart[4] << " "<< header.npart[5]  << endl;
  cout << "Mass Array: "<< header.massarr[0]<< " "<< header.massarr[1]<< " "<< header.massarr[2]<< " "<< header.massarr[3]<< " "<< header.massarr[4]<< " "<< header.massarr[5]<< " "<< endl;
  cout << "Time: "<< header.time<< endl;
  cout << "Z: "<< header.redshift<< endl;
  cout << "Flag SFR.: "<< header.flag_sfr<< endl;
  cout << "Flag Feedback: "<< header.flag_feedback<< endl;
  cout << "N. tot.: "<< header.npartTotal[0]<<" "<< header.npartTotal[1]<<" "<< header.npartTotal[2]<<" "<< header.npartTotal[3]<<" "<< header.npartTotal[4]<<" "<< header.npartTotal[5]<<" "<< endl;
  cout << "Flag cooling: "<< header.flag_cooling<< endl;
  cout << "N. files: "<< header.numfiles<< endl;
  cout << "Box size: "<< header.boxsize<< endl;
  cout << "Omega_matter: "<< header.om0<< endl;
  cout << "Omega_DE: "<< header.oml<< endl;
  cout << "h: "<< header.h<< endl;
  cout << "Flag sage: "<< header.flag_sage<< endl;
  cout << "Flag metals: "<< header.flag_metals<< endl;
  cout << "N. tot HW: "<< header.nTotalHW[0]<<" "<< header.nTotalHW[1]<<" "<< header.nTotalHW[2]<<" "<< header.nTotalHW[3]<<" "<< header.nTotalHW[4]<<" "<< header.nTotalHW[5]<<" "<< endl;
  cout << "Flag entropy: "<< header.flag_entropy<<endl;
}

void read_pos (ifstream *fin, HEADER *header, float * xx0, float * xx1, float * xx2, float * xx3, float * xx4, float * xx5){

  BLOCK block;
  int32_t blocksize;

  string str="NULL";

  while(strcmp(str.c_str(),"POS ")){

    fin->read((char *)&block, sizeof(block));
    str={block.name[0],block.name[1],block.name[2],block.name[3]};
    fin->read((char *)&blocksize, sizeof(blocksize));
    if (!strcmp(str.c_str(),"POS ")){
      fin->read((char *) xx0, 3*(header->npart[0])*sizeof(float));
      fin->read((char *) xx1, 3*(header->npart[1])*sizeof(float));
      fin->read((char *) xx2, 3*(header->npart[2])*sizeof(float));
      fin->read((char *) xx3, 3*(header->npart[3])*sizeof(float));
      fin->read((char *) xx4, 3*(header->npart[4])*sizeof(float));
      fin->read((char *) xx5, 3*(header->npart[5])*sizeof(float));
    }else
      fin->seekg(blocksize,ios_base::cur);
    fin->read((char *)&blocksize, sizeof(blocksize));
  }
}

void read_mass (ifstream *fin, HEADER *header, double * m0, double * m1, double * m2, double * m3, double * m4, double * m5){

  BLOCK block;
  int32_t blocksize;

  string str="NULL";

  while(strcmp(str.c_str(),"MASS")){
    fin->read((char *)&block, sizeof(block));
    str={block.name[0],block.name[1],block.name[2],block.name[3]};
    fin->read((char *)&blocksize, sizeof(blocksize));
    if (!strcmp(str.c_str(),"MASS")){
      cout << "Reading " << (header->npart[0])*(!(bool)(header->massarr[0])) << " masses of type 0" << endl;
      fin->read((char *) m0, (header->npart[0])*(!(bool)(header->massarr[0]))*sizeof(double));
      cout << "Reading " << (header->npart[1])*(!(bool)(header->massarr[1])) << " masses of type 1" << endl;
      fin->read((char *) m1, (header->npart[1])*(!(bool)(header->massarr[1]))*sizeof(double));
      cout << "Reading " << (header->npart[2])*(!(bool)(header->massarr[2])) << " masses of type 2" << endl;
      fin->read((char *) m2, (header->npart[2])*(!(bool)(header->massarr[2]))*sizeof(double));
      cout << "Reading " << (header->npart[3])*(!(bool)(header->massarr[3])) << " masses of type 3" << endl;
      fin->read((char *) m3, (header->npart[3])*(!(bool)(header->massarr[3]))*sizeof(double));
      cout << "Reading " << (header->npart[4])*(!(bool)(header->massarr[4])) << " masses of type 4" << endl;
      fin->read((char *) m4, (header->npart[4])*(!(bool)(header->massarr[4]))*sizeof(double));
      cout << "Reading " << (header->npart[5])*(!(bool)(header->massarr[5])) << " masses of type 5" << endl;
      fin->read((char *) m5, (header->npart[5])*(!(bool)(header->massarr[5]))*sizeof(double));
    }else
      fin->seekg(blocksize,ios_base::cur);
    fin->read((char *)&blocksize, sizeof(blocksize));
  }
}
