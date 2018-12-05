#include "util.h"

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
  p->seedcenter = std::stoi(str);//   8. Seed center
  std::getline(fin, str);
  std::getline(fin, str);
  p->seedface = std::stoi(str);//     9. Seed Face
  std::getline(fin, str);
  std::getline(fin, str);
  p->seedsign = std::stoi(str);//     10. Seed sign
  std::getline(fin, str);
  std::getline(fin, str);
  p->partinplanes = std::stoi(str);// 11. True: Each gadget particle type will have its own Map; False: One Map for All
  std::getline(fin, str);
  std::getline(fin,
       p->directory);//               12. Directory to save FITS files
  std::getline(fin, str);
  std::getline(fin,
       p->suffix);//                  13. Suffix to FITS files
  std::getline(fin, str);
  std::getline(fin, str);
  p->snopt = std::stoi(str);//        14. Shot-noise option: 0-No random Degradation; 1-Half particles degradation ...

}

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

void build_plans(double dlup, InputParams *p, int numberOfLensPerSnap, int nsnaps, vector <double> & lred, vector <double> & zl, vector <double> & dl, vector <string> & lsnap, vector <double> & ld, vector <double> & ld2, vector <int> & replication, vector <double> & zfromsnap, vector <string> & fromsnap, vector <double> & zsimlens, vector <bool> & randomize, vector <int> & pll,int myid){

  int pos, nrepi=0;
  int nrep=0;
  double ldbut;
  do{
    nrep++;
    nrepi++;
    ldbut=nrep*p->boxl/numberOfLensPerSnap;
    double dlens=ldbut-0.5*p->boxl/numberOfLensPerSnap;
    int pos_temp = getSnap(lred, dl, zl, dlens);
    if (myid==0)
      cout << " simulation snapshots = " << ldbut << "  " << getY(dl,zl,ldbut) << "  " << nrep << " from snap " << lsnap[pos_temp] << "  " << getY(dl,zl,dlens) << endl;
    ld.push_back(ldbut-p->boxl/numberOfLensPerSnap);
    ld2.push_back(ldbut);
    zfromsnap.push_back(lred[pos_temp]);
    if ( nrep != 1 && pos_temp != pos){
      for ( int i=0; i<nrepi-1; i++ ) replication.push_back(nrep-1);
      nrepi=1;
    }
    pos=pos_temp;
    zsimlens.push_back(getY(dl,zl,dlens));
    fromsnap.push_back(lsnap[pos]);
    if( nrep==1 )
      randomize.push_back(1);
    else
      randomize.push_back( !( (nrep-1)%numberOfLensPerSnap ) );

  }while(ldbut<dlup);

  for ( int i=0; i<nrepi+1; i++ ) replication.push_back(nrep); // Last plane replications
  if (myid==0){
    cout << " Comoving Distance of the last plane " << dlup << endl;
    cout << " nsnaps = " << nsnaps << endl;
  }

  ofstream planelist;
  string planes_list;
  planes_list = p->directory+"planes_list_"+p->suffix+".txt";


  if (myid==0)
    planelist.open(planes_list.c_str());

  for(int i=0;i<fromsnap.size();i++){

    if (myid==0){
      cout << zsimlens[i] << " planes = " << ld[i] << "  " << ld2[i] << "  " << replication[i] << " from snap " << fromsnap[i] << endl;
      planelist <<  i << "   " << zsimlens[i] << "   " << ld[i] << "   " << ld2[i] << "   " << replication[i] << "   " << fromsnap[i]
      << "   " << zfromsnap[i] << "  "<< randomize[i]  << endl;
    }
    pll.push_back(i);
  }

  if (myid==0)
    planelist.close();

}

void randomize_box (vector <double> & x0, vector <double> & y0, vector <double> & z0, vector<int> & face, vector<int> & sgnX, vector<int> & sgnY, vector<int> & sgnZ, vector <int> & replication, int nrandom, vector <bool> & randomize, InputParams *p, int numberOfLensPerSnap, int myid){

  for(int i=0;i<nrandom;i++){

    if ( randomize[i] ){

      srand(p->seedcenter+i/numberOfLensPerSnap*13);
      x0[i] = rand() / float(RAND_MAX);
      y0[i] = rand() / float(RAND_MAX);
      z0[i] = rand() / float(RAND_MAX);
      if(myid==0){
        cout << "  " << endl;
        cout << " random centers  for the box " << i << " = " << x0[i] << "  " << y0[i] << "  " << z0[i] << endl;
      }
      face[i] = 7;
      srand(p->seedface+i/numberOfLensPerSnap*5);
      while(face[i]>6 || face[i]<1) face[i] = int(1+rand() / float(RAND_MAX)*5.+0.5);
      if (myid==0)
        cout << " face of the dice " << face[i] << endl;
      sgnX[i] = 2;
      srand(p->seedsign+i/numberOfLensPerSnap*8);
      while(sgnX[i] > 1 || sgnX[i] < 0) sgnX[i] = int(rand() / float(RAND_MAX)+0.5);
      sgnY[i] = 2;
      while(sgnY[i] > 1 || sgnY[i] < 0) sgnY[i] = int(rand() / float(RAND_MAX)+0.5);
      sgnZ[i] = 2;
      while(sgnZ[i] > 1 || sgnZ[i] < 0) sgnZ[i] = int(rand() / float(RAND_MAX)+0.5);
      if(sgnX[i]==0) sgnX[i]=-1;
      if(sgnY[i]==0) sgnY[i]=-1;
      if(sgnZ[i]==0) sgnZ[i]=-1;
      if(myid==0)
        cout << " signs of the coordinates = " << sgnX[i] << "  " << sgnY[i] << " " << sgnZ[i] << endl;

    }
    else{

      x0[i] = x0[i-1];
      y0[i] = y0[i-1];
      z0[i] = z0[i-1];
      if(myid==0){
        cout << "  " << endl;
        cout << " random centers  for the box " << i << " = " << x0[i] << "  " << y0[i] << "  " << z0[i] << endl;
      }
      face[i] = face[i-1];
      if (myid==0)
        cout << " face of the dice " << face[i] << endl;
      sgnX[i] = sgnX[i-1];
      sgnY[i] = sgnY[i-1];
      sgnZ[i] = sgnZ[i-1];
      if(myid==0)
        cout << " signs of the coordinates = " << sgnX[i] << "  " << sgnY[i] << " " << sgnZ[i] << endl;

    }

  }

}

void map_particles (float * xx, double * m, int ptype, int sn_opt, float fovradiants, Header header, vector <double> &ld,
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

int read_header (string file_in, Header *header, ifstream & fin, bool close = false){

  /* Read the Snap Header*/
  fin.open(file_in.c_str());
  if (!fin) {cerr <<"Error in opening the file: "<<file_in<<"!\n\a"; return 1;}

  int32_t blockheader[5];
  fin.read((char *)&blockheader, sizeof(blockheader));
  fin >> *header;

  if(close)
    fin.close();
  return 0;
}

void print_header (Header header){
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
  cout << "  " << endl;
  cout << "      __________________ COSMOLOGY __________________  " << endl;
  cout << " " << endl;
  int dimmass0=0;

  for(int i=0;i<=5;i++){
     if(header.massarr[i]==0){dimmass0+=header.npart[i];}
  }

  if(dimmass0==0){

    cout << "		@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
    cout << "		@                          @" << endl;
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

  }

  cout << " .......................................................... " << endl;
  cout << "   number of particles in this snapshot: " << endl;
  cout << header.npart[0] << " " << header.npart[1] << " " << header.npart[2]
   << " " << header.npart[3] << " " << header.npart[4] << " " << header.npart[5] << endl;

  cout << "      Omegam = " << header.om0 << " " << "Omegal = " << header.oml << endl;
  cout << "           h = " << header.h   << " " << "BoxSize = " << header.boxsize << endl;

  cout << "      _______________________________________________  " << endl;
  cout << " " << endl;
  cout << "   total number of particles in the simulation: " << endl;
  cout << header.nTotalHW[0]*pow(2,32)+header.npartTotal[0] << " " << header.nTotalHW[1]*pow(2,32)+header.npartTotal[1] << " " <<
          header.nTotalHW[2]*pow(2,32)+header.npartTotal[2] << " " << header.nTotalHW[3]*pow(2,32)+header.npartTotal[3] << " " <<
          header.nTotalHW[4]*pow(2,32)+header.npartTotal[4] << " " << header.nTotalHW[5]*pow(2,32)+header.npartTotal[5] << endl;
  cout << " " << endl;
  cout << "   xparticle type mass array: " << endl;
  cout << header.massarr[0] << " " << header.massarr[1] << " " << header.massarr[2]
      << " " << header.massarr[3] << " " <<  header.massarr[4] << " " <<  header.massarr[5] << endl;


}

void fastforwardToPos (ifstream & fin, int NBLOCKS, int myid,  bool close){

  Block block;
  fin >> block;
  for (int i=0; i < NBLOCKS; i++){

    if (i>0){

      fin >> block;

    }

    if (myid==0){
      cout << "Fast Fowarding next block. Name: ";
      cout << block.name[0];
      cout << block.name[1];
      cout << block.name[2];
      cout << block.name[3] << endl;
    }

      fin.seekg(block.blocksize2,ios_base::cur);

  }

  fin >> block;
  if (myid==0){
    cout << "reading next block. Name: ";
    cout << block.name[0];
    cout << block.name[1];
    cout << block.name[2];
    cout << block.name[3] << endl;
    cout << "Should be                 POS " << endl;
  }

  if(close)
    fin.close();

}

void fastforwardToMASS (ifstream & fin, int NBLOCKS, Header *data, int myid){

   int tot = (data->npart[0]+data->npart[1]+data->npart[2]+data->npart[3]+data->npart[4]+data->npart[5]);
   Block block;

   for(int i=0;i<NBLOCKS-1;i++){

     fin >> block;
     if (myid==0){
        cout << "Size of Next Block is " << block.blocksize1 << endl;
        cout << "Should be             " << 3*sizeof(int32_t)*tot << endl;
        cout << "Fast Fowarding next block. Name: ";
        cout << block.name[0];
        cout << block.name[1];
        cout << block.name[2];
        cout << block.name[3] << endl;
     }

     fin.seekg(block.blocksize2/sizeof(int8_t),fin.cur);

   }

   fin >> block;
   if (myid==0){
     cout << "Reading next block. Name: ";
     cout << block.name[0];
     cout << block.name[1];
     cout << block.name[2];
     cout << block.name[3] << endl;
     cout << "Should be                 MASS" << endl;
   }

}

void fastforwardToBHMASS (ifstream & fin, int NBLOCKS, Header *data, int myid){

  Block block;
  fin.seekg(data->npart[5]*sizeof(int32_t)/sizeof(int8_t),fin.cur);
  if (myid==0){
    cout << "\n" <<endl;
    cout << "Fast Fowarding blocks: ";
  }
  for(int i=0;i<NBLOCKS;i++){

    fin >> block;
    if (myid==0){
        cout << block.name[0];
        cout << block.name[1];
        cout << block.name[2];
        cout << block.name[3] << ", ";
    }

    fin.seekg(block.blocksize2/sizeof(int8_t),fin.cur);
  }
  if (myid==0)
    cout<< "\n" << endl;
  fin >> block;
  if (myid==0){
      cout << "Reading BH masses from. Block: ";
      cout << block.name[0];
      cout << block.name[1];
      cout << block.name[2];
      cout << block.name[3] << endl;
  }
}



void read_pos (ifstream *fin, Header *header, float * xx0, float * xx1, float * xx2, float * xx3, float * xx4, float * xx5){

  Block block;
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

void read_mass (ifstream *fin, Header *header, double * m0, double * m1, double * m2, double * m3, double * m4, double * m5){

  Block block;
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
