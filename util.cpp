#include "util.h"
#define MAX_M 1e3    // Threshold for mass; Particler heavier than MAX_M will be attached zero mass
#define POS_U 1000.0 // Unit conversion from BoxSize unit lengh to kpc/h
#define DO_NGP false // Use NGP as the MAS

using namespace CCfits;

int readInput(struct InputParams *p, string name, string & snpix, bool *physical){
  // read fileinput
  std:: ifstream fin (name.c_str());
  if(fin.is_open());
  else{
    cerr << " Params file "<< name <<" does not exist where you are running the code " << endl;
    cerr << " I will STOP here!!! " << endl;
    exit(1);
  }
  std:: string str;
  std::getline(fin, str);
  std::getline(fin, str);
  p->npix = std::stoi(str);//         1. Number of Pixels
  std::getline(fin, str);
  std::getline(fin, str);
  p->zs = std::stof(str);//           2. Redshift Source
  std::getline(fin, str);
  std::getline(fin, str);
  p->fov = std::stof(str);//          3. Field of View
  std::getline(fin, str);
  std::getline(fin,
      p->filredshiftlist);//          4. File with the redshift list it may contain three columns: snap 1/(1+z) z
  std::getline(fin, str);
  std::getline(fin,
      p->pathsnap);//                 5. Path where the snaphosts are located
  std::getline(fin, str);
  std::getline(fin,
      p->simulation);//               6. Simulation name
  std::getline(fin, str);
  std::getline(fin, str);
  p->seedcenter = std::stoi(str);//   7. Seed center
  std::getline(fin, str);
  std::getline(fin, str);
  p->seedface = std::stoi(str);//     8. Seed Face
  std::getline(fin, str);
  std::getline(fin, str);
  p->seedsign = std::stoi(str);//     9. Seed sign
  std::getline(fin, str);
  std::getline(fin, str);
  p->partinplanes = std::stoi(str);// 10. True: Each gadget particle type will have its own Map; False: One Map for All
  std::getline(fin, str);
  std::getline(fin,
       p->directory);//               11. Directory to save FITS files
  std::getline(fin, str);
  std::getline(fin,
       p->suffix);//                  12. Suffix to FITS files
  std::getline(fin, str);
  std::getline(fin, str);
  p->snopt = std::stoi(str);//        13. Shot-noise option: 0-No random Degradation; 1-Half particles degradation ...

  if(p->npix==0)
    p->simType = "SubFind";
  else
    p->simType = "Gadget";

  *physical = (p->npix<0);
  if(! *physical)
    snpix=sconv(p->npix,fINT);
  else{
    int n = -p->npix;
    snpix=sconv(n,fINT);
    snpix+="_kpc";
    p->rgrid=n;
  }

  if(p->snopt<0){
    cerr << "Impossible value for Shot-Noise option!" << endl;
    return 1;
  }

  return 0;

}

int read_redlist(string filredshiftlist, vector <double> & snapred, vector <string> & snappath, InputParams *p){

  ifstream redlist;
  redlist.open(filredshiftlist.c_str());
  double zlast=-999.9;

  if(redlist.is_open()){
    ifstream fin; // snap file
    string name;  // snapshot directory
    double z;     // snapshot redshift
    Header header;// snapshot header
    do{
      redlist >> name;
      snappath.push_back(name);
      if( read_header( p->pathsnap+name+".0" , &header, fin, true) ){
        cerr << name << " not found!" << endl;
        return 1;
      }
      if(header.redshift<zlast){
        cerr << " Snapshots on "<< filredshiftlist << " are not sorted!" << endl;
        return 1;
      }
      else
        zlast = header.redshift;
      snapred.push_back(header.redshift);
    }while(header.redshift<p->zs);
  }else{
    cerr << " redshift list file redshift_list.txt does not " << endl;
    cerr << " exist in the Code dir ... check this out      " << endl;
    cerr << "    I will STOP here !!! " << endl;
    return 1;
  }
  return 0;
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

int getSnap (vector <double> & zsnap, gsl_spline *GetDl, gsl_interp_accel *accGetDl, double dlens){

  int pos;
  double aux=99999;

  for(int i=0;i<zsnap.size();i++){
    float test = abs( gsl_spline_eval (GetDl, zsnap[i], accGetDl)-dlens );
    if( test < aux ){
      aux = test;
      pos = i;
    }
  }

  return pos;
}

void build_planes(InputParams *p, Header *header, Lens &lens, vector <double> & snapred, vector <string> & snappath, gsl_spline *GetDl,
                 gsl_interp_accel *accGetDl, gsl_spline *GetZl, gsl_interp_accel *accGetZl, int numberOfLensPerSnap,
                                                                                                           int myid){
  size_t nsnaps = snapred.size();
  int pos, nrepi = 0;
  int nrep = 0;
  double ldbut, zdbut;
  do{
    nrep++;
    nrepi++;
    ldbut = nrep*header->boxsize/1e3/numberOfLensPerSnap;
    zdbut = gsl_spline_eval (GetZl, ldbut, accGetZl);
    double dlens = ldbut-0.5*header->boxsize/1e3/numberOfLensPerSnap;
    double zlens = gsl_spline_eval (GetZl, dlens, accGetZl);
    int pos_temp = getSnap(snapred, GetDl, accGetDl, dlens);
    if (myid==0)
      cout << " simulation snapshots = " << ldbut << "  " << zdbut << "  " << nrep << " from snap "
                                                       << snappath[pos_temp] << "  " << zlens << endl;
    lens.ld.push_back(ldbut-header->boxsize/1e3/numberOfLensPerSnap);
    lens.ld2.push_back(ldbut);
    lens.zfromsnap.push_back(snapred[pos_temp]);
    if ( nrep != 1 && pos_temp != pos){
      for ( int i=0; i<nrepi-1; i++ )
        lens.replication.push_back(nrep-1);
      nrepi=1;
    }
    pos=pos_temp;
    lens.zsimlens.push_back(zlens);
    lens.fromsnap.push_back(snappath[pos]);
    if( nrep==1 )
      lens.randomize.push_back(1);
    else
      lens.randomize.push_back( !( (nrep-1)%numberOfLensPerSnap ) );

  }while(ldbut<p->Ds);

  for ( int i=0; i<nrepi+1; i++ ) lens.replication.push_back(nrep); // Last plane replications
  if (myid==0){
    cout << " Comoving Distance of the last plane " << p->Ds << endl;
    cout << " nsnaps = " << nsnaps << "\n" << endl;
  }

  ofstream planelist;
  string planes_list;
  planes_list = p->directory+"planes_list_"+p->suffix+".txt";

  if (myid==0)
    planelist.open(planes_list.c_str());

  for(int i=0;i<lens.fromsnap.size();i++){

    if (myid==0){
      cout << lens.zsimlens[i] << " planes = " << lens.ld[i] << "  " << lens.ld2[i] << "  " << lens.replication[i] <<
                                                                     " from snap " << lens.fromsnap[i] << endl;
      planelist <<  i << "   " << lens.zsimlens[i] << "   " << lens.ld[i] << "   " << lens.ld2[i] << "   " <<
          lens.replication[i] << "   " << lens.fromsnap[i]  << "   " << lens.zfromsnap[i] << "  "
          << lens.randomize[i]  << endl;
    }

    lens.pll.push_back(i);

  }

  if (myid==0)
    planelist.close();

  lens.nplanes = lens.replication.back();
}

void randomize_box (Random & random, Lens *lens, InputParams *p, int numberOfLensPerSnap, int myid){

  size_t nrandom = lens->replication.back();
  /* Inflating Random Structure */
  random.x0.resize(nrandom); random.y0.resize(nrandom); random.z0.resize(nrandom);
  random.sgnX.resize(nrandom); random.sgnY.resize(nrandom); random.sgnZ.resize(nrandom);
  random.face.resize(nrandom);
  for(int i=0;i<nrandom;i++){

    if ( lens->randomize[i] ){

      srand(p->seedcenter+i/numberOfLensPerSnap*13);
      random.x0[i] = rand() / float(RAND_MAX);
      random.y0[i] = rand() / float(RAND_MAX);
      random.z0[i] = rand() / float(RAND_MAX);
      if(myid==0){
        cout << "  " << endl;
        cout << " random centers  for the box " << i << " = " << random.x0[i] << "  " << random.y0[i] << "  " << random.z0[i] << endl;
      }
      random.face[i] = 7;
      srand(p->seedface+i/numberOfLensPerSnap*5);
      while(random.face[i]>6 || random.face[i]<1) random.face[i] = int(1+rand() / float(RAND_MAX)*5.+0.5);
      if (myid==0)
        cout << " face of the dice " << random.face[i] << endl;
      random.sgnX[i] = 2;
      srand(p->seedsign+i/numberOfLensPerSnap*8);
      while(random.sgnX[i] > 1 || random.sgnX[i] < 0) random.sgnX[i] = int(rand() / float(RAND_MAX)+0.5);
      random.sgnY[i] = 2;
      while(random.sgnY[i] > 1 || random.sgnY[i] < 0) random.sgnY[i] = int(rand() / float(RAND_MAX)+0.5);
      random.sgnZ[i] = 2;
      while(random.sgnZ[i] > 1 || random.sgnZ[i] < 0) random.sgnZ[i] = int(rand() / float(RAND_MAX)+0.5);
      if(random.sgnX[i]==0) random.sgnX[i]=-1;
      if(random.sgnY[i]==0) random.sgnY[i]=-1;
      if(random.sgnZ[i]==0) random.sgnZ[i]=-1;
      if(myid==0)
        cout << " signs of the coordinates = " << random.sgnX[i] << "  " << random.sgnY[i] << " " << random.sgnZ[i] << endl;

    }
    else{

      random.x0[i] = random.x0[i-1];
      random.y0[i] = random.y0[i-1];
      random.z0[i] = random.z0[i-1];
      if(myid==0){
        cout << "  " << endl;
        cout << " random centers  for the box " << i << " = " << random.x0[i] << "  " << random.y0[i] << "  " << random.z0[i] << endl;
      }
      random.face[i] = random.face[i-1];
      if (myid==0)
        cout << " face of the dice " << random.face[i] << endl;
      random.sgnX[i] = random.sgnX[i-1];
      random.sgnY[i] = random.sgnY[i-1];
      random.sgnZ[i] = random.sgnZ[i-1];
      if(myid==0)
        cout << " signs of the coordinates = " << random.sgnX[i] << "  " << random.sgnY[i] << " " << random.sgnZ[i] << endl;

    }

  }

}

int test_fov(double fov, double boxl, double Ds, int myid, double *fovradiants){
  if(myid==0)
    cout << " Setting the field of view to be square in degrees " << endl;

  *fovradiants = fov/180.*M_PI;
  /* check if the field of view is too large with respect to the box size */
  if( (*fovradiants)*Ds>boxl && myid==0 ){
    cerr << " !!Field view too large!!\n !!!I will STOP here!!! " << endl;
    cerr << " Value set is = " << fov << endl;
    cerr << " Maximum value allowed " << boxl/Ds*180./M_PI << " in degrees " << endl;
    return 1;
  }
  return 0;
}

void test_hydro(InputParams *p, Header *data){

  if( p->simType.compare("Gadget")==0 ){

    int dim = data->npart[0]+data->npart[1]+data->npart[2]+data->npart[3]+data->npart[4]+data->npart[5];
    int dimmass0=0;
    for(int i=0;i<=5;i++){
       if(data->massarr[i]==0)
         dimmass0+=data->npart[i];
    }
    p->hydro=bool(dimmass0);
  }

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
  cout << "Box size: "<< header.boxsize/1e3<< endl;
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
  cout << "           h = " << header.h   << " " << "BoxSize = " << header.boxsize/1e3 << endl;

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

void fastforwardNVars (ifstream & fin, size_t size, size_t n, int myid){
  fin.seekg(n*size/sizeof(int8_t),fin.cur);
}

void fastforwardToBlock (ifstream & fin, string BLOCK, int myid){

  Block block;
  char nullchar = '\0';
  char blockStringLike [5];
  blockStringLike[4]=nullchar;
  bool iterate = true;

  while(iterate){

    fin >> block;
    for(int i=0; i<4; i++)
      blockStringLike[i] = block.name[i];

    iterate = strcmp(blockStringLike, &BLOCK[0]);
    if(iterate){
      if (myid==0)
        cout << "Fast Fowarding next block. Name: "<< blockStringLike << endl;
        fin.seekg(block.blocksize2,ios_base::cur);
      }

  }

  for(int i=0; i<4; i++)
    blockStringLike[i] = block.name[i];
  if (myid==0){
    cout << "reading next block. Name: "<< blockStringLike << endl;
    cout << "Should be                 "<< BLOCK << endl;
  }

}

void ReadPos (ifstream & fin, Header *data, InputParams *p, Random *random,
                              int isnap, float* xx[6][3], float rcase,int myid){

  float num_float1,num_float2, num_float3; // Dummy vars to read x,y, and z
  int imax = (p->simType == "Gadget") ? 6 : 2;
  if(p->simType == "Gadget"){
    fastforwardToBlock (fin, "POS ", myid);
  }else{
    fastforwardToBlock (fin, "GPOS", myid);
  }
  /* Loop on different types */
  for (int i = 0; i<imax; i++){

    if(p->simType=="SubFind" & i==1)
      fastforwardToBlock (fin, "SPOS", myid);

    for (int pp=0; pp<data->npart[i]; pp++){

      float x, y, z;
      float xb, yb, zb;

      fin.read((char *)&num_float1, sizeof(num_float1));
      fin.read((char *)&num_float2, sizeof(num_float2));
      fin.read((char *)&num_float3, sizeof(num_float3));

      xb = random->sgnX[isnap]*(((num_float1)/data->boxsize));
      yb = random->sgnY[isnap]*(((num_float2)/data->boxsize));
      zb = random->sgnZ[isnap]*(((num_float3)/data->boxsize));

      // wrapping periodic condition
      if(xb>1.) xb = xb - 1.;
      if(yb>1.) yb = yb - 1.;
      if(zb>1.) zb = zb - 1.;
      if(xb<0.) xb = 1. + xb;
      if(yb<0.) yb = 1. + yb;
      if(zb<0.) zb = 1. + zb;
      switch (random->face[isnap]){
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
      x = x - random->x0[isnap];
      y = y - random->y0[isnap];
      z = z - random->z0[isnap];
      // wrapping periodic condition again
      if(x>1.) x = x - 1.;
      if(y>1.) y = y - 1.;
      if(z>1.) z = z - 1.;
      if(x<0.) x = 1. + x;
      if(y<0.) y = 1. + y;
      if(z<0.) z = 1. + z;
      z+=double(rcase); // pile the cones
      xx[i][0][pp] = x;
      xx[i][1][pp] = y;
      xx[i][2][pp] = z;
    }
  }
}

void ReadVel (ifstream & fin, Header *data, InputParams *p, Random *random,
                              int isnap, float* vv[6][3], int myid){

  float num_float1,num_float2, num_float3; // Dummy vars to read x,y, and z velocities
  int imax = (p->simType == "Gadget") ? 6 : 2;
  int imin = (p->simType == "Gadget") ? 0 : 1;
  if(p->simType == "Gadget")
    fastforwardToBlock (fin, "VEL ", myid);
  else
    fastforwardToBlock (fin, "SVEL", myid);
  /* Loop on different types */
  for (int i = imin; i<imax; i++){

    for (int pp=0; pp<data->npart[i]; pp++){

      float x, y, z;
      float xb, yb, zb;

      fin.read((char *)&num_float1, sizeof(num_float1));
      fin.read((char *)&num_float2, sizeof(num_float2));
      fin.read((char *)&num_float3, sizeof(num_float3));

      xb = random->sgnX[isnap]*num_float1;
      yb = random->sgnY[isnap]*num_float2;
      zb = random->sgnZ[isnap]*num_float3;

      // wrapping periodic condition
      switch (random->face[isnap]){
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
      vv[i-imin][0][pp] = x;
      vv[i-imin][1][pp] = y;
      vv[i-imin][2][pp] = z;
    }
  }
}

int MapParticles(ifstream & fin, Header *data, InputParams *p, Lens *lens,
    float* xx[6][3], double fovradiants, int isnap, valarray<float>(& mapxyi)[6],
    int(& ntotxyi)[6], int myid){

  int imax = 6;
  float num_float1;
  int totPartxyi[6];
  for(int i=0; i<imax; i++){

    size_t n = data->npart[i];

    if(n>0){
      // quadrate box
      double xmin=double(*min_element(xx[i][0], &xx[i][0][n]));
      double xmax=double(*max_element(xx[i][0], &xx[i][0][n]));
      double ymin=double(*min_element(xx[i][1], &xx[i][1][n]));
      double ymax=double(*max_element(xx[i][1], &xx[i][1][n]));
      double zmin=double(*min_element(xx[i][2], &xx[i][2][n]));
      double zmax=double(*max_element(xx[i][2], &xx[i][2][n]));

      if (myid==0){
       cout << " " << endl;
       cout << " n"<<i<<" particles " << endl;
       cout << "xmin = " << xmin << endl;
       cout << "xmax = " << xmax << endl;
       cout << "ymin = " << ymin << endl;
       cout << "ymax = " << ymax << endl;
       cout << "zmin = " << zmin << endl;
       cout << "zmax = " << zmax << endl;
       cout << "  " << endl;
      }

      if(xmin<0 || ymin<0 || zmin< 0){
        cerr << "xmin = " << xmin << endl;
        cerr << "xmax = " << xmax << endl;
        cerr << "ymin = " << ymin << endl;
        cerr << "ymax = " << ymax << endl;
        cerr << "zmin = " << zmin << endl;
        cerr << "zmax = " << zmax << endl;
        cerr << "  0 type check this!!! I will STOP here!!! " << endl;
        cerr << "Aborting from Rank "<< myid << endl;
        return 1;
      }
      if (myid==0){
         cout << " Mapping type "<< i <<" particles on the grid with " << p->npix << " pixels" << endl;
         cout << "Min distance: "<< lens->ld[isnap]/data->boxsize*1.e+3/POS_U<< " "<<lens->ld2[isnap]/data->boxsize*1.e+3/POS_U << endl;
      }

      vector<float> xs(0),ys(0),ms(0);
      for(int l=0;l<data->npart[i];l++){

        if(p->hydro && data->massarr[i]==0){

          if(l==0 && i==5){
            fastforwardNVars (fin, sizeof(int32_t), data->npart[5],myid);
            fastforwardToBlock (fin, "BHMA", myid);
          }

          fin.read((char *)&num_float1, sizeof(num_float1));
          if (num_float1>MAX_M)
            num_float1=0;
        }
        else
          num_float1=data->massarr[i];

        double di = sqrt(pow(xx[i][0][l]-0.5,2) + pow(xx[i][1][l]-0.5,2) + pow(xx[i][2][l],2))*data->boxsize/1.e+3*POS_U;
        if(di>=lens->ld[isnap] && di<lens->ld2[isnap]){

          double rai,deci,dd;
          getPolar(xx[i][0][l]-0.5,xx[i][1][l]-0.5,xx[i][2][l],&rai,&deci,&dd);

          if(fabs(rai)<=fovradiants*(1.+2./p->npix)*0.5 && fabs(deci)<=fovradiants*(1.+2./p->npix)*0.5){
            xs.push_back(deci/fovradiants+0.5);
            ys.push_back(rai/fovradiants+0.5);
            if(p->snopt==0){
              ms.push_back(num_float1);
            }
            else{
              if(rand()/ float(RAND_MAX) < 1./pow(2,p->snopt)) ms.push_back(pow(2,p->snopt)*num_float1);
              else ms.push_back(0.);
            }
          }
        }
      }
      totPartxyi[i]=xs.size();
      ntotxyi[i]+=totPartxyi[i];

      if(totPartxyi[i]>0){
        mapxyi[i] = gridist_w(xs,ys,ms,p->npix,DO_NGP);
      }

      //re-normalize to the total mass!
      double mtot=0.;
      double mnorm=0.;
      for(int i=0;i<ms.size();i++)
        mnorm+=ms[i];

      if(totPartxyi[i]>0){
        for(int l=0;l<p->npix*p->npix;l++)
          mtot += mapxyi[i][l];
        if(mtot==0.)
          mtot=1.; //To avoid NaN
        for(int l=0;l<p->npix*p->npix;l++) mapxyi[i][l]*=mnorm/mtot;
      }
    }
  }

  return 0;

}

void getPolar(double x, double y, double z, double *ra, double *dec, double *d){
  *d = sqrt(x*x+y*y+z*z);
  *dec = acos(z/(*d));
  *ra = atan2(y,x);
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

void write_maps (InputParams *p, Header *data, Lens *lens, int isnap, double zsim,
                 string snappl, string snpix, valarray<float>& mapxytotrecv,
                 valarray<float>(& mapxytotirecv)[6], int(& ntotxyi)[6], int myid){

  if (myid==0){
    if(p->partinplanes==false){
     /*
      * write image array(s) to FITS files all particles in a FITS file!
      */
      long naxis = 2;
      long naxes[2]={ p->npix,p->npix };
      string fileoutput;
      fileoutput = p->directory+p->simulation+"."+snappl+".plane_"+snpix+"_"+p->suffix+".fits";
      cout << "Saving the maps on: " << fileoutput << endl;
      unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
      vector<long> naxex( 2 );
      naxex[0]=p->npix;
      naxex[1]=p->npix;
      PHDU *phxy=&ffxy->pHDU();
      phxy->write( 1, p->npix*p->npix, mapxytotrecv );
      phxy->addKey ("REDSHIFT",zsim," ");
      phxy->addKey ("PHYSICALSIZE",p->fov," ");
      phxy->addKey ("PIXELUNIT",1.e+10/data->h,"Mass unit in M_Sun");
      phxy->addKey ("DlLOW",lens->ld[isnap]/data->h,"comoving distance in Mpc");
      phxy->addKey ("DlUP",lens->ld2[isnap]/data->h,"comoving distance in Mpc");
      phxy->addKey ("nparttype0",ntotxyi[0]," ");
      phxy->addKey ("nparttype1",ntotxyi[1]," ");
      phxy->addKey ("nparttype2",ntotxyi[2]," ");
      phxy->addKey ("nparttype3",ntotxyi[3]," ");
      phxy->addKey ("nparttype4",ntotxyi[4]," ");
      phxy->addKey ("nparttype5",ntotxyi[5]," ");
      phxy->addKey ("HUBBLE",data->h," ");
      phxy->addKey ("OMEGAMATTER",data->om0," ");
      phxy->addKey ("OMEGALAMBDA",data->oml," ");
      phxy->addKey ("m0",data->massarr[0]," ");
      phxy->addKey ("m1",data->massarr[1]," ");
      phxy->addKey ("m2",data->massarr[2]," ");
      phxy->addKey ("m3",data->massarr[3]," ");
      phxy->addKey ("m4",data->massarr[4]," ");
      phxy->addKey ("m5",data->massarr[5]," ");
    }else{
    /**
    * write image array(s) to FITS files each particle type in different planes
    */
     for(int i=0; i<6; i++){

       if(ntotxyi[i]>0){
         long naxis = 2;
         long naxes[2]={ p->npix,p->npix };
         string fileoutput;
         fileoutput = p->directory+p->simulation+"."+snappl+".ptype"+sconv(i,fINT)+"_plane_"+snpix+"_"+p->suffix+".fits";
         unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
         vector<long> naxex( 2 );
         naxex[0]=p->npix;
         naxex[1]=p->npix;
         PHDU *phxy=&ffxy->pHDU();
         phxy->write( 1, p->npix*p->npix, mapxytotirecv[i] );
         phxy->addKey ("REDSHIFT",zsim," ");
         phxy->addKey ("PHYSICALSIZE",p->fov," ");
         phxy->addKey ("PIXELUNIT",1.e+10/data->h,"Mass unit in M_Sun");
         phxy->addKey ("DlLOW",lens->ld[isnap]/data->h,"comoving distance in Mpc");
         phxy->addKey ("DlUP",lens->ld2[isnap]/data->h,"comoving distance in Mpc");
         phxy->addKey ("nparttype0",ntotxyi[i]," ");
         phxy->addKey ("HUBBLE",data->h," ");
         phxy->addKey ("OMEGAMATTER",data->om0," ");
         phxy->addKey ("OMEGALAMBDA",data->oml," ");
         phxy->addKey ("m"+sconv(i,fINT),data->massarr[i]," ");
       }
     }
    }
  }
}

SubFind::SubFind(int n){
  this->id.resize(n);
  this->truez.resize(n);
  this->xx0.resize(n); this->yy0.resize(n); this->zz0.resize(n);
  this->vx0.resize(n); this->vy0.resize(n); this->vz0.resize(n);
  this->m.resize(n);
  this->theta.resize(n);
  this->phi.resize(n);
  this->vel.resize(n);
  this->obsz.resize(n);
  this->nsub.resize(n);
};

void GetGVel(SubFind &halos, SubFind *subhalos, InputParams *p, Random *random, string FILE, int ff, int isnap){

  size_t nhalos = halos.m.size();
  vector <unsigned int>::iterator it;
  vector <unsigned int> id;
  vector <float> vx0, vy0, vz0;

  for(int i=0; i<nhalos; i++){

    it = find( subhalos->id.begin(), subhalos->id.end(), halos.id[i]);
    if( it  == subhalos->id.end()){

      Header data;

      ifstream fin;
      float *vv[6][3];
      read_header (FILE+"."+sconv(ff,fINT), &data, fin, false);
      id.resize(data.npart[1]);
      vx0.resize(data.npart[1]); vy0.resize(data.npart[1]); vz0.resize(data.npart[1]);
      vv[0][0]=&vx0[0]; vv[0][1]=&vy0[0]; vv[0][2]=&vz0[0];
      for(int i=1; i<6;i++){
        for(int j=0;j<3;j++){
          vv[i][j]=nullptr;
        }
      }
      ReadVel (fin, &data, p, random, isnap, vv, 999);
      ReadBlock(fin, data.npart[1], "GRNR", &id[0], 999);

      it = find( id.begin(), id.end(), halos.id[i]);
      int isub = it-id.begin();
      halos.vx0[i] = vx0[isub];
      halos.vy0[i] = vy0[isub];
      halos.vz0[i] = vz0[isub];

    }else{
      int isub = it-subhalos->id.begin();
      halos.vx0[i] = subhalos->vx0[isub];
      halos.vy0[i] = subhalos->vy0[isub];
      halos.vz0[i] = subhalos->vz0[isub];
    }
  }
}

int GetGID(SubFind &halos, string File, int ff){

  int nhalos=0;
  for(int i=0; i<ff; i++){
    Header data;
    ifstream fin;
    if(read_header (File+"."+sconv(i,fINT), &data, fin, true))
      return 1;
    nhalos += data.npart[0];
    fin.clear();
    fin.close();
  }

  int nlocal=halos.m.size();
  for(int i=nhalos; i<nhalos+nlocal; i++){
    halos.id[i-nhalos]=i;
  }
  return 0;
}

void GetLOSVel(SubFind &halos){

  int nhalos=halos.m.size();
  double x,y,z,r;

  for(int i=0; i<nhalos; i++){

    x = halos.xx0[i]-0.5;
    y = halos.yy0[i]-0.5;
    z = halos.zz0[i];
    r = sqrt(x*x + y*y + z*z);
    x /= r; y /= r; z /= r;

    halos.vel[i] = halos.vx0[i]*x + halos.vy0[i]*y + halos.vz0[i]*z;
    halos.obsz[i] = halos.truez[i] + halos.vel[i]/speedcunit/100.0;

  }

}

void GetTrueZ(SubFind &halos, Header *data, gsl_spline *GetZl, gsl_interp_accel *accGetZl){

  int nhalos=halos.m.size();
  double x,y,z,r;

  for(int i=0; i<nhalos; i++){

    x = halos.xx0[i]-0.5;
    y = halos.yy0[i]-0.5;
    z = halos.zz0[i];
    r = sqrt(x*x + y*y + z*z);

    halos.truez[i]=gsl_spline_eval (GetZl, r*data->boxsize/1e3, accGetZl);

  }

}

void GetAngular(SubFind &halos){

  int nhalos=halos.m.size();
  double x,y,z;
  double r;
  double phi;
  double theta;

  for(int i=0; i<nhalos; i++){

    x = halos.xx0[i]-0.5;
    y = halos.yy0[i]-0.5;
    z = halos.zz0[i];

    getPolar(x, y, z, &phi, &theta , &r);
    halos.phi[i] = phi;
    halos.theta[i] = theta;

  }

}

void CreatePLC (SubFind &halos, Header *data, InputParams *p, string snappl, int ff){

  fstream fileoutput ( p->directory+p->simulation+"."+snappl+".plane_"+sconv(p->fov,fDP1)+"_"+p->suffix+".plc."+sconv(ff,fINT), ios::out | ios::binary );
  double fovradiants = p->fov/180.*M_PI;

  int nhalos = halos.m.size();

  for(int i=0; i<nhalos; i++){

    if( (abs(halos.phi[i]) <=  fovradiants/2.0) & (abs(halos.theta[i]) <=  fovradiants/2.0) & (halos.m[i]>0) ){
      /* Explicitly casting the variables as in Pinocchio PLC */
      int dummy;
      long long unsigned int id = halos.id[i];
      double xx0 = (halos.xx0[i] - 0.5)*data->boxsize/1e3;
      double yy0 = (halos.yy0[i] - 0.5)*data->boxsize/1e3;
      double zz0 = halos.zz0[i]*data->boxsize/1e3;
      double vx0 = halos.vx0[i];
      double vy0 = halos.vy0[i];
      double vz0 = halos.vz0[i];
      double m = halos.m[i];
      double obsz = halos.obsz[i];
      double truez = halos.truez[i];
      double vel = halos.vel[i];
      double theta = 90.0-halos.theta[i]*180.0/M_PI;
      double phi = halos.phi[i]*180.0/M_PI;
      if( phi<0.0 )
        phi += 360.0;

      fileoutput.write((char*)&dummy, sizeof (int));
      fileoutput.write((char*)&id, sizeof (unsigned long long int));
      fileoutput.write((char*)&truez, sizeof (double));
      fileoutput.write((char*)&xx0, sizeof (double));
      fileoutput.write((char*)&yy0, sizeof (double));
      fileoutput.write((char*)&zz0, sizeof (double));
      fileoutput.write((char*)&vx0, sizeof (double));
      fileoutput.write((char*)&vy0, sizeof (double));
      fileoutput.write((char*)&vz0, sizeof (double));
      fileoutput.write((char*)&m, sizeof (double));
      fileoutput.write((char*)&theta, sizeof (double));
      fileoutput.write((char*)&phi, sizeof (double));
      fileoutput.write((char*)&vel, sizeof (double));
      fileoutput.write((char*)&obsz, sizeof (double));
      fileoutput.write((char*)&dummy, sizeof (int));

    }

  }

  fileoutput.close();

}
