/* Routines to read Gadget2 snapshot format */
#include "gadget2io.h"

/* Reads the "file_in" snapshot and stores it on the "header" instance
   of Header. "fin" is a ifstream instance that is leaved open for further
   reading if "close" == false.
*/
int readHeader (string file_in, Header &header, ifstream & fin,
                                              bool close = false){

  /* Read the Snap Header*/
  fin.open(file_in.c_str());
  if (!fin) {cerr <<"Error in opening the file: "<<file_in<<"!\n\a"; return 1;}

  int32_t blockheader[5];
  fin.read((char *)&blockheader, sizeof(blockheader));
  fin >> header;

  if(close)
    fin.close();
  return 0;
}

/* Tests if the snapshot has Hydro particles */
void testHydro(InputParams & p, Header & data){

  if( p.simType.compare("Gadget")==0 ){

    int dim = data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]+data.npart[4]+data.npart[5];
    int dimmass0=0;
    for(int i=0;i<=5;i++){
       if(data.massarr[i]==0)
         dimmass0+=data.npart[i];
    }
    p.hydro=bool(dimmass0);
  }

}

/* Prints the Snapshot Header */
void printHeader (Header header){
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

/* Fastforward "n" variables of size "size" of the stream "fin". */
void fastforwardNVars (ifstream & fin, size_t size, size_t n){
  fin.seekg(n*size/sizeof(int8_t),fin.cur);
}

/* Fastforward the fstream "fin" to variable labled "BLOCK".
   Prints "fin" metadata if myid == 0.
*/
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

/* Reads the position block of snapshot fin according to Header "data"
   and InputParams "p". The positions are stored in the address xx
   according to the randomization plane "random" and the box replication
   rcase.

   Produces monitoration messages if myid == 0.
*/
void readPos (ifstream & fin, Header &data, InputParams &p, Random &random,
                              int isnap, float* xx[6][3], float rcase,int myid){

  float num_float1,num_float2, num_float3; // Dummy vars to read x,y, and z
  int imax = (p.simType == "Gadget") ? 6 : 2;
  if(p.simType == "Gadget"){
    fastforwardToBlock (fin, "POS ", myid);
  }else{
    fastforwardToBlock (fin, "GPOS", myid);
  }
  /* Loop on different types */
  for (int i = 0; i<imax; i++){

    if(p.simType=="SubFind" & i==1)
      fastforwardToBlock (fin, "SPOS", myid);

    for (int pp=0; pp<data.npart[i]; pp++){

      float x, y, z;
      float xb, yb, zb;

      fin.read((char *)&num_float1, sizeof(num_float1));
      fin.read((char *)&num_float2, sizeof(num_float2));
      fin.read((char *)&num_float3, sizeof(num_float3));

      xb = random.sgnX[isnap]*(((num_float1)/data.boxsize));
      yb = random.sgnY[isnap]*(((num_float2)/data.boxsize));
      zb = random.sgnZ[isnap]*(((num_float3)/data.boxsize));

      // wrapping periodic condition
      if(xb>1.) xb = xb - 1.;
      if(yb>1.) yb = yb - 1.;
      if(zb>1.) zb = zb - 1.;
      if(xb<0.) xb = 1. + xb;
      if(yb<0.) yb = 1. + yb;
      if(zb<0.) zb = 1. + zb;
      switch (random.face[isnap]){
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
      x = x - random.x0[isnap];
      y = y - random.y0[isnap];
      z = z - random.z0[isnap];
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

/* Reads the velocity block of snapshot fin according to Header "data"
   and InputParams "p". The velocities are stored in the address vv
   according to the randomization plane "random" and the box replication
   rcase.

   Produces monitoration messages if myid == 0.
*/
void readVel (ifstream & fin, Header &data, InputParams &p, Random &random,
                              int isnap, float* vv[6][3], int myid){

  float num_float1,num_float2, num_float3; // Dummy vars to read x,y, and z velocities
  int imax = (p.simType == "Gadget") ? 6 : 2;
  int imin = (p.simType == "Gadget") ? 0 : 1;
  if(p.simType == "Gadget")
    fastforwardToBlock (fin, "VEL ", myid);
  else
    fastforwardToBlock (fin, "SVEL", myid);
  /* Loop on different types */
  for (int i = imin; i<imax; i++){

    for (int pp=0; pp<data.npart[i]; pp++){

      float x, y, z;
      float xb, yb, zb;

      fin.read((char *)&num_float1, sizeof(num_float1));
      fin.read((char *)&num_float2, sizeof(num_float2));
      fin.read((char *)&num_float3, sizeof(num_float3));

      xb = random.sgnX[isnap]*num_float1;
      yb = random.sgnY[isnap]*num_float2;
      zb = random.sgnZ[isnap]*num_float3;

      // wrapping periodic condition
      switch (random.face[isnap]){
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

/* Get the "halos" group velocity from the "subhalos" group. The velocity
   is randomized according to InputParams "p" and Randomizantion Plan
   "random".

   The code searches for the subhalos.id in snapshots "FILE".{0..nfiles}
*/
int getGVel(SubFind &halos, InputParams &p, Random &random, string FILE, int isnap){

  vector<float>::iterator vx0it = halos.vx0.begin();
  vector<float>::iterator vy0it = halos.vy0.begin();
  vector<float>::iterator vz0it = halos.vz0.begin();
  float dummyv;

  int nsubhalos = 0;
  int flast = 0;
  bool finopen = false;
  int lastsub = 0;

  for(vector<uint32_t>::iterator fsubit = halos.fsub.begin(); fsubit != halos.fsub.end(); ++fsubit){

    Header data;
    if(finopen){
      ifstream fin;
      readHeader (FILE+"."+sconv(flast,fINT), data, fin, true);
      fin.clear();
    }

    ifstream fin;
    for(int f = flast; f < data.numfiles; f++){

      if(!finopen){
        readHeader (FILE+"."+sconv(f,fINT), data, fin, false);
        lastsub=0;
      }

      if(data.npart[1] + nsubhalos < *fsubit){
        fin.close();
        fin.clear();
        nsubhalos += data.npart[1];
        continue;
      }else{

        if(!lastsub)
          fastforwardToBlock(fin, "SVEL", 0);

        fastforwardNVars(fin, sizeof(double), 3*(*fsubit-nsubhalos-lastsub) );
        lastsub = 3*(*fsubit-nsubhalos-lastsub);
        fin.read((char *)&dummyv, sizeof(dummyv));
        *vx0it=dummyv; vx0it++;
        fin.read((char *)&dummyv, sizeof(dummyv));
        *vy0it=dummyv; vy0it++;
        fin.read((char *)&dummyv, sizeof(dummyv));
        *vz0it=dummyv; vz0it++;

        if(fsubit != halos.fsub.end()){
          if(data.npart[1] + nsubhalos < *(fsubit+1)){
            fin.close();
            fin.clear();
            nsubhalos += data.npart[1];
            finopen=false;
            flast++;
          }else{
          finopen=true;
          }
        }
      }
    }
  }
  if( vx0it != halos.vx0.end() || vy0it != halos.vy0.end() || vz0it != halos.vz0.end() ){

    return 1;

  }else{

    cout << "GVEL read!" << endl;

    for(vx0it = halos.vx0.begin(), vy0it = halos.vy0.begin(), vz0it = halos.vz0.begin();
        vx0it != halos.vx0.end(),  vy0it != halos.vy0.end(),  vz0it != halos.vz0.end();
        ++vx0it, ++vy0it, ++vz0it){

      float x, y, z;
      float xb = random.sgnX[isnap]*(*vx0it);
      float yb = random.sgnY[isnap]*(*vy0it);
      float zb = random.sgnZ[isnap]*(*vz0it);

      // wrapping periodic condition
      switch (random.face[isnap]){
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
      *vx0it = x;
      *vy0it = y;
      *vz0it = z;
    }

    return 0;
  }
}

/* Get the "halos" group ID corresponding according to its position on the
   snapshot "FILE.ff"
*/
int getGID(SubFind &halos, string File, int ffmin, int ffi, int &nhalos){

  for(int i=ffmin; i<ffi; i++){
    Header data;
    ifstream fin;
    if(readHeader (File+"."+sconv(i,fINT), data, fin, true))
      return 1;
    nhalos += data.npart[0];
    fin.clear();
  }

  int nlocal=halos.m.size();
  for(int i=nhalos; i<nhalos+nlocal; i++){
    halos.id[i-nhalos]=i;
  }
  nhalos += nlocal;
  return 0;
}

/* Get the "halos" LineOfSight velocity
*/
void getLOSVel(SubFind &halos){

  int nhalos=halos.m.size();
  double x,y,z,r;

  for(int i=0; i<nhalos; i++){

    x = halos.xx0[i]-0.5;
    y = halos.yy0[i]-0.5;
    z = halos.zz0[i];
    r = sqrt(x*x + y*y + z*z);
    x /= r; y /= r; z /= r;

    if(halos.truez[i]>0){

      halos.vel[i] = halos.vx0[i]*x + halos.vy0[i]*y + halos.vz0[i]*z;
      halos.obsz[i] = halos.truez[i] + halos.vel[i]/speedcunit/100.0 * (1 + halos.truez[i]);

    }else{

      halos.vel[i] =  -999.9;
      halos.obsz[i] = -999.9;
    }

  }

}

/* Get the "halos" cosmological redshift (not taking into account pec.vel)
*/
void getTrueZ(SubFind &halos, Header &data, gsl_spline *getZl, gsl_interp_accel *accGetZl, Lens lens, int isnap){

  int nhalos=halos.m.size();
  double x,y,z,r;

  for(int i=0; i<nhalos; i++){

    x = halos.xx0[i]-0.5;
    y = halos.yy0[i]-0.5;
    z = halos.zz0[i];
    r = sqrt(x*x + y*y + z*z)*data.boxsize/1.e+3*POS_U;

    if(lens.ld[isnap] <= r && r<lens.ld2[isnap]){

      halos.truez[i]=gsl_spline_eval (getZl, r, accGetZl);

    }else{

      halos.truez[i] = -999.9;

    }

  }

}

/* Get the "halos" angular positions.
   Matches Pinocchio cordinate system.
*/
void getAngular(SubFind &halos){

  int nhalos=halos.m.size();
  double x,y,z;
  double r;
  double phi;
  double theta;

  for(int i=0; i<nhalos; i++){

    x = halos.xx0[i]-0.5;
    y = halos.yy0[i]-0.5;
    z = halos.zz0[i];

    getPolar(x, y, z, phi, theta , r, true);
    halos.phi[i] = phi;
    halos.theta[i] = theta;

  }

}

/*
  Reads the snapshots available in "filredshiftlist" up to the first deeper
  than the source in InputParams "p". The path for the snapshots and its
  redshifts are stored in snappath and snapred respectively.
*/
int readRedList(string filredshiftlist, vector <double> & snapred, vector <string> & snappath, vector <double> & snapbox, InputParams &p){

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
      if( readHeader( p.pathsnap+name+".0" , header, fin, true) ){
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
      snapbox.push_back(header.boxsize);
    }while(header.redshift<p.zs);
  }else{
    cerr << " redshift list file redshift_list.txt does not " << endl;
    cerr << " exist in the Code dir ... check this out      " << endl;
    cerr << "    I will STOP here !!! " << endl;
    return 1;
  }
  return 0;
}
