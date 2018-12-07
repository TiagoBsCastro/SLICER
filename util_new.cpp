#include "util_new.h"

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

  *physical = (p->npix<0);
  if(!physical)
    snpix=sconv(p->npix,fINT);
  else{
    int n = -p->npix;
    snpix=sconv(n,fINT);
    snpix+="_kpc";
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
