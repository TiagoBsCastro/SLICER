#include "mpi.h"
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <CCfits/CCfits>
#include <util.h>
#include <cstring>
#define MAX_M 1e3 // Threshold for mass; Particler heavier than MAX_M will be attached zero mass
#define POS_U 1   // Unit conversion from BoxSize unit lengh to kpc/h
#define NBLOCKS 1 // Number of blocks to be fastforwarded
#define numberOfLensPerSnap 2 // Number of Lens to be builded from a snap

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
/*     dev. by Carlo Giocoli - cgiocoli@gmail.com                            */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*    This code has been adapted to run on Magneticum simulations:           */
/*    other functionalities removed for simplicity                           */
/*                                                                           */
/*      - Proper mass assignment to Hydro particles        		               */
/*      - Proper mass assignment to  BH   particles                          */
/*	    - Bug fixed: Proper construction of light-cones in case of few       */
/*                   snapshots                                               */
/*      - Parallelization of Map Building algorithm                          */
/*                                                                           */
/*                        adapted  by Tiago Castro tiagobscastro@gmail.com   */
/*****************************************************************************/

using namespace std;
using namespace CCfits;

istream & operator>>(istream &input, DATA &Data)
{
  input.read((char *)&Data, sizeof(Data));
  return input;
};

int main(int argc, char** argv){

  // MPI Initialization
  int myid, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if(argc!=2 && myid==0){
    cout << "No params!! Nothing to be done!" << endl;
    MPI_Abort(MPI_COMM_WORLD,-1);
  } // Should be run ./exe paramfile

  MPI_Barrier(MPI_COMM_WORLD);

  string inifile=argv[1]; // Paramfile name
  if (myid==0){
    cout << "   ------------------------------------------------------ " << endl;
    cout << "   -                                                    - " << endl;
    cout << "   -           2D Mapping Simulation Snapshot           - " << endl;
    cout << "   -                                                    - " << endl;
    cout << "   -               collapsing one dimension             - " << endl;
    cout << "   ------------------------------------------------------ " << endl;
  }
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
  string simulation; // Simulation Name
  string partinplanes; // ALL all part in one plane, NO each part type in different planes
  string directory, suffix;
  int sn_opt; //Shot-Noise option
  string FilePath;
  bool do_NGP;

  readParameters(inifile,&npix,&boxl,&zs,&fov,
                 &filredshiftlist,&filsnaplist,
                 &pathsnap,&idc,&seedcenter,
                 &seedface,&seedsign,
                 &simulation,&nfiles,
                 &partinplanes,
                 &directory,&suffix,&sn_opt,&do_NGP);  // Reading the parameters File


 if(sn_opt<0 && myid==0){
    cout << "Impossible value for Shot-Noise option!" << endl;
    MPI_Abort(MPI_COMM_WORLD,-1);;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  double Omega_matter,Omega_lambda,Omega_baryon,hubble;
  string snpix;
  bool do_as_t11=(npix<0);
  int rgrid;

  if(!do_as_t11) snpix=conv(npix,fINT);
  else{
    int n = -npix;
    snpix=conv(n,fINT);
    snpix+="_kpc";
    rgrid=n;
  }

  // ... read the redshift list and the snap_available
  // ... to build up the light-cone
  ifstream redlist;
  redlist.open(filredshiftlist.c_str());
  vector <string> lsnap; //lsnap and lred are the selected snapshots selected (z<zs + the first snapshot deeper than zs)
  vector<double> lred;

  if(redlist.is_open()){
    string buta; //snapshot number
    double butb,butc; // snapshots redshift and scale factor
    do{
      redlist >> buta >> butb >> butc;
      if (myid==0)
        cout << butb << " " << zs << endl;
	    lsnap.push_back(buta); // Selected snapshots number and redshift
	    lred.push_back(butb);
    }while(butb<zs);
  }else{
    cout << " redshift list file redshift_list.txt does not " << endl;
    cout << " exist in the Code dir ... check this out      " << endl;
    cout << "    I will STOP here !!! " << endl;
    MPI_Abort(MPI_COMM_WORLD,-1);
  }

  int nsnaps = lsnap.size();
  if (myid==0){
    cout << "  " << endl;
    cout << " opening path for snapshots: " << endl;
    cout << pathsnap << endl;
    cout << " " << endl;
    cout << " I will look for comoving distance " << endl;
    cout << "      file = " << idc << endl;
    cout << " " << endl;
  }
  ifstream infiledc;
  vector<double> zl, dl;
  infiledc.open(idc.c_str());
  if(infiledc.is_open()){
    double zi,dli;
    while(infiledc >> zi >> dli){
      zl.push_back(zi);
      dl.push_back(dli*speedcunit); // on Mpc/h
      if (myid==0)
        cout << zi << "  " << dli*speedcunit << endl;
    }
    infiledc.close();
  }
  else{
    cout << "  " << endl;
    cout << " the comoving distance file: " << idc << endl;
    cout << " does not exists " << endl;
    cout << " I will STOP here!!! " << endl;
    MPI_Abort(MPI_COMM_WORLD,-1);
  }

  if(zs>zl[zl.size()-1]){
    cout << " source redshift larger than the highest available redshift in the comoving distance file " << endl;
    cout << "  that is = " << zl[zl.size()-1] << endl;
    cout << " I will STOP here !!! " << endl;
    MPI_Abort(MPI_COMM_WORLD,-1);
  }


  Ds = getY(zl,dl,zs);         // comoving distance of the last plane
  vector <int>    replication; // Number of repetitions of the i-th snapshot box
  vector <string> fromsnap;    // From which snapshot the i-th lens plane were build
  vector <double> zsimlens;    // z of the i-th lens (z correspodent to d=1/2(ld+ld2))
  vector <double> ld;          // Start position of the i-th lens
  vector <double> ld2;         // End position of the i-th lens
  vector <double> zfromsnap;   // z from the snapshot selected to build the i-th lens
  vector <bool>   randomize;   // Bool variable to whether the positions should be re-randomized or not

  double dllow;
  double dlup;

  //Building the ligh cone
  dllow=0;
  dlup=getY(zl,dl,zs);

  double ldbut;
  int pos, nrepi=0;
  int nrep=0;
  do{
    nrep++;
    nrepi++;
    ldbut=dllow+nrep*boxl/numberOfLensPerSnap;
    double dlens=ldbut-0.5*boxl/numberOfLensPerSnap;
    int pos_temp = getSnap(lred, dl, zl, dlens);
    if (myid==0)
      cout << " simulation snapshots = " << ldbut << "  " << getY(dl,zl,ldbut) << "  " << nrep << " from snap " << lsnap[pos_temp] << "  " << getY(dl,zl,dlens) << endl;
    ld.push_back(ldbut-boxl/numberOfLensPerSnap);
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
    cout << " Comoving Distance of the last plane " << Ds << endl;
    cout << " nsnaps = " << nsnaps << endl;
  }
  //int pl=0;
  vector<int> pll;

  ofstream planelist;
  string planes_list;
  planes_list = directory+"planes_list_"+suffix+".txt";

  if (myid==0){

    planelist.open(planes_list.c_str());

    for(int i=0;i<fromsnap.size();i++){

      cout << zsimlens[i] << " planes = " << ld[i] << "  " << ld2[i] << "  " << replication[i] << " from snap " << fromsnap[i] << endl;
      planelist <<  i << "   " << zsimlens[i] << "   " << ld[i] << "   " << ld2[i] << "   " << replication[i] << "   " << fromsnap[i]
      << "   " << zfromsnap[i] << "  "<< randomize[i]  << endl;
      pll.push_back(i);

    }

    planelist.close();

  }

  // randomization of the box realizations :
  int nrandom = replication.back();
  vector<double> x0(nrandom), y0(nrandom), z0(nrandom); // ramdomizing the center of the simulation [0,1]
  vector<int> face(nrandom); // face of the dice
  vector<int> sgnX(nrandom), sgnY(nrandom),sgnZ(nrandom); // randomizing the box axis signs

  for(int i=0;i<nrandom;i++){

    if ( randomize[i] ){

      srand(seedcenter+i/numberOfLensPerSnap*13);
      x0[i] = rand() / float(RAND_MAX);
      y0[i] = rand() / float(RAND_MAX);
      z0[i] = rand() / float(RAND_MAX);
      if(myid==0){
        cout << "  " << endl;
        cout << " random centers  for the box " << i << " = " << x0[i] << "  " << y0[i] << "  " << z0[i] << endl;
      }
      face[i] = 7;
      srand(seedface+i/numberOfLensPerSnap*5);
      while(face[i]>6 || face[i]<1) face[i] = int(1+rand() / float(RAND_MAX)*5.+0.5);
      if (myid==0)
        cout << " face of the dice " << face[i] << std:: endl;
      sgnX[i] = 2;
      srand(seedsign+i/numberOfLensPerSnap*8);
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
  if(myid==0)
    cout << " set the field of view to be square in degrees " << endl;
  double h0,fovradiants;
  double om0, omL0;
  fovradiants = fov/180.*M_PI;
  // check if the field of view is too large with respect to the box size
  if(fovradiants*Ds>boxl && myid==0){
    cout << " field view too large ... I will STOP here!!! " << endl;
    cout << " value set is = " << fov << endl;
    cout << " maximum value allowed " << boxl/Ds*180./M_PI << " in degrees " << endl;
    MPI_Abort(MPI_COMM_WORLD,-1);;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // loop on snapshots ****************************************
  if(myid==0){
    cout << " now loop on " << replication.back() << " planes " << endl;
    cout << "  " << endl;
  }
  for(int nsnap=0;nsnap<replication.back();nsnap++){

    if(ld2[nsnap]-ld[nsnap] < 0 && myid==0){
      cout << " comoving distance of the starting point " << ld[nsnap] << endl;
      cout << " comoving distance of the final    point " << ld2[nsnap] << endl;
      cout << " please check this out! I will STOP here!!! " << endl;
      MPI_Abort(MPI_COMM_WORLD,-1);;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //float rcase = ld[nsnap]/boxl*numberOfLensPerSnap; //Number of piled boxes
    float rcase = floor(ld[nsnap]/boxl); //Number of piled boxes

    if (do_as_t11) npix=int((ld2[nsnap]+ld[nsnap])/2*fovradiants/rgrid*1e3)+1; //Override npix if do_as_t11 is True

    valarray<float> mapxytot( npix*npix );
    // type 0
    int ntotxy0;
    ntotxy0 = 0;
    valarray<float> mapxytot0( npix*npix );
    // type 1
    int ntotxy1;
    ntotxy1 = 0;
    valarray<float> mapxytot1( npix*npix );
    // type 2
    int ntotxy2;
    ntotxy2 = 0;
    valarray<float> mapxytot2( npix*npix );
    // type 3
    int ntotxy3;
    ntotxy3 = 0;
    valarray<float> mapxytot3( npix*npix );
    // type 4
    int ntotxy4;
    ntotxy4 = 0;
    valarray<float> mapxytot4( npix*npix );
    // type 5
    int ntotxy5;
    ntotxy5 = 0;
    valarray<float> mapxytot5( npix*npix );

    string File = pathsnap+fromsnap[nsnap];

    string snappl;
    if( pll[nsnap]<10) snappl = "00"+conv(pll[nsnap],fINT);
    else if(pll[nsnap]>=10 && pll[nsnap]<100 ) snappl = "0"+conv(pll[nsnap],fINT);
    else snappl = conv(pll[nsnap],fINT);
    string Filesub;

    if(ifstream(directory+simulation+"."+snappl+".plane_"+snpix+"_"+suffix+".fits") && partinplanes == "ALL"){
      if (myid==0)
        cout << directory+simulation+"."+snappl+".plane_"+snpix+"_"+suffix+".fits" << " "<< "Already exists" <<endl;
      continue; // If this lens plane was already created go to the next one
    }

    float num_float1, num_float2, num_float3; // aux variables to read x,y, and z positions

    // redshift and dl of the simulation
    double zsim, dlsim;
    int intdiv, remaindiv,ffmin,ffmax;
    //Computing Working Balance
    intdiv=nfiles/numprocs;
    remaindiv=nfiles%numprocs;

    if(myid!=numprocs-1){
      ffmin=myid*intdiv;
      ffmax=(myid+1)*intdiv;
    }else{
      ffmin=myid*intdiv;
      ffmax=(myid+1)*intdiv+remaindiv;
    }

    for (unsigned int ff=ffmin; ff<ffmax; ff++){
      // map for each mass type
      valarray<float> mapxy0( npix*npix );
      valarray<float> mapxy1( npix*npix );
      valarray<float> mapxy2( npix*npix );
      valarray<float> mapxy3( npix*npix );
      valarray<float> mapxy4( npix*npix );
      valarray<float> mapxy5( npix*npix );

      // GADGET has 6 different particle type
      vector<float> xx0(0), yy0(0), zz0(0);
      vector<float> xx1(0), yy1(0), zz1(0);
      vector<float> xx2(0), yy2(0), zz2(0);
      vector<float> xx3(0), yy3(0), zz3(0);
      vector<float> xx4(0), yy4(0), zz4(0);
      vector<float> xx5(0), yy5(0), zz5(0);

      string file_in = File+"."+conv(ff,fINT);

      ifstream fin(file_in.c_str());
      if (!fin) {cerr <<"Error in opening the file: "<<file_in<<"!\n\a"; MPI_Abort(MPI_COMM_WORLD,-1);}

      if (myid==0)
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
        cout << "Box size: "<< data.boxsize*POS_U<< endl;
        cout << "Omega_matter: "<< data.om0<< endl;
        cout << "Omega_DE: "<< data.oml<< endl;
        cout << "h: "<< data.h<< endl;
        cout << "Flag sage: "<< data.flag_sage<< endl;
        cout << "Flag metals: "<< data.flag_metals<< endl;
        cout << "N. tot HW: "<< data.nTotalHW[0]<<" "<< data.nTotalHW[1]<<" "<< data.nTotalHW[2]<<" "<< data.nTotalHW[3]<<" "<< data.nTotalHW[4]<<" "<< data.nTotalHW[5]<<" "<< endl;
        cout << "Flag entropy: "<< data.flag_entropy<<endl;


        if(data.nTotalHW[0]+data.nTotalHW[1]+data.nTotalHW[2]+data.nTotalHW[3]+data.nTotalHW[4]+data.nTotalHW[5]!=0){

	      cout << "More than 2^32 particles!!" << endl;
	      cout << "Total Number of Particles is actually: "<<endl;
              cout << (long long) (data.nTotalHW[0]*pow(2,32)+data.npartTotal[0]) <<" "<<
                      (long long) (data.nTotalHW[1]*pow(2,32)+data.npartTotal[1]) <<" "<<
                      (long long) (data.nTotalHW[2]*pow(2,32)+data.npartTotal[2]) <<" "<<
                      (long long) (data.nTotalHW[3]*pow(2,32)+data.npartTotal[3]) <<" "<<
                      (long long) (data.nTotalHW[4]*pow(2,32)+data.npartTotal[4]) <<" "<<
                      (long long) (data.nTotalHW[5]*pow(2,32)+data.npartTotal[5]) << endl;

        }

      }

      BLOCK block;
      if(NBLOCKS>0){
        fin >> block;
      }

      for (int i=0; i < NBLOCKS; i++){

        if (i==0 && myid==0){
          cout << "Size of Header is " << sizeof(data) << endl;
          cout << "Should be         " << block.blocksize1 << endl;
        }else{
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

      // total number of particles as the sum of all of them
      int dim = data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]+data.npart[4]+data.npart[5];
      int dimmass0=0;
      int hydro=0;

      for(int i=0;i<=5;i++){
	       if(data.massarr[i]==0){dimmass0+=data.npart[i];}
      }


      if(dimmass0==0){
           if (myid==0){
	            cout << "		@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
              cout << "		@                          @" << endl;
	            cout << "		@  !!DM only simulation!!  @" << endl;
	            cout << "		@                          @" << endl;
	            cout << "		@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
          }
      }
      else{
         if (myid==0){
	          cout << "		@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
	          cout << "		@                          @" << endl;
	          cout << "		@  !!Hydro   simulation!!  @" << endl;
	          cout << "		@                          @" << endl;
	          cout << "		@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
         }
	       hydro=1;
      }

      if (myid==0){
        cout << " .......................................................... " << endl;
        cout << "   number of particles in this snapshot: " << endl;
        cout << data.npart[0] << " " << data.npart[1] << " " << data.npart[2]
	       << " " << data.npart[3] << " " << data.npart[4] << " " << data.npart[5] << endl;
      }

      if(ff==ffmin){
	        // compute the comoving angular diameter distance at simulation redshift
	        zsim = data.redshift;
	        dlsim = getY(zl,dl,zsim);
	        h0 = data.h;
	        cout << "  " << endl;
	        cout << "      __________________ COSMOLOGY __________________  " << endl;
	        cout << " " << endl;
	        om0 = data.om0;
	        omL0 = data.oml;
          if (myid==0){
	           cout << "      Omegam = " << data.om0 << " " << "Omegal = " << data.oml << endl;
	           cout << "           h = " << data.h   << " " << "BoxSize = " << data.boxsize*POS_U << endl;
	           cout << "      redshift = " << zsim <<   " " << "Dl (comoving) = " << dlsim << endl;
          }
          if(abs(boxl - data.boxsize*POS_U/1.e+3)>1.e-2 && myid==0){
	           cout << " set boxl and data.size differ ... check it! " << std:: endl;
	           cout << "  boxl = " << boxl << "  " << " data.boxsize = " << data.boxsize/1.e+3*POS_U << endl;
	           MPI_Abort(MPI_COMM_WORLD,-1);
	        }

          MPI_Barrier(MPI_COMM_WORLD);

          if (myid==0){

	           cout << "      _______________________________________________  " << endl;
	           cout << " " << endl;
	           cout << "   total number of particles in the simulation: " << endl;
	           cout << data.nTotalHW[0]*pow(2,32)+data.npartTotal[0] << " " << data.nTotalHW[1]*pow(2,32)+data.npartTotal[1] << " " <<
                        data.nTotalHW[2]*pow(2,32)+data.npartTotal[2] << " " << data.nTotalHW[3]*pow(2,32)+data.npartTotal[3] << " " <<
                        data.nTotalHW[4]*pow(2,32)+data.npartTotal[4] << " " << data.nTotalHW[5]*pow(2,32)+data.npartTotal[5] << endl;
	           cout << " " << endl;
	           cout << "   xparticle type mass array: " << endl;
	           cout << data.massarr[0] << " " << data.massarr[1] << " " << data.massarr[2]
	            << " " << data.massarr[3] << " " <<  data.massarr[4] << " " <<  data.massarr[5] << endl;
          }
	        m0 = data.massarr[0];
	        m1 = data.massarr[1];
	        m2 = data.massarr[2];
	        m3 = data.massarr[3];
	        m4 = data.massarr[4];
	        m5 = data.massarr[5];
      }

	    // type0
	    for (int pp=0; pp<data.npart[0]; pp++){
	      fin.read((char *)&num_float1, sizeof(num_float1));
	      fin.read((char *)&num_float2, sizeof(num_float2));
	      fin.read((char *)&num_float3, sizeof(num_float3));

	      float x, y, z;

	      float xb, yb, zb;

	      xb = sgnX[nsnap]*(((num_float1)/data.boxsize));
	      yb = sgnY[nsnap]*(((num_float2)/data.boxsize));
	      zb = sgnZ[nsnap]*(((num_float3)/data.boxsize));

	      // wrapping periodic condition
	      if(xb>1.) xb = xb - 1.;
	      if(yb>1.) yb = yb - 1.;
	      if(zb>1.) zb = zb - 1.;
	      if(xb<0.) xb = 1. + xb;
	      if(yb<0.) yb = 1. + yb;
	      if(zb<0.) zb = 1. + zb;
	      switch (face[nsnap]){
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
	      x = x - x0[nsnap];
	      y = y - y0[nsnap];
	      z = z - z0[nsnap];
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

	      xb = sgnX[nsnap]*(((num_float1)/data.boxsize));
	      yb = sgnY[nsnap]*(((num_float2)/data.boxsize));
	      zb = sgnZ[nsnap]*(((num_float3)/data.boxsize));

	      // wrapping periodic condition
	      if(xb>1.) xb = xb - 1.;
	      if(yb>1.) yb = yb - 1.;
	      if(zb>1.) zb = zb - 1.;
	      if(xb<0.) xb = 1. + xb;
	      if(yb<0.) yb = 1. + yb;
	      if(zb<0.) zb = 1. + zb;
	      switch (face[nsnap]){
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
	      x = x - x0[nsnap];
	      y = y - y0[nsnap];
	      z = z - z0[nsnap];
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

	      xb = sgnX[nsnap]*(((num_float1)/data.boxsize));
	      yb = sgnY[nsnap]*(((num_float2)/data.boxsize));
	      zb = sgnZ[nsnap]*(((num_float3)/data.boxsize));

	      // wrapping periodic condition
	      if(xb>1.) xb = xb - 1.;
	      if(yb>1.) yb = yb - 1.;
	      if(zb>1.) zb = zb - 1.;
	      if(xb<0.) xb = 1. + xb;
	      if(yb<0.) yb = 1. + yb;
	      if(zb<0.) zb = 1. + zb;
	      switch (face[nsnap]){
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
	      x = x - x0[nsnap];
	      y = y - y0[nsnap];
	      z = z - z0[nsnap];
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
	    for (int pp=data.npart[0]+data.npart[1]+data.npart[2];pp<data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]; pp++) {
	      fin.read((char *)&num_float1, sizeof(num_float1));
	      fin.read((char *)&num_float2, sizeof(num_float2));
	      fin.read((char *)&num_float3, sizeof(num_float3));
	      float x, y, z;
	      float xb, yb, zb;

	      xb = sgnX[nsnap]*(((num_float1)/data.boxsize));
	      yb = sgnY[nsnap]*(((num_float2)/data.boxsize));
	      zb = sgnZ[nsnap]*(((num_float3)/data.boxsize));

	      // wrapping periodic condition
	      if(xb>1.) xb = xb - 1.;
	      if(yb>1.) yb = yb - 1.;
	      if(zb>1.) zb = zb - 1.;
	      if(xb<0.) xb = 1. + xb;
	      if(yb<0.) yb = 1. + yb;
	      if(zb<0.) zb = 1. + zb;
	      switch (face[nsnap]){
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
	      x = x - x0[nsnap];
	      y = y - y0[nsnap];
	      z = z - z0[nsnap];
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

	      xb = sgnX[nsnap]*(((num_float1)/data.boxsize));
	      yb = sgnY[nsnap]*(((num_float2)/data.boxsize));
	      zb = sgnZ[nsnap]*(((num_float3)/data.boxsize));

	      // wrapping periodic condition
	      if(xb>1.) xb = xb - 1.;
	      if(yb>1.) yb = yb - 1.;
	      if(zb>1.) zb = zb - 1.;
	      if(xb<0.) xb = 1. + xb;
	      if(yb<0.) yb = 1. + yb;
	      if(zb<0.) zb = 1. + zb;
	      switch (face[nsnap]){
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
	      x = x - x0[nsnap];
	      y = y - y0[nsnap];
	      z = z - z0[nsnap];
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

	      xb = sgnX[nsnap]*(((num_float1)/data.boxsize));
	      yb = sgnY[nsnap]*(((num_float2)/data.boxsize));
	      zb = sgnZ[nsnap]*(((num_float3)/data.boxsize));

	      // wrapping periodic condition
	      if(xb>1.) xb = xb - 1.;
	      if(yb>1.) yb = yb - 1.;
	      if(zb>1.) zb = zb - 1.;
	      if(xb<0.) xb = 1. + xb;
	      if(yb<0.) yb = 1. + yb;
	      if(zb<0.) zb = 1. + zb;
	      switch (face[nsnap]){
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
	      x = x - x0[nsnap];
	      y = y - y0[nsnap];
	      z = z - z0[nsnap];
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

      if(hydro){

	       int tot = (data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]+data.npart[4]+data.npart[5]);

      	 BLOCK block; fin >> block;
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

         fin >> block;
         if (myid==0){
      	    cout << "Size of Header is " << block.blocksize1 << endl;
	          cout << "Should be         " << 12*tot*sizeof(int8_t) << endl;
            cout << "Fast Fowarding next block. Name: ";
            cout << block.name[0];
            cout << block.name[1];
            cout << block.name[2];
            cout << block.name[3] << endl;
         }

         fin.seekg(block.blocksize2/sizeof(int8_t),fin.cur);

         BLOCK block2; fin >> block2;
         if (myid==0){
           cout << "Reading next block. Name: ";
           cout << block2.name[0];
           cout << block2.name[1];
           cout << block2.name[2];
           cout << block2.name[3] << endl;
           cout << "Should be                 MASS" << endl;
         }

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
        if (myid==0){
	       cout << " " << endl;
	       cout << " n0 particles " << endl;
	       cout << "xmin = " << xmin << endl;
	       cout << "xmax = " << xmax << endl;
	       cout << "ymin = " << ymin << endl;
	       cout << "ymax = " << ymax << endl;
	       cout << "zmin = " << zmin << endl;
	       cout << "zmax = " << zmax << endl;
	       cout << "  " << endl;
        }

	      if(xmin<0 || ymin<0 || zmin< 0){
	        cout << "xmin = " << xmin << endl;
	        cout << "xmax = " << xmax << endl;
	        cout << "ymin = " << ymin << endl;
	        cout << "ymax = " << ymax << endl;
	        cout << "zmin = " << zmin << endl;
	        cout << "zmax = " << zmax << endl;
	        cout << "  0 type check this!!! I will STOP here!!! " << endl;
	        MPI_Abort(MPI_COMM_WORLD,-1);
	      }
        if (myid==0){
	         cout << " ... mapping type 0 particles on the grid with " << npix << " pixels" << endl;
           cout << "Min distance: "<< ld[nsnap]/data.boxsize*1.e+3/POS_U<< " "<<ld2[nsnap]/data.boxsize*1.e+3/POS_U << endl;
           cout << "Rcase       : "<< rcase << endl;
	          // 2Dgrid
        }
	      vector<float> xs(0),ys(0),ms(0);

	      for(int l=0;l<n0;l++){

          if(hydro && data.massarr[0]==0){
 			      fin.read((char *)&num_float1, sizeof(num_float1));
            if (num_float1>MAX_M)
              num_float1=0;
          }
		      else{num_float1=m0;}

	        double di = sqrt(pow(xx0[l]-0.5,2) + pow(yy0[l]-0.5,2) + pow(zz0[l],2))*data.boxsize/1.e+3*POS_U;;
	        if(di>=ld[nsnap] && di<ld2[nsnap]){
	          double rai,deci,dd;
	          getPolar(xx0[l]-0.5,yy0[l]-0.5,zz0[l],&rai,&deci,&dd);

	          if(fabs(rai)<=fovradiants*(1.+2./npix)*0.5 && fabs(deci)<=fovradiants*(1.+2./npix)*0.5){
              xs.push_back(deci/fovradiants+0.5);
	            ys.push_back(rai/fovradiants+0.5);
              if(sn_opt==0){
                ms.push_back(num_float1);
              }
              else{
                if(rand()/ float(RAND_MAX) < 1./pow(2,sn_opt)) ms.push_back(pow(2,sn_opt)*num_float1);
                else ms.push_back(0.);
              }
	          }
	        }
	      }
	      totPartxy0=xs.size();
	      ntotxy0+=totPartxy0;

	      if(totPartxy0>0){
          mapxy0 = gridist_w(xs,ys,ms,npix,do_NGP);
        }

	      //re-normalize to the total mass!
	      double mtot0=0.;
        double mnorm=0.;
        for(int i=0;i<ms.size();i++)
          mnorm+=ms[i];

	      if(totPartxy0>0){
	        for(int l=0;l<npix*npix;l++)
            mtot0 += mapxy0[l];
	        if(mtot0==0.)
            mtot0=1.; //To avoid NaN
	        for(int l=0;l<npix*npix;l++) mapxy0[l]*=mnorm/mtot0;
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
        if (myid==0){
	       cout << " " << endl;
	       cout << " n1 particles " << endl;
	       cout << "xmin = " << xmin << endl;
	       cout << "xmax = " << xmax << endl;
	       cout << "ymin = " << ymin << endl;
	       cout << "ymax = " << ymax << endl;
	       cout << "zmin = " << zmin << endl;
	       cout << "zmax = " << zmax << endl;
	       cout << "  " << endl;
        }
	      if(xmin<0 || ymin<0 || zmin< 0){
	        cout << "xmin = " << xmin << endl;
	        cout << "xmax = " << xmax << endl;
	        cout << "ymin = " << ymin << endl;
	        cout << "ymax = " << ymax << endl;
	        cout << "zmin = " << zmin << endl;
	        cout << "zmax = " << zmax << endl;
	        cout << "  1 type check this!!! I will STOP here!!! " << endl;
	        MPI_Abort(MPI_COMM_WORLD,-1);
	      }
        if (myid==0){
	        cout << " ... mapping type 1 particles on the grid with " << npix << " pixels" << endl;
          cout << "n1:" << n1 << endl;
        }
	      // 2Dgrid
	      vector<float> xs(0),ys(0),ms(0);

	      for(int l=0;l<n1;l++){
		      if(hydro && data.massarr[1]==0){
            fin.read((char *)&num_float1, sizeof(num_float1));
            if (num_float1>MAX_M)
              num_float1=0;
          }
		      else num_float1=m1;

	        double di = sqrt(pow(xx1[l]-0.5,2) + pow(yy1[l]-0.5,2) + pow(zz1[l],2))*data.boxsize*POS_U/1.e+3;
	        if(di>=ld[nsnap] && di<ld2[nsnap]){
	          double rai,deci,dd;
	          getPolar(xx1[l]-0.5,yy1[l]-0.5,zz1[l],&rai,&deci,&dd);
	          if(fabs(rai)<=fovradiants*(1.+2./npix)*0.5 && fabs(deci)<=fovradiants*(1.+2./npix)*0.5){
              xs.push_back(deci/fovradiants+0.5);
	            ys.push_back(rai/fovradiants+0.5);
              if(sn_opt==0){
                ms.push_back(num_float1);
              }
              else{
                if(rand()/ float(RAND_MAX) < 1./pow(2,sn_opt)) ms.push_back(pow(2,sn_opt)*num_float1);
                else ms.push_back(0.);
              }
	          }
	        }
	      }
	      totPartxy1=xs.size();
	      ntotxy1+=totPartxy1;

	      if(totPartxy1>0) 	mapxy1 = gridist_w(xs,ys,ms,npix,do_NGP);

	      // re-normalize to the total mass!
	      double mtot1=0;
	      double mnorm=accumulate(ms.begin(), ms.end(), 0.);

	      if(totPartxy1>0){
	        for(int l=0;l<npix*npix;l++) mtot1 += mapxy1[l];
          if(mtot1==0.) mtot1=1.; //To avoid NaN
	        for(int l=0;l<npix*npix;l++) mapxy1[l]=mapxy1[l]/mtot1*mnorm;
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
        if (myid==0){
	         cout << " " << endl;
	         cout << " n2 particles " << endl;
	         cout << "xmin = " << xmin << endl;
	         cout << "xmax = " << xmax << endl;
	         cout << "ymin = " << ymin << endl;
	         cout << "ymax = " << ymax << endl;
	         cout << "zmin = " << zmin << endl;
	         cout << "zmax = " << zmax << endl;
	         cout << "  " << endl;
        }
	      if(xmin<0 || ymin<0 || zmin< 0){
	        cout << "xmin = " << xmin << endl;
	        cout << "xmax = " << xmax << endl;
	        cout << "ymin = " << ymin << endl;
	        cout << "ymax = " << ymax << endl;
	        cout << "zmin = " << zmin << endl;
	        cout << "zmax = " << zmax << endl;
	        cout << "  2 type check this!!! I will STOP here!!! " << endl;
	        MPI_Abort(MPI_COMM_WORLD,-1);
	      }
        if (myid==0)
	       cout << " ... mapping type 2 particles on the grid with " << npix << " pixels" << endl;
	       // 2Dgrid
	      vector<float> xs(0),ys(0),ms(0);
	      for(int l=0;l<n2;l++){
		      if(hydro && data.massarr[2]==0){
            fin.read((char *)&num_float1, sizeof(num_float1));
            if (num_float1>MAX_M)
              num_float1=0;
          }
          else num_float1=m2;

	        double di = sqrt(pow(xx2[l]-0.5,2) + pow(yy2[l]-0.5,2) + pow(zz2[l],2))*data.boxsize*POS_U/1.e+3;
	        if(di>=ld[nsnap] && di<ld2[nsnap]){
	          double rai,deci,dd;
	          getPolar(xx2[l]-0.5,yy2[l]-0.5,zz2[l],&rai,&deci,&dd);
	          if(fabs(rai)<=fovradiants*(1.+2./npix)*0.5 && fabs(deci)<=fovradiants*(1.+2./npix)*0.5){
	            xs.push_back(deci/fovradiants+0.5);
	            ys.push_back(rai/fovradiants+0.5);
              if(sn_opt==0){
                ms.push_back(num_float1);
              }
              else{
                if(rand()/ float(RAND_MAX) < 1./pow(2,sn_opt)) ms.push_back(pow(2,sn_opt)*num_float1);
                else ms.push_back(0.);
              }
	          }
	        }
	      }
	      totPartxy2=xs.size();
	      // cout << " n2: totPartxy2 " << totPartxy2 << endl;
	      ntotxy2+=totPartxy2;

	      if(totPartxy2>0) mapxy2 = gridist_w(xs,ys,ms,npix,do_NGP);
	      // re-normalize to the total mass!
	      double mtot2=0;
	      double mnorm=accumulate(ms.begin(),ms.end(),0.);
	      if(totPartxy2>0){
	        for(int l=0;l<npix*npix;l++) mtot2 += mapxy2[l];
	        if(mtot2==0.) mtot2=1.; //To avoid NaN
          for(int l=0;l<npix*npix;l++) mapxy2[l]=mapxy2[l]/mtot2*mnorm;

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
        if (myid==0){
	         cout << " " << endl;
	         cout << " n3 particles " << endl;
	         cout << "xmin = " << xmin << endl;
	         cout << "xmax = " << xmax << endl;
	         cout << "ymin = " << ymin << endl;
	         cout << "ymax = " << ymax << endl;
	         cout << "zmin = " << zmin << endl;
	         cout << "zmax = " << zmax << endl;
	         cout << "  " << endl;
        }
	      if(xmin<0 || ymin<0 || zmin< 0){
	        cout << "xmin = " << xmin << endl;
	        cout << "xmax = " << xmax << endl;
	        cout << "ymin = " << ymin << endl;
	        cout << "ymax = " << ymax << endl;
	        cout << "zmin = " << zmin << endl;
	        cout << "zmax = " << zmax << endl;
	        cout << "  3 type check this!!! I will STOP here!!! " << endl;
	        MPI_Abort(MPI_COMM_WORLD,-1);
	      }
        if (myid==0)
	       cout << " ... mapping type 3 particles on the grid with " << npix << " pixels" << endl;
	       // 2Dgrid
	      vector<float> xs(0),ys(0),ms(0);
	      for(int l=0;l<n3;l++){

		      if(hydro && data.massarr[3]==0){
            fin.read((char *)&num_float1, sizeof(num_float1));
            if (num_float1>MAX_M)
              num_float1=0;
          }
          else num_float1=m3;

	        double di = sqrt(pow(xx3[l]-0.5,2) + pow(yy3[l]-0.5,2) + pow(zz3[l],2))*data.boxsize*POS_U/1.e+3;
	        if(di>=ld[nsnap] && di<ld2[nsnap]){
	          double rai,deci,dd;
	          getPolar(xx3[l]-0.5,yy3[l]-0.5,zz3[l],&rai,&deci,&dd);
	          if(fabs(rai)<=fovradiants*(1.+2./npix)*0.5 && fabs(deci)<=fovradiants*(1.+2./npix)*0.5){
	            xs.push_back(deci/fovradiants+0.5);
	            ys.push_back(rai/fovradiants+0.5);
              if(sn_opt==0){
                ms.push_back(num_float1);
              }
              else{
                if(rand()/ float(RAND_MAX) < 1./pow(2,sn_opt)) ms.push_back(pow(2,sn_opt)*num_float1);
                else ms.push_back(0.);
              }
	          }
	        }
	      }

        totPartxy3=xs.size();
	      // cout << " n3: totPartxy3 " << totPartxy3 << endl;
	      ntotxy3+=totPartxy3;

	      if(totPartxy3>0) mapxy3 = gridist_w(xs,ys,ms,npix,do_NGP);

	      // re-normalize to the total mass!
	      double mtot3=0;
	      double mnorm=accumulate(ms.begin(),ms.end(),0.);
	      if(totPartxy3>0){
	        for(int l=0;l<npix*npix;l++) mtot3 += mapxy3[l];
          if(mtot3==0.) mtot3=1.; //To avoid NaN
	        for(int l=0;l<npix*npix;l++) mapxy3[l]=mapxy3[l]/mtot3*mnorm;
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
        if (myid==0){
	         cout << " " << endl;
	         cout << " n4 particles " << endl;
	         cout << "xmin = " << xmin << endl;
	         cout << "xmax = " << xmax << endl;
	         cout << "ymin = " << ymin << endl;
	         cout << "ymax = " << ymax << endl;
	         cout << "zmin = " << zmin << endl;
	         cout << "zmax = " << zmax << endl;
	         cout << "  " << endl;
        }
	      if(xmin<0 || ymin<0 || zmin< 0){
	        cout << "xmin = " << xmin << endl;
	        cout << "xmax = " << xmax << endl;
	        cout << "ymin = " << ymin << endl;
	        cout << "ymax = " << ymax << endl;
	        cout << "zmin = " << zmin << endl;
	        cout << "zmax = " << zmax << endl;
	        cout << "  4 type check this!!! I will STOP here!!! " << endl;
	        MPI_Abort(MPI_COMM_WORLD,-1);
	      }
	      cout << " ... mapping type 4 particles on the grid with " << npix << " pixels" << endl;
	      // 2Dgrid
	      vector<float> xs(0),ys(0),ms(0);
	      // vector<double> ra(0),dec(0);
	      for(int l=0;l<n4;l++){

		      if(hydro && data.massarr[4]==0){
            fin.read((char *)&num_float1, sizeof(num_float1));
            if (num_float1>MAX_M)
              num_float1=0;
          }
		      else num_float1=m4;

	        double di = sqrt(pow(xx4[l]-0.5,2) + pow(yy4[l]-0.5,2) + pow(zz4[l],2))*data.boxsize*POS_U/1.e+3;
	        if(di>=ld[nsnap] && di<ld2[nsnap]){
	          double rai,deci,dd;
	          getPolar(xx4[l]-0.5,yy4[l]-0.5,zz4[l],&rai,&deci,&dd);
	          if(fabs(rai)<=fovradiants*(1.+2./npix)*0.5 && fabs(deci)<=fovradiants*(1.+2./npix)*0.5){
	            xs.push_back(deci/fovradiants+0.5);
	            ys.push_back(rai/fovradiants+0.5);
              if(sn_opt==0){
                ms.push_back(num_float1);
              }
              else{
                if(rand()/ float(RAND_MAX) < 1./pow(2,sn_opt)) ms.push_back(pow(2,sn_opt)*num_float1);
                else ms.push_back(0.);
              }
	          }
	        }
	      }
	      totPartxy4=xs.size();
	      // cout << " n4: totPartxy4 " << totPartxy4 << endl;
	      ntotxy4+=totPartxy4;

	      if(totPartxy4>0) mapxy4 = gridist_w(xs,ys,ms,npix,do_NGP);

	      // re-normalize to the total mass!
	      double mtot4=0;
        double mnorm = accumulate(ms.begin(),ms.end(),0.);
	      if(totPartxy4>0){
	        for(int l=0;l<npix*npix;l++) mtot4 += mapxy4[l];
	        if(mtot4==0.) mtot4=1.; //To avoid NaN
          for(int l=0;l<npix*npix;l++) mapxy4[l]=mapxy4[l]/mtot4*mnorm;
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
        if (myid==0){
	         cout << " " << endl;
	         cout << " n5 particles " << endl;
	         cout << "xmin = " << xmin << endl;
	         cout << "xmax = " << xmax << endl;
	         cout << "ymin = " << ymin << endl;
	         cout << "ymax = " << ymax << endl;
	         cout << "zmin = " << zmin << endl;
	         cout << "zmax = " << zmax << endl;
	         cout << "  " << endl;
        }
	      if(xmin<0 || ymin<0 || zmin< 0){
	        cout << "xmin = " << xmin << endl;
	        cout << "xmax = " << xmax << endl;
	        cout << "ymin = " << ymin << endl;
	        cout << "ymax = " << ymax << endl;
	        cout << "zmin = " << zmin << endl;
	        cout << "zmax = " << zmax << endl;
	        cout << "  5 type check this!!! I will STOP here!!! " << endl;
	        MPI_Abort(MPI_COMM_WORLD,-1);
	      }
        if (myid==0)
	       cout << " ... mapping type 5 particles on the grid with " << npix << " pixels" << endl;
	      // 2Dgrid
	      vector<float> xs(0),ys(0),ms(0);

	      for(int l=0;l<n5;l++){

		      if(hydro && data.massarr[5]==0){
            if(l==0){
              fin.seekg(n5*sizeof(int32_t)/sizeof(int8_t),fin.cur);
              if (myid==0){
                cout << "\n" <<endl;
                cout << "Fast Fowarding blocks: ";
              }
		          for(int i=0;i<9;i++){

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

 		        fin.read((char *)&num_float1, sizeof(num_float1));
            if (num_float1>MAX_M)
              num_float1=0;
          }
		      else num_float1=m5;

	        double di = sqrt(pow(xx5[l]-0.5,2) + pow(yy5[l]-0.5,2) + pow(zz5[l],2))*data.boxsize*POS_U/1.e+3;
	        if(di>=ld[nsnap] && di<ld2[nsnap]){
	          double rai,deci,dd;
	          getPolar(xx5[l]-0.5,yy5[l]-0.5,zz5[l],&rai,&deci,&dd);
	          if(fabs(rai)<=fovradiants*(1.+2./npix)*0.5 && fabs(deci)<=fovradiants*(1.+2./npix)*0.5){
	            xs.push_back(deci/fovradiants+0.5);
	            ys.push_back(rai/fovradiants+0.5);
              if(sn_opt==0){
                ms.push_back(num_float1);
              }
              else{
                if(rand()/ float(RAND_MAX) < 1./pow(2,sn_opt)) ms.push_back(pow(2,sn_opt)*num_float1);
                else ms.push_back(0.);
              }
	          }
	        }
	      }
	      totPartxy5=xs.size();
	      // cout << " n5: totPartxy5 " << totPartxy5 << endl;
	      ntotxy5+=totPartxy5;

	      if(totPartxy5>0) mapxy5 = gridist_w(xs,ys,ms,npix,do_NGP);

	      if(hydro){

		      fin >> block;
          if (myid==0){
            cout << endl << "If masses were read correctly next box should be BHMD " << endl;
	          cout << "Name:                                            ";
	          cout << block.name[0];
	          cout << block.name[1];
	          cout << block.name[2];
	          cout << block.name[3] << endl;
          }
	      }

	      // re-normalize to the total mass!
	      double mtot5=0;
	      double mnorm = accumulate(ms.begin(),ms.end(),0.);
	      if(totPartxy5>0){
	        for(int l=0;l<npix*npix;l++) mtot5 += mapxy5[l];
          if(mtot5==0.) mtot5=1.; //To avoid NaN
	          for(int l=0;l<npix*npix;l++) mapxy5[l]=mapxy5[l]/mtot5*mnorm;
	        }
      }
      cout << " maps done! from Rank:" << myid << endl;

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

      MPI_Reduce( MPI_IN_PLACE, &mapxytot[0],  npix*npix, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce( MPI_IN_PLACE, &mapxytot0[0], npix*npix, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce( MPI_IN_PLACE, &mapxytot1[0], npix*npix, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce( MPI_IN_PLACE, &mapxytot2[0], npix*npix, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce( MPI_IN_PLACE, &mapxytot3[0], npix*npix, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce( MPI_IN_PLACE, &mapxytot4[0], npix*npix, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce( MPI_IN_PLACE, &mapxytot5[0], npix*npix, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      if (myid==0)
        cout << " done map*tot " << endl;
    }
    if (myid==0){
      if(partinplanes=="ALL"){
	       /*
	        * write image array(s) to FITS files all particles in a FITS file!
	        */
	        long naxis = 2;
	        long naxes[2]={ npix,npix };
	        string fileoutput;
	        // fileoutput = simulation+"."+snapnum+".plane_"+snpix+".fits";
	        fileoutput = directory+simulation+"."+snappl+".plane_"+snpix+"_"+suffix+".fits";
          cout << "Saving the maps on: " << fileoutput << endl;
	        unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	        vector<long> naxex( 2 );
	        naxex[0]=npix;
	        naxex[1]=npix;
	        PHDU *phxy=&ffxy->pHDU();
	        phxy->write( 1, npix*npix, mapxytot );
	        // phxy->addKey ("x0",x0," unit of the boxsize");
	        // phxy->addKey ("y0",y0," unit of the boxsize");
	        phxy->addKey ("REDSHIFT",zsim," ");
	        phxy->addKey ("PHYSICALSIZE",fov," ");
	        phxy->addKey ("PIXELUNIT",1.e+10/h0,"Mass unit in M_Sun");
	        phxy->addKey ("DlLOW",ld[nsnap]/h0,"comoving distance in Mpc");
	        phxy->addKey ("DlUP",ld2[nsnap]/h0,"comoving distance in Mpc");
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
	         unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	         vector<long> naxex( 2 );
	         naxex[0]=npix;
	         naxex[1]=npix;
	         PHDU *phxy=&ffxy->pHDU();
	         phxy->write( 1, npix*npix, mapxytot0 );
	         phxy->addKey ("REDSHIFT",zsim," ");
	         phxy->addKey ("PHYSICALSIZE",fov," ");
	         phxy->addKey ("PIXELUNIT",1.e+10/h0,"Mass unit in M_Sun");
	         phxy->addKey ("DlLOW",ld[nsnap]/h0,"comoving distance in Mpc");
	         phxy->addKey ("DlUP",ld2[nsnap]/h0,"comoving distance in Mpc");
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
	        unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	        vector<long> naxex( 2 );
	        naxex[0]=npix;
	        naxex[1]=npix;
	        PHDU *phxy=&ffxy->pHDU();
	        phxy->write( 1, npix*npix, mapxytot1 );
	        phxy->addKey ("REDSHIFT",zsim," ");
	        phxy->addKey ("PHYSICALSIZE",fov," ");
	        phxy->addKey ("PIXELUNIT",1.e+10/h0,"Mass unit in M_Sun");
	        phxy->addKey ("DlLOW",ld[nsnap]/h0,"comoving distance in Mpc");
	        phxy->addKey ("DlUP",ld2[nsnap]/h0,"comoving distance in Mpc");
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
	        unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	        vector<long> naxex( 2 );
	        naxex[0]=npix;
	        naxex[1]=npix;
	        PHDU *phxy=&ffxy->pHDU();
	        phxy->write( 1, npix*npix, mapxytot2 );
	        phxy->addKey ("REDSHIFT",zsim," ");
	        phxy->addKey ("PHYSICALSIZE",fov," ");
	        phxy->addKey ("PIXELUNIT",1.e+10/h0,"Mass unit in M_Sun");
	        phxy->addKey ("DlLOW",ld[nsnap]/h0,"comoving distance in Mpc");
	        phxy->addKey ("DlUP",ld2[nsnap]/h0,"comoving distance in Mpc");
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
	        unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	        vector<long> naxex( 2 );
	        naxex[0]=npix;
	        naxex[1]=npix;
	        PHDU *phxy=&ffxy->pHDU();
          phxy->write( 1, npix*npix, mapxytot3 );
	        phxy->addKey ("REDSHIFT",zsim," ");
	        phxy->addKey ("PHYSICALSIZE",fov," ");
	        phxy->addKey ("PIXELUNIT",1.e+10/h0,"Mass unit in M_Sun");
	        phxy->addKey ("DlLOW",ld[nsnap]/h0,"comoving distance in Mpc");
	        phxy->addKey ("DlUP",ld2[nsnap]/h0,"comoving distance in Mpc");
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
	         unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	         vector<long> naxex( 2 );
	         naxex[0]=npix;
	         naxex[1]=npix;
	         PHDU *phxy=&ffxy->pHDU();
	         phxy->write( 1, npix*npix, mapxytot4 );
	         phxy->addKey ("REDSHIFT",zsim," ");
	         phxy->addKey ("PHYSICALSIZE",fov," ");
	         phxy->addKey ("PIXELUNIT",1.e+10/h0,"Mass unit in M_Sun");
	         phxy->addKey ("DlLOW",ld[nsnap]/h0,"comoving distance in Mpc");
	         phxy->addKey ("DlUP",ld2[nsnap]/h0,"comoving distance in Mpc");
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
	         unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	         vector<long> naxex( 2 );
	         naxex[0]=npix;
	         naxex[1]=npix;
	         PHDU *phxy=&ffxy->pHDU();
	         phxy->write( 1, npix*npix, mapxytot5 );
	         phxy->addKey ("REDSHIFT",zsim," ");
	         phxy->addKey ("PHYSICALSIZE",fov," ");
	         phxy->addKey ("PIXELUNIT",1.e+10/h0,"Mass unit in M_Sun");
	         phxy->addKey ("DlLOW",ld[nsnap]/h0,"comoving distance in Mpc");
	         phxy->addKey ("DlUP",ld2[nsnap]/h0,"comoving distance in Mpc");
	         phxy->addKey ("nparttype5",ntotxy5," ");
	         phxy->addKey ("HUBBLE",h0," ");
	         phxy->addKey ("OMEGAMATTER",om0," ");
	         phxy->addKey ("OMEGALAMBDA",omL0," ");
	         phxy->addKey ("m5",m5," ");
	       }
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  cout << " end of work ... ;-)  " << endl;
  return 0;
}
