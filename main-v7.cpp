#include "mpi.h"
#include "util_new.h"
#define MAX_M 1e3 // Threshold for mass; Particler heavier than MAX_M will be attached zero mass
#define POS_U 1.0 // Unit conversion from BoxSize unit lengh to kpc/h
#define NBLOCKSTOPOS 1 // Number of blocks to be fastforwarded to get to POS
#define NBLOCKSTOMASS 3 // Number of blocks to be fastforwarded from POS to get to MASS
#define NBLOCKSTOBHMASS 9 // Number of blocks to be fastforwarded from MASS to get to BHMASS
#define DO_NGP false // Use NGP as the MAS
#define numberOfLensPerSnap  1 // Number of Lens to be builded from a snap
#define neval 1000

/*****************************************************************************/
/*                                                                           */
/*            This code creates PLC from snapshots written as:               */
/*             - Density maps for Gadget2 snapshots                          */
/*             - File with Halo positions and Halo                           */
/*               properties for Halos inside the PLC                         */
/*                                                                           */
/*          The PLC builder algorithm, although written from scratch,        */
/*          resembles that presented in the original MapSim                  */
/*                         (Contact Carlo Giocoli - cgiocoli@gmail.com)      */
/*                                                                           */
/*                      developed  by Tiago Castro tiagobscastro@gmail.com   */
/*****************************************************************************/

using namespace std;

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

  string inifile=argv[1]; // Paramfile name
  if (myid==0){
    cout << "   ------------------------------------------------------ " << endl;
    cout << "   -        SLICER - Simulation LIght ConE BuildeR      - " << endl;
    cout << "   -                       v1.0                         - " << endl;
    cout << "   ------------------------------------------------------ " << endl;
  }

  // Reading InputFiles
  struct InputParams p;
  string snpix;
  bool physical;
  double fovradiants;
  size_t nsnap;

  if(readInput(&p, inifile, snpix, &physical))
    MPI_Abort(MPI_COMM_WORLD,-1);

  // Reading the Snapshots available
  vector <string> snappath; // List of SnapShots paths
  vector <double> snapred; // List of SnapShots redshift
  if( read_redlist(p.filredshiftlist, snapred, snappath, &p) )
    MPI_Abort(MPI_COMM_WORLD,-1);
  nsnap = snapred.size(); // Number of Snapshots

  /* Read the Snap Header and get Cosmological Values*/
  Header simdata;
  ifstream fin;
  if( read_header (p.pathsnap+snappath[0]+".0", &simdata, fin, true) )
    MPI_Abort(MPI_COMM_WORLD,-1);
  test_hydro(&p, &simdata);


  /* Creating an Instance of the cosmology class to compute distances (!!h=1!!) */
  cosmology cosmo(simdata.om0,simdata.oml,1.0,-1.0);
  /* Creating a table with redshifts and comovingdistances to be interpolated*/
  vector <double> zl(neval),dl(neval);
  for(int i=0;i<neval;i++){
    zl[i] = i * (p.zs+1.0)/(neval-1);
    dl[i] = cosmo.comovDist(zl[i])*speedcunit;
  }
  /* Initializing the auxiliary functions to get Dc given z and vice-versa */
  gsl_interp_accel *accGetDl = gsl_interp_accel_alloc ();
  gsl_interp_accel *accGetZl = gsl_interp_accel_alloc ();
  gsl_spline *GetDl = gsl_spline_alloc (gsl_interp_cspline, neval);
  gsl_spline *GetZl = gsl_spline_alloc (gsl_interp_cspline, neval);
  gsl_spline_init (GetDl, &zl[0], &dl[0], neval);
  gsl_spline_init (GetZl, &dl[0], &zl[0], neval);
  /* Comoving distance of the last plane*/
  p.Ds = gsl_spline_eval (GetDl, p.zs, accGetDl);
  if(test_fov(p.fov, simdata.boxsize/1e3, p.Ds, myid, &fovradiants))
    MPI_Abort(MPI_COMM_WORLD,-1);

  /* Creating an Instance of the Lens and Building the Planes */
  Lens lens;
  build_planes(&p, &simdata, lens, snapred, snappath, GetDl, accGetDl, GetZl, accGetZl, numberOfLensPerSnap, myid);
  /* Creating an Instance of the Randomization plan */
  Random random;
  randomize_box (random, &lens, &p, numberOfLensPerSnap, myid);

  /* Looping on the Snapshots */
  if(myid==0){
    cout << " Now loop on " << lens.nplanes << " planes " << endl;
    cout << "  " << endl;
  }
  for(int isnap=0; isnap < lens.nplanes; isnap++){

    float rcase = floor(lens.ld[nsnap]/simdata.boxsize*1e3);
    /* Override p.npix if physical is True */
    if (physical)
      p.npix=int( (lens.ld2[nsnap]+lens.ld[nsnap])/2*fovradiants/p.rgrid*1e3 )+1;
    /* Get the Snapshot Name */
    string File = p.pathsnap+lens.fromsnap[isnap];
    string snappl;
    if( lens.pll[isnap]<10)
      snappl = "00"+sconv(lens.pll[isnap],fINT);
    else if(lens.pll[isnap]>=10 && lens.pll[isnap]<100 )
      snappl = "0"+sconv(lens.pll[isnap],fINT);
    else
      snappl = sconv(lens.pll[isnap],fINT);

    /* If this lens plane was already created go to the next one */
    if(ifstream(p.directory+p.simulation+"."+snappl+".plane_"+snpix+"_"+p.suffix+".fits") && p.partinplanes == false){
      if (myid==0)
        cout << p.directory+p.simulation+"."+snappl+".plane_"+snpix+"_"+p.suffix+".fits" << " "<< "Already exists" <<endl;
      continue;
    }

    /* Computing Working Balance  */
    int intdiv, remaindiv,ffmin,ffmax;
    intdiv=simdata.numfiles/numprocs;
    remaindiv=simdata.numfiles%numprocs;

    if(myid!=numprocs-1){
      ffmin=myid*intdiv;
      ffmax=(myid+1)*intdiv;
    }else{
      ffmin=myid*intdiv;
      ffmax=(myid+1)*intdiv+remaindiv;
    }

    /* Starting the loop on different Snapshot subfiles */
    valarray<float> mapxytot( p.npix*p.npix );
    int ntotxyi[6];
    valarray<float> mapxytoti[6];
    for(int i=0;i<6;i++){
      ntotxyi[i]=0;
      mapxytoti[i].resize( p.npix*p.npix );
    }
    for (unsigned int ff=ffmin; ff<ffmax; ff++){

      Header data;
      string file_in = File+"."+sconv(ff,fINT);
      ifstream fin;
      if (read_header (file_in, &data, fin, false))
        MPI_Abort(MPI_COMM_WORLD,-1);
      if(ff==0)
        print_header(data);
      /* Creating the pointers for the Data structures */
      Gadget *gadget;
      SubFind *subfind;
      float *xx[6][3];

      if( p.simType.compare("SubFind") == 0 ){
        subfind = new SubFind;
        subfind->xx0.resize(data.npart[0]); subfind->yy0.resize(data.npart[0]); subfind->zz0.resize(data.npart[0]);
        xx[0][0]=&subfind->xx0[0]; xx[0][1]=&subfind->yy0[0]; xx[0][2]=&subfind->zz0[0];
        xx[1][0]=nullptr; xx[1][1]=nullptr; xx[1][2]=nullptr;
        xx[2][0]=nullptr; xx[2][1]=nullptr; xx[2][2]=nullptr;
        xx[3][0]=nullptr; xx[3][1]=nullptr; xx[3][2]=nullptr;
        xx[4][0]=nullptr; xx[4][1]=nullptr; xx[4][2]=nullptr;
        xx[5][0]=nullptr; xx[5][1]=nullptr; xx[5][2]=nullptr;
      }else{
        gadget = new Gadget;
        gadget->xx0.resize(data.npart[0]); gadget->yy0.resize(data.npart[0]); gadget->zz0.resize(data.npart[0]);
        gadget->xx1.resize(data.npart[1]); gadget->yy1.resize(data.npart[1]); gadget->zz1.resize(data.npart[1]);
        gadget->xx2.resize(data.npart[2]); gadget->yy2.resize(data.npart[2]); gadget->zz2.resize(data.npart[2]);
        gadget->xx3.resize(data.npart[3]); gadget->yy3.resize(data.npart[3]); gadget->zz3.resize(data.npart[3]);
        gadget->xx4.resize(data.npart[4]); gadget->yy4.resize(data.npart[4]); gadget->zz4.resize(data.npart[4]);
        gadget->xx5.resize(data.npart[5]); gadget->yy5.resize(data.npart[5]); gadget->zz5.resize(data.npart[5]);
        xx[0][0]=&gadget->xx0[0]; xx[0][1]=&gadget->yy0[0]; xx[0][2]=&gadget->zz0[0];
        xx[1][0]=&gadget->xx1[0]; xx[1][1]=&gadget->yy1[0]; xx[1][2]=&gadget->zz1[0];
        xx[2][0]=&gadget->xx2[0]; xx[2][1]=&gadget->yy2[0]; xx[2][2]=&gadget->zz2[0];
        xx[3][0]=&gadget->xx3[0]; xx[3][1]=&gadget->yy3[0]; xx[3][2]=&gadget->zz3[0];
        xx[4][0]=&gadget->xx4[0]; xx[4][1]=&gadget->yy4[0]; xx[4][2]=&gadget->zz4[0];
        xx[5][0]=&gadget->xx5[0]; xx[5][1]=&gadget->yy5[0]; xx[5][2]=&gadget->zz5[0];
      }

      ReadPos (fin,  &data, &p, &random, isnap, xx, rcase, myid);

      /* If Hydro run I have to read the masses, otherwise close the snapshot*/
      if(p.hydro)
         fastforwardToBlock (fin, "MASS", myid);
      else{
        fin.clear();
        fin.close();
      }

      // map for each mass type
      valarray<float> mapxyi[6];
      int ntotxyi[6];
      for(int i=0; i<6; i++)
        mapxyi[i].resize(p.npix*p.npix);

      if(MapParticles(fin, &data, &p, &lens, xx, fovradiants, isnap, mapxyi,
                                                     ntotxyi, myid))
        MPI_Abort(MPI_COMM_WORLD,-1);

      if(p.hydro){
        fin.clear();
        fin.close();
      }

      mapxytot+=mapxyi[0]+mapxyi[1]+mapxyi[2]+mapxyi[3]+mapxyi[4]+mapxyi[5];
      for(int i=0; i<6; i++)
        mapxytoti[i]+=mapxyi[i];

      if (myid==0)
        cout << " done map*tot " << endl;

      if( p.simType.compare("SubFind") == 0 )
        delete subfind;
      else
        delete gadget;

    }

    valarray<float> mapxytotrecv( p.npix*p.npix );
    valarray<float> mapxytotirecv[6];
    for(int i=0; i<6; i++)
      mapxytotirecv[i].resize( p.npix*p.npix );

    cout << " maps done! from Rank:" << myid << endl;

    MPI_Reduce( &mapxytot[0],  &mapxytotrecv[0],  p.npix*p.npix, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    for(int i=0; i<6; i++)
      MPI_Reduce( &mapxytoti[i][0],  &mapxytotirecv[i][0],  p.npix*p.npix, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  }

  gsl_spline_free (GetDl);gsl_spline_free (GetZl);
  gsl_interp_accel_free (accGetDl);gsl_interp_accel_free (accGetZl);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  cout << " end of work ... ;)  " << endl;

  return 0;
}
