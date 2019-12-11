#include "mpi.h"
#include "cosmology.h"
#include "densitymaps.h"
#include "writeplc.h"
#define neval 1000             // Number of Points to interpolate the comoving distance

/*****************************************************************************/
/*                                                                           */
/*            This code creates PLC from snapshots written as:               */
/*             - Density maps for Gadget2 snapshots                          */
/*             - File with Halo positions and Halo                           */
/*               properties for Halos inside the PLC                         */
/*                                                                           */
/*          The PLC builder algorithm, although re-written from scratch,     */
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
    cout << "   -                       v2.0                         - " << endl;
    cout << "   -                                                    - " << endl;
    cout << "   -               Running on "<<numprocs<<" processes              - " << endl;
    cout << "   -                                                    - " << endl;
    cout << "   ------------------------------------------------------ " << endl;
  }

  // Reading InputFiles
  struct InputParams p;
  double fovradiants;
  size_t nsnap;

  if(readInput(p, inifile))
    MPI_Abort(MPI_COMM_WORLD,-1);

  if(p.simType == "SubFind" && numprocs>1){
    cerr << "!! Some of the routines for SubFind PLC creates concurrency between the process !!" << endl;
    cerr << "!!                       and this is not sorted out: sorry :(                   !!" << endl;
    MPI_Abort(MPI_COMM_WORLD,-1);
  }

  // Reading the Snapshots available
  vector <string> snappath; // List of SnapShots paths
  vector <double> snapred; // List of SnapShots redshift
  vector <double> snapbox; // List of different boxes sizes
  if( readRedList(p.filredshiftlist, snapred, snappath, snapbox, p) )
    MPI_Abort(MPI_COMM_WORLD,-1);
  nsnap = snapred.size(); // Number of Snapshots

  /* Read the Snap Header and get Cosmological Values*/
  Header simdata;
  ifstream fin;
  if( readHeader (p.pathsnap+snappath[0]+".0", simdata, fin, true) )
    MPI_Abort(MPI_COMM_WORLD,-1);
  testHydro(p, simdata);

  /* Creating an Instance of the cosmology class to compute distances (!!h=1!!) */
  cosmology cosmo(simdata.om0,simdata.oml,1.0,-1.0);
  /* Creating a table with redshifts and comoving distances to be interpolated*/
  vector <double> zl(neval),dl(neval);
  for(int i=0;i<neval;i++){
    zl[i] = i * (p.zs+1.0)/(neval-1);
    dl[i] = cosmo.comovDist(zl[i])*speedcunit;
  }
  /* Initializing the auxiliary functions to get Dc given z and vice-versa */
  gsl_interp_accel *accGetDl = gsl_interp_accel_alloc ();
  gsl_interp_accel *accGetZl = gsl_interp_accel_alloc ();
  gsl_spline *getDl = gsl_spline_alloc (gsl_interp_cspline, neval);
  gsl_spline *getZl = gsl_spline_alloc (gsl_interp_cspline, neval);
  gsl_spline_init (getDl, &zl[0], &dl[0], neval);
  gsl_spline_init (getZl, &dl[0], &zl[0], neval);
  /* Comoving distance of the last plane*/
  p.Ds = gsl_spline_eval (getDl, p.zs, accGetDl);

  /* Creating an Instance of the Lens and Building the Planes */
  Lens lens;
  buildPlanes(p, lens, snapred, snappath, snapbox, getDl, accGetDl, getZl, accGetZl, numberOfLensPerSnap, myid);
  /* Testing if the PLC fits inside the piled boxes*/
  for(int i = 0; i< lens.ld.size(); i++){
#ifndef ReplicationOnPerpendicularPlane
    /*Check if the FOV does not require repetitions on the perpendicular plane*/
    if(i==0)
      lens.nrepperp.resize(lens.ld.size(), 0);
    if(testFov(p.fov, snapbox[lens.fromsnapi[i]]/1e3*POS_U, lens.ld2[i], myid, fovradiants))
      MPI_Abort(MPI_COMM_WORLD,-1);
#else
    /*Compute the number of repetitions required on the perpendicular plane*/
    if(i==0)
      lens.nrepperp.resize(lens.ld.size());
    computeReplications(p.fov, snapbox[lens.fromsnapi[i]]/1e3*POS_U, lens.ld2[i], myid, fovradiants, lens.nrepperp[i]);
    if(myid==0){
      if(i==0)
          cout << "\n"<< " Computing the Replications on the Perpendicular Plane" << endl;
      cout << " Plane "<< i <<" Dlow = " << lens.ld[i] << " Dlup = " << lens.ld2[i] << endl;
      cout << " Repetitions on the perpendicular plane = "<< pow(lens.nrepperp[i] + 1, 2) - 1 << endl;
    }
#endif
  }
  /* Creating an Instance of the Randomization plan */
  Random random;
  randomizeBox (random, lens, p, numberOfLensPerSnap, myid);

  /* Looping on the Snapshots */
  if(myid==0){
    cout << " Now loop on " << lens.nplanes << " planes " << endl;
    cout << "  " << endl;
  }
  for(int isnap=0; isnap < lens.nplanes; isnap++){

    /* Override p.npix if physical is True */
    if (p.physical)
      p.npix=int( (lens.ld2[isnap]+lens.ld[isnap])/2*fovradiants/p.rgrid*1e3/POS_U )+1;

    /* Get the Snapshot Name */
    string File = p.pathsnap+lens.fromsnap[isnap];
    string snappl;
    /*Override simdata*/
    Header simdata;
    ifstream fin;
    if( readHeader (File+".0", simdata, fin, true) )
      MPI_Abort(MPI_COMM_WORLD,-1);

    if( lens.pll[isnap]<10)
      snappl = "00"+sconv(lens.pll[isnap],fINT);
    else if(lens.pll[isnap]>=10 && lens.pll[isnap]<100 )
      snappl = "0"+sconv(lens.pll[isnap],fINT);
    else
      snappl = sconv(lens.pll[isnap],fINT);

    /* Computing Working Balance */
    int intdiv, remaindiv, ffmin, ffmax;
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
    if(p.simType == "Gadget"){

      /* If this lens plane was already created go to the next one */
      if(p.partinplanes == false){
        if( ifstream( fileOutput(p, snappl) ) ){
          if (myid==0)
            cout << fileOutput(p, snappl) << " Already exists" <<endl;
          continue;
        }
      }else{
        if (myid==0)
          cout << "!It is not possible to resume a Gadget run with partinplanes == true!" << endl <<
                  "!!               Files on Output folder will be overwritten        !!" << endl;
      }

      valarray<float> mapxytot;
      int ntotxyi[6];
      valarray<float> mapxytoti[6];
      float rcase;
      /*Computing the minimum distance for the lens in units of the current box size*/
      rcase = floor(lens.ld[isnap]/snapbox[lens.fromsnapi[isnap]]*1e3/POS_U);

      if( createDensityMaps (p, lens, random, isnap, ffmin, ffmax, File, fovradiants,
                            rcase, getDl,accGetDl, getZl, accGetZl, mapxytot, mapxytoti,
                            ntotxyi, myid))
        MPI_Abort(MPI_COMM_WORLD,-1);

      valarray <float> mapxytotrecv( p.npix*p.npix );
      valarray <float> mapxytotirecv[6];
      for(int i=0; i<6; i++)
        mapxytotirecv[i].resize( p.npix*p.npix );

      MPI_Reduce( &mapxytot[0],  &mapxytotrecv[0],  p.npix*p.npix, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

      for(int i=0; i<6; i++)
        MPI_Reduce( &mapxytoti[i][0],  &mapxytotirecv[i][0],  p.npix*p.npix, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

      double zsim = gsl_spline_eval (getZl, (lens.ld2[isnap]+lens.ld[isnap])/2.0, accGetZl);
      writeMaps (p, simdata, lens, isnap, zsim, snappl, p.snpix, mapxytotrecv,
                       mapxytotirecv, ntotxyi,  myid);

    }else{

      /* number of halos from snap.0 up to snap.ff. It will be modified by getGID */
      int nhalos = 0;
      for (unsigned int ff=ffmin; ff<ffmax; ff++){

        if( ifstream( fileOutput(p, "groups."+snappl, ff) ) && ifstream( fileOutput(p, "subgroups."+snappl, ff) ) ){
          if (myid==0)
            cout << fileOutput(p, snappl, ff) << " Already exists" <<endl;
          continue;
        }

        Header data;
        string file_in = File+"."+sconv(ff,fINT);
        ifstream fin;
        if (readHeader (file_in, data, fin, false))
          MPI_Abort(MPI_COMM_WORLD,-1);
        int fastforwardheader = fin.tellg();
        if(ff==0)
          printHeader(data);
        /* Creating the pointers for the Data structures */
        SubFind *halos;
        SubFind *subhalos;
        float *xx[6][3];
        float *vv[6][3];

        halos = new SubFind(data.npart[0], 1 );
        subhalos = new SubFind(data.npart[1], 0);
        xx[0][0]=&halos->xx0[0]; xx[0][1]=&halos->yy0[0]; xx[0][2]=&halos->zz0[0];
        xx[1][0]=&subhalos->xx0[0]; xx[1][1]=&subhalos->yy0[0]; xx[1][2]=&subhalos->zz0[0];
        vv[0][0]=&subhalos->vx0[0]; vv[0][1]=&subhalos->vy0[0]; vv[0][2]=&subhalos->vz0[0];
        for(int i=1; i<6;i++){
          for(int j=0;j<3;j++){
            if(i>1)
              xx[i][j]=nullptr;
            vv[i][j]=nullptr;
          }
        }

        readPos (fin,  data, p, random, isnap, xx, floor(isnap*1.0/numberOfLensPerSnap), myid);
        fin.seekg(fastforwardheader);
        readBlock(fin, data.npart[0], "MCRI", &halos->m[0], myid);
        readBlock(fin, data.npart[0], "NSUB", &halos->nsub[0], myid);
        readBlock(fin, data.npart[0], "FSUB", &halos->fsub[0], myid);
        readBlock(fin, data.npart[1], "MSUB", &subhalos->m[0], myid);
        readVel (fin, data, p, random, isnap, vv, myid);
        readBlock(fin, data.npart[1], "GRNR", &subhalos->id[0], myid);
        if(ff==ffmin){
          if(getGID(*halos, File, 0, ff, nhalos))
            MPI_Abort(MPI_COMM_WORLD,-1);
        }
        else{
          if(getGID(*halos, File, ff, ff, nhalos))
            MPI_Abort(MPI_COMM_WORLD,-1);
        }

        getGVel(*halos, p, random, File, isnap);
        getTrueZ(*halos, data, getZl, accGetZl, lens, isnap);
        getTrueZ(*subhalos, data, getZl, accGetZl, lens, isnap);
        getLOSVel(*halos);
        getLOSVel(*subhalos);
        getAngular(*halos);
        getAngular(*subhalos);

        writePLC (*halos,    data, p,    "groups."+snappl, ff);
        writePLC (*subhalos, data, p, "subgroups."+snappl, ff);

        fin.clear();
        fin.close();
        delete halos;
        delete subhalos;

      }
    }
  }

  gsl_spline_free (getDl);gsl_spline_free (getZl);
  gsl_interp_accel_free (accGetDl);gsl_interp_accel_free (accGetZl);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  cout << " end of work ... ;)  " << endl;

  return 0;
}
