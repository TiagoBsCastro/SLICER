#include "mpi.h"
#include "densitymaps.h"
#define numberOfLensPerSnap  1 // Number of Lens to be builded from a snap
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
    cout << "   -                       v1.0                         - " << endl;
    cout << "   ------------------------------------------------------ " << endl;
  }

  // Reading InputFiles
  struct InputParams p;
  string snpix;
  bool physical;
  double fovradiants;
  size_t nsnap;

  if(readInput(p, inifile, snpix, physical))
    MPI_Abort(MPI_COMM_WORLD,-1);

  // Reading the Snapshots available
  vector <string> snappath; // List of SnapShots paths
  vector <double> snapred; // List of SnapShots redshift
  if( read_redlist(p.filredshiftlist, snapred, snappath, p) )
    MPI_Abort(MPI_COMM_WORLD,-1);
  nsnap = snapred.size(); // Number of Snapshots

  /* Read the Snap Header and get Cosmological Values*/
  Header simdata;
  ifstream fin;
  if( read_header (p.pathsnap+snappath[0]+".0", simdata, fin, true) )
    MPI_Abort(MPI_COMM_WORLD,-1);
  test_hydro(p, simdata);

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
  gsl_spline *GetDl = gsl_spline_alloc (gsl_interp_cspline, neval);
  gsl_spline *GetZl = gsl_spline_alloc (gsl_interp_cspline, neval);
  gsl_spline_init (GetDl, &zl[0], &dl[0], neval);
  gsl_spline_init (GetZl, &dl[0], &zl[0], neval);
  /* Comoving distance of the last plane*/
  p.Ds = gsl_spline_eval (GetDl, p.zs, accGetDl);
  if(test_fov(p.fov, simdata.boxsize/1e3, p.Ds, myid, fovradiants))
    MPI_Abort(MPI_COMM_WORLD,-1);

  /* Creating an Instance of the Lens and Building the Planes */
  Lens lens;
  build_planes(p, simdata, lens, snapred, snappath, GetDl, accGetDl, GetZl, accGetZl, numberOfLensPerSnap, myid);
  /* Creating an Instance of the Randomization plan */
  Random random;
  randomize_box (random, lens, p, numberOfLensPerSnap, myid);

  /* Looping on the Snapshots */
  if(myid==0){
    cout << " Now loop on " << lens.nplanes << " planes " << endl;
    cout << "  " << endl;
  }
  for(int isnap=0; isnap < lens.nplanes; isnap++){

    float rcase = floor(lens.ld[isnap]/simdata.boxsize*1e3);

    /* Override p.npix if physical is True */
    if (physical)
      p.npix=int( (lens.ld2[isnap]+lens.ld[isnap])/2*fovradiants/p.rgrid*1e3 )+1;

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
    if(ifstream(p.directory+p.simulation+"."+snappl+".plane_"+snpix+"_"+p.suffix+".fits") && p.partinplanes == false && p.simType == "Gadget"){
      if (myid==0)
        cout << p.directory+p.simulation+"."+snappl+".plane_"+snpix+"_"+p.suffix+".fits" << " "<< "Already exists" <<endl;
      continue;
    }

    /* Computing Working Balance */
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
    if(p.simType == "Gadget"){
      valarray<float> mapxytot;
      int ntotxyi[6];
      valarray<float> mapxytoti[6];

      if( CreateDensityMaps (p, lens, random, isnap, ffmin, ffmax, File, fovradiants,
                            rcase, GetDl,accGetDl, GetZl, accGetZl, mapxytot, mapxytoti,
                            ntotxyi, myid))
        MPI_Abort(MPI_COMM_WORLD,-1);

      valarray <float> mapxytotrecv( p.npix*p.npix );
      valarray <float> mapxytotirecv[6];
      for(int i=0; i<6; i++)
        mapxytotirecv[i].resize( p.npix*p.npix );

      MPI_Reduce( &mapxytot[0],  &mapxytotrecv[0],  p.npix*p.npix, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

      for(int i=0; i<6; i++)
        MPI_Reduce( &mapxytoti[i][0],  &mapxytotirecv[i][0],  p.npix*p.npix, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

      double zsim = gsl_spline_eval (GetZl, (lens.ld2[isnap]+lens.ld[isnap])/2.0, accGetZl);
      write_maps (p, simdata, lens, isnap, zsim, snappl, snpix, mapxytotrecv,
                       mapxytotirecv, ntotxyi,  myid);

    }else{

      for (unsigned int ff=ffmin; ff<ffmax; ff++){

        Header data;
        string file_in = File+"."+sconv(ff,fINT);
        ifstream fin;
        if (read_header (file_in, data, fin, false))
          MPI_Abort(MPI_COMM_WORLD,-1);
        int fastforwardheader = fin.tellg();
        if(ff==0)
          print_header(data);
        /* Creating the pointers for the Data structures */
        SubFind *halos;
        SubFind *subhalos;
        float *xx[6][3];
        float *vv[6][3];

        halos = new SubFind(data.npart[0]);
        subhalos = new SubFind(data.npart[1]);
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

        ReadPos (fin,  data, p, random, isnap, xx, rcase, myid);
        fin.seekg(fastforwardheader);
        ReadBlock(fin, data.npart[0], "MCRI", &halos->m[0], myid);
        ReadBlock(fin, data.npart[0], "NSUB", &halos->nsub[0], myid);
        ReadBlock(fin, data.npart[1], "MSUB", &subhalos->m[0], myid);
        ReadVel (fin, data, p, random, isnap, vv, myid);
        ReadBlock(fin, data.npart[1], "GRNR", &subhalos->id[0], myid);

        if(GetGID(*halos, File, ff))
          MPI_Abort(MPI_COMM_WORLD,-1);
        GetGVel(*halos, *subhalos, p, random, File, isnap);
        GetTrueZ(*halos, data, GetZl, accGetZl);
        GetTrueZ(*subhalos, data, GetZl, accGetZl);
        GetLOSVel(*halos);
        GetLOSVel(*subhalos);
        GetAngular(*halos);
        GetAngular(*subhalos);

        CreatePLC (*halos,    data, p,    "groups."+snappl, ff);
        CreatePLC (*subhalos, data, p, "subgroups."+snappl, ff);

        fin.clear();
        fin.close();
        delete halos;
        delete subhalos;

      }
    }
  }

  gsl_spline_free (GetDl);gsl_spline_free (GetZl);
  gsl_interp_accel_free (accGetDl);gsl_interp_accel_free (accGetZl);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  cout << " end of work ... ;)  " << endl;

  return 0;
}
