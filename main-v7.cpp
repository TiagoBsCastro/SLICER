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
    cout << "   -           SLICE - Simulation LIght ConE            - " << endl;
    cout << "   -                       v1.0                         - " << endl;
    cout << "   ------------------------------------------------------ " << endl;
  }

  // Reading InputFiles
  struct InputParams p;
  string snpix;
  bool physical;

  if(readInput(&p, inifile, snpix, &physical))
    MPI_Abort(MPI_COMM_WORLD,-1);

  // Reading the Snapshots available
  vector <string> snappath; // List of SnapShots paths
  vector <double> snapred; // List of SnapShots redshift
  size_t nsnap;
  if( read_redlist(p.filredshiftlist, snapred, snappath, &p) )
    MPI_Abort(MPI_COMM_WORLD,-1);
  nsnap = snapred.size(); // Number of Snapshots

  /* Read the Snap Header and get Cosmological Values*/
  Header simdata;
  ifstream fin;
  if( read_header (p.pathsnap+snappath[0]+".0", &simdata, fin, true) )
    MPI_Abort(MPI_COMM_WORLD,-1);

  /* Creating an Instance of the cosmology class to compute distances (!!h=1!!) */
  cosmology cosmo(simdata.om0,simdata.oml,1.0,-1.0);
  /* Creating a table with redshifts and comovingdistances to be interpolated*/
  vector <double> zl(neval),dl(neval);
  for(int i=0;i<neval;i++){
    zl[i] = i * (p.zs+1.0)/(neval-1);
    dl[i] = cosmo.comovDist(zl[i])*speedcunit;
  }
  /*Initializing the auxiliary functions to get Dc given z and vice-versa*/
  gsl_interp_accel *accGetDl = gsl_interp_accel_alloc ();
  gsl_interp_accel *accGetZl = gsl_interp_accel_alloc ();
  gsl_spline *GetDl = gsl_spline_alloc (gsl_interp_cspline, neval);
  gsl_spline *GetZl = gsl_spline_alloc (gsl_interp_cspline, neval);
  gsl_spline_init (GetDl, &zl[0], &dl[0], neval);
  gsl_spline_init (GetZl, &dl[0], &zl[0], neval);
  /* Comoving distance of the last plane*/
  p.Ds = gsl_spline_eval (GetDl, p.zs, accGetDl);

  /* Creating an Instance of the Lens and Building the Planes */
  Lens lens;
  build_planes(&p, &simdata, lens, snapred, snappath, GetDl, accGetDl, GetZl, accGetZl, numberOfLensPerSnap, myid);
  /* Creating an Instance of the Randomization plan */
  Random random;
  randomize_box (random, &lens, &p, numberOfLensPerSnap, myid);

  gsl_spline_free (GetDl);gsl_spline_free (GetZl);
  gsl_interp_accel_free (accGetDl);gsl_interp_accel_free (accGetZl);

  return 0;
}
