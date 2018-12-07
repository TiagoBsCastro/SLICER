#include "mpi.h"
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <CCfits/CCfits>
#include <util.h>
#include <utilities.h>
#include <cstring>
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
/*      - Re-writeen for better maintenance                                  */
/*                                                                           */
/*                        adapted  by Tiago Castro tiagobscastro@gmail.com   */
/*****************************************************************************/

using namespace std;
using namespace CCfits;

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
  // ******************** to be read in the INPUT file ********************
  // ... project - number of pixels
  double Ds;

  struct InputParams p;
  readInput(&p, inifile);

  if(p.snopt<0 && myid==0){
    cout << "Impossible value for Shot-Noise option!" << endl;
    MPI_Abort(MPI_COMM_WORLD,-1);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  string snpix;
  bool do_as_t11=(p.npix<0);
  int rgrid;
  if(!do_as_t11) snpix=sconv(p.npix,fINT);
  else{
    int n = -p.npix;
    snpix=sconv(n,fINT);
    snpix+="_kpc";
    rgrid=n;
  }

  // ... read the redshift list and the snap_available
  // ... to build up the light-cone
  vector <string> lsnap; //lsnap and lred are the selected snapshots selected (z<zs + the first snapshot deeper than zs)
  vector<double> lred;
  if( read_redlist(p.filredshiftlist, lred, lsnap, p.zs) )
    MPI_Abort(MPI_COMM_WORLD,-1);

  int nsnaps = lsnap.size();
  if (myid==0){
    cout << "  " << endl;
    cout << " opening path for snapshots: " << endl;
    cout << p.pathsnap << endl;
  }

  /* Read the Snap Header and get Cosmological Values*/
  string file_in = p.pathsnap+lsnap[0]+"."+sconv(myid,fINT);
  cout << file_in << endl;
  Header header;
  ifstream fin;
  if( read_header (file_in, &header, fin, true) )
    MPI_Abort(MPI_COMM_WORLD,-1);
  /* Creating an Instance of the cosmology class to compute distances (!!h=1!!) */
  double h0 = header.h;
  double om0 = header.om0;
  double omL0 = header.oml;
  cosmology cosmo(header.om0,header.oml,1.0,-1.0);

  vector <double> zl(neval),dl(neval);   // Redshift and ComovDistance array to be interpolated
  vector <int>    replication;  // Number of repetitions of the i-th snapshot box
  vector <string> fromsnap;     // From which snapshot the i-th lens plane were build
  vector <double> zsimlens;     // z of the i-th lens (z correspodent to d=1/2(ld+ld2))
  vector <double> ld;           // Start position of the i-th lens
  vector <double> ld2;          // End position of the i-th lens
  vector <double> zfromsnap;    // z from the snapshot selected to build the i-th lens
  vector <bool> randomize;      // Bool variable to whether the positions should be re-randomized or not

  /* Creating a table with redshifts and comovingdistances to be interpolated*/
  for(int i=0;i<neval;i++){

    zl[i] = i * (p.zs+1.0)/(neval-1);
    dl[i] = cosmo.comovDist(zl[i])*speedcunit;

  }
  Ds = getY(zl,dl,p.zs);          // comoving distance of the last plane
  vector<int> pll;

  build_plans(Ds, &p, numberOfLensPerSnap, nsnaps, lred, zl, dl, lsnap, ld, ld2, replication, zfromsnap, fromsnap, zsimlens,randomize, pll, myid);

  int nrandom = replication.back();
  vector<double> x0(nrandom), y0(nrandom), z0(nrandom); // ramdomizing the center of the simulation [0,1]
  vector<int> face(nrandom); // face of the dice
  vector<int> sgnX(nrandom), sgnY(nrandom),sgnZ(nrandom); // randomizing the box axis signs

  // randomization of the box realizations :
  randomize_box (x0, y0, z0, face, sgnX, sgnY, sgnZ, replication, nrandom, randomize, &p, numberOfLensPerSnap, myid);

  if(myid==0)
    cout << " set the field of view to be square in degrees " << endl;
  double fovradiants;
  fovradiants = p.fov/180.*M_PI;
  // check if the field of view is too large with respect to the box size
  if(fovradiants*Ds>p.boxl && myid==0){
    cout << " field view too large ... I will STOP here!!! " << endl;
    cout << " value set is = " << p.fov << endl;
    cout << " maximum value allowed " << p.boxl/Ds*180./M_PI << " in degrees " << endl;
    MPI_Abort(MPI_COMM_WORLD,-1);
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
      MPI_Abort(MPI_COMM_WORLD,-1);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    float rcase = floor(ld[nsnap]/p.boxl);
    if (do_as_t11)
      p.npix=int((ld2[nsnap]+ld[nsnap])/2*fovradiants/rgrid*1e3)+1; //Override p.npix if do_as_t11 is True
    valarray<float> mapxytot( p.npix*p.npix );
    int ntotxyi[6];
    valarray<float> mapxytoti[6];
    for(int i=0;i<6;i++){
      ntotxyi[i]=0;
      mapxytoti[i].resize( p.npix*p.npix );
    }
    string File = p.pathsnap+fromsnap[nsnap];
    string snappl;
    if( pll[nsnap]<10)
      snappl = "00"+sconv(pll[nsnap],fINT);
    else if(pll[nsnap]>=10 && pll[nsnap]<100 )
      snappl = "0"+sconv(pll[nsnap],fINT);
    else
      snappl = sconv(pll[nsnap],fINT);
    string Filesub;

    if(ifstream(p.directory+p.simulation+"."+snappl+".plane_"+snpix+"_"+p.suffix+".fits") && p.partinplanes == false){
      if (myid==0)
        cout << p.directory+p.simulation+"."+snappl+".plane_"+snpix+"_"+p.suffix+".fits" << " "<< "Already exists" <<endl;
      continue; // If this lens plane was already created go to the next one
    }

    float num_float1, num_float2, num_float3; // aux variables to read x,y, and z positions

    // redshift and dl of the simulation
    double zsim, dlsim;

    //Computing Working Balance (!!THIS IS A VERY STUPID WORKBALANCE!!)
    int intdiv, remaindiv,ffmin,ffmax;
    intdiv=header.numfiles/numprocs;
    remaindiv=header.numfiles%numprocs;

    if(myid!=numprocs-1){
      ffmin=myid*intdiv;
      ffmax=(myid+1)*intdiv;
    }else{
      ffmin=myid*intdiv;
      ffmax=(myid+1)*intdiv+remaindiv;
    }
    Header data;
    for (unsigned int ff=ffmin; ff<ffmax; ff++){

      string file_in = File+"."+sconv(ff,fINT);
      ifstream fin;
      if (read_header (file_in, &data, fin, false))
        MPI_Abort(MPI_COMM_WORLD,-1);

      // map for each mass type
      valarray<float> mapxyi[6];
      for(int i=0; i<6; i++)
        mapxyi[i].resize(p.npix*p.npix);

      // GADGET has 6 different particle type
      vector<float> xx0(data.npart[0]), yy0(data.npart[0]), zz0(data.npart[0]);
      vector<float> xx1(data.npart[1]), yy1(data.npart[1]), zz1(data.npart[1]);
      vector<float> xx2(data.npart[2]), yy2(data.npart[2]), zz2(data.npart[2]);
      vector<float> xx3(data.npart[3]), yy3(data.npart[3]), zz3(data.npart[3]);
      vector<float> xx4(data.npart[4]), yy4(data.npart[4]), zz4(data.npart[4]);
      vector<float> xx5(data.npart[5]), yy5(data.npart[5]), zz5(data.npart[5]);
      float *xx[6][3];
      xx[0][0]=&xx0[0]; xx[0][1]=&yy0[0]; xx[0][2]=&zz0[0];
      xx[1][0]=&xx1[0]; xx[1][1]=&yy1[0]; xx[1][2]=&zz1[0];
      xx[2][0]=&xx2[0]; xx[2][1]=&yy2[0]; xx[2][2]=&zz2[0];
      xx[3][0]=&xx3[0]; xx[3][1]=&yy3[0]; xx[3][2]=&zz3[0];
      xx[4][0]=&xx4[0]; xx[4][1]=&yy4[0]; xx[4][2]=&zz4[0];
      xx[5][0]=&xx5[0]; xx[5][1]=&yy5[0]; xx[5][2]=&zz5[0];

      if(ff==0)
        print_header(data);
      fastforwardToPos(fin, NBLOCKSTOPOS, myid, false);

      // total number of particles as the sum of all of them
      int dim = data.npart[0]+data.npart[1]+data.npart[2]+data.npart[3]+data.npart[4]+data.npart[5];
      int dimmass0=0;
      for(int i=0;i<=5;i++){
	       if(data.massarr[i]==0)
           dimmass0+=data.npart[i];
      }
      bool hydro=bool(dimmass0);

      if(ff==ffmin){
	        // compute the comoving angular diameter distance at simulation redshift
	        zsim = data.redshift;
	        dlsim = getY(zl,dl,zsim);

          if(myid == 0)
            cout << "      redshift = " << zsim <<   " " << "Dl (comoving) = " << dlsim << endl;

          if(abs(p.boxl - data.boxsize*POS_U/1.e+3)>1.e-2 && myid==0){
	           cout << " set boxl and data.size differ ... check it! " << std:: endl;
	           cout << "  boxl = " << p.boxl << "  " << " data.boxsize = " << data.boxsize/1.e+3*POS_U << endl;
	           MPI_Abort(MPI_COMM_WORLD,-1);
	        }

          MPI_Barrier(MPI_COMM_WORLD);

      }
      // Loop on different types
      for (int i = 0; i<6; i++){

	      for (int pp=0; pp<data.npart[i]; pp++){
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
           xx[i][0][pp] = x;
           xx[i][1][pp] = y;
           xx[i][2][pp] = z;
	      }
      }
      /* If Hydro run I have to read the masses, otherwise close the snapshot*/
      if(hydro)
         fastforwardToMASS (fin, NBLOCKSTOMASS, &data, myid);
      else{
        fin.clear();
        fin.close();
      }

      if (myid==0){

	      cout << "  " << endl;
	      cout << data.npart[0] <<"   type (0) particles selected until now"<<endl;
	      cout << data.npart[1] <<"   type (1) particles selected until now"<<endl;
	      cout << data.npart[2] <<"   type (2) particles selected until now"<<endl;
	      cout << data.npart[3] <<"   type (3) particles selected until now"<<endl;
	      cout << data.npart[4] <<"   type (4) particles selected until now"<<endl;
	      cout << data.npart[5] <<"   type (5) particles selected until now"<<endl;
	      cout << "  " << endl;

      }

      int totPartxyi[6];
      for(int i=0; i<6; i++){

        size_t n = data.npart[i];

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
  	        cout << "xmin = " << xmin << endl;
	          cout << "xmax = " << xmax << endl;
	          cout << "ymin = " << ymin << endl;
	          cout << "ymax = " << ymax << endl;
	          cout << "zmin = " << zmin << endl;
	          cout << "zmax = " << zmax << endl;
	          cout << "  0 type check this!!! I will STOP here!!! " << endl;
	          cout << "Aborting from Rank "<< myid << endl;
		        MPI_Abort(MPI_COMM_WORLD,-1);
	        }
          if (myid==0){
	           cout << " ... mapping type "<< i <<" particles on the grid with " << p.npix << " pixels" << endl;
             cout << "Min distance: "<< ld[nsnap]/data.boxsize*1.e+3/POS_U<< " "<<ld2[nsnap]/data.boxsize*1.e+3/POS_U << endl;
             cout << "Rcase       : "<< rcase << endl;
          }

          vector<float> xs(0),ys(0),ms(0);
	        for(int l=0;l<data.npart[i];l++){

            if(hydro && data.massarr[i]==0){

              if(l==0 && i==5)
                fastforwardToBHMASS (fin, NBLOCKSTOBHMASS, &data, myid);

 			        fin.read((char *)&num_float1, sizeof(num_float1));
              if (num_float1>MAX_M)
                num_float1=0;
            }
		        else
              num_float1=data.massarr[i];

	          double di = sqrt(pow(xx[i][0][l]-0.5,2) + pow(xx[i][1][l]-0.5,2) + pow(xx[i][2][l],2))*data.boxsize/1.e+3*POS_U;
	          if(di>=ld[nsnap] && di<ld2[nsnap]){

  	          double rai,deci,dd;
	            getPolar(xx[i][0][l]-0.5,xx[i][1][l]-0.5,xx[i][2][l],&rai,&deci,&dd);

	            if(fabs(rai)<=fovradiants*(1.+2./p.npix)*0.5 && fabs(deci)<=fovradiants*(1.+2./p.npix)*0.5){
                xs.push_back(deci/fovradiants+0.5);
	              ys.push_back(rai/fovradiants+0.5);
                if(p.snopt==0){
                  ms.push_back(num_float1);
                }
                else{
                  if(rand()/ float(RAND_MAX) < 1./pow(2,p.snopt)) ms.push_back(pow(2,p.snopt)*num_float1);
                  else ms.push_back(0.);
                }
	            }
	          }
	        }
	        totPartxyi[i]=xs.size();
	        ntotxyi[i]+=totPartxyi[i];

	        if(totPartxyi[i]>0){
            mapxyi[i] = gridist_w(xs,ys,ms,p.npix,DO_NGP);
          }

	        //re-normalize to the total mass!
	        double mtot=0.;
          double mnorm=0.;
          for(int i=0;i<ms.size();i++)
            mnorm+=ms[i];

  	      if(totPartxyi[i]>0){
	          for(int l=0;l<p.npix*p.npix;l++)
              mtot += mapxyi[i][l];
	          if(mtot==0.)
              mtot=1.; //To avoid NaN
	          for(int l=0;l<p.npix*p.npix;l++) mapxyi[i][l]*=mnorm/mtot;
	        }
        }
      }

      if(hydro){
        Block block;
        fin >> block;
        if (myid==0){
          cout << endl << "If masses were read correctly next box should be BHMD " << endl;
          cout << "Name:                                            ";
          cout << block.name[0];
          cout << block.name[1];
          cout << block.name[2];
          cout << block.name[3] << endl;
        }
	      fin.clear();
	      fin.close(); // Close file here
      }
      // sum up all the maps
      mapxytot+=mapxyi[0]+mapxyi[1]+mapxyi[2]+mapxyi[3]+mapxyi[4]+mapxyi[5];
      for(int i=0; i<6; i++)
        mapxytoti[i]+=mapxyi[i];

      if (myid==0)
        cout << " done map*tot " << endl;
    }

    valarray<float> mapxytotrecv( p.npix*p.npix );
    valarray<float> mapxytotirecv[6];
    for(int i=0; i<6; i++)
      mapxytotirecv[i].resize( p.npix*p.npix );

    cout << " maps done! from Rank:" << myid << endl;

    MPI_Reduce( &mapxytot[0],  &mapxytotrecv[0],  p.npix*p.npix, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    for(int i=0; i<6; i++)
      MPI_Reduce( &mapxytoti[i][0],  &mapxytotirecv[i][0],  p.npix*p.npix, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (myid==0){
      if(p.partinplanes==false){
	       /*
	        * write image array(s) to FITS files all particles in a FITS file!
	        */
	        long naxis = 2;
	        long naxes[2]={ p.npix,p.npix };
	        string fileoutput;
	        fileoutput = p.directory+p.simulation+"."+snappl+".plane_"+snpix+"_"+p.suffix+".fits";
          cout << "Saving the maps on: " << fileoutput << endl;
	        unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
	        vector<long> naxex( 2 );
	        naxex[0]=p.npix;
	        naxex[1]=p.npix;
	        PHDU *phxy=&ffxy->pHDU();
	        phxy->write( 1, p.npix*p.npix, mapxytotrecv );
	        phxy->addKey ("REDSHIFT",zsim," ");
	        phxy->addKey ("PHYSICALSIZE",p.fov," ");
	        phxy->addKey ("PIXELUNIT",1.e+10/h0,"Mass unit in M_Sun");
	        phxy->addKey ("DlLOW",ld[nsnap]/h0,"comoving distance in Mpc");
	        phxy->addKey ("DlUP",ld2[nsnap]/h0,"comoving distance in Mpc");
	        phxy->addKey ("nparttype0",ntotxyi[0]," ");
	        phxy->addKey ("nparttype1",ntotxyi[1]," ");
	        phxy->addKey ("nparttype2",ntotxyi[2]," ");
	        phxy->addKey ("nparttype3",ntotxyi[3]," ");
	        phxy->addKey ("nparttype4",ntotxyi[4]," ");
	        phxy->addKey ("nparttype5",ntotxyi[5]," ");
	        phxy->addKey ("HUBBLE",h0," ");
	        phxy->addKey ("OMEGAMATTER",om0," ");
	        phxy->addKey ("OMEGALAMBDA",omL0," ");
	        phxy->addKey ("m0",data.massarr[0]," ");
	        phxy->addKey ("m1",data.massarr[1]," ");
	        phxy->addKey ("m2",data.massarr[2]," ");
	        phxy->addKey ("m3",data.massarr[3]," ");
	        phxy->addKey ("m4",data.massarr[4]," ");
	        phxy->addKey ("m5",data.massarr[5]," ");
      }else{
	    /**
	     * write image array(s) to FITS files each particle type in different planes
	     */
       for(int i=0; i<6; i++){

         if(ntotxyi[i]>0){
  	       long naxis = 2;
  	       long naxes[2]={ p.npix,p.npix };
  	       string fileoutput;
           fileoutput = p.directory+p.simulation+"."+snappl+".ptype"+sconv(i,fINT)+"_plane_"+snpix+"_"+p.suffix+".fits";
  	       unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
  	       vector<long> naxex( 2 );
  	       naxex[0]=p.npix;
  	       naxex[1]=p.npix;
  	       PHDU *phxy=&ffxy->pHDU();
  	       phxy->write( 1, p.npix*p.npix, mapxytotirecv[i] );
   	       phxy->addKey ("REDSHIFT",zsim," ");
  	       phxy->addKey ("PHYSICALSIZE",p.fov," ");
  	       phxy->addKey ("PIXELUNIT",1.e+10/h0,"Mass unit in M_Sun");
  	       phxy->addKey ("DlLOW",ld[nsnap]/h0,"comoving distance in Mpc");
  	       phxy->addKey ("DlUP",ld2[nsnap]/h0,"comoving distance in Mpc");
  	       phxy->addKey ("nparttype0",ntotxyi[i]," ");
  	       phxy->addKey ("HUBBLE",h0," ");
  	       phxy->addKey ("OMEGAMATTER",om0," ");
  	       phxy->addKey ("OMEGALAMBDA",omL0," ");
  	       phxy->addKey ("m0",data.massarr[i]," ");
  	     }
       }
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  cout << " end of work ... ;-)  " << endl;
  return 0;
}
