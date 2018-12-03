#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <util.h>
#include <cstring>

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
/*                                                                           */
/*                              Co-dev.  by Tiago Castro castro@if.ufrj.br   */
/*****************************************************************************/

using namespace std;
using namespace CCfits;

int main(int argc, char** argv){

  if(argc!=2){cout << "No params!! Nothing to be done!" << endl; return -1;} // Should be run ./exe paramfile

  string inifile=argv[1]; // Paramfile name

  cout << "   ------------------------------------------------------ " << endl;
  cout << "   -                                                    - " << endl;
  cout << "   -           2D Mapping Simulation Snapshot           - " << endl;
  cout << "   -                                                    - " << endl;
  cout << "   -               collapsing one dimension             - " << endl;
  cout << "   ------------------------------------------------------ " << endl;

  // ... masses of the different type of particles
  double m [6];
  HEADER header;
  // ******************** to be read in the INPUT file ********************
  int npix;// ... project - number of pixels
  // ... set by hand the redshift of the source and than loop only up to
  // ... the needed snaphost when creating the light cone!
  double zs,Ds;
  string filredshiftlist,filsnaplist;
  string pathsnap;
  double boxl; // Mpc/h
  string idc; // comoving distance file
  int seedcenter, seedface, seedsign,nfiles;
  double fov;
  string simulation; // Simulation Name
  string partinplanes; // ALL all part in one plane, NO each part type in different planes
  string directory, suffix;
  int sn_opt; //Shot-Noise option
  bool do_NGP;

  readParameters(inifile,&npix,&boxl,&zs,&fov,
                 &filredshiftlist,&filsnaplist,
                 &pathsnap,&idc,&seedcenter,
                 &seedface,&seedsign,
                 &simulation,&nfiles,
                 &partinplanes,
                 &directory,&suffix,&sn_opt,&do_NGP);  // Reading the parameters File

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
  vector <int> lsnap; /*lsnap and lred are the selected snapshots */
  vector <double> lred;

  read_redlist(filredshiftlist, lred, lsnap, zs);
  redlist.close();

  int nsnaps = lsnap.size();

  cout << "  " << endl;
  cout << " opening path for snapshots: " << endl;
  cout << pathsnap << endl;
  cout << " " << endl;
  cout << " I will look for comoving distance " << endl;
  cout << "      file = " << idc << endl;
  cout << " " << endl;

  vector<double> zl, dl;
  read_dl(idc,zl,dl,zs);

  Ds = getY(zl,dl,zs);  // comoving distance of the last plane
  vector<int> replication; // Number of repetitions of the i-th snapshot box
  vector<int> fromsnap; // From which snapshot the i-th lens plane were build
  vector<double> zsimlens;//z of the i-th lens (z correspodent to d=1/2(ld+ld2))
  vector<double> ld; //Start position of the i-th lens
  vector<double> ld2;//End position of the i-th lens

  double h0,fovradiants;
  fovradiants = fov/180.*M_PI;
  // check if the field of view is too large with respect to the box size
  if(fovradiants*Ds>boxl){
    cout << " field view too large ... I will STOP here!!! " << endl;
    cout << " value set is = " << fov << endl;
    cout << " maximum value allowed " << boxl/Ds*180./M_PI << " in degrees " << endl;
    exit(1);
  }

  //Building the ligh cone
  plans_builder (zl, dl, boxl, zs, lsnap, lred, ld, ld2, zsimlens, fromsnap, replication);
  cout << " Comoving Distance of the last plane " << Ds << endl;
  cout << " nsnaps = " << nsnaps << endl;

  vector<int> pll;

  ofstream planelist;
  string planes_list;
  planes_list = directory+"planes_list_"+suffix+".txt";
  planelist.open(planes_list.c_str());
  for(int i=0;i<fromsnap.size();i++){

    cout << zsimlens[i] << " planes = " << ld[i] << "  " << ld2[i] << "  " << replication[i] << " from snap " << fromsnap[i] << endl;
    planelist <<  i << "   " << zsimlens[i] << "   " << ld[i] << "   " << ld2[i] << "   " << replication[i] << "   " << fromsnap[i] << endl;
    pll.push_back(i);
  }
  planelist.close();

  // randomization of the box realizations :
  int nrandom = replication.back();
  vector<double> x0(nrandom), y0(nrandom), z0(nrandom); // randomizing the center of the simulation [0,1]
  vector<int> face(nrandom); // face of the dice
  vector<int> sgnX(nrandom), sgnY(nrandom),sgnZ(nrandom); // randomizing the box axis signs

  box_randomize(x0, y0, z0, sgnX, sgnY, sgnZ, face, nrandom, seedsign, seedface, seedcenter);

  // loop on snapshots ****************************************
  cout << " now loop on " << replication.back() << " planes " << endl;
  cout << "  " << endl;
  for(int nsnap=0;nsnap<replication.back();nsnap++){

    RANDOMIZATION rand_plan (face[nsnap], sgnX[nsnap], sgnY[nsnap], sgnZ[nsnap], x0[nsnap], y0[nsnap], z0[nsnap]);

    if(ld2[nsnap]-ld[nsnap] < 0){
      cout << " comoving distance of the starting point " << ld[nsnap] << endl;
      cout << " comoving distance of the final    point " << ld2[nsnap] << endl;
      cout << " please check this out! I will STOP here!!! " << endl;
      exit(1);
    }
    float rcase = ld[nsnap]/boxl; //Number of piled boxes

    if (do_as_t11)
      npix=int((ld2[nsnap]+ld[nsnap])/2*fovradiants/rgrid*1e3)+1; //Override npix if do_as_t11 is True

    valarray<float> mapxytot( npix*npix );
    int ntotxyptype[]={0,0,0,0,0,0};
    valarray<float> mapxytotptype [6];
    for(int i=0; i<6; i++)
      mapxytotptype[i].resize(npix*npix);

    string snapnum;
    if( fromsnap[nsnap]<10) snapnum = "00"+conv(fromsnap[nsnap],fINT);
    else if(fromsnap[nsnap]>=10 && fromsnap[nsnap]<100 ) snapnum = "0"+conv(fromsnap[nsnap],fINT);
    else snapnum = conv(fromsnap[nsnap],fINT);

    string File = pathsnap+"/snapdir_"+snapnum+"/"+simulation+"_"+snapnum;
    string snappl;
    if( pll[nsnap]<10) snappl = "00"+conv(pll[nsnap],fINT);
    else if(pll[nsnap]>=10 && pll[nsnap]<100 ) snappl = "0"+conv(pll[nsnap],fINT);
    else snappl = conv(pll[nsnap],fINT);

    if(ifstream(directory+simulation+"."+snappl+".plane_"+snpix+"_"+suffix+".fits") && partinplanes == "ALL"){
      cout << directory+simulation+"."+snappl+".plane_"+snpix+"_"+suffix+".fits" << " "<< "Already exists" <<endl;
      continue; // If this lens plane was already created go to the next one
    }
    // redshift and dl of the simulation
    double zsim, dlsim;

    for (unsigned int ff=0; ff<nfiles; ff++){

      string file_in = File+"."+conv(ff,fINT);
      ifstream fin(file_in.c_str());
      cout <<"    reading the input file: "<<file_in<<endl;
      if (!fin) {cerr <<"Error in opening the file: "<<file_in<<"!\n\a"; exit(1);}

      read_header(&fin, header);
      if(header.nTotalHW[0]+header.nTotalHW[1]+header.nTotalHW[2]+header.nTotalHW[3]+header.nTotalHW[4]+header.nTotalHW[5]!=0){
	      cout << "More than 2^32 particles!! The code has to be changed! Exiting now!" << endl;
	      return -1;
      }

      bool hydro [] = {false,false,false,false,false,false};
      for(int i=0;i<=5;i++)
	      if(header.massarr[i]==0 && header.npart[i]>0){
          hydro[i]=true;
        }

      // maps for each mass type
      valarray <float> mapxy0( npix*npix );
      valarray <float> mapxy1( npix*npix );
      valarray <float> mapxy2( npix*npix );
      valarray <float> mapxy3( npix*npix );
      valarray <float> mapxy4( npix*npix );
      valarray <float> mapxy5( npix*npix );

      valarray <float> xx0 (3*header.npart[0]);
      valarray <float> xx1 (3*header.npart[1]);
      valarray <float> xx2 (3*header.npart[2]);
      valarray <float> xx3 (3*header.npart[3]);
      valarray <float> xx4 (3*header.npart[4]);
      valarray <float> xx5 (3*header.npart[5]);

      double *m0, *m1, *m2, *m3, *m4, *m5;

      if(ff==0){
        print_header(header);
  	    // compute the comoving angular diameter distance at simulation redshift
  	    zsim = header.redshift;
  	    dlsim = getY(zl,dl,zsim);
  	    h0 = header.h;
  	    cout << "  " << endl;
  	    cout << "      __________________ COSMOLOGY __________________  " << endl;
  	    cout << " " << endl;
  	    cout << "      Omegam = " << header.om0 << " " << "Omegal = " << header.oml << endl;
  	    cout << "           h = " << header.h   << " " << "BoxSize = " << header.boxsize << endl;
  	    cout << "      redshift = " << zsim <<   " " << "Dl (comoving) = " << dlsim << endl;
        cout << "      _______________________________________________  " << endl;
  	    cout << " " << endl;
        if(abs(boxl - header.boxsize/1.e+3)>1.e-2 ){
  	      cout << " set boxl and header.size differ ... check it! " << std:: endl;
  	      cout << "  boxl = " << boxl << "  " << " header.boxsize = " << header.boxsize/1.e+3 << endl;
  	      exit(1);
  	    }
  	    m[0] = header.massarr[0];
        m0=&m[0];
  	    m[1] = header.massarr[1];
        m1=&m[1];
  	    m[2] = header.massarr[2];
        m2=&m[2];
  	    m[3] = header.massarr[3];
        m3=&m[3];
  	    m[4] = header.massarr[4];
        m4=&m[4];
  	    m[5] = header.massarr[5];
        m5=&m[5];
      }

      read_pos(&fin,& header, &xx0[0],&xx1[0],&xx2[0],&xx3[0],&xx4[0],&xx5[0]);
      xx0/ (float) header.boxsize; xx1/ (float) header.boxsize;
      xx2/ (float) header.boxsize; xx3/ (float) header.boxsize;
      xx4/ (float) header.boxsize; xx5/ (float) header.boxsize;

      rand_pos (&xx0[0], header, 0, nsnap, rcase, rand_plan);
      rand_pos (&xx1[0], header, 1, nsnap, rcase, rand_plan);
      rand_pos (&xx2[0], header, 2, nsnap, rcase, rand_plan);
      rand_pos (&xx3[0], header, 3, nsnap, rcase, rand_plan);
      rand_pos (&xx4[0], header, 4, nsnap, rcase, rand_plan);
      rand_pos (&xx5[0], header, 5, nsnap, rcase, rand_plan);

      if(!max_element(hydro,hydro+6)){
           cout << endl;
	         cout << "		@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
           cout << "		@                          @" << endl;
	         cout << "		@  !!DM only simulation!!  @" << endl;
	         cout << "		@                          @" << endl;
	         cout << "		@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
           cout << endl;
      }
      else{
        cout << endl;
        cout << "		@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
	      cout << "		@                          @" << endl;
	      cout << "		@  !!Hydro   simulation!!  @" << endl;
	      cout << "		@                          @" << endl;
	      cout << "		@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
        cout << endl;

        if (hydro[0] && header.npart[0])
          m0 = new double [header.npart[0]];
        if (hydro[1] && header.npart[1])
          m1 = new double [header.npart[1]];
        if (hydro[2] && header.npart[2])
          m2 = new double [header.npart[2]];
        if (hydro[3] && header.npart[3])
          m3 = new double [header.npart[3]];
        if (hydro[4] && header.npart[4])
          m4 = new double [header.npart[4]];
        if (hydro[5] && header.npart[5])
          m5 = new double [header.npart[5]];

        read_mass(&fin, &header, &m0[0], &m1[0], &m2[0], &m3[0], &m4[0], &m5[0]);
        read_block(&fin, &m5[0], "BHMA");
      }

      int n0 = header.npart[0];
      int n1 = header.npart[1];
      int n2 = header.npart[2];
      int n3 = header.npart[3];
      int n4 = header.npart[4];
      int n5 = header.npart[5];

      cout << "  " << endl;
      cout << n0 <<"   type (0) particles read"<<endl;
      cout << n1 <<"   type (1) particles read"<<endl;
      cout << n2 <<"   type (2) particles read"<<endl;
      cout << n3 <<"   type (3) particles read"<<endl;
      cout << n4 <<"   type (4) particles read"<<endl;
      cout << n5 <<"   type (5) particles read"<<endl;
      cout << "  " << endl;

      if(n0>0)
        map_particles (&xx0[0], &m[0], 0, sn_opt, fovradiants, header, ld,
          ld2, nsnap, mapxy0, &ntotxyptype[0], do_NGP, hydro[0]);

      if(n1>0)
        map_particles (&xx1[0], &m[1], 1, sn_opt, fovradiants, header, ld,
          ld2, nsnap, mapxy1, &ntotxyptype[1], do_NGP, hydro[1]);

      if(n2>0)
        map_particles (&xx2[0], &m[2], 2, sn_opt, fovradiants, header, ld,
          ld2, nsnap, mapxy2, &ntotxyptype[2], do_NGP, hydro[2]);

      if(n3>0)
        map_particles (&xx3[0], &m[3], 3, sn_opt, fovradiants, header, ld,
          ld2, nsnap, mapxy3, &ntotxyptype[3], do_NGP, hydro[3]);

      if(n4>0)
        map_particles (&xx4[0], &m[4], 4, sn_opt, fovradiants, header, ld,
          ld2, nsnap, mapxy4, &ntotxyptype[4], do_NGP, hydro[4]);

      if(n5>0)
        map_particles (&xx5[0], &m[5], 5, sn_opt, fovradiants, header, ld,
          ld2, nsnap, mapxy5, &ntotxyptype[5], do_NGP, hydro[5]);

      cout << " maps done! " << endl;

      fin.clear();
	    fin.close(); // Close file here

      if (hydro[0] && header.npart[0])
        delete[] m0;
      if (hydro[1] && header.npart[1])
        delete[] m1;
      if (hydro[2] && header.npart[2])
        delete[] m2;
      if (hydro[3] && header.npart[3])
        delete[] m3;
      if (hydro[4] && header.npart[4])
        delete[] m4;
      if (hydro[5] && header.npart[5])
        delete[] m5;

      // sum up all the maps
      mapxytot+=mapxy0+mapxy1+mapxy2+mapxy3+mapxy4+mapxy5;
      mapxytotptype[0]+=mapxy0;
      mapxytotptype[1]+=mapxy1;
      mapxytotptype[2]+=mapxy2;
      mapxytotptype[3]+=mapxy3;
      mapxytotptype[4]+=mapxy4;
      mapxytotptype[5]+=mapxy5;

      cout << " done map*tot " << endl;
    }

    write_map (partinplanes, directory, simulation, snappl, snpix, suffix, nsnap, ntotxyptype,
          mapxytot, mapxytotptype, fov, zsim, ld, ld2, h0, header.om0, header.oml, m);
  }

  cout << " end of work ... ;-)  " << endl;
  return 0;
}
