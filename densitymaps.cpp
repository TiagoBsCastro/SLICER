#include "densitymaps.h"
#include "gadget2io.h"

/*
 * Reads the vector with snapshots redshifts "zsnap" and returns the position
 * of the snapshot that the comoving distance computed on its redshift is the
 * closest to "dlens". Comoving distance is Computed with the GetDl interpolated
 * function.
*/
int getSnap (vector <double> & zsnap, gsl_spline *GetDl,
                    gsl_interp_accel *accGetDl, double dlens){

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

/*
 * Creates the lens Planes plan according to the InputParams "p", the sim.
 * specification contained in Header "header". The plan is stored in the
 * "lens" Structure
 *
 * - snapred is the vector with the snapshots redshifts
 * - snappath is the vector with the snapshots paths
 * - GetDl, accGetDl, GetZl, and accGetZl are auxiliary interp. func. to compute
 *   the comoving distance given z and vice-versa
 * - Each lens will be BoxSize/numberOfLensPerSnap thick
 * - if myid == 0: monitoring messages are produced
 */
void buildPlanes(InputParams &p, Lens &lens,
  vector <double> & snapred, vector <string> & snappath, vector <double> & snapbox,
  gsl_spline *GetDl, gsl_interp_accel *accGetDl, gsl_spline *GetZl, gsl_interp_accel *accGetZl,
  int numOfLensPerSnap, int myid){

  size_t nsnaps = snapred.size();
  int pos=0;
  int nrepi = 0;
  int nrep = 0;
  double zdbut, ldbut = 0.0;

  do{
    nrep++;
    nrepi++;
    double ztest=9999;
    int pos_temp = pos;
    // checking which snapshot is the closest to the next lens plane
    for(int i=pos_temp; i<nsnaps; i++){
      double dtest = ldbut + snapbox[i]/1e3/numOfLensPerSnap;
      int    itest = getSnap(snapred, GetDl, accGetDl, dtest);
      double dz    = fabs( snapred[itest]-gsl_spline_eval (GetZl, dtest, accGetZl) );
      if(dz < ztest){
        if( nrep==1 || ( !bool((nrep-1)%numOfLensPerSnap) || snapbox[itest] == snapbox[pos] ) ){
          pos_temp = itest;
          ztest    = dz;
        }
      }
    }
    ldbut += snapbox[pos_temp]/1e3/numOfLensPerSnap;
    zdbut  = gsl_spline_eval (GetZl, ldbut, accGetZl);
    double dlens = ldbut-0.5*snapbox[pos_temp]/1e3/numOfLensPerSnap;
    double zlens = gsl_spline_eval (GetZl, dlens, accGetZl);
    pos_temp = getSnap(snapred, GetDl, accGetDl, dlens);
    if (myid==0)
      cout << " simulation snapshots = " << ldbut << "  " << zdbut << "  " << nrep << " from snap "
                                                       << snappath[pos_temp] << "  " << zlens << endl;
    lens.ld.push_back(ldbut-snapbox[pos_temp]/1e3/numOfLensPerSnap);
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
    lens.fromsnapi.push_back(pos);
    if( nrep==1 )
      lens.randomize.push_back(1);
    else
      lens.randomize.push_back( !( (nrep-1)%numOfLensPerSnap ) );

  }while(ldbut<p.Ds);

  for ( int i=0; i<nrepi+1; i++ ) lens.replication.push_back(nrep); // Last plane replications
  if (myid==0){
    cout << " Comoving Distance of the last plane " << p.Ds << endl;
    cout << " nsnaps = " << nsnaps << "\n" << endl;
  }

  ofstream planelist;
  string planes_list;
  planes_list = p.directory+"planes_list_"+p.suffix+".txt";

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

/*
 * Creates the Randomization plan according to the "lens" planes
 * and the params inside the InputParams "p". The plan is stored
 * it in the "random" structure.
 * - Each lens will be BoxSize/numberOfLensPerSnap thick
 * - if myid == 0: monitoring messages are produced
 *
 */
void randomizeBox (Random & random, Lens & lens, InputParams & p,
                                int numOfLensPerSnap, int myid){

  size_t nrandom = lens.replication.back();
  /* Inflating Random Structure */
  random.x0.resize(nrandom); random.y0.resize(nrandom); random.z0.resize(nrandom);
  random.sgnX.resize(nrandom); random.sgnY.resize(nrandom); random.sgnZ.resize(nrandom);
  random.face.resize(nrandom);
  for(int i=0;i<nrandom;i++){

    if ( lens.randomize[i] ){

      srand(p.seedcenter+i/numOfLensPerSnap*13);
      random.x0[i] = rand() / float(RAND_MAX);
      random.y0[i] = rand() / float(RAND_MAX);
      random.z0[i] = rand() / float(RAND_MAX);
      if(myid==0){
        cout << "  " << endl;
        cout << " random centers  for the box " << i << " = " << random.x0[i] << "  " << random.y0[i] << "  " << random.z0[i] << endl;
      }
      random.face[i] = 7;
      srand(p.seedface+i/numOfLensPerSnap*5);
      while(random.face[i]>6 || random.face[i]<1) random.face[i] = int(1+rand() / float(RAND_MAX)*5.+0.5);
      if (myid==0)
        cout << " face of the dice " << random.face[i] << endl;
      random.sgnX[i] = 2;
      srand(p.seedsign+i/numOfLensPerSnap*8);
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


/*
 * Test if the chosen angular aperture is allowed.
 * !! SLICER do not allow repetitions of the Box in the perpendicular plane!!
 * !! to the PLC axis unless THE DIRECTIVE ReplicationOnPerpendicularPlane  !!
 */
int testFov(double fov, double boxl, double Ds, int myid, double & fovradiants){

  fovradiants = fov/180.*M_PI;
  /* check if the field of view is too large with respect to the box size */
  if( (fovradiants)*Ds>boxl && myid==0 ){
    cerr << " !!Field view too large!!\n !!!I will STOP here!!! " << endl;
    cerr << " Value set is = " << fov << endl;
    cerr << " Maximum value allowed " << boxl/Ds*180./M_PI << " in degrees " << endl;
    cerr << " For the lens at "<< Ds << endl;
    return 1;
  }
  return 0;
}

/*
 * Compute the number of replications on the perpendicular plane are necessary
 * !!!! ONLY USED IF THE DIRECTIVE ReplicationOnPerpendicularPlane is defined !!!!
 */
 void computeReplications(double fov, double boxl, double Ds, int myid, double & fovradiants, int & nrepperp){

  fovradiants = fov/180.*M_PI;
  if(Ds>0)
    nrepperp = ceil( Ds * tan(fovradiants/2) / boxl * 2 ) - 1;
  else
    nrepperp = 0;

}

/*
 * Maps the Particles on the grids "mapxyi[6]" using the TSC MAS
 * (unless USE_DGP is True) according to the InputParams "p", the
 * sim. specifications read from Header "data", the lens plan "lens"
 * the particle positions at "xx".
 *
 * - fovradiants is the PLC angular aperture in radians
 * - isnap is the number of boxes repetitions so far
 * - fin is the fstream snapshot file. Used to read masses for hydro part.
 * - if myid == 0: monitoring messages are produced
 *
 */
int mapParticles(ifstream & fin, Header &data, InputParams &p, Lens &lens,
    float* xx[6][3], double fovradiants, int isnap, valarray<float>(& mapxyi)[6],
    int(& ntotxyi)[6], int myid){

  int imax = 6;
  float num_float1;
  int totPartxyi[6];
  for(int i=0; i<imax; i++){

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
      double minDist = isnap*1.0/numberOfLensPerSnap;
      double maxDist = isnap/numberOfLensPerSnap + (isnap%numberOfLensPerSnap + 1) * (lens.ld2[isnap]-lens.ld[isnap])/data.boxsize*1.e+3/POS_U;
      if (myid==0){
         cout << " Mapping type "<< i <<" particles on the grid with " << p.npix << " pixels" << endl;
         cout << "Distance Range: "<< minDist << " " << maxDist << endl;
      }

      vector<float> xs(0),ys(0),ms(0);
      for(int l=0;l<data.npart[i];l++){

        if(p.hydro && data.massarr[i]==0){

          if(l==0 && i==5){
            fastforwardNVars (fin, sizeof(int32_t), data.npart[5]);
            fastforwardToBlock (fin, "BHMA", myid);
          }

          fin.read((char *)&num_float1, sizeof(num_float1));
          if (num_float1>MAX_M)
            num_float1=0;
        }
        else
          num_float1=data.massarr[i];

        if(xx[i][2][l]>=minDist && xx[i][2][l]<maxDist){
          for(int ni = -lens.nrepperp[isnap]; ni<=lens.nrepperp[isnap]; ni++)
            for(int nj = -lens.nrepperp[isnap]; nj<=lens.nrepperp[isnap]; nj++){

              double di = sqrt(pow(xx[i][0][l]+ni-0.5,2) + pow(xx[i][1][l]+nj-0.5,2) + pow(xx[i][2][l],2));
              double rai,deci,dd;
              getPolar(xx[i][0][l]+ni-0.5,xx[i][1][l]+nj-0.5,xx[i][2][l],rai,deci,dd, true);
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
      }
      totPartxyi[i]=xs.size();
      ntotxyi[i]+=totPartxyi[i];

      if(totPartxyi[i]>0){
        mapxyi[i] = gridist_w(xs,ys,ms,p.npix,DO_NGP);
      }

      //re-normalize to the total mass!
      double mtot=0.;
      double mnorm=0.;
      for(int ii=0;ii<ms.size();ii++){
        mnorm+=ms[ii];
      }

      if(totPartxyi[i]>0){
        for(int l=0;l<p.npix*p.npix;l++)
          mtot += mapxyi[i][l];
        if(mtot==0.)
          mtot=1.; //To avoid NaN
        for(int l=0;l<p.npix*p.npix;l++)
          mapxyi[i][l]*=mnorm/mtot;
      }
    }
  }

  return 0;

}

/*
 * Creates the density maps mapping the particles in the planes mapxy.
 * See the description of the mapParticles routine
 */
int createDensityMaps (InputParams &p, Lens &lens, Random &random, int isnap,
  int ffmin, int ffmax, string File, double fovradiants, double rcase,
  gsl_spline *GetDl, gsl_interp_accel *accGetDl, gsl_spline *GetZl,
  gsl_interp_accel *accGetZl, valarray<float> &mapxytot,
  valarray<float> (&mapxytoti)[6], int (&ntotxyi)[6], int myid){

  mapxytot.resize( p.npix*p.npix );
  for(int i=0;i<6;i++){
    ntotxyi[i]=0;
    mapxytoti[i].resize( p.npix*p.npix );
  }
  for (unsigned int ff=ffmin; ff<ffmax; ff++){

    Header data;
    string file_in = File+"."+sconv(ff,fINT);
    ifstream fin;
    if (readHeader (file_in, data, fin, false))
      return 1;
    if(ff==0)
      printHeader(data);
    /* Creating the pointers for the Data structures */
    Gadget *gadget;
    float *xx[6][3];

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

    readPos (fin,  data, p, random, isnap, xx, rcase, myid);

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

    if(mapParticles(fin, data, p, lens, xx, fovradiants, isnap, mapxyi,
                                                 ntotxyi, myid)) return 1;


    if(p.hydro){
      fin.clear();
      fin.close();
    }

    mapxytot+=mapxyi[0]+mapxyi[1]+mapxyi[2]+mapxyi[3]+mapxyi[4]+mapxyi[5];
    for(int i=0; i<6; i++)
      mapxytoti[i]+=mapxyi[i];

    if (myid==0)
      cout << " done map*tot " << endl;

    delete gadget;

  }

  if(myid == 0)
    cout << " maps done! from Rank:" << myid << endl;
  return 0;

}

/*
 * Writes down the density maps.
 * See the description of the mapParticles routine
 */
void writeMaps (InputParams &p, Header &data, Lens &lens, int isnap, double zsim,
                 string snappl, string snpix, valarray<float>& mapxytotrecv,
                 valarray<float>(& mapxytotirecv)[6], int(& ntotxyi)[6], int myid){

  if (myid==0){
    if(p.partinplanes==false){
     /*
      * write image array(s) to FITS files all particles in a FITS file!
      */
      long naxis = 2;
      long naxes[2]={ p.npix,p.npix };
      string fileoutput;
      fileoutput = fileOutput(p, snappl);
      cout << "Saving the maps on: " << fileoutput << endl;
      unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
      vector<long> naxex( 2 );
      naxex[0]=p.npix;
      naxex[1]=p.npix;
      PHDU *phxy=&ffxy->pHDU();
      phxy->write( 1, p.npix*p.npix, mapxytotrecv );
      phxy->addKey ("REDSHIFT",zsim," ");
      phxy->addKey ("PHYSICALSIZE",p.fov," ");
      phxy->addKey ("PIXELUNIT",1.e+10/data.h,"Mass unit in M_Sun");
      phxy->addKey ("DlLOW",lens.ld[isnap]/data.h,"comoving distance in Mpc");
      phxy->addKey ("DlUP",lens.ld2[isnap]/data.h,"comoving distance in Mpc");
      phxy->addKey ("nparttype0",ntotxyi[0]," ");
      phxy->addKey ("nparttype1",ntotxyi[1]," ");
      phxy->addKey ("nparttype2",ntotxyi[2]," ");
      phxy->addKey ("nparttype3",ntotxyi[3]," ");
      phxy->addKey ("nparttype4",ntotxyi[4]," ");
      phxy->addKey ("nparttype5",ntotxyi[5]," ");
      phxy->addKey ("HUBBLE",data.h," ");
      phxy->addKey ("OMEGAMATTER",data.om0," ");
      phxy->addKey ("OMEGALAMBDA",data.oml," ");
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
         fileoutput = fileOutput(p, snappl, i);
         unique_ptr<FITS> ffxy( new FITS( fileoutput, FLOAT_IMG, naxis, naxes ) );
         vector<long> naxex( 2 );
         naxex[0]=p.npix;
         naxex[1]=p.npix;
         PHDU *phxy=&ffxy->pHDU();
         phxy->write( 1, p.npix*p.npix, mapxytotirecv[i] );
         phxy->addKey ("REDSHIFT",zsim," ");
         phxy->addKey ("PHYSICALSIZE",p.fov," ");
         phxy->addKey ("PIXELUNIT",1.e+10/data.h,"Mass unit in M_Sun");
         phxy->addKey ("DlLOW",lens.ld[isnap]/data.h,"comoving distance in Mpc");
         phxy->addKey ("DlUP",lens.ld2[isnap]/data.h,"comoving distance in Mpc");
         phxy->addKey ("nparttype0",ntotxyi[i]," ");
         phxy->addKey ("HUBBLE",data.h," ");
         phxy->addKey ("OMEGAMATTER",data.om0," ");
         phxy->addKey ("OMEGALAMBDA",data.oml," ");
         phxy->addKey ("m"+sconv(i,fINT),data.massarr[i]," ");
       }
     }
    }
  }
}

/*
 * Output file name
 * - label is the file suffix for SubFind files or the particle Type for partinplanes == true
 */
string fileOutput (InputParams p, string snappl, int label){

  if(p.simType == "Gadget" && p.partinplanes == false)
    return p.directory+p.simulation+"."+snappl+".plane_"+p.snpix+"_"+p.suffix+".fits";
  else if(p.simType == "Gadget" && p.partinplanes == true)
    return p.directory+p.simulation+"."+snappl+".ptype"+sconv(label,fINT)+"_plane_"+p.snpix+"_"+p.suffix+".fits";
  else if(p.simType == "SubFind")
    return p.directory+p.simulation+"."+snappl+".plane_"+sconv(p.fov,fDP1)+"_"+p.suffix+".plc."+sconv(label,fINT);
  else{
    throw invalid_argument("Output name format not recognized");
  }

}
