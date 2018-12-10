#include "densitymaps.h"

int CreateDensityMaps (InputParams *p, Lens *lens, Random *random, int isnap, int ffmin, int ffmax, string File, double fovradiants, double rcase, gsl_spline *GetDl, gsl_interp_accel *accGetDl,
                        gsl_spline *GetZl, gsl_interp_accel *accGetZl, valarray<float> &mapxytot, valarray<float> (&mapxytoti)[6], int (&ntotxyi)[6], int myid){

  mapxytot.resize( p->npix*p->npix );
  for(int i=0;i<6;i++){
    ntotxyi[i]=0;
    mapxytoti[i].resize( p->npix*p->npix );
  }
  for (unsigned int ff=ffmin; ff<ffmax; ff++){

    Header data;
    string file_in = File+"."+sconv(ff,fINT);
    ifstream fin;
    if (read_header (file_in, &data, fin, false))
      return 1;
    if(ff==0)
      print_header(data);
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

    ReadPos (fin,  &data, p, random, isnap, xx, rcase, myid);

    /* If Hydro run I have to read the masses, otherwise close the snapshot*/
    if(p->hydro)
      fastforwardToBlock (fin, "MASS", myid);
    else{
      fin.clear();
      fin.close();
    }

    // map for each mass type
    valarray<float> mapxyi[6];
    int ntotxyi[6];
    for(int i=0; i<6; i++)
      mapxyi[i].resize(p->npix*p->npix);

    if(MapParticles(fin, &data, p, lens, xx, fovradiants, isnap, mapxyi,
                                                 ntotxyi, myid))
      return 1;

    if(p->hydro){
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
