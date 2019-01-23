#include "util.h"
#define MAX_M 1e3    // Threshold for mass; Particler heavier than MAX_M will be attached zero mass
#define POS_U 1.0    // Unit conversion from BoxSize unit lengh to kpc/h
#define DO_NGP false // Use NGP as the MAS

using namespace CCfits;

int readInput(struct InputParams &p, string name, string & snpix, bool &physical){
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
  p.npix = std::stoi(str);//         1. Number of Pixels
  std::getline(fin, str);
  std::getline(fin, str);
  p.zs = std::stof(str);//           2. Redshift Source
  std::getline(fin, str);
  std::getline(fin, str);
  p.fov = std::stof(str);//          3. Field of View
  std::getline(fin, str);
  std::getline(fin,
      p.filredshiftlist);//          4. File with the redshift list it may contain three columns: snap 1/(1+z) z
  std::getline(fin, str);
  std::getline(fin,
      p.pathsnap);//                 5. Path where the snaphosts are located
  std::getline(fin, str);
  std::getline(fin,
      p.simulation);//               6. Simulation name
  std::getline(fin, str);
  std::getline(fin, str);
  p.seedcenter = std::stoi(str);//   7. Seed center
  std::getline(fin, str);
  std::getline(fin, str);
  p.seedface = std::stoi(str);//     8. Seed Face
  std::getline(fin, str);
  std::getline(fin, str);
  p.seedsign = std::stoi(str);//     9. Seed sign
  std::getline(fin, str);
  std::getline(fin, str);
  p.partinplanes = std::stoi(str);// 10. True: Each gadget particle type will have its own Map; False: One Map for All
  std::getline(fin, str);
  std::getline(fin,
       p.directory);//               11. Directory to save FITS files
  std::getline(fin, str);
  std::getline(fin,
       p.suffix);//                  12. Suffix to FITS files
  std::getline(fin, str);
  std::getline(fin, str);
  p.snopt = std::stoi(str);//        13. Shot-noise option: 0-No random Degradation; 1-Half particles degradation ...

  if(p.npix==0)
    p.simType = "SubFind";
  else
    p.simType = "Gadget";

  physical = (p.npix<0);
  if(! physical)
    snpix=sconv(p.npix,fINT);
  else{
    int n = -p.npix;
    snpix=sconv(n,fINT);
    snpix+="_kpc";
    p.rgrid=n;
  }

  if(p.snopt<0){
    cerr << "Impossible value for Shot-Noise option!" << endl;
    return 1;
  }

  return 0;

}

int read_redlist(string filredshiftlist, vector <double> & snapred, vector <string> & snappath, InputParams &p){

  ifstream redlist;
  redlist.open(filredshiftlist.c_str());
  double zlast=-999.9;

  if(redlist.is_open()){
    ifstream fin; // snap file
    string name;  // snapshot directory
    double z;     // snapshot redshift
    Header header;// snapshot header
    do{
      redlist >> name;
      snappath.push_back(name);
      if( read_header( p.pathsnap+name+".0" , header, fin, true) ){
        cerr << name << " not found!" << endl;
        return 1;
      }
      if(header.redshift<zlast){
        cerr << " Snapshots on "<< filredshiftlist << " are not sorted!" << endl;
        return 1;
      }
      else
        zlast = header.redshift;
      snapred.push_back(header.redshift);
    }while(header.redshift<p.zs);
  }else{
    cerr << " redshift list file redshift_list.txt does not " << endl;
    cerr << " exist in the Code dir ... check this out      " << endl;
    cerr << "    I will STOP here !!! " << endl;
    return 1;
  }
  return 0;
}

SubFind::SubFind(int n){
  this->id.resize(n);
  this->truez.resize(n);
  this->xx0.resize(n); this->yy0.resize(n); this->zz0.resize(n);
  this->vx0.resize(n); this->vy0.resize(n); this->vz0.resize(n);
  this->m.resize(n);
  this->theta.resize(n);
  this->phi.resize(n);
  this->vel.resize(n);
  this->obsz.resize(n);
  this->nsub.resize(n);
};

void CreatePLC (SubFind &halos, Header &data, InputParams &p, string snappl, int ff){

  fstream fileoutput ( p.directory+p.simulation+"."+snappl+".plane_"+sconv(p.fov,fDP1)+"_"+p.suffix+".plc."+sconv(ff,fINT), ios::out | ios::binary );
  double fovradiants = p.fov/180.*M_PI;

  int nhalos = halos.m.size();

  for(int i=0; i<nhalos; i++){

    if( (abs(halos.phi[i]) <=  fovradiants/2.0) & (abs(halos.theta[i]) <=  fovradiants/2.0) & (halos.m[i]>0) ){
      /* Explicitly casting the variables as in Pinocchio PLC */
      int dummy;
      long long unsigned int id = halos.id[i];
      double xx0 = (halos.xx0[i] - 0.5)*data.boxsize/1e3;
      double yy0 = (halos.yy0[i] - 0.5)*data.boxsize/1e3;
      double zz0 = halos.zz0[i]*data.boxsize/1e3;
      double vx0 = halos.vx0[i];
      double vy0 = halos.vy0[i];
      double vz0 = halos.vz0[i];
      double m = halos.m[i];
      double obsz = halos.obsz[i];
      double truez = halos.truez[i];
      double vel = halos.vel[i];
      double theta, phi, r;

      getPolar(xx0, yy0, zz0, theta, phi, r, false);
      theta = 90.0-theta*180.0/M_PI;
      phi = phi*180.0/M_PI;
      if( phi<0.0 )
        phi += 360.0;

      fileoutput.write((char*)&dummy, sizeof (int));
      fileoutput.write((char*)&id, sizeof (unsigned long long int));
      fileoutput.write((char*)&truez, sizeof (double));
      fileoutput.write((char*)&xx0, sizeof (double));
      fileoutput.write((char*)&yy0, sizeof (double));
      fileoutput.write((char*)&zz0, sizeof (double));
      fileoutput.write((char*)&vx0, sizeof (double));
      fileoutput.write((char*)&vy0, sizeof (double));
      fileoutput.write((char*)&vz0, sizeof (double));
      fileoutput.write((char*)&m, sizeof (double));
      fileoutput.write((char*)&theta, sizeof (double));
      fileoutput.write((char*)&phi, sizeof (double));
      fileoutput.write((char*)&vel, sizeof (double));
      fileoutput.write((char*)&obsz, sizeof (double));
      fileoutput.write((char*)&dummy, sizeof (int));

    }

  }

  fileoutput.close();

}
