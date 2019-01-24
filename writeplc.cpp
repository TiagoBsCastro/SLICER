#include "writeplc.h"

void writePLC (SubFind &halos, Header &data, InputParams &p, string snappl, int ff){

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
