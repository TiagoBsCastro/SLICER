/* Routines to read InputFiles */
#include "data.h"

/**
 * @brief Reads the input parameters file and stores the values in the provided InputParams structure.
 *
 * This function opens and reads the contents of a named input parameters file. It parses the file
 * and populates the provided InputParams structure with the values read from the file. If the file
 * cannot be opened, the program will print an error message and terminate.
 *
 * @param p Reference to an InputParams structure where read parameters will be stored.
 * @param name Name of the file to read. The file should be in the current working directory.
 * @return int Returns 0 on success, or exits with 1 if the file cannot be opened.
 */
int readInput(struct InputParams &p, string name){
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
  std::getline(fin, str);
  std::getline(fin, str);
  p.w = std::stof(str);    //        14. Dark-Energy EOS w


  if(p.npix==0)
    p.simType = "SubFind";
  else
    p.simType = "Gadget";

  p.physical = (p.npix<0);
  if(! p.physical)
    p.snpix=sconv(p.npix,fINT);
  else{
    int n = -p.npix;
    p.snpix=sconv(n,fINT);
    p.snpix+="_kpc";
    p.rgrid=n;
  }

  if(p.snopt<0){
    cerr << "Impossible value for Shot-Noise option!" << endl;
    return 1;
  }

  return 0;

}

/**
 * @brief Constructor for the SubFind object, initializing its internal state and preparing it for use.
 *
 * This constructor sets up the initial configuration of the SubFind object, based on the provided
 * parameters. It ensures that all necessary properties are set to default values if not specified
 * and prepares any internal resources needed for the object's operations.
 *
 * @param param1 Description of what param1 represents and how it's used in initialization.
 * @param param2 Description of what param2 does and its significance to the SubFind object.
 */
SubFind::SubFind(int n, bool halos){
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
  if(halos)
    this->fsub.resize(n);
};
