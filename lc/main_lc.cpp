#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define PY_SSIZE_T_CLEAN
#include <cmath>
#include <chrono> 
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <functional>
#include <stdio.h>     
#include <stdlib.h>     
#include <ctime>
#include <bits/stdc++.h> 
#include <typeinfo>
#include <cassert>
#include <dirent.h>
//#include <readSUBFIND.h>
#include "/usr/include/hdf5/serial/H5Cpp.h"
#include </usr/include/eigen3/Eigen/Dense>
#include "readTNGParticle.h"
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include "functions.h"

/*****************************************************************************/
/*                                                                           */
/*             FORECAST - lightcone construction module                      */
/*                                                                           */
/*  original dark matter-only code by cgiocoli@gmail.com                     */
/*  updated to its final form by flaminia.fortuni@inaf.it                    */
/*  if you use it or do any mods, please cite Fortuni et al. (2023).         */
/*                                                                           */
/*                                                                           */  
/*  this is the first module of the FORECAST code;                           */
/*  it builds the structure of the lightcone and arrange particles in fov.   */
/*  it runs on one snapshot per time.                                        */
/*  - input: snapshot files; groups catalogs (if needed for subhalo IDs)     */
/*  - output: particle in fov catalog needed for next module                 */
/*                                                                           */
/*                                                                           */
/*  for a comprehensive guide, visit                                         */
/*                          https://github.com/flaminiafortuni/FORECAST      */
/*  for a full description of the software                                   */     
/*       https://ui.adsabs.harvard.edu/abs/2023arXiv230519166F/abstract      */
/*                                                                           */
/*                                                                           */           
/*****************************************************************************/

using namespace std;
using namespace arma;
using namespace std::chrono; 


const int bleft = 24;
const double speedcunit = 2.99792458e+3;
const double speedcunitas = 2.9979e+18; //AA/s
const double ergsa = 3.839e+33;
const double mpM = 8.4089382e-58; //proton mass in Msun
const double mpg = 1.6726219e-24; //proton mass in g

const double h0 = 0.6774;



int main(int argc, char** argv){
  
  auto start = high_resolution_clock::now();
  cout << "----------------------------------------------------------------------" << endl;
  cout << " " << endl; 
  cout << "   ------------------------------------------------------ " << endl;
  cout << "   -                                                    - " << endl;
  cout << "   -           2D Mapping Simulation Snapshot           - " << endl;
  cout << "   -                                                    - " << endl;
  cout << "   -               building the light-cone              - " << endl;
  cout << "   ------------------------------------------------------ " << endl;

  // check if the file restart exsits ... 
  std:: string fileplstart = std::string(argv[2])+".d";
  int iplrestart=0;
  std:: ifstream infileplstart;
  infileplstart.open(fileplstart.c_str());
  if(infileplstart.is_open()){
    std::cout << " " << std:: endl;
    std:: cout << " I will read the restart file >> " << fileplstart << std:: endl;
    infileplstart >> iplrestart;
    std:: cout << " iplrestart = " << iplrestart << std:: endl;
    std:: cout << " " << std:: endl;
    infileplstart.close();
  }

  
  // ******************** to be read in the INPUT file ********************
  // ... simulation box_length in [Mpc/h], highest redshift & distance reached by the cone, dimension of fov in degrees, resolution of the image in arcsecs
  double boxl,zs,Ds,fov,res;
  // ... files
  string filredshiftlist,filsnaplist, filtimelist, filfilters,idc;
  // ... paths
  string pathsnap, rdir;
  // ... seeds for randomization of the simulation box
  long seedcenter, seedface, seedsign;
  // ... sim
  string sim;

  readParameters(&boxl,&zs,&fov,&res,
		 &filredshiftlist,&filsnaplist,&filtimelist,&idc,
		 &pathsnap,&rdir,
		 &seedcenter,&seedface,&seedsign,&sim);

  
  //pixels of the image
  int truenpix=int(fov*3600/res);
  int bufferpix=int(ceil((truenpix + 1)*20 / 14142));  // add bufferpix/2 in each side!
  int npix=truenpix+bufferpix;

  cout << "N. pixels: " << truenpix << "; buffer pixels: " << bufferpix << endl;
  cout << endl;
  
  
  // ... read the redshift list and the snap_available in redshift_list
  ifstream redlist;
  redlist.open(filredshiftlist.c_str());
  vector <int> tsnaplist;
  vector <double> dtsnaplist;
  vector <double> tredlist;
  int nmax = 1024;
  vector<double> snapToRedshift(nmax); // no way to have more than 1024 snaphosts
  for(int i=0;i<nmax;i++){
    snapToRedshift[i] = -1;
  }
  if(redlist.is_open()){
    int buta;
    double butb,butc;
    while(redlist >> buta >> butb >> butc){
      tsnaplist.push_back(buta);
      tredlist.push_back(butc);
      dtsnaplist.push_back(buta);
      if(buta>nmax){
	cout << " check nmax variable and increase it! " << endl;
	cout << " for now I will STOP here!!! " << endl;
	exit(1);
      }
      snapToRedshift[buta] = butc;
    }
    redlist.close();
  }else{   
    cout << " redshift list file " << filredshiftlist << " does not " << endl;
    cout << " exist in the Code dir ... check this out      " << endl;
    cout << "    I will STOP here !!! " << endl;
    exit(1);
  }

  
  // open and read timelist 
  ifstream timelist;
  timelist.open(filtimelist.c_str());
  vector <int> snap_tlist;
  vector <double> z_tlist;
  vector <double> age_tlist;
  if(timelist.is_open()){
    int buta;
    double butb,butc;
    while(timelist >> buta >> butb >> butc){
      snap_tlist.push_back(buta);
      z_tlist.push_back(butb);
      age_tlist.push_back(butc);
    }
    timelist.close();
  }
  		 
  else{
    cout << " time list file " << filtimelist << "  does not " << endl;
    cout << " exist in the Code dir ... check this out      " << endl;
    cout << "    I will STOP here !!! " << endl;
    exit(1);
  }

  // open and read snaplist
  ifstream snaplist;
  snaplist.open(filsnaplist.c_str());
  vector<int> vsnaps;
  vector <int> lsnap;
  vector<double> lred;
  vector<double> lD;
  if(snaplist.is_open()){
    int s;
    double zn;
    int check=0;
    while(snaplist >> s){
      zn = getY(dtsnaplist,tredlist,double(s));
      vsnaps.push_back(s);
      if(zn>zs) check=1;
      if(check==0){
	lsnap.push_back(s);
	lred.push_back(zn);
      }
    }
    snaplist.close();
  }else{
    cout << filsnaplist << " does not exist in the code dir " << endl;
    cout << "   I will STOP here!!! " << endl;
    exit(1);
  }
  
  int nsnaps = lsnap.size();

  cout << "  " << endl;
  cout << " Opening path for snapshots >>  " << pathsnap << endl;
  cout << " " << endl;
  cout << " I will look for comoving distance file >>  " << idc << endl;
  cout << " " << endl;
  
  //open and read LCDM-comovingdistance - this file is created with astropy 
  ifstream infiledc;
  vector<double> zl, dl, dlum;  
  infiledc.open(idc.c_str());
  if(infiledc.is_open()){
    double zi,dci,dli;
    while(infiledc >> zi >> dci >> dli){
      zl.push_back(zi);
      dl.push_back(dci*speedcunit); //comoving Mpc/h. dl is a comoving distance, not a luminosity distance!
      dlum.push_back(dli*speedcunit);//comoving Mpc/h. Luminosity distance.
    }
    infiledc.close();
  }
  else{
    cout << "  " << endl;
    cout << " the comoving distance file: " << idc << endl;
    cout << " does not exists " << endl;
    cout << " I will STOP here!!! " << endl;
    exit(1);
  }
  
  if(zs>zl[zl.size()-1]){
    cout << " source redshift larger than the highest available redshift in the comoving distance file " << endl;
    cout << "  that is = " << zl[zl.size()-1] << endl;
    cout << " I will STOP here !!! " << endl;
    exit(1);
  }
  
  Ds = getY(zl,dl,zs);  // comoving distance of the last plane;
  int nreplications = int(Ds/boxl)+1;
  
  cout << " nreplications = " << nreplications << ";   Ds (cMpc/h) = " << Ds << "  " << std:: endl;
  vector<int> replication;
  vector<int> fromsnap;
  vector<double> lD2;

  
  for(int j=0;j<nreplications;j++){
    for(int i=0;i<nsnaps;i++){
      double ldbut = getY(zl,dl,lred[i]);
      if(ldbut>=j*boxl && ldbut<=(j+1.)*boxl){
	replication.push_back(j);
	fromsnap.push_back(lsnap[i]);
	lD.push_back(ldbut); //Mpc/h
      }
    }
  }
  cout << " " << endl;
  cout << "  " << endl;
  cout << " ... reorganazing the planes ... " << std:: endl;
  cout << " " << endl;
  

  // creating the lightcone
  for(int i=0;i<lD.size();i++){
    if(i<(lD.size()-1)) lD2.push_back(lD[i+1]);
    else lD2.push_back(Ds);
  }
  for(int i=0;i<lD.size();i++){
    for(int k=1;k<=512;k++){
      if(lD[i]<double(k)*boxl && lD2[i]>double(k)*boxl){
	lD[i+1] = lD2[i]; 
	lD2[i]=double(k)*boxl;
      }
    }
    if(lD[i]<513*boxl && lD2[i]>513*boxl){
      cout << " exiting ... increase the number of replications by hand in the file it is now 512 !!!! " << endl;
      exit(1);
    }
  } 
  
  std:: cout << "  " << endl;
  vector<double> zsimlens(nsnaps);
  cout << " N_snaphot (including replications) = " << nsnaps << endl;
  cout << " " << endl;
  

  for(int i=0;i<nsnaps;i++){
    if(i<nsnaps-1){
      if(lD[i+1]-lD2[i]>boxl*1e-9){
	fromsnap[i] = -fromsnap[i];
      }
    }
    double dlbut = (lD[i] + lD2[i])*0.5;    
    // half distance between the two!
    zsimlens[i] = getY(dl,zl,dlbut);
  }
  
  vector<double> bfromsnap,blD,blD2,bzsimlens,blred;
  vector<int> breplication, blsnap;

  vector<double> bbfromsnap,bblD,bblD2,bbzsimlens,bblred;
  vector<int> bbreplication, bblsnap;

  vector<double> bbbfromsnap,bbblD,bbblD2,bbbzsimlens,bbblred;
  vector<int> bbbreplication, bbblsnap;
  
  vector<double> Bfromsnap,BlD,BlD2,Bzsimlens,Blred;
  vector<int> Breplication, Blsnap;

  int pl=0;
  vector<int> pll;
  
  for(int i=0;i<nsnaps;i++){
    Bfromsnap.push_back(fabs(fromsnap[i]));
    BlD.push_back(lD[i]);
    BlD2.push_back(lD2[i]);
    Bzsimlens.push_back(zsimlens[i]);
    Breplication.push_back(replication[i]);
    Blsnap.push_back(lsnap[i]);
    Blred.push_back(lred[i]);
    if(fromsnap[i]<0){
      Bfromsnap.push_back(-fromsnap[i]);
      BlD.push_back(lD2[i]); 
      BlD2.push_back(lD[i+1]);
      double dlbut = (lD[i+1] + lD2[i])*0.5;    
      // half distance between the two!
      Bzsimlens.push_back(getY(dl,zl,dlbut)); 
      Breplication.push_back(replication[i+1]);
      Blsnap.push_back(lsnap[i]);
      Blred.push_back(lred[i]);
    }   
  }
  cout << "  " << endl;
  cout << " ... re-reorganazing the planes ..." << std:: endl;
  cout << " " << endl;
  cout << " nsnaps (including replications) = " << nsnaps << endl;
  cout << " " << endl;

   for (auto i = 0; i < Blred.size(); i++) {
    float delta = BlD2[i] - BlD[i];
    float delta1=BlD2[i+1]-BlD[i+1];
    int n = static_cast<int>(std::ceil(delta / boxl));
    if (delta > boxl) {
        double add = delta / static_cast<double>(n);

        for (int k = 1; k <= n; k++) {
            bblD.push_back(BlD[i] + (k - 1) * add);
            bblD2.push_back(BlD[i] + k * add);
            double dlbut =BlD[i] + (2 * k - 1) * 0.5 * add;// bblD.back() + (k - 0.5) * add;
            bbzsimlens.push_back(getY(dl, zl, dlbut));
            bbreplication.push_back(Breplication[i - 1] + k);
            bblsnap.push_back(Blsnap[i]);
            bblred.push_back(Blred[i]);
            bbfromsnap.push_back(Bfromsnap[i]);
        }
    } else {
        bblD.push_back(BlD[i]);
        bblD2.push_back(BlD[i] + delta);
        double dlbut = BlD[i] + 0.5 * delta;  
        bbzsimlens.push_back(getY(dl, zl, dlbut));
        bbreplication.push_back(Breplication[i - 1] + 1);
        bblsnap.push_back(Blsnap[i]);
        bblred.push_back(Blred[i]);
        bbfromsnap.push_back(Bfromsnap[i]);
    }
    
    if((delta+delta1)/2.<=boxl && Bfromsnap[i]==Bfromsnap[i+1] ){ //set equal delta
      //firts half
      bblD.push_back(BlD[i]);
      bblD2.push_back((BlD[i]+BlD2[i+1])/2.);
      double dlbut=(BlD[i]+(BlD2[i+1]+BlD[i])/2.)/2.;
      bbzsimlens.push_back(getY(dl,zl,dlbut));
      bbreplication.push_back(Breplication[i]);
      bblsnap.push_back(Blsnap[i]);
      bblred.push_back(Blred[i]);
      bbfromsnap.push_back(Bfromsnap[i]);
      //second half
      bblD.push_back((BlD2[i+1]+BlD[i])/2.);
      bblD2.push_back(BlD2[i+1]);
      double dlbut2=(BlD2[i+1]+(BlD2[i+1]+BlD[i])/2.)/2.;
      bbzsimlens.push_back(getY(dl,zl,dlbut2));
      bbreplication.push_back(Breplication[i]+1);
      bblsnap.push_back(Blsnap[i]);
      bblred.push_back(Blred[i]);
      bbfromsnap.push_back(Bfromsnap[i]);
      i=i+1;
    }
    
    
  }
  
  // merge planes if they are too small
  for (auto i=0;i<bblred.size();i++){ 
    float delta=bblD2[i]-bblD[i];
    float delta1=bblD2[i+1]-bblD[i+1];
    if(delta+delta1<=boxl && bbfromsnap[i]==bbfromsnap[i+1] ){
      bbbfromsnap.push_back(bbfromsnap[i]);
      bbblD.push_back(bblD[i]);
      bbblD2.push_back(bblD2[i+1]);
      double dlbut=(bblD[i]+bblD2[i+1])/2.;
      bbbzsimlens.push_back(getY(dl,zl,dlbut));
      bbbreplication.push_back(bbreplication[i]+1);
      bbblsnap.push_back(bblsnap[i]);
      bbblred.push_back(bblred[i]);
      i=i+1;
    }
    else{
      bbbfromsnap.push_back(bbfromsnap[i]);
      bbblD.push_back(bblD[i]);
      bbblD2.push_back(bblD2[i]);
      bbbzsimlens.push_back(bbzsimlens[i]);
      bbbreplication.push_back(bbreplication[i]);
      bbblsnap.push_back(bblsnap[i]);
      bbblred.push_back(bblred[i]);
    }
  }
  //again
  for (auto i=0;i<bbblred.size();i++){ 
    float delta=bbblD2[i]-bbblD[i];
    float delta1=bbblD2[i+1]-bbblD[i+1];
    if(delta+delta1<=boxl && bbbfromsnap[i]==bbbfromsnap[i+1] ){
      bfromsnap.push_back(bbbfromsnap[i]);
      blD.push_back(bbblD[i]);
      blD2.push_back(bbblD2[i+1]);
      double dlbut=(bbblD[i]+bbblD2[i+1])/2.;
      bzsimlens.push_back(getY(dl,zl,dlbut));
      breplication.push_back(bbbreplication[i]+1);
      blsnap.push_back(bbblsnap[i]);
      blred.push_back(bbblred[i]);
      i = i + 1;
    }
    else{
      bfromsnap.push_back(bbbfromsnap[i]);
      blD.push_back(bbblD[i]);
      blD2.push_back(bbblD2[i]);
      bzsimlens.push_back(bbbzsimlens[i]);
      breplication.push_back(bbbreplication[i]);
      blsnap.push_back(bbblsnap[i]);
      blred.push_back(bbblred[i]);
    }
  }

  //erase repeated rows
  vector<int> er_inx(0);
  for (auto i = 0; i < blD.size() - 1; i++) {
    if (blD[i] == blD[i + 1] && blD2[i] != blD2[i + 1]) {
      er_inx.push_back(i);
    }
    else if (blD[i] == blD[i + 1] && blD2[i] != blD2[i + 1]) { //INTERROTTO QUI
      er_inx.push_back(i);
    }
  }

  for (auto it = er_inx.rbegin(); it != er_inx.rend(); ++it) {
    bfromsnap.erase(bfromsnap.begin() + *it);
    blD.erase(blD.begin() + *it);
    blD2.erase(blD2.begin() + *it);
    bzsimlens.erase(bzsimlens.begin() + *it);
    breplication.erase(breplication.begin() + *it);
    blsnap.erase(blsnap.begin() + *it);
    blred.erase(blred.begin() + *it);
  }
    
  Bfromsnap.clear();
  Bfromsnap.shrink_to_fit();
  BlD.clear();
  BlD.shrink_to_fit();
  BlD2.clear();
  BlD2.shrink_to_fit();
  Bzsimlens.clear();
  Bzsimlens.shrink_to_fit();
  Breplication.clear();
  Breplication.shrink_to_fit();
  Blsnap.clear();
  Blsnap.shrink_to_fit();
  Blred.clear();
  Blred.shrink_to_fit();
  
  bbfromsnap.clear();
  bbfromsnap.shrink_to_fit();
  bblD.clear();
  bblD.shrink_to_fit();
  bblD2.clear();
  bblD2.shrink_to_fit();
  bbzsimlens.clear();
  bbzsimlens.shrink_to_fit();
  bbreplication.clear();
  bbreplication.shrink_to_fit();
  bblsnap.clear();
  bblsnap.shrink_to_fit();
  bblred.clear();
  bblred.shrink_to_fit();
  breplication.clear();
  breplication.shrink_to_fit();
  

  // creating planes_list.txt file: structure of the lighcone
  ofstream planelist;
  planelist.open("planes_list.txt"); 
  nsnaps = bfromsnap.size();  
  for(int i=0;i<nsnaps;i++){
    breplication.push_back(i);      
    pl++;
    planelist <<  pl << "   " <<  bzsimlens[i] << "   " << blD[i] << "   " << blD2[i] << "   " <<  breplication[i] << "   " << bfromsnap[i] << "   " << blred[i] << std:: endl;
    pll.push_back(pl);    
  }

  
  
  if(blD2[nsnaps-1]<Ds-1e-6){
    // we need to add one more snaphost
    snaplist.open(filsnaplist.c_str());
    double s;
    while(snaplist >> s){
      if(s<bfromsnap[nsnaps-1]){
	bfromsnap.push_back(s);
	double dlbut = (Ds+blD2[nsnaps-1])*0.5;
	double zbut = getY(dl,zl,dlbut);
	blD.push_back(blD2[nsnaps-1]);
	blD2.push_back(Ds);
	bzsimlens.push_back(zbut);
	breplication.push_back(breplication[nsnaps-1]+1);
	double zn = getY(dtsnaplist,tredlist,double(s));
	blsnap.push_back(s);
	lred.push_back(zn);
	nsnaps++;
	snaplist.close();	
      }
    }
    pl++;
    planelist << pl << "   " << bzsimlens[nsnaps-1] << "   " << blD[nsnaps-1] << "   " << blD2[nsnaps-1] << "   " << breplication[nsnaps-1] << "   " << bfromsnap[nsnaps-1] << "   " << blred[nsnaps-1] << std:: endl;
    pll.push_back(pl);	 
  }
  
  planelist.close();

  std:: cout << " " << std:: endl;
  std:: cout << " ... re-reorganazing the planes ..." << std:: endl;
  std:: cout << " " << std:: endl;
  std:: cout << " nsnaps (including replications) at the end = " << nsnaps << std::endl;
  std:: cout << "  " << endl;

  
  // randomizzation of the box realizations :
  int nrandom = breplication[nsnaps-1]+1;
  vector<double> x0(nrandom), y0(nrandom), z0(nrandom); // ramdomizing the center of the simulation [0,1]
  vector<int> face(nrandom); // face of the dice
  vector<int> sgnX(nrandom), sgnY(nrandom),sgnZ(nrandom);
  
  for(int i=0;i<nrandom;i++){
    if(seedcenter>0){
      srand(seedcenter+i*13);
      x0[i] = rand() / float(RAND_MAX);
      y0[i] = rand() / float(RAND_MAX);
      z0[i] = rand() / float(RAND_MAX);
    }else{
      x0[i] = 0.;
      y0[i] = 0.;
      z0[i] = 0.;
    }
    face[i] = 7;
    if(seedface>0){
      srand(seedface+i*5);
      while(face[i]>6 || face[i]<1) face[i] = int(1+rand() / float(RAND_MAX)*5.+0.5);
    }else{
      face[i]=1;
    }
    sgnX[i] = 2;
    if(seedsign>0){
      srand(seedsign+i*8);
      while(sgnX[i] > 1 || sgnX[i] < 0) sgnX[i] = int(rand() / float(RAND_MAX)+0.5);
      sgnY[i] = 2;
      while(sgnY[i] > 1 || sgnY[i] < 0) sgnY[i] = int(rand() / float(RAND_MAX)+0.5);
      sgnZ[i] = 2;
      while(sgnZ[i] > 1 || sgnZ[i] < 0) sgnZ[i] = int(rand() / float(RAND_MAX)+0.5);
      if(sgnX[i]==0) sgnX[i]=-1;
      if(sgnY[i]==0) sgnY[i]=-1;
      if(sgnZ[i]==0) sgnZ[i]=-1;
    }else{
      sgnX[i]=1;
      sgnY[i]=1;
      sgnZ[i]=1;
    }
  }
  
  double truefov = fov;
  std:: cout << "  " << endl;
  cout << "  " << endl;
  cout << " set the field of view to be square in degrees. " << endl;
  cout << " fov's side  is " << fov << " in degrees and " << fov*3600. << " in arcsec. " << endl;
  double fovradiants;
  double om0, omL0;
  fovradiants = fov/180.*M_PI;
  double truefovradiants = fovradiants;

  // check of the field of view is too large with respect to the box size
  std:: cout << " [ maximum fov side's value allowed " << boxl/Ds*180./M_PI << " in degrees] " << std:: endl;
  if(fovradiants*Ds>boxl){
    std:: cout << " field view too large ... I will STOP here!!! " << std:: endl;
    std:: cout << " value set is = " << fov << std:: endl;
    std:: cout << " maximum value allowed " << boxl/Ds*180./M_PI << " in degrees " << std:: endl;
    std:: cout << " Check it out ... I will STOP here!!! " << endl;

    exit(1);
  }
  
  //for selection of the particles
  fov = truefov*double(npix)/double(truenpix);
  fovradiants = fov/180.*M_PI;
  cout << " " << endl;
  cout << " FORECAST uses " << nsnaps << " snapshots (including replications) for this image simulations" << endl;
  
  int nsnap=iplrestart;
  
  if(blD2[nsnap]-blD[nsnap] < 0){
    cout << " comoving distance of the starting point " << blD[nsnap] << endl;
    cout << " comoving distance of the final    point " << blD2[nsnap] << endl;
    cout << " please check this out! I will STOP here!!! " << endl;
    exit(1);
  }

  int rcase = breplication[nsnap];
  // get current snapshot number
  string snappl=conv(blsnap[nsnap],fINT);       
  int sourceID=blsnap[nsnap];
 
  cout << "" << endl;
  cout << "----------------------------------------------------------------------" << endl;
  cout << " " << endl;
  cout << "... Starting to read TNG for snapshot " << snappl << "..." << endl;
  cout << " " << endl;
 
  // ... masses of the different type of particles
  double m0,m1,m2,m3,m4,m5;
  // ... star particles
  vector<int> ID4_(0), ID4(0);                                                            
  vector<float> X4_(0), Y4_(0), Z4_(0), IM4_(0), M4_(0), Met4_(0), FTIME4_(0), ZF4_(0);             
  vector<float> X4(0), Y4(0), Z4(0), IM4(0), M4(0), Met4(0), ZF4(0), FTIME4(0),AGE4(0);           
  vector<float> xx4(0), yy4(0), zz4(0), zzs(0), zzred4(0), org_z(0);
  // ... 2Dgrid stars
  vector<int> ids(0), idSH4(0), idshs(0);
  vector<float> xs(0),ys(0),zstar(0),zreds(0),ms4(0),ims4(0),ages4(0),dtmp(0),orgZ(0);
  vector<long double> mets4(0);
  vector<double> fiub(0);    
  // ... other
  vector<int> gcLenTypeS(0), gcOffsetsTypeS(0);
  vector<long double> met_bc03 { 0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05};
  double zsim, dlsim;
  double bs,time;  
  // ... constants
  float convm=(1.e10/h0);
 
  // ... readind simulation 
  // initialize the reading particles class
  readTNGParticle tngParticle;
  tngParticle.Initialize(pathsnap, sourceID); 
  
  // read header
  tngParticle.readHeader(0);
  // get header info
  bs = tngParticle.getBoxSize(); //[ckpc/h]
  om0 = tngParticle.getOmegaZero();
  omL0 = tngParticle.getOmegaLambda();
  time = tngParticle.getTime();
  zsim = tngParticle.getRedshift();
  dlsim = getY(zl,dl,zsim); //[Mpc/h];
  std::vector<double> mass = tngParticle.getMassTable(); //[1e10 Msun/h]
  std::vector<int> npart = tngParticle.getNumPartTotal();
  float bsh0=bs*1./h0; //bs (boxsize) is in kpc/h
  float h0bs=h0*1./bs;
  
  // reserve space for star particles
  X4_.reserve(npart[4]);
  Y4_.reserve(npart[4]);
  Z4_.reserve(npart[4]);
  IM4_.reserve(npart[4]);
  M4_.reserve(npart[4]);
  Met4_.reserve(npart[4]);
  FTIME4_.reserve(npart[4]);
  ZF4_.reserve(npart[4]);
  ID4_.reserve(npart[4]);
    
  // Count files in snapdir and groups directories
  std::string workdir_sn = pathsnap + "/snapdir_0"+snappl+"/"; // snapshot path; each snapshot stored in its own folder
  std::string workdir_g = pathsnap + "/groups_0"+snappl+"/"; // groups path; each group in its folder. These files need to be separated from snapshot and offset files
  std::string workdir_os = pathsnap + "/offsets_0"+snappl+".hdf5"; // offset path - all together
  int nf_sn=countHDF5Files(workdir_sn, "hdf5");
  int nf_g=countHDF5Files(workdir_g, "hdf5");
  
  // read stellar particles (PartType4) looping through nf_sn *.hdf5 files
  for (int i=0; i<nf_sn; i++){
    tngParticle.readStars(i);
    std::vector<double> metal = tngParticle.getMetallicity();
    std::vector<double> inMass = tngParticle.getInitialMass();
    std::vector<double> sTime = tngParticle.getStellarFormationTime();
    std::vector<double> Mass = tngParticle.getMasses();
    std::vector<double> X_ = tngParticle.getX();
    std::vector<double> Y_ = tngParticle.getY();
    std::vector<double> Z_ = tngParticle.getZ();	
    Met4_.insert(Met4_.end(), metal.begin(), metal.end());
    FTIME4_.insert(FTIME4_.end(), sTime.begin(), sTime.end());
    IM4_.insert(IM4_.end(), inMass.begin(), inMass.end());
    M4_.insert(M4_.end(), Mass.begin(), Mass.end());
    X4_.insert(X4_.end(), X_.begin(), X_.end());
    Y4_.insert(Y4_.end(), Y_.begin(), Y_.end());
    Z4_.insert(Z4_.end(), Z_.begin(), Z_.end());
  }
   
  // matching ids between stars and subhalos - this is specific for TNG simulation
  for (int i=0; i<nf_g; i++){
    tngParticle.readFof(i);
    std::vector<int> sLTS = tngParticle.getStarsLenType();
    gcLenTypeS.insert(gcLenTypeS.end(), sLTS.begin(), sLTS.end());
  }  
  tngParticle.readOffset();
  gcOffsetsTypeS = tngParticle.getStarsByType();
  
  // stars: filtering out spurious particles, computing age etc - this is specific for TNG simulation
  int datalen= Met4_.size();
  for (int i=0; i<datalen; i++){
    ID4_.push_back(i);
  }
  
  for (int i=0; i<datalen; i++){
    if (FTIME4_[i]>0.){
      Met4.push_back(Met4_[i]);
      FTIME4.push_back(FTIME4_[i]);
      IM4.push_back(IM4_[i]*convm); //[Msun/h]
      M4.push_back(M4_[i]*convm); //[Msun/h]
      X4.push_back(X4_[i]); //[kpc/h]
      Y4.push_back(Y4_[i]);
      Z4.push_back(Z4_[i]);
      ZF4.push_back(1./FTIME4_[i]-1.);
      ID4.push_back(ID4_[i]);
    }
  }
  idSH4 = inverseMap_sh_idv3(gcLenTypeS,gcOffsetsTypeS,ID4, true);
      
  cout << "TNG FULL snapshot - N particles: " <<ID4.size() << endl;
           
  Met4_.clear();
  Met4_.shrink_to_fit();
  IM4_.clear();
  IM4_.shrink_to_fit();
  M4_.clear();
  M4_.shrink_to_fit();
  X4_.clear();
  X4_.shrink_to_fit();
  Y4_.clear();
  Y4_.shrink_to_fit();
  Z4_.clear();
  Z4_.shrink_to_fit();
  FTIME4_.clear();
  FTIME4_.shrink_to_fit();
  ID4_.clear();
  ID4_.shrink_to_fit();

  // compute age of the particle, knowing the formation redshift of the particle - this is specific for TNG
  int finaldatalen=Y4.size();      
  for (int i=0; i<finaldatalen; i++){
    float tf=getY(z_tlist,age_tlist,ZF4[i]); // this interpolation holds true because z- Dc the spacing is thick
    float t0=getY(z_tlist, age_tlist,zsim );
    float tmp=0.;
    tmp=t0-tf;
    AGE4.push_back(tmp);
  }
  
     
  
  // print simulation header info
  cout << "  " << endl;
  cout << "      __________________ COSMOLOGY __________________  " << endl;
  cout << " " << endl;
  cout << "      Omegam = " << om0 << "      " << "Omegal = " << omL0 << endl;
  cout << "           h = " << h0   << "      " << "BoxSize (cMpc/h)= " << bs/(1.e+3) << endl;
  cout << "      redshift = " << zsim <<   "   " << "Dl (cMpc/h) = " << dlsim << endl;
  if(abs(boxl - bs/1.e+3)>1.e-2 ){
    cout << " set boxl and data.size differ ... check it! " << std:: endl;
    cout << "  boxl = " << boxl << "  " << " data.boxsize = " << bs/1.e+3 << endl;
    exit(1);
  }
  
  cout << "      _______________________________________________  " << endl;
  cout << " " << endl;
  cout << "   gas (0); dm (1); tracers (3); stars (4); bh (5)   //  (2) unused. " << endl;   
  cout << "   total number of particles in the simulation: " << endl; 
  cout << npart[0] << " " << npart[1]  << " " << npart[3] << " " << npart[4]
       << " " << npart[5] <<  endl;
  cout << " " << endl;
  cout << "   particle type mass array: " << endl; 
  cout << mass[0] << " " << mass[1] << " " << mass[3] 
       << " " << mass[4] << " " <<  mass[5]  << endl;
  cout << " " << endl;
  
  m0 = mass[0]; //gas, empty
  m1 = mass[1]; //dm
  m3 = mass[3]; //tracers
  m4 = mass[4]; //stars, empty
  m5 = mass[5]; //bh, empty
  
  cout << " " << endl;
  cout << npart[4] << "   type (4)   - STAR particles firstly selected in the snapshot."<<endl;
  cout << "" <<endl;
  

  // manipulating the coordinates in the box
  for (int pp=0; pp<X4.size(); pp++) {
    float x, y, z, zred, zred_int, dlz;
    float xb, yb, zb, zr, orgz;
    
    xb = sgnX[rcase]*(((X4[pp])/bs));
    yb = sgnY[rcase]*(((Y4[pp])/bs));
    zb = sgnZ[rcase]*(((Z4[pp])/bs));
    zr = sgnZ[rcase]*(((Z4[pp])/bs));
    orgz=Z4[pp];
    // wrapping periodic condition 
    if(xb>1.) xb = xb - 1.;
    if(yb>1.) yb = yb - 1.;
    if(zb>1.) zb = zb - 1.;
    if(zr>1.) zr = zr - 1.;
    if(xb<0.) xb = 1. + xb;
    if(yb<0.) yb = 1. + yb;
    if(zb<0.) zb = 1. + zb;
    if(zr<0.) zr = 1. + zr;
    switch (face[rcase]){
    case(1):
      x = xb;
      y = yb;
      z = zb;
      zred = zr;
      break;
    case(2):
      x = xb;
      y = zb;
      z = yb;
      zred = yb;
      break;
    case(3):
      x = yb;
      y = zb;
      z = xb;
      zred = xb;
      break;
    case(4):
      x = yb;
      y = xb;
      z = zb;
      zred = zr;
      break;
    case(5):
      x = zb;
      y = xb;
      z = yb;
      zred = yb;
      break;
    case(6):
      x = zb;
      y = yb;
      z = xb;
      zred = xb;
      break;
    }
    // recenter
    x = x - x0[rcase];
    y = y - y0[rcase];
    z = z - z0[rcase];
    zred = zred - z0[rcase];
    // wrapping periodic condition again
    if(x>1.) x = x - 1.;
    if(y>1.) y = y - 1.;
    if(z>1.) z = z - 1.;
    if(zred>1.) zred = zred - 1.;
    if(x<0.) x = 1. + x;
    if(y<0.) y = 1. + y;
    if(z<0.) z = 1. + z;
    if(zred<0.) zred = 1. + zred;
    zzs.push_back(z); // z in [0,1] normalized to bs 
    z+=double(rcase); // pile the cones double(blred[nsnap])
    xx4.push_back(x);
    yy4.push_back(y);
    zz4.push_back(z); // N replica
    org_z.push_back(orgz);
  }
  
  X4.clear();
  X4.shrink_to_fit();
  Y4.clear();
  Y4.shrink_to_fit();
  Z4.clear();
  Z4.shrink_to_fit();
  
  
  
  int n4 = xx4.size();
  int totPartxy4;
  
  cout << n4 <<"   type (4) - STAR particles re-arranged in the snapshot."<<endl;
  cout << "" << endl;
  
  if(n4>0){
    // print some infos of the quadrate box
    double xmin=double(*min_element(xx4.begin(), xx4.end()));
    double xmax=double(*max_element(xx4.begin(), xx4.end()));  
    double ymin=double(*min_element(yy4.begin(), yy4.end()));
    double ymax=double(*max_element(yy4.begin(), yy4.end()));  
    double zmin=double(*min_element(zz4.begin(), zz4.end()));
    double zmax=double(*max_element(zz4.begin(), zz4.end()));
    double zstarmin=double(*min_element(zzs.begin(), zzs.end()));
    double zstarmax=double(*max_element(zzs.begin(), zzs.end()));
    cout << " " << endl;
    cout << " n4 particles " << endl;
    cout << "xmin = " << xmin << endl;
    cout << "xmax = " << xmax << endl;
    cout << "ymin = " << ymin << endl;
    cout << "ymax = " << ymax << endl;
    cout << "zmin = " << zstarmin << endl;
    cout << "zmax = " << zstarmax << endl;
    cout << "boxmin = " << zmin << endl;
    cout << "boxmax = " << zmax << endl;
    cout << "  " << endl;
    if(xmin<0 || ymin<0 || zmin< 0){
      cout << "xmin = " << xmin << endl;
      cout << "xmax = " << xmax << endl;
      cout << "ymin = " << ymin << endl;
      cout << "ymax = " << ymax << endl;
      cout << "boxmin = " << zmin << endl;
      cout << "boxmax = " << zmax << endl;
      cout << "  4 type check this!!! I will STOP here!!! " << endl;
      exit(1);
    }
    cout << " ... selecting only particles in fov ..." << endl;
    cout << "" << endl;
    
    
    // select particles in the field of view
    vector<double> dtmp(0), ra(0), decl(0), xtmp(0), ytmp(0), ztmp(0), redtmp(0), infv(0);
    for(int l=0;l<n4;l++){
      
      double Ztmp= (zzs[l]*(bs/1.e+3)+blD[nsnap])/(bs/1.e+3); // [Mpc/h]
      double di = sqrt(pow(xx4[l]-0.5,2) + pow(yy4[l]-0.5,2) + pow(Ztmp,2))*bs/1.e+3; // distance between particle and observer at (.5,.5,0)*bs
      
      if(di>=blD[nsnap] && di<blD2[nsnap]){	    
	double rai,deci,dd;
	getPolar(xx4[l]-0.5,yy4[l]-0.5,Ztmp,&rai,&deci,&dd);
	if(fabs(rai)<=fovradiants*0.5 && fabs(deci)<=fovradiants*0.5){	  
	  double fovinunitbox = fovradiants*di/(bs/1.e+3);		
	  ids.push_back(ID4[l]);
	  fiub.push_back(fovinunitbox);
	  xs.push_back((xx4[l]-0.5)/fovinunitbox+0.5);
	  ys.push_back((yy4[l]-0.5)/fovinunitbox+0.5);
	  zstar.push_back(zzs[l]/fovinunitbox);
	  zreds.push_back(getY(dl,zl,Ztmp*(bs/1.e+3)));
	  ms4.push_back((M4[l]));
	  ims4.push_back(IM4[l]);
	  mets4.push_back(Met4[l]);
	  ages4.push_back(AGE4[l]);
	  dtmp.push_back(di);		
	  xtmp.push_back(xx4[l]*(bs/1.e+3)); //Mpc/h [0,75]
	  ytmp.push_back(yy4[l]*(bs/1.e+3)); //Mpc/h [0,75]
	  orgZ.push_back(org_z[l]);
	}
      }	    
    }
    
    Met4.clear();
    Met4.shrink_to_fit();
    IM4.clear();
    IM4.shrink_to_fit();
    M4.clear();
    M4.shrink_to_fit();
    AGE4.clear();
    AGE4.shrink_to_fit();
    FTIME4.clear();
    FTIME4.shrink_to_fit();
    ZF4.clear();
    ZF4.shrink_to_fit();
    ID4.clear();
    ID4.shrink_to_fit();

    cout << " " << endl;
    cout << " ... Coordinates & redshift limits ..." << endl;
    cout << "xmin, xmax in Mpc/h " << double(*min_element(xtmp.begin(), xtmp.end())) << ", " << double(*max_element(xtmp.begin(), xtmp.end())) << endl;
    cout << "ymin, ymax in Mpc/h " << double(*min_element(ytmp.begin(), ytmp.end())) << ", " << double(*max_element(ytmp.begin(), ytmp.end())) << endl;
    cout << "zmin, zmax in Mpc/h " << double(*min_element(zreds.begin(), zreds.end())) << ", " << double(*max_element(zreds.begin(), zreds.end())) << endl;
    cout << " " << endl;
    
    // find SHid of only the particles included in the field of view
    idshs =inverseMap_sh_idv3(gcLenTypeS,gcOffsetsTypeS,ids, true);
    
    //filtering out idshs==-1, which are fuzzes
    vector<int>idshtmp;
    for(auto i=0;i<idshs.size();i++){
      if (idshs[i]>-1)
	idshtmp.push_back(idshs[i]);
    }
    vector<int> idsh;
    std::unique_copy(idshtmp.begin(), idshtmp.end(), std::back_inserter(idsh));
    idshtmp.clear();
    idshtmp.shrink_to_fit();
	
    // CMz for stars
    vector<double> cmzsh(0);
    for (auto k=0;k<idsh.size();k++){
      double mzs=0.;
      double mtot=0.;
      for(auto l=0;l<idshs.size();l++){
	if (idshs[l]==idsh[k]){
	  mzs= mzs + (orgZ[l]*ms4[l]);
	  mtot+=ms4[l];
	}
      }
      if (mtot != 0.0) { 
	cmzsh.push_back(1.0 * mzs / mtot);
      } else {
	cmzsh.push_back(0.0);
      }
    }	  
   
    // counting number of particles in each SH from TNG cat
    std::vector<int> NpTNG(0);
    for (int i = 0; i < idsh.size(); i++) {
      int countshTNG=0;
      for (auto l=0;l<idSH4.size();l++){
	if (idSH4[l]==idsh[i]){
	  countshTNG++;
	}
      }
      NpTNG.push_back(countshTNG);
    }

    vector<double> CMzsh(0), NpshTNG;
    for (auto k=0;k<idsh.size();k++){
      for(auto l=0;l<idshs.size();l++){
	if (idshs[l]==idsh[k]){
	  CMzsh.push_back(cmzsh[k]);
	  NpshTNG.push_back(NpTNG[k]);
	}
      }
    }

    idSH4.clear();
    idSH4.shrink_to_fit();
    NpTNG.clear();
    cmzsh.clear();
    NpTNG.shrink_to_fit();
    cmzsh.shrink_to_fit();
    
    
    totPartxy4=xs.size();
    cout << totPartxy4 <<"   type (4) - STAR particles selected in the FoV."<<endl;
    
    
    if (totPartxy4>0)
      cout <<  "first subhalo ID: " << idshs[0] << "  |  last subhalo ID: " << idshs[totPartxy4-1] << endl;
    cout << "" << endl;
    cout << " ... Writing the output file (coords.lc.snap_plane.txt)... " << endl;
    cout << " .. content: #(ids,x,y,redshift,mass,intial mass,metallicity, -, age, -, CM of sh with ids, Nparticles of sh with ids in sim) for stellar particles;" << endl;
    cout << endl;
    
    // write particles in fov on a txt file:
    string coord_path_="coords.lc."+snappl+"_"+conv(iplrestart,fINT)+".txt";
    ofstream myfile2_;
    myfile2_.open(coord_path_);
    for (int i=0;i<totPartxy4; i++){
      if(idshs[i]>-1)
	myfile2_ << idshs[i] << " " << xs[i] << " " << ys[i] << " " << zreds[i] << " " << ms4[i] << " " << ims4[i] << " " << mets4[i] << " " << fiub[i] << " " << ages4[i] << " " << zstar[i]  << " " << CMzsh[i] << " " << NpshTNG[i] <<  endl;
    }
    myfile2_.close();
    
    
	cout << ".. Coord file written .." << endl;
        cout << "" << endl;
	
  }
  
       
  cout << " First step is concluded. " << endl;
  cout << " ... I'm trying to free your mind, Neo ... " << endl;
  cout << endl;
  //clock
  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start);
  cout << "execution time in [h] " <<  2.77778e-10*duration.count() << endl;
  std::cout << "----------------------------------------------------------------------" << endl; 
  exit(1);
  
}
