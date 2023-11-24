#include <cmath> 
#include <chrono> 
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <iterator>
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <stdio.h>     
#include <stdlib.h>     
#include <ctime>
#include <bits/stdc++.h> 
#include <typeinfo>
#include <cassert>
#include "/usr/include/hdf5/serial/H5Cpp.h"
#include </usr/include/eigen3/Eigen/Dense>
#include <iostream>
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include <typeinfo>

#include "readTNGParticle.h"
#include "SED.h"
#include "functions.h"


/*****************************************************************************/
/*                                                                           */
/*             FORECAST - dust-corrected fluxes calculations module          */
/*                                                                           */
/*  original dark matter-only code by cgiocoli@gmail.com                     */
/*  updated to its final form by flaminia.fortuni@inaf.it                    */
/*  if you use it or do any mods, please cite Fortuni et al. (2023).         */
/*                                                                           */
/*                                                                           */ 
/*  this is the third module of the FORECAST code;                           */
/*  it computes dust-corrected fluxes, including k-correction (no IGm).      */
/*  it runs on one snapshot per time.                                        */
/*  - input: dust-free particle catalog from previous module df              */
/*  - output: dust-corrected particle catalog in N filters                   */ 
/*            needed for the next module                                     */
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
using namespace std::chrono;
using namespace arma;

const double speedcunit = 2.99792458e+3;
const double mpM = 8.4089382e-58; //proton mass in Msun
const double mpg = 1.6726219e-24; //proton mass in g
const double mpc2tocm2 = 9.523e+42; //Mpc^2 to cm^2; used to convert NHI (division)
const double h0 = 0.6774;
//4th grade fit; these are from snap 33 and xH=0.76 in T and 1000000 points
const float a1=2.31105169e-02;
const float a2=-8.91893982e-01;
const float a3=1.14885315e+01;
const float a4=-6.27073890e+01;
const float a5= 1.18956375e+02;

int main (int argc, char** argv){
  
  auto start = high_resolution_clock::now();

  cout << "----------------------------------------------" << endl;
  cout << "-                                            -" << endl;
  cout << "::::DUST POST-PROCESSING for the LIGHTCONE::::" << endl;
  cout << "-                                            -" << endl;
  cout << "----------------------------------------------" << endl;


  // checking if the file restart exsits ... 
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

  // ******************** from INPUT file ********************
  
  // ... box_length of the simulation in [Mpc/h]
  double boxl;
  float maglim;
  // ... files
  string filfilters,filsnaplist, filtimelist;
  // ... directories
  string pathsnap, bc03dir, rdir, idc;
  //... user can choose to read an existing dust catalog (read==0) or to compute dust from scratch from the simulation (read==1)
  int read;
  // ... SSP model & imf
  string imf, model;
  string pixunits="uJy";
  
  // reading parameters from MapSimDust.ini
  readParameters(&boxl,&maglim,
		 &filfilters,&filsnaplist,&filtimelist,&idc,
		 &pathsnap,&bc03dir,&rdir,
		 &read,
		 &imf,&model);
  
  // reading files
  int sourceID=std::atoi(argv[1]);
  int n_pl;
  float blD, blD2, zsim;
  string snappl;
  string fplane= "../lc/planes_list.txt";
  ifstream oplane;  
  oplane.open(fplane.c_str());
  if(oplane.is_open()){
    int npl,reppl,blsnappl;
    float zpl, blDpl, blD2pl,zpltrue;
    
    while(oplane >> npl >> zpl >> blDpl >> blD2pl >> reppl >> blsnappl >> zpltrue){
      if(npl-1==iplrestart && blsnappl==sourceID){
	n_pl=npl;
	blD=blDpl;
	blD2=blD2pl;
	snappl=conv(blsnappl,fINT);
	zsim=zpltrue;
      }
      
    }
    oplane.close();
  }
  else{
    cout << "  " << endl;
    cout << " planes list  " << fplane << endl;
    cout << " does not exists for this snapshot. " << endl;
    cout << " I will STOP here!!! " << endl;
    exit(1);
  }
  
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
    cout << " time list file " << filtimelist << " does not " << endl;
    cout << " exist in the Code dir ... check this out      " << endl;
    cout << "    I will STOP here !!! " << endl;
    exit(1);
  }

  //open age_bc03.txt or age_cb16.txt
  ifstream agelist;
  string filagelist="../files/age_"+model+".txt";
  agelist.open(filagelist.c_str());
  vector <double> nage;
  vector <long double> age_bc03;
  if(agelist.is_open()){
    double buta,butb;
    while(agelist >> buta >> butb){
      nage.push_back(buta);
      age_bc03.push_back(butb);
    }
    agelist.close();
  }
  
  else{
    cout << " time list file " << filagelist << " does not " << endl;
    cout << " exist in the Code dir ... check this out      " << endl;
    cout << "    I will STOP here !!! " << endl;
    exit(1);
  }
  
  //open extinction curve coefficients and compute A(lambda)/Av=a+b/Rv; lambda in [um] - CCM+89
  ifstream extcurve;
  string filextc=rdir+"extinctionCCM89.dat";
  extcurve.open(filextc.c_str());
  vector <long double> ex_l;
  vector <long double> ex_curve;
  double Rv=3.1;  
  if(extcurve.is_open()){
    std::string  str1;    
    for (int i=0; i<1; i++)      
      std::getline(extcurve, str1); // skip 1 lines    
    float buta,butb,butc,butd,bute;    
    while(extcurve >> buta >> butb >> butc >> butd >> bute){      
      ex_l.push_back(buta*1.e4); //[um] to [AA]      
      ex_curve.push_back(butb+butc/Rv);     
    }
    extcurve.close();    
  }  
  else{
    cout << " Extinction curve file " << filextc << " does not " << endl;    
    cout << " exist in the Code dir ... check this out      " << endl;    
    cout << "    I will STOP here !!! " << endl;    
    exit(1);
  }

  // read Cosmology file (consistent with simulation's Cosmology)
  ifstream infiledc;
  vector<double> zl, dl, dlum;
  infiledc.open(idc.c_str());
  if(infiledc.is_open()){
    double zi,dci,dli;
    while(infiledc >> zi >> dci >> dli){
      zl.push_back(zi);
      dl.push_back(dci*speedcunit);   //comoving distance [cMpc/h]
      dlum.push_back(dli*speedcunit); //luminosity distance [cMpc/h]
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

  // filter list in filters.dat
  ifstream listfilter;
  listfilter.open(filfilters.c_str());
  vector <string> nfilter;
  if(listfilter.is_open()){
    string buta;
    while(listfilter >> buta){
      nfilter.push_back(buta);
    }
    listfilter.close();
  }  		 
  else{
    cout << " filter list file does not " << endl;					     
    cout << " exist in the Code dir ... check this out      " << endl;
    cout << "    I will STOP here !!! " << endl;
    exit(1);
  }

  // reading filters  
  int Nfilters= nfilter.size();
  std::vector <vector<long double> > fresp(Nfilters,vector<long double>(0));
  std::vector <vector<long double> > fwaves(Nfilters,vector<long double>(0));
  
  for(auto N=0;N<Nfilters;N++){
    string filterin="../files/filters/" + nfilter[N] + ".dat";
    ifstream filterlist;
    filterlist.open(filterin.c_str());
    if(filterlist.is_open()){
      long double  sa, sb;	
      while((filterlist >> sa >> sb )){
	fwaves[N].push_back(sa);
	fresp[N].push_back(sb);
      }
      filterlist.close();
    }
    else{
      cout << "This filter " << nfilter[N] << " is missing in filters folder." << endl;
      cout << " Check it out ... I will STOP here!!! " << endl;
      exit(1);
    }
  }
    
  //SSP EXTRACTION
  std::vector<long double> waves;//[AA]
  vector<vector<double> > time_grid;
  vector<vector<vector<double> > > full_table;
  if(model=="bc03"){
    time_grid.resize(6);
    full_table.resize(6);
    for (int i = 0; i < 6; ++i) {
      readSSPTables(bc03dir, model, imf, i, waves, time_grid[i], full_table[i]);
    }
    
  }
  else if(model=="cb16"){
    time_grid.resize(14);
    full_table.resize(14);
    for (int i = 0; i < 14; ++i) {
      readSSPTables(bc03dir, model, imf, i, waves, time_grid[i], full_table[i]);
    }
  }

 
  // variable, vector and constant declarations
  double bs, om0, omL0, dlsim;
  //gas cell                                                                            
  vector<int> ID0(0), idshgas(0), idsh0(0);
  //star particles
  vector<int> idSH4(0),ID4_(0), ID4(0);
  vector<float> Z4_(0), IM4_(0), M4_(0), Met4_(0), FTIME4_(0);             
  vector<float> Z4(0), IM4(0), M4(0), Met4(0), ZF4(0), FTIME4(0),AGE4(0);
  // 2Dgrid gas
  vector<int> idg(0),idshg(0);
  // 2Dgrid stars
  vector<int> ids(0), idshs(0), idshsmap(0), idsmap(0), NpshTNG(0);
  vector<float> ims4(0), zredsmap(0);
  vector<float> ages4(0);
  vector<long double> mets4(0);
  vector <double> NHIm(0), Zm(0);
  vector<float> xmap(0),ymap(0), ms4map(0), fluxgaldusty(0), deltaflux(0), z4map(0), ims4map(0), ages4map(0), mets4map(0), cmzshmap(0);
  vector<vector<double>> fluxmap;
  //other
  vector<long double> met_bc03 { 0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05};
  vector<long double> met_cb16 { 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, 0.010, 0.014, 0.017, 0.020, 0.030, 0.040 };
  SED sedy; 
  float convrho = (1.e10/h0)/(pow(1./h0,3.));
  float convm= 1.e10/h0;
  float convZsun=1./0.0127;
  double invmpM=1./mpM;   //[Msun-1]
  double kB=1.38e-16;     //[erg/K]
  double gamma=5./3.;     //adiabatic index
  double xH=0.76;         // H fraction
  double NHI_0=2.e64;     //[kpc-2]
  double invNHI_0=1./NHI_0;

  cout << " " << endl;
  cout << " " << endl;
  cout << "... Reading dust-free particle-based catalog produced in the previous module (flux.df.snap_plane.txt) ..." << endl;
  cout << " .. content: #(ids,x,y,redshift,mass,intial mass,metallicity, CM_sh, age, N_sim) for stellar particles;" << endl;
  cout << " ..          #(dust-free flux in Nfilters) for stellar particles." << endl;
  cout << endl;

  // reading df output catalogue (dust-free fluxes)
  string filoutcat="../df/flux.df."+snappl+"_"+conv(iplrestart,fINT)+".txt";
  ifstream ocf;
  ocf.open(filoutcat.c_str());
  if(ocf.is_open()){
    int shid, NPTNG;
    float xi,yi,zi,zri,mi,imi,ai,fi, sni, cmzsh;
    long double meti, mabi;
    while (ocf >> shid >> xi >> yi  >> zri >> mi >> imi >> meti >> cmzsh >> ai >> NPTNG) { 
    idshsmap.push_back(shid);
    xmap.push_back(xi);
    ymap.push_back(yi);
    zredsmap.push_back(zri);
    ms4map.push_back(mi);
    ims4map.push_back(imi);
    mets4map.push_back(meti);
    ages4map.push_back(ai);
    cmzshmap.push_back(cmzsh);        
    NpshTNG.push_back(NPTNG);
    std::vector<double> temp_flux;
    double temp_val;
    int cf = 0;
    while (cf < Nfilters) {
        if (ocf >> temp_val) {
            temp_flux.push_back(temp_val);
            cf += 1;
        } else {
            break;
        }
    }
    fluxmap.push_back(temp_flux);
    ocf.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }    
    ocf.close();
  }
  else{
    cout << "  " << endl;
    cout << " OUT_CAT file: " << filoutcat << endl;
    cout << " does not exists for this snapshot. " << endl;
    cout << " I will STOP here!!! " << endl;
    exit(1);
  }
  
  int totPartMapxy4=idshsmap.size();
  cout << idshsmap.size() << "   type (4)   - STAR particles in dust-free catalogue " << endl;

  // printing some parameters for check    
  cout << " " << endl; 
  cout << "... Some parameters ..." << endl;
  cout << " (D_l, D_l2) = (" << blD << ", " << blD2 << ") cMpc; z between (" << (getY(dl,zl,blD)) << ", " << (getY(dl,zl,blD2))  << "); snapshot " << sourceID << ", plane " << n_pl <<  endl;
  cout << " filters: " << endl;
  for (auto i=0;i<nfilter.size();i++){
	cout << nfilter[i] << "  -  " ;
  }
 
  cout << "" << endl;
  cout << "....................... " << endl;
  cout << " " << endl;
  cout << " " << endl;
  
  cout << "... Creating SH cat from input dust-free star particles cat ..." << endl;
  // subhalo catalogs
  std::vector<int> idsh, Npsh;
  std::vector<float> zrsh, CMzsh;
  std::vector<std::vector<float>> fluxsh;   
  std::unique_copy(idshsmap.begin(), idshsmap.end(), std::back_inserter(idsh));  
  for (int i = 0; i < idsh.size(); i++) {
    int countsh = 0;
    float redsh = 0., cmzsh=0.;
    std::vector<float> fsh(fluxmap[0].size(), 0.0);
    for (auto l = 0; l < idshsmap.size(); l++) {
        if (idshsmap[l] == idsh[i]) {	 
	  countsh++;
	  redsh += zredsmap[l];
	  cmzsh=cmzshmap[l];	  
	  for (size_t col = 0; col < fluxmap[l].size(); col++) {
	    fsh[col] += fluxmap[l][col];
	  }
        }
    }
    CMzsh.push_back(cmzsh);
    Npsh.push_back(countsh);
    zrsh.push_back(redsh * 1.0 / countsh);
    fluxsh.push_back(fsh);
  }


  // mass of galaxies in input catalog
  vector<float>  ms4sh(0);
  for (auto k=0;k<idsh.size();k++){
    float tmas=0.;
    for(int l=0; l<idshsmap.size(); l++){ 
      if (idsh[k]==idshsmap[l]){
	tmas+=ms4map[l];
      }
    }
    ms4sh.push_back(tmas);
  }

  cout << " " << endl;
  cout << " " << endl;

  cout << "... Starting to read TNG cat ..." << endl;
  cout << " " << endl;
 
  if (read==1){
  
    // TNG particle class
    readTNGParticle tngParticle;
    tngParticle.Initialize(pathsnap, sourceID);
    tngParticle.readHeader(0); //chunk file id:0   
    bs = tngParticle.getBoxSize(); //[ckpc/h]
    om0 = tngParticle.getOmegaZero();
    omL0 = tngParticle.getOmegaLambda();
    zsim = tngParticle.getRedshift();
    dlsim = getY(zl,dl,zsim); //[Mpc/h]
    std::vector<int> npart = tngParticle.getNumPartTotal();
    cout << endl;
    
    // count files in snapdir and groups directories
    std::string workdir_sn = pathsnap + "/snapdir_0"+conv(sourceID,fINT)+"/"; //snapshot path
    std::string workdir_g = pathsnap + "/groups_0"+conv(sourceID,fINT)+"/"; //groups path. These files need to be separated from snapshot and offset files
    std::string workdir_os = pathsnap + "/offsets_0"+conv(sourceID,fINT)+".hdf5"; //offset path
    int nf_sn=countHDF5Files(workdir_sn, "hdf5");
    int nf_g=countHDF5Files(workdir_g, "hdf5");
    
    cout << " .. Extracting gas from TNG files .." << endl;

    // reading subhalos
    vector<int> gcLenTypeG(0), gcOffsetsTypeG(0);
    for (int i=0; i<nf_g; i++){ 
      tngParticle.readFof(i);
      std::vector<int> sLTG = tngParticle.getGasLenType();
      gcLenTypeG.insert(gcLenTypeG.end(), sLTG.begin(), sLTG.end());
    }  
    tngParticle.readOffset();
    gcOffsetsTypeG = tngParticle.getGasByType();

    vector<float> zgas(0), massgas(0),rhogas(0),  metgas(0), ugas(0), xegas(0);
    xegas.reserve(npart[0]);
    ugas.reserve(npart[0]);
    massgas.reserve(npart[0]);
    metgas.reserve(npart[0]);
    rhogas.reserve(npart[0]);
    zgas.reserve(npart[0]);
    idshgas.reserve(npart[0]);
    vector<double>  metal0(0), Mass0(0), Rho_0(0), Z_0(0), u_0(0), x_e0(0);
    int minxm = *min_element(idshsmap.begin(), idshsmap.end());
    int maxxm = *max_element(idshsmap.begin(), idshsmap.end());
    std::vector<int> Indexes;
    int current_index = 0;

    // reading simulation files
    for (auto i = 0; i < nf_sn; i++) {
      tngParticle.readGAS(i);
      metal0 = tngParticle.getGASMetallicity();
      Mass0 = tngParticle.getGASMasses(); // [1e10 Msun/h]
      Rho_0 = tngParticle.getGASDensity(); // [(1e10 Msun/h)/(kpc/h)^3]
      Z_0 = tngParticle.getGASZ(); 	
      u_0 = tngParticle.getGASIntEnergy(); // internal energy u [(km/s)^2]; to convert it in internal temperature, check IllustrisTNG doc
      x_e0 = tngParticle.getGASeAbundance(); // abundace of electron in H n_e= x_e0 * n_H
      
      int metal_length = metal0.size();
      
      std::vector<int> current_indexes(metal_length);
      std::iota(current_indexes.begin(), current_indexes.end(), current_index);
      
      idsh0 = inverseMap_sh_idv3(gcLenTypeG,gcOffsetsTypeG,current_indexes,true);
      
      for(auto k=0;k<idsh0.size();k++){
	if(idsh0[k]>=minxm and idsh0[k]<=maxxm){
	  idshgas.push_back(idsh0[k]);
	  metgas.push_back(metal0[k]);
	  massgas.push_back(Mass0[k]);
	  rhogas.push_back(Rho_0[k]);
	  zgas.push_back(Z_0[k]);
	  ugas.push_back(u_0[k]);
	  xegas.push_back(x_e0[k]);  
	}	
      }
      current_index += metal_length;
      
      metal0.clear();
      metal0.shrink_to_fit();
      idsh0.clear();
      idsh0.shrink_to_fit();
      x_e0.clear();
      u_0.clear();
      Z_0.clear();
      Rho_0.clear();
      Mass0.clear();    
      x_e0.resize(0);
      u_0.resize(0);
      Z_0.resize(0);
      Rho_0.resize(0);
      Mass0.resize(0);      
    }

 
    cout << idshgas.size()     << "   type (0)   - GAS  cells     firstly selected in the snapshot." << endl;
    cout << "" << endl;
        
    //mag of galaxies in input catalog
    std::vector<std::vector<float>> magsh; 
    for (int j = 0; j < fluxsh.size(); j++) {
      std::vector<float> magsh_j;    
      for (int i = 0; i < fluxsh[0].size(); i++) {
	magsh_j.push_back(2.5 * (29 - log10(fluxsh[j][i])) - 48.6);
      }    
      magsh.push_back(magsh_j);
    }
    
    cout << ".. Picking gas cells .." << endl;
    //pick only gas cells with same idSH as star particles
    vector<float> CMzsh_mlim(0);
    vector<float> Z0(0), M0(0), Rho0(0), Met0(0), U0(0), eab(0);
    vector<int> idSH0(0);
    eab.reserve(idshgas.size());
    U0.reserve(idshgas.size());
    M0.reserve(idshgas.size());
    Z0.reserve(idshgas.size());
    Rho0.reserve(idshgas.size());
    Met0.reserve(idshgas.size());
    idSH0.reserve(idshgas.size());
    int gasc=0;
    for(auto i=0;i<idshgas.size();i++){
      for(auto j=0;j<idsh.size();j++){
	if (idshgas[i]==idsh[j]){
	  idSH0.push_back(idshgas[i]);
	  eab.push_back(xegas[i]);
	  U0.push_back(ugas[i]);
	  Z0.push_back(zgas[i]);
	  Rho0.push_back(rhogas[i]);
	  M0.push_back(massgas[i]);
	  Met0.push_back(metgas[i]);
	  if (magsh[0][j] <= maglim){ //MAG LIM on filter n=0
	    CMzsh_mlim.push_back(CMzsh[j]);
	  }
	  else{
	    CMzsh_mlim.push_back(-9); 
	  }
	  gasc++;
	}
      }
    }
  
    CMzsh.clear();
    CMzsh.resize(0);
    magsh.clear();
    magsh.resize(0);
    
    eab.resize(gasc);
    U0.resize(gasc);
    M0.resize(gasc);
    Met0.resize(gasc);
    Rho0.resize(gasc);
    Z0.resize(gasc);
    idSH0.resize(gasc);
    
    xegas.clear();
    ugas.clear();
    zgas.clear();
    rhogas.clear();
    massgas.clear();
    metgas.clear();
    idshgas.clear();
    ID0.clear();
    ID0.resize(0);
    xegas.resize(0);
    ugas.resize(0);
    zgas.resize(0);
    rhogas.resize(0);
    massgas.resize(0);
    metgas.resize(0);
    idshgas.resize(0);
    
    // selecting gas that redden star light: only for cells with z> zsh_CM
    cout << "... Selecting TNG gas cell reddening the flux..." << endl;
    vector<long double> metg0(0), NHIgInt(0) , Mratiog(0);
    vector<float> mg0(0);
    vector<int> idshg(0);
    
    idshg.reserve(idSH0.size());
    mg0.reserve(idSH0.size());
    metg0.reserve(idSH0.size());
    NHIgInt.reserve(idSH0.size());
    Mratiog.reserve(idSH0.size());
    int validCount=0;
 
    for(auto l=0;l<idSH0.size();l++){
      if (Z0[l]>=CMzsh_mlim[l]){
	idshg.push_back(idSH0[l]);
	mg0.push_back(M0[l]*convm); //Msun
	double V0;      
	if (Rho0[l]!=0)
	  V0 = ((M0[l]*convm)/(Rho0[l]*convrho));
	else if (Rho0[l]==0)
	  V0=0.;
	double L0 = pow(V0,1./3.); //[ckpc]
	double mu0= (4.*mpg)/(1+3.*0.76+4.*0.76*eab[l]); //[g] and fixed xH=0.76
	double T0=(gamma-1)*(U0[l]/kB)*1.e10*mu0; //[K]
	double logT0=log10(T0);
	double logMratio=a1*pow(logT0,4.) + a2*pow(logT0,3.) + a3*pow(logT0,2.) + a4*logT0 + a5;
	double Mratio=pow(10.,logMratio); 
	double nHIgInt=Mratio*Rho0[l]*convrho*invmpM; //[kpc-3]
	metg0.push_back(Met0[l]);
	NHIgInt.push_back(nHIgInt*L0); //[kpc-2]
	Mratiog.push_back(Mratio);
	validCount++;
      }
      else if(((Z0[l]<CMzsh_mlim[l]) || CMzsh_mlim[l]==-9 )){
	idshg.push_back(idSH0[l]);
	mg0.push_back(M0[l]*convm); //Msun
	metg0.push_back(Met0[l]);
	NHIgInt.push_back(0.); //[kpc-2]
	Mratiog.push_back(0.);
	validCount++;
      }    
    }
    
    Met0.clear();
    Met0.resize(0);
    M0.clear();
    M0.resize(0);
    Rho0.clear();
    Rho0.resize(0);
    Z0.clear();
    Z0.resize(0);
    U0.clear();
    U0.resize(0);
    eab.clear();
    eab.resize(0);
    CMzsh_mlim.clear();
    CMzsh_mlim.resize(0);
    idSH0.clear();
    idSH0.resize(0);

    
    auto stop2 = high_resolution_clock::now(); 
    auto duration2 = duration_cast<microseconds>(stop2 - start);
    cout << "execution time in [h] " <<  2.77778e-10*duration2.count() << endl;

    cout << " " << endl;
    cout << totPartMapxy4 << "   type (4)  - STAR particles in dust-free catalogue."<<endl;
    cout << " " << endl;    
    int totPartxy0=idshg.size();
    cout << totPartxy0 <<"        type (0) - GAS cells with same SH ids."<<endl; 
    cout << "" << endl; 
    cout << "... Selecting TNG gas cells in front of MapSim star particles ..." << endl;
    
    //geometric selection of gas cells in front of star particles - deriving dust extinction from gas properties - model C - (not manipulated coords)
    NHIm.reserve(idsh.size());
    Zm.reserve(idsh.size());
          
    //compute NHImean and ZHImean for each subhaloID in idsh[l]
    for (auto l=0;l<idsh.size();l++){ //galaxy-based
      double NHIsum=0.;
      double ZHIsum=0.;
      double MHIsum=0.;
      for (int k=0; k<totPartxy0; k++){
	if (idshg[k]==idsh[l]){
	  NHIsum=NHIsum+(mg0[k]*Mratiog[k]*NHIgInt[k]);
	  ZHIsum=ZHIsum+(mg0[k]*Mratiog[k]*metg0[k]);
	  MHIsum=MHIsum+(mg0[k]*Mratiog[k]);
	}
      }
      //mean NHI and Zg (weighted with MHI mass) for each subhalo with idsh ID
      if (MHIsum!=0.){
	NHIm.push_back((1.*NHIsum / MHIsum) * invNHI_0); //idsh dimension
	Zm.push_back((1.*ZHIsum / MHIsum) * convZsun);
      }
      else{ //no gas, no NHI, no met
	NHIm.push_back(0.); //idsh dimension
	Zm.push_back(0.);
      }
    }

    mg0.clear();
    metg0.clear();
    NHIgInt.clear();
    Mratiog.clear();
    idshg.clear();

    mg0.resize(0);
    metg0.resize(0);
    NHIgInt.resize(0);
    Mratiog.resize(0);
    idshg.resize(0);
    
    cout << "::: gas done ::: " << endl;
    cout << "" << endl;
    cout << " ... Writing gas catalog (coords.dc.snap_plane.txt) ..." << endl;
    cout << " .. content: #(ids,x,y,redshift,mass,intial mass,metallicity, age) for stellar particles;" << endl;
    cout << " ..          #(Zgas_w, NHIgas_w) computed on galaxy-basis, for galaxies with id==is." << endl;
    cout << " ..          #(dust-free flux in Nfilters) for stellar particles." << endl;
    cout << endl;
 

    // Write coord catalogue - Nsh and NshTNG included
    string gcoord_path =rdir + "coords.dc." + snappl + "_" + conv(iplrestart, fINT) + ".txt";
    ofstream myfile2g;
    myfile2g.open(gcoord_path);
  
    for (auto k = 0; k < idsh.size(); k++) {
      for (int i = 0; i < idshsmap.size(); i++) {
	if (idsh[k] == idshsmap[i]) {
	  myfile2g << idshsmap[i] << " " << xmap[i] << " " << ymap[i] << " " << zredsmap[i] << " "
		   << ms4map[i] << " " << ims4map[i] << " " << mets4map[i] << " " << ages4map[i] << " "
		   << Zm[k] << " " << NHIm[k] << " " ;       
	  for (int j = 0; j < fluxmap[i].size(); j++) {
	    myfile2g << fluxmap[i][j] << " ";
	  }
	  
	  myfile2g << Npsh[k] << " " << NpshTNG[i] << endl;
	}
      }
    }
    myfile2g.close();
    
  }

  else if (read==0){

    // Cosmology TNG100
    bs = boxl*1000.; // [ckpc/h]
    om0 = 0.3089;
    omL0 = 0.6911;
    dlsim=getY(zl,dl,zsim);
    //
    cout << " .. Extracting gas properties from gas catalog (coords.dc.snap_plane.txt) .." << endl;
    cout << " .. ! This is an old version of coords.xx files ! .." << endl;
    cout << " .. content: #(ids,x,y,redshift,mass,intial mass,metallicity, age) for stellar particles;" << endl;
    cout << " ..          #(Zgas_w, NHIgas_w) computed on galaxy-basis, for galaxies with id==is." << endl;
    cout << " ..          #(dust-free flux in Nfilters) for stellar particles." << endl;
    cout << " ..          #(Np_sh, N_sim): number of particles with ids in the lightcone, number of particle with ids in the original simulation." << endl;
    
    cout << endl;

    // Read gas properties PARTICLE-BASED -> to be transformed into galaxy-based
    vector<int> idsh0p(0), Npshp(0), NpshTNGp(0);
    vector<double> Zmp(0), NHImp(0);
    string filoutcatgas = rdir + "/coords.dc." + snappl + "_" + conv(iplrestart, fINT) + ".txt";
    ifstream ocfgas;
    ocfgas.open(filoutcatgas.c_str());
    if (ocfgas.is_open()) {
      string line;
      while (getline(ocfgas, line)) {
        stringstream ss(line);
        int shid, npsh, nptng;
        double ZM, NHIM;
	float xi,yi,zri,mi,imi,ai;
	long double meti;
        ss >> shid >> xi >> yi >> zri >> mi >> imi >> meti >> ai >> ZM >>  NHIM;
        for (int i = 0; i < Nfilters; i++) {
	  double temp;
            ss >> temp;
        }
        ss >> npsh >> nptng;
        idsh0p.push_back(shid);
        Zmp.push_back(ZM);
        NHImp.push_back(NHIM);
        NpshTNGp.push_back(nptng);
    }
      ocfgas.close();
    } else {
      cout << "  " << endl;
      cout << " OUT_CAT file: " << filoutcatgas << endl;
      cout << " does not exist for this snapshot. " << endl;
      cout << " I will STOP here!!! " << endl;
      exit(1);
    }
    
    //gas (galaxy based)
    std::unique_copy(idsh0p.begin(), idsh0p.end(), std::back_inserter(idsh0));
    
    for(auto i=0;i<idsh0.size();i++){
      int k=-1;
      for (auto j=0;j<idsh0p.size();j++){
	if(idsh0[i]==idsh0p[j]){
	  if (k==-1){
	    NpshTNG.push_back(NpshTNGp[j]);
	    Zm.push_back(Zmp[j]);
	    NHIm.push_back(NHImp[j]);
	  }
	  k=k+1;
	}
      }
    }

    cout << "::: gas read ::: " << endl;
    cout << "" << endl;

  }
 
  cout << "... Now assigning fluxes ..." << endl;
  cout << " " << endl;
  std::vector <vector<double> > flux(0), df(0);  

  // assigning  fluxes to galaxies in the FoV - as a sum of particles fluxes and applying mean dust extinction
  for (auto k=0;k<idsh.size();k++){ //galaxy-based cycle
    vector<long double> tmpSpe;
    vector <vector<long double> > SPE;
    float zr=0.;    
    std::vector<long double> modC;
    std::vector<long double> wavesc;
    if (model=="bc03"){
      wavesc.resize(1221);
      std::copy ( waves.begin(), waves.end(), wavesc.begin() );
    }
    else if (model=="cb16"){  
      wavesc.resize(13391);
      std::copy ( waves.begin(), waves.end(), wavesc.begin() );
    }
    
    std::vector<float> ex_ls(15);
    
    for(int l=0; l<totPartMapxy4; l++){ //particle-based cycle
      if (idshsmap[l]==idsh[k]){

	zr=zrsh[k];  
	std::vector<long double> spe;	
	std::vector<long double> modB(0);
	std::vector<long double> nspe(0);

        // approximating ages4 to age_bc03 (we need it when opening bc03 vec)
        int age_inx=index_closest(age_bc03.begin(), age_bc03.end(), ages4map[l]);

	//approximatin mets to met_model and extract SED
	int met_inx;
	if (model=="bc03"){
	  met_inx = index_closest(met_bc03.begin(), met_bc03.end(), mets4map[l]);      
	  SEDbc03_interp_2spec(full_table[met_inx], time_grid[met_inx], age_inx, ages4map[l], spe);
	}
	else if (model=="cb16"){     
	  met_inx = index_closest(met_cb16.begin(), met_cb16.end(), mets4map[l]);
	  SEDcb16_extract_spec(full_table[met_inx], time_grid[met_inx], age_inx, ages4map[l], spe);
	}
       
	// reddened spectra: model B here; model C in outer loop with mean quantities
	//modB for unresolved dust Nelson+19 (rest-frame values)	
	for (auto i=0;i<wavesc.size();i++){
	  if (ages4map[l]<=0.01){ //Gyr 
	    modB.push_back(expl(-1.*pow((wavesc[i]/5500.),-0.7))); 
	  }
	  else if(ages4map[l]>0.01){
	    modB.push_back(expl(-0.3*pow((wavesc[i]/5500.),-0.7)));
	  }
	}
	
	for (auto i=0;i<wavesc.size();i++){
	  nspe.push_back(spe[i]*modB[i]); 
	} 
	spe.clear();
	spe.resize(0);
	
	modB.clear();
	modB.resize(0);
	// end modB

	for(auto i=0;i<nspe.size();i++){
	  tmpSpe.push_back(nspe[i]*ims4map[l]);
	}
	
	nspe.clear();
	nspe.resize(0);

	SPE.push_back(tmpSpe);

	tmpSpe.clear();
	tmpSpe.resize(0);
      }
    }
       
    //galaxy spectrum
    auto galSpe = m_col_add(SPE);   //with modB

    //compute modC Nelson+19 with mean values   
    if (NHIm[k]!=0.){
      //model C Nelson+19
      sedy.modelC(zr, ex_curve,ex_l,ex_ls, modC, NHIm[k],Zm[k]); //obs-frame
    }
    else{
      for (auto k=0;k<ex_curve.size();k++){
	modC.push_back(1.);
      }
    }
    
    // z evol
    sedy.z_evol(zr,wavesc,galSpe,zl,dlum ); //obs-frame apparent mag (z and Dl^2): model B applied
   
    // dust attenuation (fix that it goes through this also if modC==0)
    sedy.dust_attenuation(ex_ls,modC,wavesc,galSpe); //obs-frame apparent mag: model B and C now applied
    
    // computing apparent magnitude in chosen filter; flux in [uJ]
    vector<double> f(0);
    vector<double> dfn(Nfilters, 0.0); // , m(0);
    for(auto N=0; N<Nfilters; N++){
      double m_ = sedy.compute_mab(zr,wavesc,galSpe,fwaves[N],fresp[N]);    //obs-frame apparent mag in chosen filter; zr is useless parameter, B+C
      double f_= pow(10,(29-(m_+48.6)/2.5)); //uJ ;
      f.push_back(f_);

      if (fluxsh[k][N] > f_) {
	dfn[N] = f_ / fluxsh[k][N];
      } else if(fluxsh[k][N]<=f_){
	dfn[N] = 1.0;
      }
    }
    
    flux.push_back(f);
    df.push_back(dfn);
    
    SPE.clear();
    galSpe.clear();
    modC.clear();
    ex_ls.clear();
    wavesc.clear();
    f.clear();

    SPE.resize(0);
    galSpe.resize(0);
    modC.resize(0);
    ex_ls.resize(0);
    f.resize(0);
    wavesc.resize(0);
    
  }

  auto stop3 = high_resolution_clock::now(); 
  auto duration3 = duration_cast<microseconds>(stop3 - start);
  cout << "execution time in [h] " <<  2.77778e-10*duration3.count() << endl;
  		
  cout << " " << endl;
  cout << "... Reddening the flux particle by particle, for the image ..." << endl;

  std::vector<std::vector<float>> finalflux; 
  finalflux.resize(idshsmap.size()); 
  
  for (int l = 0; l < idshsmap.size(); l++) {
    for (int k = 0; k < idsh.size(); k++) {
      if (idsh[k] == idshsmap[l]) {
	for (int j = 0; j < fluxmap[l].size(); j++) {
	  finalflux[l].push_back(fluxmap[l][j] * df[k][j]);
	}
      }
    }
  }

  cout << "... Writing output file (flux.dc.snap_plane.txt) from this module ..." << endl;
  cout << " .. content: #(ids,x,y,redshift,mass,intial mass,metallicity, age) for stellar particles;" << endl;
  cout << " ..          #(Zgas_w, NHIgas_w) computed on galaxy-basis, for galaxies with id==is." << endl;
  cout << " ..          #(Np_sh, N_sim): number of particles with ids in the lightcone, number of particle with ids in the original simulation." << endl;
  cout << " ..          #(reddened flux in Nfilters) for stellar particles." << endl;
  cout << endl;
  
  
  // print particles in fov on a txt file -- physical chars + unreddened fluxes followed by reddened fluxes
  string ncoord_path =rdir + "flux.dc." + snappl + "_" + conv(iplrestart, fINT) + ".txt";
 
  ofstream myfile2p;
  myfile2p.open(ncoord_path);
  for (auto k = 0; k < idsh.size(); k++) {
    for (int i = 0; i < idshsmap.size(); i++) {
      if (idsh[k] == idshsmap[i]) {
	myfile2p << idshsmap[i] << " " << xmap[i] << " " << ymap[i] << " " << zredsmap[i] << " "
		 << ms4map[i] << " " << ims4map[i] << " " << mets4map[i] << " " << ages4map[i] << " "
		 << Zm[k] << " " << NHIm[k] << " " << Npsh[k] << " " << NpshTNG[k] << " ";
   
	for (int f = 0; f < Nfilters; f++) {
	  myfile2p << finalflux[i][f] << " ";
	}
	
	myfile2p << endl;
      }
    }
  }
  myfile2p.close();



  cout << " " << endl;
  cout << " ... Flux file written ... " << endl;
  cout << "" << endl;
  cout << " " << endl;
  
  int ntotxy4=finalflux.size();
  cout << ntotxy4 << "        - (4) particles in final dust-corrected fov. " << endl;
  cout << " " << endl;
  cout << " " << endl;
  cout << "... But I can only show you the door. You're the one that has to walk through it ..." << endl;
  //clock
  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start);
  cout << "execution time in [h] " <<  2.77778e-10*duration.count() << endl;
  cout << " ------------------------------" << endl;


}
