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
#include <CCfits/CCfits>
#include <typeinfo>

#include "readTNGParticle.h"
#include "SED.h"
#include "IGM_Ino14.h"
#include "functions.h"

/*****************************************************************************/
/*                                                                           */
/*             FORECAST - IGM-corrected fluxes calculations module           */
/*                                                                           */
/*  original dark matter-only code by cgiocoli@gmail.com                     */
/*  updated to its final form by flaminia.fortuni@inaf.it                    */
/*  if you use it or do any mods, please cite Fortuni et al. (2023).         */
/*                                                                           */
/*                                                                           */ 
/*  this is the fourth (and last)  module of the FORECAST code;              */
/*  it assumes IGM recipe from Inoue et al. (2014).                          */
/*  it computes IGM-corrected fluxes, including k-correction.                */
/*  it runs on one snapshot per time.                                        */
/*  - input: dust-orrected particle catalog from previous module df          */
/*  - output: IGM-corrected particle catalog in N filters                    */ 
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
using namespace CCfits;

const double speedcunit = 2.99792458e+3;
const double mpM = 8.4089382e-58; //proton mass in Msun
const double mpg = 1.6726219e-24; //proton mass in g
const double mpc2tocm2 = 9.523e+42; //Mpc^2 to cm^2; used to convert NHI (division)
const double h0 = 0.6774;


int main (int argc, char** argv){
  
  auto start = high_resolution_clock::now();

  cout << "----------------------------------------------" << endl;
  cout << "-                                            -" << endl;
  cout << "::::IGM POST-PROCESSING for the LIGHTCONE::::" << endl;
  cout << "-                                            -" << endl;
  cout << "----------------------------------------------" << endl;


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

  // ******************** from INPUT file ******************** 
  // ... box_length of the simulation [Mpc/h]
  double boxl;
  // ... files
  string filfilters,filsnaplist, filtimelist,idc;
  // ... directories
  string pathsnap, bc03dir, rdir;
  // ... SSP model & IMF
  string model,imf;
  string pixunits="uJy";
 
  readParameters(&boxl,
		 &filfilters,&filsnaplist,&filtimelist,&idc,
		 &pathsnap,&bc03dir,&rdir,
		 &model,&imf);

  // reading files
  
  // lightcone
  int sourceID=std::atoi(argv[1]);
  int n_pl=0;
  float blD, blD2;
  double zsim;
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
    cout << " planes_list.txt: " << fplane << endl;
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
    cout << " time list file timelist_TNG.txt does not " << endl;
    cout << " exist in the Code dir ... check this out      " << endl;
    cout << "    I will STOP here !!! " << endl;
    exit(1);
  }
  
  // filters list in filters.dat
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
    string filterin=rdir + "filters/" + nfilter[N] + ".dat";
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

  // open age_bc03.txt/age_cb16.txt
  ifstream agelist;
  string filagelist=rdir+"age_"+model+".txt";
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
    cout << " time list file " << filagelist << ".txt does not " << endl;
    cout << " exist in the Code dir ... check this out      " << endl;
    cout << "    I will STOP here !!! " << endl;
    exit(1);
  }
  
  // read Cosmology file: z, comovindg_D, luminosity_D (consistent with simulation's Cosmology)
  ifstream infiledc;
  vector<double> zl, dl, dlum;
  infiledc.open(idc.c_str());
  if(infiledc.is_open()){
    double zi,dci,dli;
    while(infiledc >> zi >> dci >> dli){
      zl.push_back(zi);
      dl.push_back(dci*speedcunit); 
      dlum.push_back(dli*speedcunit);
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
  
  //SSP EXTRACTION
  std::vector<long double> waves; //[AA]
  vector<vector<double> > time_grid;
  vector<vector<vector<double> > > full_table;
  if(model=="bc03"){
    waves.resize(1221);
    time_grid.resize(6);
    full_table.resize(6);
    for (int i = 0; i < 6; ++i) {
      readSSPTables(bc03dir, model, imf, i, waves, time_grid[i], full_table[i]);
    }
    
  }
  else if(model=="cb16"){
    waves.resize(13391);
    time_grid.resize(14);
    full_table.resize(14);
    for (int i = 0; i < 14; ++i) {
      readSSPTables(bc03dir, model, imf, i, waves, time_grid[i], full_table[i]);
    }
  }
    
  //Variables declaration
  // 2Dgrid stars
  vector<int> ids(0), idshs(0), idshsmap(0), idsmap(0), Npshmap(0), NpshTNGmap(0);
  vector<float> xs(0),ys(0),zreds(0),ms4(0),ims4(0), zredsmap(0),  mbc, ages4(0);
  vector <double> NHImap(0), Zmap(0);
  vector<long double> mets4(0);
  vector<float> xmap(0),ymap(0), ms4map(0), fluxgaldusty(0), deltaflux(0), deltamag(0), z4map(0), ims4map(0), ages4map(0), mets4map(0);
  vector<vector<double>> fluxmap;
  // other
  vector<long double> met_bc03 { 0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05};
  vector<long double> met_cb16 { 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, 0.010, 0.014, 0.017, 0.020, 0.030, 0.040 }; 
  // classes
  SED sedy;
  IGM igmy;
 
  
  cout << " " << endl;
  cout << "... Reading input file (flux.xxDUST_yy.txt) from previous module ..." << endl;
  cout << " .. content: #(ids,x,y,redshift,mass,intial mass,metallicity, age) for stellar particles;" << endl;
  cout << " ..          #(Zgas_w, NHIgas_w) computed on galaxy-basis, for galaxies with id==is." << endl;
  cout << " ..          #(Np_sh, N_sim): number of particles with ids in the lightcone, number of particle with ids in the original simulation." << endl;
  cout << " ..          #(reddened flux in Nfilters) for stellar particles." << endl;
  cout << endl;
  

  // reading dc16 output catalogue (dust-corrected fluxes)
  string filoutcat="../dc/flux.jwst16."+snappl+"DUST_"+conv(iplrestart,fINT)+".txt";
  ifstream ocf;
  ocf.open(filoutcat.c_str());  
  if(ocf.is_open()){
    int shid,npsh,npshtng;
    float xi,yi,zi,zri,mi,imi,ai,Zmi, NHImi;
    long double meti;
    while (ocf >> shid >> xi >> yi  >> zri >> mi >> imi >> meti >> ai >>  Zmi >> NHImi >> npsh >> npshtng) { 
    idshsmap.push_back(shid);
    xmap.push_back(xi);
    ymap.push_back(yi);
    zredsmap.push_back(zri);
    ms4map.push_back(mi);      
    ims4map.push_back(imi);
    mets4map.push_back(meti);
    ages4map.push_back(ai);
    NHImap.push_back(NHImi);
    Zmap.push_back(Zmi);
    Npshmap.push_back(npsh);
    NpshTNGmap.push_back(npshtng);
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
						

  double bs = boxl*1000.; // [ckpc/h]
  double dlsim=getY(zl,dl,zsim);
  cout << " " << endl;
  cout << " " << endl;
  
  cout << "... Creating SH cat from input dusty star particles cat ..." << endl;
  // subhalo catalogs
  std::vector<int> idsh, Npsh;
  std::vector<float> zrsh;
  std::vector<std::vector<float>> fluxsh;   
  std::unique_copy(idshsmap.begin(), idshsmap.end(), std::back_inserter(idsh));  
  for (int i = 0; i < idsh.size(); i++) {
    int countsh = 0;
    float redsh = 0.;
    std::vector<float> fsh(fluxmap[0].size(), 0.0);
    for (auto l = 0; l < idshsmap.size(); l++) {
        if (idshsmap[l] == idsh[i]) {	 
	  countsh++;
	  redsh += zredsmap[l];	  
	  for (size_t col = 0; col < fluxmap[l].size(); col++) {
	    fsh[col] += fluxmap[l][col];
	  }
        }
    }
    Npsh.push_back(countsh);
    zrsh.push_back(redsh * 1.0 / countsh);
    fluxsh.push_back(fsh);
  }

  
  //stellar particles to be affected by IGM
  for (auto k=0;k<idshsmap.size();k++){
    zreds.push_back(zredsmap[k]);
    idshs.push_back(idshsmap[k]);
    xs.push_back(xmap[k]);
    ys.push_back(ymap[k]);
    ms4.push_back((ms4map[k]));
    ims4.push_back(ims4map[k]);
    mets4.push_back(mets4map[k]);
    ages4.push_back(ages4map[k]); //[Gyr]
    //z4.push_back(z4map[k]);
  }
  

  cout << ""<< endl;
  cout << "... Now assigning fluxes ..." << endl;
  cout << " " << endl;
  int totPartxy4=xmap.size();
  cout << "N star particles: " << totPartxy4 << endl;
  cout << ""<< endl;

  
  std::vector <vector<double> > flux(0), df(0);//(totPartxy4,vector<double>(Nfilters));  
  
  //assigning  fluxes to galaxies in the FoV - as a sum of particles fluxes and applying IGM extinction
  for (auto k=0;k<idsh.size();k++){
    vector<long double> tmpSpe;
    vector <vector<long double> > SPE;
    float zr=0.;
  
    std::vector<long double> wavesc;
    if (model=="bc03"){
      wavesc.resize(1221);
      std::copy ( waves.begin(), waves.end(), wavesc.begin() );
    }
    else if (model=="cb16"){  
      wavesc.resize(13391);
      std::copy ( waves.begin(), waves.end(), wavesc.begin() );
    }
 
    for(int l=0; l<totPartxy4; l++){ //totPartxy4
      if (idshsmap[l]==idsh[k]){
			  
	zr=zrsh[k];  
	std::vector<long double> spe;	
       
        //APPROXIMATING ages4 to age_bc03 (we need it when opening bc03 vec)
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


	for(auto i=0;i<spe.size();i++){
	  tmpSpe.push_back(spe[i]*ims4map[l]);
	}
	
	spe.clear();
	spe.shrink_to_fit();

	SPE.push_back(tmpSpe);

	tmpSpe.clear();
	tmpSpe.shrink_to_fit();


      }
    }
       
    //galaxy spectrum
    auto galSpe = m_col_add(SPE);   //with mass normalization
    
    //z evol
    sedy.z_evol(zr,wavesc,galSpe,zl,dlum ); //obs-frame apparent mag (z and Dl^2); dust already applied
   
    //IGM attenuation
    igmy.igm_absorption(zr,wavesc,galSpe); //obs-frame igm

    // computing apparent magnitude in chosen filter; flux in [uJ]
    vector<double> f(0);
    vector<double> dfn(Nfilters, 0.0);
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
    wavesc.clear();
    f.clear();

    SPE.shrink_to_fit();
    galSpe.shrink_to_fit();    
    wavesc.shrink_to_fit();
    f.resize(0);
    
  }
  		
  cout << " " << endl;
  cout << "... IGM effect on the particle flux ..." << endl;

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

  cout << " " << endl;
  cout << "... Writing output file (flux.xxigm_yy.txt) from this module ..." << endl;
  cout << " .. content: #(ids,x,y,redshift,mass,intial mass,metallicity, age) for stellar particles;" << endl;
  cout << " ..          #(Zgas_w, NHIgas_w) computed on galaxy-basis, for galaxies with id==is." << endl;
  cout << " ..          #(Np_sh, N_sim): number of particles with ids in the lightcone, number of particle with ids in the original simulation." << endl;
  cout << " ..          #(reddened+igm flux in Nfilters) for stellar particles." << endl;
  cout << endl;
  // print particles in fov on a txt file -- physical chars + unreddened fluxes followed by reddened fluxes
  string ncoord_path=rdir+"flux."+snappl+"igm_"+conv(iplrestart,fINT)+".txt";
  ofstream myfile2p;
  myfile2p.open(ncoord_path);
  for (auto k = 0; k < idsh.size(); k++) {
    for (int i = 0; i < idshsmap.size(); i++) {
      if (idsh[k] == idshsmap[i]) {
	myfile2p << idshsmap[i] << " " << xmap[i] << " " << ymap[i] << " " << zredsmap[i] << " "
		 << ms4map[i] << " " << ims4map[i] << " " << mets4map[i] << " " << ages4map[i] << " "
		 << Zmap[i] << " " << NHImap[i] << " " << Npshmap[i] << " " << NpshTNGmap[i] << " ";
   
	for (int f = 0; f < Nfilters; f++) {
	  myfile2p << finalflux[i][f] << " ";
	}
	
	myfile2p << endl;
      }
    }
  }
  myfile2p.close();

  
  cout << " IGM file written." << endl;
  cout << " " << endl;
  cout << "... You are done, Neo ..." << endl;
  //clock
  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start);
  cout << "execution time in [h] " <<  2.77778e-10*duration.count() << endl;


}
