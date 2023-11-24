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
#ifdef COSMOLIB
#include <cosmo.h>
#endif
#include "/usr/include/hdf5/serial/H5Cpp.h"
#include </usr/include/eigen3/Eigen/Dense>
#include "readTNGParticle.h"
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include "SED.h"
#include "functions.h"


/*****************************************************************************/
/*                                                                           */
/*             FORECAST - dust-free fluxes calculations module               */
/*                                                                           */
/*  original dark matter-only code by cgiocoli@gmail.com                     */
/*  updated to its final form by flaminia.fortuni@inaf.it                    */
/*  if you use it or do any mods, please cite Fortuni et al. (2023).         */
/*                                                                           */
/*                                                                           */  
/*  this is the second module of the FORECAST code;                          */
/*  it computes dust-free fluxes, including k-correction (no dust, no IGm).  */
/*  it runs on one snapshot per time.                                        */
/*  - input: lightcone built with previous module lc                         */
/*  - output: dust-free particle catalog in N filters                        */ 
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
using namespace arma;
using namespace std::chrono; 


const int bleft = 24;
const double speedcunit = 2.99792458e+3;
const double speedcunitas = 2.9979e+18; //AA/s
const double ergsa = 3.839e+33;
const double mpM = 8.4089382e-58; //proton mass in Msun
const double mpg = 1.6726219e-24; //proton mass in g

// ... TNG Cosmology
const double h0 = 0.6774;
const double om0 = 0.3089;
const double omL0 = 0.6911;




int main(int argc, char** argv){
  auto start = high_resolution_clock::now();
  cout << "----------------------------------------------------------------------" << endl;
  cout << " " << endl; 
  cout << "   ------------------------------------------------------ " << endl;
  cout << "   -                                                    - " << endl;
  cout << "   -           2D Mapping Simulation Snapshot           - " << endl;
  cout << "   -                                                    - " << endl;
  cout << "   -                 Let there be light!                - " << endl;
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
  // ... box_length in [Mpc/h], highest redshift and distance reached by the cone
  double boxl,zs;
  // ... input files
  string filredshiftlist,filsnaplist, filtimelist, filfilters, idc;
  // ... directories
  string pathsnap, bc03dir, rdir; 
  // ! SED units: L_lambda; sed resolution for bc03: lr.
  string model,imf;
  string pixunits="uJy";
  
  int sourceID=std::atoi(argv[1]);
  float blD, blD2;
  string snappl;
  double zsim;
  int n_pl=0;
  
  readParameters(&boxl,&zs,
		 &filfilters,&filredshiftlist,&filsnaplist,&filtimelist,&idc,
		 &pathsnap,&bc03dir,&rdir,
		 &model,&imf);

  
  // planes list created by lc module
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
    cout << " planes list file " << fplane << endl;
    cout << " does not exists for this snapshot. " << endl;
    cout << " I will STOP here!!! " << endl;
    exit(1);
  }
  cout << "... Reading " << fplane << " file ... " << endl;
  cout << "( " << blD << ", " << blD2 << " ), " << snappl << ", " << zsim << endl;
  cout << " " << endl;
  

  // redshift list redshift_list
  ifstream redlist;
  redlist.open(filredshiftlist.c_str());
  vector <int> tsnaplist;
  vector <double> dtsnaplist;
  vector <double> tredlist;
  int nmax=1024;
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
    }
    redlist.close();
  }else{
    cout << " redshift list file " << filredshiftlist << " does not " << endl;
    cout << " exist in the Code dir ... check this out      " << endl;
    cout << "    I will STOP here !!! " << endl;
    exit(1);
  }
 
  cout << "... Reading filter list ..." << endl;
  // read filter list
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
    cout << " filter list file " << filfilterlist << "does not " << endl;   
    cout << " exist in the Code dir ... check this out      " << endl;
    cout << "    I will STOP here !!! " << endl;
    exit(1);
  }

  // read filter files  
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
 
  // read timelist
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
    cout << " time list file " << filtimelist  <<" does not " << endl;
    cout << " exist in the Code dir ... check this out      " << endl;
    cout << "    I will STOP here !!! " << endl;
    exit(1);
  }
  
  // open and read LCDM-comovingdistance
  ifstream infiledc;
  vector<double> zl, dl, dlum;
  infiledc.open(idc.c_str());
  if(infiledc.is_open()){
    double zi,dci,dli;
    while(infiledc >> zi >> dci >> dli){
      zl.push_back(zi);
      dl.push_back(dci*speedcunit); // comoving Mpc/h. dl is a comoving distance, not a luminosity distance!
      dlum.push_back(dli*speedcunit); // comoving Mpc/h. Luminosity distance.
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
  
 
  
  cout << "" << endl;
  cout << "-------------------------------------------" << endl;
  cout << " " << endl;
  cout << " " ;
  cout << "... Starting to read TNG for snapshot " << snappl << " in filters " << endl ;
  for (int i = 0; i < nfilter.size(); i++) cout << nfilter[i] << "  ||  " ;
  cout << " ..." << endl;
  cout << " " << endl;
  
 
  double bs = boxl*1000.; // [ckpc/h];	          
  // .. Star particles                                                            
  vector<int> ids(0), idSH4(0), idshs(0), NpshTNG(0);
  vector<float> xs(0),ys(0),zstar(0),zreds(0),ms4(0),ims4(0),ages4(0),CMzsh(0);
  vector<long double> mets4(0);
  vector<double> fiub(0);
  // other
  vector<long double> met_bc03 { 0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05};
  vector<long double> met_cb16  { 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, 0.010, 0.014, 0.017, 0.020, 0.030, 0.040 };
  // .. Class to manipulate the SED
  SED sedy;
  
  
  //open age_bc03/age_cb16
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

  // SSP EXTRACTION
  std::vector<long double> waves; //[AA]
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
 
    
  cout << ".. Reading input cone particle-based catalog from lc module.." << endl;
  cout << " .. content: #(ids,x,y,redshift,mass,intial mass,metallicity, -, age, -, CM of sh with ids, Nparticles of sh with ids in sim) for stellar particles;" << endl;
  cout << endl;
  // Reading particle catalog created by the previous module lc
  string filoutcat="../lc/coords.lc."+snappl+"_"+conv(iplrestart,fINT)+".txt";
  ifstream ocf;
  ocf.open(filoutcat.c_str());
  if(ocf.is_open()){
    int shid,idsi, NPTNG;
    float xi,yi,zi,zri,mi,imi,ai,disttmp, sni, fiubb, cmzsh;
    long double meti, mabi;
    while(ocf >> shid >> xi >> yi  >> zri >> mi >> imi >> meti >> fiubb >> ai >> zi >> cmzsh >>  NPTNG ){
      if (shid>-1){
	idshs.push_back(shid);
	xs.push_back(xi);
	ys.push_back(yi);
	zreds.push_back(zri);
	ms4.push_back(mi);      
	zstar.push_back(zi);
	ims4.push_back(imi);
	mets4.push_back(meti);
	ages4.push_back(ai);     
	fiub.push_back(fiubb);
	CMzsh.push_back(cmzsh);
	NpshTNG.push_back(NPTNG);
      }  
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

  
  // Check print
  cout << " " << endl; 
  cout << "... Some parameters ..." << endl;
  cout << " (D_l, D_l2) = (" << blD << ", " << blD2 << ") cMpc; z between (" << (getY(dl,zl,blD)) << ", " << (getY(dl,zl,blD2))  << "); snapshot " << sourceID << ", plane " << n_pl <<  endl;
  cout << "" << endl;
  cout << "....................... " << endl;
  cout << " " << endl;
  cout << " " << endl;
  cout << "... Starting to read TNG cat ..." << endl;
  cout << " " << endl;
    
  
  int n4 = xs.size();  
  int totPartxy4;
  
  cout << "... Min and max coordinates of star particles in fov ..."<<endl;
  cout << "" << endl;
  
  if(n4>0){
    // quadrate box
    double xmin=double(*min_element(xs.begin(), xs.end()));
    double xmax=double(*max_element(xs.begin(), xs.end()));  
    double ymin=double(*min_element(ys.begin(), ys.end()));
    double ymax=double(*max_element(ys.begin(), ys.end()));  
    double zmin=double(*min_element(zstar.begin(), zstar.end()));
    double zmax=double(*max_element(zstar.begin(), zstar.end()));
    double redmin=double(*min_element(zreds.begin(), zreds.end()));
    double redmax=double(*max_element(zreds.begin(), zreds.end()));
    cout << " " << endl;
    cout << " n4 particles >>>> " << n4 << endl;
    cout << "xmin = " << xmin << endl;
    cout << "xmax = " << xmax << endl;
    cout << "ymin = " << ymin << endl;
    cout << "ymax = " << ymax << endl;
    cout << "zmin = " << zmin << endl;
    cout << "zmax = " << zmax << endl;
    cout << "min redshift = " << redmin << endl;
    cout << "max redshift = " << redmax << endl;
    
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
    
		
    totPartxy4=xs.size();
    std::vector <vector<double> > flux(0);	
    cout << totPartxy4 <<"   type (4) - STAR particles selected in fov."<<endl;
    if (totPartxy4>0)
      cout <<  "first subhalo ID: " << idshs[0] << "  |  last subhalo ID: " << idshs[totPartxy4-1] << endl;
    cout << "" << endl;   
    cout << "      +++... Dust is not included in this run ...+++" << endl;	
    cout << "" << endl;
    cout << "... Now assigning fluxes ..." << endl;
    cout << "" << endl;
    
    // assigning  photometry to particles in the FoV
    for(int l=0; l<totPartxy4; l++){
      
      float zr=zreds[l];	                      
      
      std::vector<long double> spe;
      std::vector<long double> wavesc;
      std::vector<float> ex_ls(15);
                  
      // approximating ages4 to age_bc03 (we need it when opening bc03/cb16 vec)
      int age_inx=index_closest(age_bc03.begin(), age_bc03.end(), ages4[l]);         
	  
      // taking the SSP's SED from BC03 or CB16 tables 
      int met_inx;
      if (model=="bc03"){
	wavesc.resize(1221);
	std::copy ( waves.begin(), waves.end(), wavesc.begin() );
	met_inx = index_closest(met_bc03.begin(), met_bc03.end(), mets4[l]);      
	SEDbc03_interp_2spec(full_table[met_inx], time_grid[met_inx], age_inx, ages4[l], spe);
      }
      else if (model=="cb16"){
	wavesc.resize(13391);
	std::copy ( waves.begin(), waves.end(), wavesc.begin() );
	met_inx = index_closest(met_cb16.begin(), met_cb16.end(), mets4[l]);
	SEDcb16_extract_spec(full_table[met_inx], time_grid[met_inx], age_inx, ages4[l], spe);
      }

      // z evolution	 
      sedy.z_evol(zr,wavesc,spe,zl,dlum );

      // computing apparent magnitude in chosen filter; flux in [uJ]
      vector<double> f(0); // , m(0);
      for(auto N=0; N<Nfilters; N++){
	double m_ = sedy.compute_mab(zr,wavesc,spe,fwaves[N],fresp[N]);    //obs-frame apparent mag in chosen filter; zr is useless parameter, B+C
	m_ = m_ - 2.5*log10(ims4[l]);
	double f_= pow(10,(29-(m_+48.6)/2.5)); //uJ ;
	f.push_back(f_);	
      }
      flux.push_back(f);

      f.clear();
      wavesc.clear();
      spe.clear();
      ex_ls.clear();
      
      wavesc.shrink_to_fit();
      spe.shrink_to_fit();
      ex_ls.shrink_to_fit();
      f.shrink_to_fit();
      
    }
   
    if (flux.size()==totPartxy4){
      cout << " ****-fluxes computed-****" << endl;
      cout << "" << endl;
    }
    else{
      cout << "" << endl;
      cout << " ####-fluxes not computed, please check-####" << endl;
      cout << "" << endl;
    }

    cout << " .. Writing the particle-based, dust-free flux catalogue .." << endl;
    cout << " .. content: #(ids,x,y,redshift,mass,intial mass,metallicity, CM of sh with ids, age, Nparticles of sh with ids in sim) for stellar particles;" << endl;
    cout << " ..          #(dust-free flux in Nfilters) for stellar particles." << endl;
    cout << endl;
    
    //printing particles within fov on a txt file    
    string coord_path="flux.df."+snappl+"_"+conv(iplrestart,fINT)+".txt";
    ofstream myfile2;
    myfile2.open(coord_path);
    for (int i=0;i<totPartxy4; i++){ //	   
      myfile2 << idshs[i] << " " << xs[i] << " " << ys[i] << " " << zreds[i] << " "
	      << ms4[i] << " " << ims4[i] << " " << mets4[i] << " " << CMzsh[i] << " "
	      << ages4[i] << " " << NpshTNG[i]  << " " ;
      for(auto N=0; N<flux[0].size(); N++){
	myfile2 << flux[i][N] << " " ;
      }
      myfile2 << endl;
    }
    myfile2.close();

    cout << ".. Flux file written .." << endl;

    // making M_tot, F_tot sums over particles and print them 
    vector<long double> F_tot(Nfilters,0.);
    long double M_tot=0.;
    
    for(auto N=0; N<Nfilters; N++){
      for (int i=0;i<totPartxy4; i++)  {
	F_tot[N] += flux[i][N];
      }
    }
    for (int i=0;i<totPartxy4; i++)
	M_tot += ms4[i];
    
    for(auto N=0;N<Nfilters;N++)
      cout << " F_tot (total flux in fov [uJy]) for  " << nfilter[N] << " filter: "  << F_tot[N] << endl;
    cout << "" << endl;
    cout << " M_tot (total mass in fov [M_sun]): " << M_tot << endl;
    cout << "" << endl;	
    cout << " N4_tot (total stellar particles in fov):  " << totPartxy4 << endl;
    cout << "" << endl;
      
  }//If n4>0

  cout << " " << endl;
  cout << ""<< endl;
   		
  // end of the work
  std::cout << " ... I'm trying to free your mind, Neo ... once again ..." << endl;
  std::cout << endl;
  auto stop = high_resolution_clock::now(); 
  auto duration = duration_cast<microseconds>(stop - start);
  cout << "execution time in [h] " <<  2.77778e-10*duration.count() << endl;
  std::cout << "----------------------------------------------------------------------" << endl;             

}
