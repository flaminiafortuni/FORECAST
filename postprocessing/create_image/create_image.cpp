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
#include <CCfits/CCfits>
#include <FitsError.h>
#include "functions.h"


/***************************************************************************/
/*                                                                         */
/*            postprocessing FORECAST - creating the raw image             */
/*                                                                         */
/*  if you use it or do any mods, please cite Fortuni et al. (2023).       */
/*                                         flaminia.fortuni@inaf.it        */
/*                                                                         */
/*                                                                         */  
/*  this is a postprocessing module for the FORECAST code;                 */
/*  it creates the raw image from the particle catalog                     */
/*  it runs on one snapshot/plane per time.                                */
/*  - input: particle catalog created from either df,dc or igm module      */
/*  - output: raw image in the chosen filter, one snapshot/plane per time  */ 
/*                                                                         */
/*                                                                         */
/*  for a comprehensive guide of FORECAST, visit                           */
/*                          https://github.com/flaminiafortuni/FORECAST    */
/*  for a full description of the software                                 */     
/*       https://ui.adsabs.harvard.edu/abs/2023arXiv230519166F/abstract    */
/*                                                                         */   
/***************************************************************************/


using namespace std;
using namespace CCfits;

const double h0=0.6774;


int main(int argc, char** argv){

  cout << "----------------------------------------------------------------------" << endl;
  cout << " " << endl; 
  cout << "   ------------------------------------------------------ " << endl;
  cout << "   -                                                    - " << endl;
  cout << "   -           2D Mapping Simulation Snapshot           - " << endl;
  cout << "   -                                                    - " << endl;
  cout << "   -               creating the final image             - " << endl;
  cout << "   ------------------------------------------------------ " << endl;
  cout << " " << endl;
  cout << " " << endl; 
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

  // reading files
  int sourceID=std::atoi(argv[1]);
  int n_pl;
  float blD, blD2, zsim, invh0=1./h0;
  string fplane= "../../lc/planes_list.txt";
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
	zsim=zpl;
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

  string snappl= std::string(argv[1]);
  string plane= std::string(argv[2]);
  // ******************** to be read in the INPUT file ********************
  // ... fov in degrees
  float fov,res;
  // ... filter & module to read
  string filter, module, rdir;
  // .. filter position in filter list
  int pos;

  readParameters(&fov,&res,
		 &filter,
		 &module,
		 &rdir);
  
  //pixels of the image
  unsigned long int truenpix=int(fov*3600/res);
  int bufferpix=int(ceil((truenpix + 1)*20 / 14142));  // add bufferpix/2 in each side!
  unsigned long int npix=truenpix+bufferpix;

  cout << "N. pixels: " << truenpix << "; buffer pixels: " << bufferpix << endl;
  cout << endl;

  vector<float> idshs(0), xs(0),ys(0),fH(0);
  vector<vector<double>> fluxmap;

  cout << "... Reading filter list ..." << endl;

  // read filter list
  string filfilters="../../"+module+"/filters.dat";
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
    cout << filfilters << " filter list file does not " << endl;   
    cout << " exist in the Code dir ... check this out      " << endl;
    cout << "    I will STOP here !!! " << endl;
    exit(1);
  }

  //find filter position in filter list
  int Nfilters= nfilter.size();
  for(auto i=0;i<Nfilters;i++){
    if (nfilter[i]==filter){
      pos=i;
    }
  }


  if(module=="dc" || module=="igm"){

    cout << " " << endl;
    cout << "... Reading input file (flux."<< module<<".snap_plane.txt) from " << module << "  module ..." << endl;
    cout << " .. content: #(ids,x,y,redshift,mass,intial mass,metallicity, age) for stellar particles;" << endl;
    cout << " ..          #(Zgas_w, NHIgas_w) computed on galaxy-basis, for galaxies with id==ids." << endl;
    cout << " ..          #(Np_sh, N_sim): number of particles with ids in the lightcone, number of particle with ids in the original simulation." << endl;
    cout << " ..          #(reddened flux in Nfilters) for stellar particles." << endl;
    cout << endl;
    
    // reading input catalogue (dust-corrected or dc+igm-correctedfluxes)
    string filoutcat="../../"+module+"/flux."+module+"."+snappl+"_"+plane+".txt";         
    ifstream ocf;
    ocf.open(filoutcat.c_str());  
    if(ocf.is_open()){
      int shid,npsh,npshtng;
      float xi,yi,zi,zri,mi,imi,ai,Zmi, NHImi;
      long double meti;
      while (ocf >> shid >> xi >> yi  >> zri >> mi >> imi >> meti >> ai >>  Zmi >> NHImi >> npsh >> npshtng) {
	idshs.push_back(shid);
	xs.push_back(xi);
	ys.push_back(yi);
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
      cout << filoutcat << "  file " << endl;
      cout << " does not exists for this snapshot. " << endl;
      cout << " I will STOP here!!! " << endl;
      exit(1);
    }

  }
  else if(module=="df"){

    cout << "... Reading dust-free particle-based catalog produced in the " << module << " module ..." << endl;
    cout << " .. content: #(ids,x,y,redshift,mass,intial mass,metallicity, CM_sh, age, N_sim) for stellar particles;" << endl;
    cout << " ..          #(dust-free flux in Nfilters) for stellar particles." << endl;
    cout << endl;
    
    // reading df output catalogue (dust-free fluxes)
    string filoutcat="../../"+module+"/flux."+module+"."+snappl+"_"+plane+".txt";
    ifstream ocf;
    ocf.open(filoutcat.c_str());
    if(ocf.is_open()){
      int shid, NPTNG;
      float xi,yi,zi,zri,mi,imi,ai,fi, sni, cmzsh;
      long double meti, mabi;
      while (ocf >> shid >> xi >> yi  >> zri >> mi >> imi >> meti >> cmzsh >> ai >> NPTNG) { 
	idshs.push_back(shid);
	xs.push_back(xi);
	ys.push_back(yi);
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
      cout << filoutcat << " file: " << filoutcat << endl;
      cout << " does not exists for this snapshot. " << endl;
      cout << " I will STOP here!!! " << endl;
      exit(1);
    }
  }

  // collecting fluxes for the final image
  for (int j = 0; j < fluxmap.size(); j++) {
    std::vector<float> magsh_j;    
      fH.push_back(fluxmap[j][pos]);
  }

  
  std::valarray<float> mapxy4( npix*npix );
  string pixu="flux [uJ]";
  string fileoutput=rdir+"/"+filter+"."+module+"."+snappl+"_"+plane+".fits";
  
  //make the map
  mapxy4=gridist_nok(xs,ys,fH,npix);
  int ntotxy4= xs.size();
  cout << "N. particles: " << xs.size() << " at snapshot " << snappl << " and plane " << plane << endl;
  cout << endl;
  if(ntotxy4>0){
    
          long naxis = 2;
	  long naxes[2]={ truenpix,truenpix };
	  string count="1";
	  
	  std::unique_ptr<FITS> ffxy( new FITS(fileoutput, FLOAT_IMG, naxis, naxes ) ); 
	  std::vector<long> naxex( 2 );
	  naxex[0]=truenpix;
	  naxex[1]=truenpix;
	  PHDU *phxy=&ffxy->pHDU();
	  // phxy->write( 1, truenpix*truenpix, mapxytot4 );
	  valarray<float> pmap(truenpix*truenpix);
	  pmap = rescalemap(mapxy4,npix,truenpix);
	  cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
	  cout << "@  ... Writing .fits file " << fileoutput << " ... @" << endl;
	  cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
	  cout << " " << endl;
	  long  fpixel(1);
	  phxy->write( fpixel, truenpix*truenpix , pmap );
	  phxy->addKey ("N_PIXEL_BY_SIDE",truenpix," ");
	  phxy->addKey ("BUFFER_PIXELS",bufferpix," ");
	  phxy->addKey ("PIXELUNIT",pixu," ");
	  phxy->addKey ("N_STELLAR_PARTICLES",ntotxy4," ");
	  phxy->addKey ("FILTER",filter," ");	  
	  phxy->addKey ("REDSHIFT",zsim," "); 
	  phxy->addKey ("PHYSICALSIZE_BY_SIDE",fov,"degrees "); 
	  phxy->addKey ("Dl_LOW",blD*invh0,"comoving distance in Mpc"); 
	  phxy->addKey ("Dl_UP",blD2*invh0,"comoving distance in Mpc");
	  phxy->addKey ("HUBBLE",h0," ");
  }
  cout << "--------------------------" << endl;
}


