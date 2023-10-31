#include <string>
#include <vector>
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <dirent.h>
#include <map>
#include "functions.h"
using namespace std;


// id coupling function - this is built for data structure of TNG
vector<int> inverseMap_sh_idv3(vector<int> &gcLenType, vector<int> &gcOffsetsType, vector<int> &part_id, bool flagfuzz){  
  vector<int> val = searchsorted(gcOffsetsType, part_id);
  if (flagfuzz == true){
    vector<int> gcOffsetsMax= vec_sum(gcOffsetsType,gcLenType);    
    std::vector<int> ww;
    ww.reserve(val.size());
    for(size_t i=0 ; i<val.size() ; ++i)
      if(part_id[i] > gcOffsetsMax[val[i]])
        ww.push_back(i);    
    if (!ww.empty()){
      for (auto i=0;i<ww.size();i++){	
	val[ww[i]]= -1;
      }
    }
  }
  return val;
}

vector<int> searchsorted(vector<int> &a,vector<int> &par_ind){
  vector <int> diff=vec_diff(a);
  vector<int> val;  
  for (int i=0; i< par_ind.size(); i++){
    auto lower = std::lower_bound(diff.begin(), diff.end(), par_ind[i]);
    // check that value has been found
    const bool found = lower != diff.end() && *lower == par_ind[i];
    auto idx = std::distance(diff.begin(), lower);
    val.push_back(idx-1);
  }
  return val;
}

vector<int> vec_diff(vector<int> &a){
  vector<int> diff;
  for (int i=0; i<a.size(); i++){
    diff.push_back(a[i]-1);
  }
  return diff;
}

vector<int> vec_sum(vector<int> &a,vector<int> &b){
  vector<int> summy;
  for (int i=0; i<a.size(); i++){
    summy.push_back(a[i]+b[i]-1);
  }
  return summy;
}

// interpolating the corresponding yi, giving xi and comparing with vectors (x,y)
double getY(std:: vector<double> x, std:: vector<double> y,double xi){
  int nn = x.size();
  if(x[0]<x[nn-1]){         
    if(xi>x[nn-1]) return y[nn-1];
    if(xi<x[0]) return y[0];
  }  
  else{
    if(xi<x[nn-1]) return y[nn-1];
    if(xi>x[0]) return y[0];  
  }
  int i = locate (x,xi);
  i = std::min (std::max (i,0), int (nn)-2);
  double f=(xi-x[i])/(x[i+1]-x[i]);
  if(i>1 && i<nn-2){
    double a0,a1,a2,a3,f2;                                     
    f2 = f*f;
    a0 = y[i+2] - y[i+1] - y[i-1] + y[i];
    a1 = y[i-1] - y[i] - a0;
    a2 = y[i+1] - y[i-1];
    a3 = y[i];
    return a0*f*f2+a1*f2+a2*f+a3;
  }                                                                      
  else return f*y[i+1]+(1-f)*y[i];           
}

// polar coords
void getPolar(double x, double y, double z, double *ra, double *dec, double *d){
  *d = sqrt(x*x+y*y+z*z);
  *dec = asin(x/(*d));
  *ra = atan2(y,z);
}

// locating x in vector v
template <class T>
int locate (const std::vector<T> &v, const T x){
  size_t n = v.size ();
  int jl = -1;
  int ju = n;
  bool as = (v[n-1] >= v[0]);
  while (ju-jl > 1){
    int jm = (ju+jl)/2;
    if ((x >= v[jm]) == as)
      jl=jm;
    else
      ju=jm;
  }
  if (x == v[0])
    return 0;
  else if (x == v[n-1])
    return n-2;
  else
    return jl;
}

int countHDF5Files(const std::string& path, const std::string& extension) {
    int count = 0;
    DIR* dir = opendir(path.c_str());
    if (dir) {
        struct dirent* entry;
        while ((entry = readdir(dir))) {
            std::string filename = entry->d_name;
            if (filename.find(extension) != std::string::npos) {
                count++;
            }
        }
        closedir(dir);
    } else {
        std::cerr << "Error opening directory: " << path << std::endl;
    }
    return count;
}

// read .ini file parameters
void readParameters(int *npix, int *bufferpix, double *boxl,
		    double *zs, double *fov,
		    string *filredshiftlist,string *filsnaplist, string *filtimelist, string *idc,
		    string *pathsnap, string *rdir, 
		    long *seedcenter, long *seedface, long *seedsign){ 

  string butstr;
  ifstream inputf;
  inputf.open("MapSim.ini");
  if(inputf.is_open()){
    inputf >> butstr; // number of pixels by side
    inputf >> *npix;
    inputf >> butstr; // number of buffer pixels in total
    inputf >> *bufferpix;
    inputf >> butstr; // box_length of the simulation in [Mpc/h]
    inputf >> *boxl;
    inputf >> butstr; // source redshift
    inputf >> *zs;
    inputf >> butstr; // field of view in degrees by side
    inputf >> *fov;
    inputf >> butstr; // file with the redshift list; it must contain three columns: [snap, 1/(1+z), z]
    inputf >> *filredshiftlist;
    inputf >> butstr; // file with the snapshot list available 
    inputf >> *filsnaplist;
    inputf >> butstr; // file with the snapshot vs time list available 
    inputf >> *filtimelist;
    inputf >> butstr; // path and file name of the comoving distance file (if not available create with astropy and sim cosmology)
    inputf >> *idc;
    inputf >> butstr; // path where the snaphosts are located
    inputf >> *pathsnap;
    inputf >> butstr; // path where outputs are located
    inputf >> *rdir;    
    inputf >> butstr; // seed for the random location of the center
    inputf >> *seedcenter;
    inputf >> butstr; // seed for the random selection of the dice face
    inputf >> *seedface;
    inputf >> butstr; // seed for the selection of the sign of the coordinates
    inputf >> *seedsign;  
    inputf.close();
  }else{
    cout << " INPUT file does not exist ... I will stop here!!! " << endl;
    exit(1);
  }
};

