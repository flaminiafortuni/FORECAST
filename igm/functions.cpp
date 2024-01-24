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

//polar coords
void getPolar(double x, double y, double z, double *ra, double *dec, double *d){
  *d = sqrt(x*x+y*y+z*z);
  *dec = asin(x/(*d));
  *ra = atan2(y,z);
}

// extract from cb16 ised file
FileStructIterator::FileStructIterator(const std::string& filename):
    m_filename(filename), m_state(STATE_START) 
{}

std::string FileStructIterator::getNext() {
    switch (m_state) {
    case STATE_START:
        m_state = STATE_LOOP1;
        if (m_filename.find("cb2016") != std::string::npos or m_filename.find("bc2003") != std::string::npos) {
            m_badLineNumber = 5;
        } else {
            m_badLineNumber = 6;
        }
        m_i = 0;
        return "data";

    case STATE_LOOP1:
        ++m_i;
        if (m_i >= m_badLineNumber) {
            m_state = STATE_LOOP2;
        }
        return "unused";

    case STATE_LOOP2:
        return "data";

    default:
        throw std::logic_error("bad state");
    }
}

void read_bc03_ssp(string infile, std::vector<double> &time_grid,  std::vector <vector<double> > &full_table, std::vector<long double> &waves ){

  std::ifstream ifile(infile.c_str());
  if (!ifile) {
    cerr <<"Error in opening the file: "<<infile<<". EXITING now.\n\a";
    exit(1);
  }
  else{
    full_table.clear();
    std::string line;
    int counter=0;
    std::vector<double> tmp_table;
    auto it = FileStructIterator(infile.c_str());
    std::string what_line = it.getNext();
  
    for( string line; getline( ifile, line ); ){
      std::istringstream iss(line);
      std::string token;
      if (what_line=="data"){
	for (std::string token; iss >> token; ){
	  if (counter==0){
	    counter=std::stod(token);
	  }
	  else{
	    tmp_table.push_back(stod(token.c_str()));
	    counter -= 1;
	    if (counter==0){
	      if (time_grid.empty()){
		for (int k=0; k<tmp_table.size(); k++)
		  time_grid.push_back(tmp_table[k]*1e-9);
	      }
	      else if (waves.empty()){
		for (int k=0; k<tmp_table.size(); k++)
		  waves.push_back(tmp_table[k]);

	      }
	      else if (tmp_table.size()>250){
		
		full_table.push_back(tmp_table);
	      }

	      tmp_table.clear();    
	    }
	  }
	  
	}
      }
      
      if (counter==0){
	what_line=it.getNext();
      }
      
    }
    ifile.close();
  }
  //time_grid.erase(time_grid.begin()); //if you delete it, full_table[indx] and not full_tale[indx+1] when reading bc03 tables
  
}

void read_bc03_mass(std::string infile, std::vector <double> &m){

  std::ifstream ifile(infile.c_str());
  if (!ifile) {
    cerr <<"Error in opening the file: "<<infile<<". EXITING now.\n\a";
    exit(1);
  }
      
  std::string line, str;
  
  for (int i=0; i<29; i++)
    std::getline(ifile, str); // skip 29 lines
  for( string line; getline( ifile, line );){
    
    std::istringstream iss(line);
    std::string token;
    int counter=0;
    for (std::string token; iss >> token; ){
      double temp=std::stod(token);
      counter++;
      if (counter==7){
	m.push_back(temp);
      }
    }
  }
  m.insert(m.begin(),1); //insert M=1 Msun in the first position (corresponding to the SED with t=0 Gyr;

}

int index_closest(std::vector<long double>::iterator begin, std::vector<long double>::iterator end, long double value) {   
  auto it = std::lower_bound(begin, end, value);//, std::less<long double>());  
  if ( it == begin )
    return 0; 
  if (it == end)
    return std::distance(begin, end) - 1;  
  double aft = *it - value;
  double bef = value - *(it-1);
  if (aft < bef){
    return std::distance(begin, it);
  }  
  else if(aft > bef){
    return std::distance(begin,it) - 1;
  }    
}

int index_closest(std::vector<double>::iterator begin, std::vector<double>::iterator end, double value) {   
  auto it = std::lower_bound(begin, end, value);//, std::less<long double>());  
  if ( it == begin )
    return 0; 
  if (it == end)
    return std::distance(begin, end) - 1;  
  double aft = *it - value;
  double bef = value - *(it-1);
  if (aft < bef){
    return std::distance(begin, it);
  }  
  else if(aft > bef){
    return std::distance(begin,it) - 1;
  }    
}

int index_closest(std::vector<float>::iterator begin, std::vector<float>::iterator end, float value) {   
  auto it = std::lower_bound(begin, end, value);//, std::less<long double>());  
  if ( it == begin )
    return 0; 
  if (it == end)
    return std::distance(begin, end) - 1;  
  double aft = *it - value;
  double bef = value - *(it-1);
  if (aft < bef){
    return std::distance(begin, it);
  }  
  else if(aft > bef){
    return std::distance(begin,it) - 1;
  }    
}

vector<int> vec_diff(vector<int> &a){
  vector<int> diff;
  for (int i=0; i<a.size(); i++){
    diff.push_back(a[i]-1);
  }
  return diff;
}

vector<double> vec_diff(vector<double> &a){
  vector<double> diff;
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

//read .ini file parameters
void readParameters(double *boxl,
		    string *filfilters, string *filsnaplist, string *filtimelist,string *idc,
		    string *pathsnap,string *bc03dir, string *rdir,  
		    string *model, string *imf){ 

  string butstr;
  ifstream inputf;
  inputf.open("igm.ini");
  if(inputf.is_open()){
    inputf >> butstr; // box_length of the sim box [Mpc/h]
    inputf >> *boxl;
    inputf >> butstr; // filters list
    inputf >> *filfilters;
    inputf >> butstr; // file with the snapshot list available 
    inputf >> *filsnaplist;
    inputf >> butstr; // file with the snapshot list available 
    inputf >> *filtimelist;
    inputf >> butstr; // path and file name of the comoving distance file (if not available use astropy consistently with sim cosmology)
    inputf >> *idc;
    inputf >> butstr; // path where the snaphosts are located
    inputf >> *pathsnap;
    inputf >> butstr; // path where bc03 software is  located
    inputf >> *bc03dir;
    inputf >> butstr; // path where outputs are located
    inputf >> *rdir; 
    inputf >> butstr; // SSP model bc03 OR cb16
    inputf >> *model;
    inputf >> butstr; // IMF for SSP chabrier OR salpeter
    inputf >> *imf;
    inputf.close();
  }else{
    cout << " INPUT file does not exsit ... I will stop here!!! " << endl;
    exit(1);
  }
};


std::string getBC03MetallicityCode(int metallicityIndex) {
    std::map<int, std::string> bc03MetallicityCodes = {
        {0, "m22"},
        {1, "m32"},
        {2, "m42"},
        {3, "m52"},
        {4, "m62"},
        {5, "m72"},
        {6, "m82"}
    };
    if (metallicityIndex < 0 || metallicityIndex > 6) {
        std::cerr << "Incorrect metallicity index for bc03 model." << std::endl;
        return "";
    }
    return bc03MetallicityCodes[metallicityIndex];
}

std::string getCB16MetallicityCode(int metallicityIndex) {
    // Mappa che associa gli indici di metallicità ai codici CB16
    std::map<int, std::string> cb16MetallicityCodes = {
        {0, "z0001"},
        {1, "z0002"},
        {2, "z0005"},
        {3, "z001"},
        {4, "z002"},
        {5, "z004"},
        {6, "z006"},
        {7, "z008"},
        {8, "z010"},
        {9, "z014"},
        {10, "z017"},
        {11, "z020"},
        {12, "z030"},
        {13, "z040"}
    };

    // Verifica se l'indice di metallicità è valido
    if (metallicityIndex < 0 || metallicityIndex > 13) {
        std::cerr << "Incorrect metallicity index for cb16 model." << std::endl;
        return "";
    }

    // Restituisci il codice CB16 di metallicità corrispondente
    return cb16MetallicityCodes[metallicityIndex];
}


void readSSPTables(
    std::string bc03dir,
    const std::string& model,          // "bc03" o "cb16"
    const std::string& imf,            // "salpeter" o "chab"
    int metallicityIndex,              // 1 per "z0001", 2 per "z0002", ecc.
    std::vector<long double>& waves,   // Vettore per le lunghezze d'onda
    std::vector<double>& timeGrid,     // Vettore per la griglia temporale
    std::vector<std::vector<double>>& fullTable // Tabella dati
) {

  waves.clear();
  fullTable.clear();
  timeGrid.clear();
    if (model == "bc03") {
        if (imf != "salpeter" && imf != "chabrier") {
            std::cerr << "Invalid IMF. Choose chab or salp." << std::endl;
            return;
        }

	std::string imf_dir;
	if(imf.compare("chabrier")==0){
	  imf_dir="chab";
	}
	else{
	  imf_dir="salp";
	}
	
	waves.reserve(1221);
	timeGrid.reserve(221);
	fullTable.resize(221, std::vector<double>(1221));
        std::string bc03MetallicityCode = getBC03MetallicityCode(metallicityIndex);
        if (bc03MetallicityCode.empty() || (metallicityIndex < 0 || metallicityIndex > 6)) {
	  cout << "Incorrect metallicity index. Choose between 1 and 7." << endl;
	  return;
        }
        std::string inputFile = bc03dir + "/models/Padova1994/" + imf + "/bc2003_lr_" + bc03MetallicityCode + "_" + imf_dir + "_ssp.ised_ASCII";
        read_bc03_ssp(inputFile, timeGrid, fullTable, waves);



    } else if (model == "cb16") {
      waves.reserve(13391);
      timeGrid.reserve(221);
      fullTable.resize(221, std::vector<double>(13391));
      std::string cb16MetallicityCode = getCB16MetallicityCode(metallicityIndex);

      if (cb16MetallicityCode.empty() || metallicityIndex < 0 || metallicityIndex > 13) {   ///Uno dei due IF
        cout << "Incorrect metallicity index. Choose between 1 and 14." << endl;
        return;
      }
      
      
        std::string cb16Files[] = {
            "_u1p0",
            "_u1p5",
            "_u2p0",
            "_u2p5",
            "_u3p0",
            "_u3p0",
            "_u3p5",
            "_u3p5",
            "_u3p5",
            "_u3p5",
            "_u3p5",
            "_u3p5",
            "_u4p0",
            "_u4p0"
        };

	std::string inputFile = bc03dir + "/models/cb16/star_gas/cb2016_" + cb16MetallicityCode + cb16Files[metallicityIndex] +"_xi3_n2_mup100_C100_noLya.ised_ASCII";
	//cout << inputFile << " " <<  cb16MetallicityCode  <<endl;
	
        read_bc03_ssp(inputFile, timeGrid, fullTable, waves);
    } else {
        std::cerr << "Error: unvalid model. Please use 'bc03' or 'cb16'." << std::endl;
        return;
    }
}







///SED
void SEDbc03_interp_2spec(std::vector <vector<double> > &full_table, std::vector<double> &time_grid,  int a_indx, float ages4,  std::vector<long double> &spe){

  //Interpolating between two spectra of contiguous ages.
  double ergsa = 3.9e+33;
  double t_1,t_2,a_1,a_2;
  vector<double>f_1,f_2;

  if(ages4<=0.){
    for (int i=0; i<1221; i++){
      spe.push_back(full_table[0][i]);//*ergsa);
    }
  }
  

  else if(ages4>=20.){
    for (int i=0; i<1221; i++){
      spe.push_back(full_table[220][i]);//*ergsa);
    }
  }
  
  
  else if (ages4>0. & ages4<20.){
    if ((ages4<=time_grid[a_indx]) & (ages4>time_grid[a_indx-1])){
      cout << "cae i-1 and i" << endl;
      t_1 = time_grid[a_indx];
      t_2 = time_grid[a_indx-1];
      a_1 = (ages4-t_2)/(t_1-t_2);
      a_2 = 1-a_1;
      f_1 = full_table[a_indx];
      f_2 = full_table[a_indx-1];
    }
    else if ((ages4<time_grid[a_indx+1]) & (ages4>=time_grid[a_indx])){
      cout << "case i and i+1" << endl;
      t_1 = time_grid[a_indx+1];
      t_2 = time_grid[a_indx];
      a_1 = (ages4-t_2)/(t_1-t_2);
      a_2 = 1-a_1;
      f_1 = full_table[a_indx+1];
      f_2 = full_table[a_indx];
    }
    
    for (int i=0; i<1221; i++){
      spe.push_back((a_1*f_1[i])+ (a_2*f_2[i]));//*ergsa)) su entrambi;
    }
  }
  
}




void SEDcb16_extract_spec(std::vector <vector<double> > &full_table, std::vector<double> &time_grid, int a_indx, float ages4, std::vector<long double> &spe){
  
//Interpolating between two spectra of contiguous ages not needed with lines

  double ergsa = 3.9e+33;
  
  if(ages4<=0.){
    for (int i=0; i<13391; i++){
      spe.push_back(full_table[0][i]*ergsa);
    }
  } else if(ages4>=20.){
    for (int i=0; i<13391; i++){
      spe.push_back(full_table[220][i]*ergsa);
    }
  }
  else if (ages4>0. & ages4<20.){
    for (int i=0; i<13391; i++){
      spe.push_back(full_table[a_indx][i]*ergsa);
    }
  }
}

