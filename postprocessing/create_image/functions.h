#include <string>
#include <vector>
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <dirent.h>
#include <map>
#include <valarray>
using namespace std;

// function declarations
vector<int> inverseMap_sh_idv3(vector<int> &,vector<int> &,  vector<int> &, bool );

vector<int> searchsorted(vector<int> &,vector<int> &);

vector<int>vec_diff(vector<int>&);

vector<int>vec_sum(vector<int> &,vector<int> &);

vector<double>vec_diff(vector<double>&);

double getY(std:: vector<double> , std:: vector<double> ,double );

template <class T> int locate(const std::vector<T> &, const T );

void getPolar(double, double, double, double *, double *, double *);

class FileStructIterator {
public:
    explicit FileStructIterator(const std::string& filename);

    std::string getNext();

private:
    std::string m_filename;

    enum State { STATE_START, STATE_LOOP1, STATE_LOOP2 };
    State m_state;

    size_t m_badLineNumber;
    size_t m_i;
};

void read_bc03_ssp(string, std::vector<double> &,  std::vector <vector<double> > &, std::vector<long double> &);

void read_bc03_mass(std::string, std::vector <double> &);

//void read_bc03_ssp(string, std::vector<double> &,  std::vector <vector<double> > &, std::vector<long double> &);

int index_closest(std::vector<long double>::iterator, std::vector<long double>::iterator, long double);

int index_closest(std::vector<double>::iterator, std::vector<double>::iterator, double);

int index_closest(std::vector<float>::iterator, std::vector<float>::iterator, float);

void readParameters(float* fov,float* res,
                    std::string* filter,
                    std::string* module, std::string* rdir);

bool allNegOne(const std::vector<float>& vec);

// conversion: double or int -> string                                   
static const char fINT[] = "%i";
static const char fLONG[] = "%lli";
static const char fDP0[] = "%1.0f";
static const char fDP1[] = "%2.1f";
static const char fDP2[] = "%3.2f";
static const char fDP3[] = "%4.3f";
static const char fDP4[] = "%5.4f";
static const char fDP5[] = "%6.5f";
static const char ee3[] = "%4.3e";

template <class T> string conv (T &val, const char *fact)
{
  char VAL[20]; sprintf (VAL, fact, val);
  return string(VAL);
}

template <typename T>
std::vector<T> m_col_add(std::vector<std::vector<T>> const& mat) {
    std::vector<T> res;
    const auto column_size = mat[0].size();
    for (size_t x = 0; x < column_size; ++x)
        res.push_back(std::accumulate(mat.begin(), mat.end(), T{}, [x](T const& a, std::vector<T> const& row) {
            return a + row[x];
        }));
    return res;
}


int countHDF5Files(const std::string& path, const std::string& extension);


void readSSPTables(
    std::string ,
    const std::string& ,        
    const std::string& ,          
    int ,           
    std::vector<long double>&, 
    std::vector<double>& ,    
    std::vector<std::vector<double>>&  );

std::string getCB16MetallicityCode(int);

std::string getBC03MetallicityCode(int);



///SED
void SEDbc03_interp_2spec(std::vector<std::vector<double> > &full_table, std::vector<double> &time_grid,  int a_indx, float ages4, std::vector<long double> &spe);

void SEDcb16_extract_spec(std::vector<std::vector<double> > &full_table, std::vector<double> &time_grid, int a_indx, float ages4, std::vector<long double> &spe);

//make image functions
valarray<float> gridist_nok(vector<float>, vector<float>, vector<float>,  unsigned long int);
valarray<float> ncounts(vector<float>, vector<float>, int);
valarray<float> rescalemap(std::valarray<float>,unsigned long int,unsigned long int);
