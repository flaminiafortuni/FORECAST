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
using namespace std;

vector<int> inverseMap_sh_idv3(vector<int> &,vector<int> &,  vector<int> &, bool);

vector<int> searchsorted(vector<int> &,vector<int> &);

vector<int>vec_diff(vector<int>&);

vector<int>vec_sum(vector<int> &,vector<int> &);

double getY(std::vector<double> , std::vector<double>,double);

void getPolar(double, double, double, double *, double *, double *);

template <class T> int locate(const std::vector<T> &, const T );

int countHDF5Files(const std::string& path, const std::string& extension);

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

void readParameters(double *, double *, double *,double *, 
                    string *,string *,string *,string *,
                    string *, string *,
                    long *, long *, long *, string *);
